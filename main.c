#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <sys/times.h>
#include <string.h>
#include "cubic.h"
#include "nausparse.h"

/* Nauty worksize */
#define WORKSIZE 50 * MAXM

/* Variables for nauty */
int lab[MAXN], ptn[MAXN], orbits[MAXN];
static DEFAULTOPTIONS_SPARSEGRAPH(options);
statsblk stats;
setword workspace[WORKSIZE];

//For the timing
#define time_factor sysconf(_SC_CLK_TCK)

sparsegraph sg; /* Sparse graph datastructure for nauty */
sparsegraph sg_canon; /* Sparse graph datastructure for nauty */

int cubicIndex = 0;
Graph **cubicGraphs = NULL;

// Nauty and Sparse Graph Functions

void init_nauty_options() {
    //options.getcanon = TRUE;
    //options.userautomproc = save_generators;

    /* Init the nauty datastructures */
    SG_INIT(sg);
    SG_ALLOC(sg, MAX_VERTICES, 3 * MAX_VERTICES, "malloc");

    //sg.v and sg.d only have to be set once
    int i;
    for (i = 0; i < MAX_VERTICES; i++) {
        sg.v[i] = i * MAX_DEGREE;
        sg.d[i] = 3;
    }

    SG_INIT(sg_canon);
    SG_ALLOC(sg_canon, MAX_VERTICES, 3 * MAX_VERTICES, "malloc");
}

/**
 * Copy a graph represented by an AdjList to a sparsegraph struct for nauty representation.
 *
 * @param graph The input graph represented by an AdjList.
 */
void copy_sparse_graph(Graph *graph) {

    sg.nv = graph->current_nodes;
    sg.nde = 3 * graph->current_nodes;

    int i, j;
    for (i = 0; i < graph->V; i++) {
        //These values were already set in init_nauty_options()
        //sg.v[i] = i * REG;
        //sg.d[i] = 3;
        struct AdjListNode *adjListNode = graph->nodes[i].head;
        for (j = 0; j < 3 && adjListNode != NULL; j++) {
            sg.e[i * 3 + j] = adjListNode->label;
            adjListNode = adjListNode->next;
        }
    }
}

/**
 * Print a sparse graph in the nauty format.
 *
 * @param sparse_graph The sparsegraph struct representing the graph to be printed.
 */
void print_sparse_graph_nauty(sparsegraph sparse_graph) {
    int i, j;
    fprintf(stderr, "Printing sparse graph nauty:\n");
    for (i = 0; i < sparse_graph.nv; i++) {
        fprintf(stderr, "%d :", i);
        for (j = 0; j < sparse_graph.d[i]; j++) {
            //fprintf(stderr, " %d", sparse_graph.e[i * REG + j]);
            fprintf(stderr, " %d", sparse_graph.e[sparse_graph.v[i] + j]);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "Number of directed edges: %lu\n", (unsigned long) sparse_graph.nde);
}

/* Basic Graph Operations */

/**
 * Counts the degree of a specific vertex in a cubic graph.
 *
 * This function calculates the degree (number of adjacent vertices) of a given vertex in a cubic graph.
 * A cubic graph is a graph in which every vertex has exactly three neighbors.
 *
 * @param graph A pointer to the Graph structure representing the cubic graph.
 * @param label The label of the vertex whose degree needs to be counted.
 * @return An integer representing the degree of the specified vertex. If the vertex does not exist in the graph,
 *         or if the graph parameter is NULL, the function returns 0.
 */
int countVertexDegree(Graph *graph, int label) {
    int count = 0;
    struct AdjListNode *nodeDest = graph->nodes[label].head;
    while (nodeDest != NULL) {
        count++;
        nodeDest = nodeDest->next;
    }
    return count;
}

/**
 * Removes a specific node from the adjacency list of a vertex in a cubic graph, if present.
 *
 * This function removes a specified node with the label 'dest' from the adjacency list of the vertex 'src'
 * in a cubic graph. If the node is found and removed, its memory is freed. The function assumes that both
 * 'src' and 'dest' are valid vertices in the graph.
 *
 * @param graph A pointer to the Graph structure representing the cubic graph.
 * @param src The label of the source vertex from which to remove the node.
 * @param dest The label of the node to be removed from the adjacency list.
 */
void removeNodeFromAdjList(Graph *graph, int src, int dest) {
    struct AdjListNode *nodeDest = graph->nodes[src].head;
    struct AdjListNode *placeHolder = NULL;
    if (nodeDest->label == dest) {
        placeHolder = nodeDest;
        graph->nodes[src].head = nodeDest->next;
        free(placeHolder);
        return; //It wasn't there before
    }
    while (nodeDest->next != NULL) {
        if (nodeDest->next->label == dest) {
            placeHolder = nodeDest->next;
            nodeDest->next = nodeDest->next->next;
            free(placeHolder);
            return; //There was a break before
        }
        nodeDest = nodeDest->next;
    }
}

/**
 * Adds an undirected edge between two vertices in a cubic graph, if possible.
 *
 * This function attempts to add an undirected edge between two vertices, 'src' and 'dest', in a cubic graph.
 * It checks if the graph has not reached its maximum degree constraint (MAX_DEGREE) for either 'src' or 'dest'.
 * If either vertex has reached the maximum degree or if the edge already exists, the addition is skipped.
 *
 * @param graph A pointer to the Graph structure representing the cubic graph.
 * @param src The label of the source vertex.
 * @param dest The label of the destination vertex.
 */
void addEdge(Graph *graph, int src, int dest) {
    // Check
    if (countVertexDegree(graph, src) >= MAX_DEGREE || countVertexDegree(graph, dest) >= MAX_DEGREE)
        return;
    if (graph->nodes[src].head == NULL)
        graph->current_nodes++;
    if (graph->nodes[dest].head == NULL)
        graph->current_nodes++;

    //Add Edge to AdjList
    struct AdjListNode *nodeDest = graph->nodes[src].head;
    while (nodeDest) {
        if (nodeDest->label == dest)
            return;
        nodeDest = nodeDest->next;
    }

    struct AdjListNode *newNode = newAdjListNode(dest);
    newNode->next = graph->nodes[src].head;
    graph->nodes[src].head = newNode;

    newNode = newAdjListNode(src);
    newNode->next = graph->nodes[dest].head;
    graph->nodes[dest].head = newNode;

    //Add Edge to EdgeList
    EdgeListNode *edge = newEdgeListNode(src, dest);
    edge->next = graph->edges->head;
    graph->edges->head = edge;

    //Update total edges
    graph->current_edges++;
}

/**
 * Removes an undirected edge between two vertices in a cubic graph, if present.
 *
 * This function removes an undirected edge between 'src' and 'dest' vertices in a cubic graph, if the edge exists.
 * It removes both occurrences of the edge from the adjacency lists of both vertices and the edge list.
 * The function assumes that 'src' and 'dest' are valid vertices in the graph.
 *
 * @param graph A pointer to the Graph structure representing the cubic graph.
 * @param src The label of one endpoint of the edge.
 * @param dest The label of the other endpoint of the edge.
 */
void removeEdge(Graph *graph, int src, int dest) {

    //Remove Edge from AdjList
    removeNodeFromAdjList(graph, src, dest);
    removeNodeFromAdjList(graph, dest, src);

    //Remove Edge from edgeList
    removeEdgeFromEdgeList(graph->edges, src, dest);
    removeEdgeFromEdgeList(graph->edges, dest, src);

    //Update current edges
    graph->current_edges--;
}

/**
 * Removes an edge node from the edge list, if present.
 *
 * This function removes an edge node with the given source 'src' and destination 'dest' from the edge list
 * represented by 'edgeList', if the edge node exists in the list. The memory occupied by the removed edge node
 * is freed. The function assumes that the 'edgeList' is a valid EdgeList structure.
 *
 * @param edgeList A pointer to the EdgeList structure representing the list of edges.
 * @param src The label of the source vertex of the edge.
 * @param dest The label of the destination vertex of the edge.
 */
void removeEdgeFromEdgeList(EdgeList *edgeList, int src, int dest) {
    EdgeListNode *edgeNode = edgeList->head;
    EdgeListNode *placeHolder = NULL;
    if (edgeNode->edge->src == src && edgeNode->edge->dest == dest) {
        placeHolder = edgeNode;
        edgeList->head = edgeNode->next;
        //free( placeHolder->edge);
        free(placeHolder);
        return; //It wasn't there before
    }
    while (edgeNode->next != NULL) {
        if (edgeNode->next->edge->src == src && edgeNode->next->edge->dest == dest) {
            EdgeListNode *placeHolder2 = edgeNode->next;
            edgeNode->next = edgeNode->next->next;
            //free(placeHolder2->edge);
            free(placeHolder2);
            return; // There was a break before
        }
        edgeNode = edgeNode->next;
    }
}

/**
 * Inserts a new edge into the edge list if it does not already exist.
 *
 * This function inserts a new edge, with the source 'src' and destination 'dest', into the edge list
 * represented by 'edgeList', but only if the edge does not already exist in the list. The function first
 * verifies the non-existence of the edge using the 'verifyEdgeExistence' function. If the edge is not found,
 * a new edge node is created and added to the front of the edge list.
 *
 * @param edgeList A pointer to the EdgeList structure representing the list of edges.
 * @param src The label of the source vertex of the new edge.
 * @param dest The label of the destination vertex of the new edge.
 */
void insertInEdgeList(EdgeList *edgeList, int src, int dest) {
    if (!verifyEdgeExistence(src, dest, edgeList)) {
        EdgeListNode *edge = newEdgeListNode(src, dest);
        edge->next = edgeList->head;
        edgeList->head = edge;
    }
}

/**
 * Removes the first edge node from the edge list, if the list is not empty.
 *
 * This function removes the first edge node from the edge list represented by 'edgeList',
 * if the list is not empty. The memory occupied by the removed edge node is freed.
 * The function assumes that 'edgeList' is a valid EdgeList structure.
 *
 * @param edgeList A pointer to the EdgeList structure representing the list of edges.
 */
void popFromEdgeList(EdgeList *edgeList) {
    EdgeListNode *headEdgeListNode = edgeList->head;
    if (headEdgeListNode != NULL) {
        edgeList->head = headEdgeListNode->next;
        //free(headEdgeListNode->edge);
        free(headEdgeListNode);
        return;
    }
    return;
}

/**
 * Compares two edges for equality.
 *
 * This function compares two edges 'e1' and 'e2' to check if they are equal or represent the same edge,
 * regardless of the directionality. It returns 1 if the edges are equal or represent the same edge,
 * and 0 otherwise.
 *
 * @param e1 A pointer to the first Edge structure to be compared.
 * @param e2 A pointer to the second Edge structure to be compared.
 * @return 1 if the edges are equal or represent the same edge, 0 otherwise.
 */
int compareEdge(Edge *e1, Edge *e2) {
    return ((e1->src == e2->src && e1->dest == e2->dest) || (e1->src == e2->dest && e1->dest == e2->src));
}

/**
 * Verify if an edge exists in the given edge list.
 *
 * This function checks if an edge with the given source 'src' and destination 'dest' exists in the specified edge list.
 * It traverses the edge list and compares each edge's source and destination nodes with the provided values.
 *
 * @param src The source node of the edge to be checked.
 * @param dest The destination node of the edge to be checked.
 * @param edgeList A pointer to the EdgeList where the edge existence will be verified.
 * @return 1 if the edge (src, dest) or (dest, src) is found in the edge list, 0 otherwise.
 */
int verifyEdgeExistence(int src, int dest, EdgeList *edgeList) {
    struct EdgeListNode *edgeNode = edgeList->head;
    while (edgeNode != NULL) {
        if ((edgeNode->edge->src == src && edgeNode->edge->dest == dest) ||
            (edgeNode->edge->src == dest && edgeNode->edge->dest == src))
            return 1;
        edgeNode = edgeNode->next;
    }
    return 0;
}

/* Graph Manipulation */

/**
 * Generates a new edge and two nodes to connect to existing vertices in the graph.
 *
 * This function generates a new edge and two new nodes to define it in the graph 'graph'.
 * It also connects one of the new nodes to vertices 'src1' and 'dest1', and the other to vertices 'src2' and 'dest2'.
 * Before connecting the new nodes, it removes the edges (src1, dest1) and (src2, dest2) from the graph.
 * The new nodes are created using spare vertices from the graph, which must have at least 2 spare vertices available.
 * If there are not enough spare vertices (less than 2), the function returns 0 to indicate the operation failed.
 * Otherwise, it successfully adds the new edge and nodes and returns 1.
 *
 * @param graph A pointer to the Graph structure representing the graph.
 * @param src1 The source vertex of the first edge to be removed.
 * @param dest1 The destination vertex of the first edge to be removed.
 * @param src2 The source vertex of the second edge to be removed.
 * @param dest2 The destination vertex of the second edge to be removed.
 * @return 1 if the operation is successful, 0 if there are not enough spare vertices in the graph.
 */
int edgeInsertion(Graph *graph, int src1, int dest1, int src2, int dest2) {
    int size = 2;
    int spareVertices[size];
    int retrVertices = getSpareVertices(graph, spareVertices, size);
    if (retrVertices < 2) {
        return 0;
    }

    removeEdge(graph, src1, dest1);
    removeEdge(graph, src2, dest2);
    addEdge(graph, spareVertices[0], spareVertices[1]);
    addEdge(graph, spareVertices[0], src1);
    addEdge(graph, spareVertices[0], dest1);
    addEdge(graph, spareVertices[1], src2);
    addEdge(graph, spareVertices[1], dest2);

    return 1;
}

/**
 * Generates a K4- and connects it to existing vertices in the graph.
 *
 * This function generates a K4- (K4 with one edge removed) using six spare vertices from the graph 'graph'.
 * It then connects the two degree-2 nodes of the K4- structure to two new nodes.
 * The new nodes are connected to vertices 'src1' and 'dest1', and 'src2' and 'dest2' of the given graph.
 * Before connecting the new nodes, it removes the edges (src1, dest1) and (src2, dest2) from the graph.
 * The operation requires at least six spare vertices in the graph.
 * If there are not enough spare vertices, the function returns 0 to indicate the operation failed.
 * Otherwise, it successfully adds the K4- structure and connects it to the specified vertices, and returns 1.
 *
 * @param graph A pointer to the Graph structure representing the graph.
 * @param src1 The source vertex of the first edge to be removed.
 * @param dest1 The destination vertex of the first edge to be removed.
 * @param src2 The source vertex of the second edge to be removed.
 * @param dest2 The destination vertex of the second edge to be removed.
 * @return 1 if the operation is successful, 0 if there are not enough spare vertices in the graph.
 */
int nonAdjDiamondInsertion(Graph *graph, int src1, int dest1, int src2, int dest2) {
    int size = 6;
    int spareVertices[size];
    int retrVertices = getSpareVertices(graph, spareVertices, size);
    if (retrVertices < size) {
        return 0;
    }
    // Create K4-
    addEdge(graph, spareVertices[0], spareVertices[1]);
    addEdge(graph, spareVertices[1], spareVertices[2]);
    addEdge(graph, spareVertices[2], spareVertices[3]);
    addEdge(graph, spareVertices[3], spareVertices[0]);
    addEdge(graph, spareVertices[0], spareVertices[2]);

    addEdge(graph, spareVertices[1], spareVertices[4]);
    addEdge(graph, spareVertices[3], spareVertices[5]);

    removeEdge(graph, src1, dest1);
    removeEdge(graph, src2, dest2);
    addEdge(graph, src1, spareVertices[4]);
    addEdge(graph, spareVertices[4], dest1);
    addEdge(graph, src2, spareVertices[5]);
    addEdge(graph, spareVertices[5], dest2);

    return 1;
}

/**
* Generates a K4- and connects it to existing vertices in the graph.
*
* This function generates a K4- (K4 with one edge removed) using four spare vertices from the graph 'graph'.
* It then connects the two degree-2 nodes of the K4- structure to the nodes 'src' and 'dest' of the given graph.
* Before connecting the new nodes, it removes the edge (src, dest) from the graph.
* The operation requires at least four spare vertices in the graph.
* If there are not enough spare vertices, the function returns 0 to indicate the operation failed.
* Otherwise, it successfully adds the K4- structure and connects it to the specified vertices, and returns 1.
*
* @param graph A pointer to the Graph structure representing the graph.
* @param src The source vertex of the edge to be removed.
* @param dest The destination vertex of the edge to be removed.
* @return 1 if the operation is successful, 0 if there are not enough spare vertices in the graph.
*/
int adjDiamondInsertion(Graph *graph, int src, int dest) {
    int size = 4;
    int spareVertices[size];
    int retrVertices = getSpareVertices(graph, spareVertices, size);
    if (retrVertices < size) {
        return 0;
    }
    // Create K4-
    addEdge(graph, spareVertices[0], spareVertices[1]);
    addEdge(graph, spareVertices[1], spareVertices[2]);
    addEdge(graph, spareVertices[2], spareVertices[3]);
    addEdge(graph, spareVertices[3], spareVertices[0]);
    addEdge(graph, spareVertices[0], spareVertices[2]);

    removeEdge(graph, src, dest);
    addEdge(graph, src, spareVertices[3]);
    addEdge(graph, spareVertices[1], dest);

    return 1;
}

/**
 * Generates a K4+ and connects it to existing vertices in the graph.
 *
 * This function generates a K4+ (K4 with one edge removed and a new node connected to the endpoints of the removed edge)
 * using six spare vertices from the graph 'graph'.
 * It then connects the degree-2 node of the K4+ structure to a new node, which is then connected to vertices 'src' and 'dest'
 * of the given graph.
 * Before connecting the new nodes, it removes the edge (src, dest) from the graph.
 * The operation requires at least six spare vertices in the graph.
 * If there are not enough spare vertices, the function returns 0 to indicate the operation failed.
 * Otherwise, it successfully adds the K4+ structure and connects it to the specified vertices, and returns 1.
 *
 * @param graph A pointer to the Graph structure representing the graph.
 * @param src The source vertex of the edge to be removed.
 * @param dest The destination vertex of the edge to be removed.
 * @return 1 if the operation is successful, 0 if there are not enough spare vertices in the graph.
 */
int lollipopInsertion(Graph *graph, int src, int dest) {
    int size = 6;
    int spareVertices[size];
    int retrVertices = getSpareVertices(graph, spareVertices, size);
    if (retrVertices < size) {
        return 0;
    }

    // Create K4+
    addEdge(graph, spareVertices[0], spareVertices[1]);
    addEdge(graph, spareVertices[1], spareVertices[2]);
    addEdge(graph, spareVertices[3], spareVertices[0]);
    addEdge(graph, spareVertices[0], spareVertices[2]);
    addEdge(graph, spareVertices[3], spareVertices[1]);
    addEdge(graph, spareVertices[3], spareVertices[4]);
    addEdge(graph, spareVertices[2], spareVertices[4]);

    addEdge(graph, spareVertices[4], spareVertices[5]);
    removeEdge(graph, src, dest);
    addEdge(graph, src, spareVertices[5]);
    addEdge(graph, spareVertices[5], dest);

    return 1;
}

/* Irreducibility Checks */

/**
 * Performs Depth-First Search (DFS) traversal for finding bridges in a graph.
 *
 * This function performs a Depth-First Search (DFS) traversal starting from the vertex 'start' in the given graph 'graph'.
 * It uses an array 'visited' to mark visited vertices during the traversal.
 * The DFS traversal explores the edges of the graph to detect bridges.
 * Bridges are the edges whose removal would increase the number of connected components in the graph.
 * During the DFS, if an unvisited edge leads to an unvisited vertex, it recursively calls itself for that vertex.
 *
 * @param graph A pointer to the Graph structure representing the graph.
 * @param start The label of the starting vertex for DFS traversal.
 * @param visited An array to mark visited vertices during DFS traversal.
 */
void dfsForBridge(Graph *graph, int start, int *visited) {
    visited[start] = 1;
    EdgeListNode *edgeListNode = graph->edges->head;
    while (edgeListNode) {
        if (edgeListNode->edge->src == start && !visited[edgeListNode->edge->dest])
            dfsForBridge(graph, edgeListNode->edge->dest, visited);
        if (edgeListNode->edge->dest == start && !visited[edgeListNode->edge->src])
            dfsForBridge(graph, edgeListNode->edge->src, visited);
        edgeListNode = edgeListNode->next;
    }
    /*struct AdjListNode* nodeDest = graph->nodes[start].head;
    while(nodeDest != NULL){
        if(!visited[nodeDest->label])
            dfsForBridge(graph, nodeDest->label, visited);
        nodeDest = nodeDest->next;
    }*/
}

/**
 * Checks if an Edge is a bridge in the graph.
 *
 * This function checks if the given Edge 'edge' is a bridge in the graph 'graph'.
 * It temporarily removes the edge from the graph, performs Depth-First Search (DFS) traversal using 'dfsForBridge' function
 * to find bridges, and then restores the edge in the graph.
 * A bridge is an edge whose removal would increase the number of connected components in the graph.
 * If the destination vertex of the edge is not visited after the temporary removal, the function returns 1,
 * indicating that the edge is a bridge. Otherwise, it returns 0, indicating that it is not a bridge.
 *
 * @param graph A pointer to the Graph structure representing the graph.
 * @param edge A pointer to the Edge structure to be checked for being a bridge.
 * @return 1 if 'edge' is a bridge in the graph, 0 otherwise.
 */
int isABridge(Graph *graph, Edge *edge) {
    removeEdge(graph, edge->src, edge->dest);
    int visited[graph->V];
    for (int i = 0; i < (graph->V); i++) {
        visited[i] = 0;
    }
    dfsForBridge(graph, edge->src, visited);

    addEdge(graph, edge->src, edge->dest);
    return !visited[edge->dest];
}

// Unused function
int isABridgeList(Graph *graph) {
    struct EdgeListNode *edgeNode = graph->edges->head;
    while (edgeNode != NULL) {
        if (!isABridge(graph, edgeNode->edge))
            return 0;
        edgeNode = edgeNode->next;
    }
    return 1;
}

/**
 * Checks if one of the edge endpoints is part of a triangle, and the edge isn't.
 *
 * This function checks if one of the endpoints of the given 'edge' is part of a triangle in the graph,
 * while the 'edge' itself is not part of any triangle. It returns 1 if the condition is met,
 * indicating that the edge is reducible, and 0 otherwise, indicating that the edge is irreducible.
 *
 * @param graph A pointer to the Graph structure representing the graph.
 * @param edge A pointer to the Edge structure representing the edge to be checked.
 * @return 1 if the edge is reducible, 0 if the edge is irreducible.
 */
int irreducibilityCondition2(Graph *graph, Edge *edge) {
    TriangleList *triangleList = graph->triangles;
    TriangleListNode *triangleNode = triangleList->head;
    int src = edge->src;
    int dest = edge->dest;

    while (triangleNode) {
        if (!compareEdge(edge, triangleNode->triangle->e1) && !compareEdge(edge, triangleNode->triangle->e2) &&
            !compareEdge(edge, triangleNode->triangle->e3))
            if (src == triangleNode->triangle->e1->src || src == triangleNode->triangle->e1->dest ||
                src == triangleNode->triangle->e2->src || src == triangleNode->triangle->e2->dest ||
                src == triangleNode->triangle->e3->src || src == triangleNode->triangle->e3->dest ||
                dest == triangleNode->triangle->e1->src || dest == triangleNode->triangle->e1->dest ||
                dest == triangleNode->triangle->e2->src || dest == triangleNode->triangle->e2->dest ||
                dest == triangleNode->triangle->e3->src || dest == triangleNode->triangle->e3->dest)
                return 1;

        triangleNode = triangleNode->next;

    }
    return 0;

}

/**
 * Checks if both endpoints of the edge are part of the same 4-gon, and the edge isn't.
 *
 * This function checks if both endpoints of the given 'edge' are part of the same 4-gon (tetragon) in the graph,
 * while the 'edge' itself is not part of the 4-gon. It returns 1 if the condition is met,
 * indicating that the edge is reducible, and 0 otherwise, indicating that the edge is irreducible.
 *
 * @param graph A pointer to the Graph structure representing the graph.
 * @param edge A pointer to the Edge structure representing the edge to be checked.
 * @return 1 if the edge is reducible, 0 if the edge is irreducible.
 */
int irreducibilityCondition3(Graph *graph, Edge *edge) {
    TetraList *tList = graph->tetragons;
    int src = edge->src;
    int dest = edge->dest;
    TetraListNode *tNode = tList->head;
    while (tNode) {
        if (!compareEdge(edge, tNode->tetragon->e1) && !compareEdge(edge, tNode->tetragon->e2) &&
            !compareEdge(edge, tNode->tetragon->e3) && !compareEdge(edge, tNode->tetragon->e4))
            if ((src == tNode->tetragon->e1->src || src == tNode->tetragon->e1->dest ||
                 src == tNode->tetragon->e2->src || src == tNode->tetragon->e2->dest ||
                 src == tNode->tetragon->e3->src || src == tNode->tetragon->e3->dest ||
                 src == tNode->tetragon->e4->src || src == tNode->tetragon->e4->dest) &&
                (dest == tNode->tetragon->e1->src || dest == tNode->tetragon->e1->dest ||
                 dest == tNode->tetragon->e2->src || dest == tNode->tetragon->e2->dest ||
                 dest == tNode->tetragon->e3->src || dest == tNode->tetragon->e3->dest ||
                 dest == tNode->tetragon->e4->src || dest == tNode->tetragon->e4->dest))
                return 1;

        tNode = tNode->next;
    }
    return 0;

}

/**
 * Checks if the given graph is an irreducible graph.
 *
 * This function checks if the given graph 'g' is an irreducible graph by analyzing its edges and triangle/tetragon lists.
 * It checks the three irreducibility conditions: whether any edge is a bridge ('isABridge'), whether there exist edges that satisfy
 * 'irreducibilityCondition2', and whether there exist edges that satisfy 'irreducibilityCondition3'.
 *
 * @param g A pointer to the Graph structure representing the graph to be checked.
 * @return 1 if the graph is irreducible, 0 if the graph is reducible.
 */
int isAIrreducibleGraph(Graph *g) {
    TriangleList *trianglelist = (TriangleList *) malloc(sizeof(TriangleList));

    trianglelist->head = NULL;
    findTriangles(g, 0, trianglelist);
    insertTriangles(g, trianglelist);
    TetraList *tetralist = (TetraList *) malloc(sizeof(TetraList));

    tetralist->head = NULL;
    findTetragons(g, 0, tetralist);
    insertTetragons(g, tetralist);

    EdgeListNode *edgeListNode = g->edges->head;
    EdgeListNode *edgeListNodeSucc = NULL;

    while (edgeListNode) {
        edgeListNodeSucc = edgeListNode->next;
        Edge *edge = edgeListNode->edge;

        int b = isABridge(g, edge);
        //edge = edgeListNode->edge;
        int irr2 = irreducibilityCondition2(g, edge);
        int irr3 = irreducibilityCondition3(g, edge);

        if (!(b || irr2 || irr3))
            return 0;

        edgeListNode = edgeListNodeSucc;
    }
    return 1;
}

/**
 * Initializes the given graph as an irreducible graph.
 *
 * This function initializes the given graph 'graph' as an irreducible graph by creating a K4 structure using spare vertices.
 *
 * @param graph A pointer to the Graph structure to be initialized as an irreducible graph.
 */
void initIrreducibleGraph(Graph *graph) {

    int size = 4;
    int spareVertices[size];
    int retrVertices = getSpareVertices(graph, spareVertices, size);
    if (retrVertices < size) {
        return;
    }

    addEdge(graph, spareVertices[0], spareVertices[1]);
    addEdge(graph, spareVertices[1], spareVertices[2]);
    addEdge(graph, spareVertices[2], spareVertices[3]);
    addEdge(graph, spareVertices[3], spareVertices[0]);
    addEdge(graph, spareVertices[0], spareVertices[2]);
    addEdge(graph, spareVertices[1], spareVertices[3]);

}

/* Triangle and Tetragon Operations */

/**
 * Perform Depth-First Search (DFS) to find cycles of length 'n' starting from the specified vertex 'start' in the graph 'g'.
 *
 * @param g The Graph on which the DFS is performed.
 * @param trianglelist The TriangleList to store found triangles.
 * @param tetraList The TetraList to store found tetragons.
 * @param pathEdgeList The current path (sequence of edges) being explored during the DFS.
 * @param visited An array that keeps track of visited vertices during the DFS.
 * @param n The desired length of the cycle to be found.
 * @param vert The current vertex being explored.
 * @param start The starting vertex of the cycle.
 * @param count A pointer to an integer that keeps track of the number of cycles found.
 * @param triangles A flag (1 or 0) indicating whether to look for triangles (1) or tetragons (0) in the graph.
 */
void
DFSforCycle(Graph *g, TriangleList *trianglelist, TetraList *tetraList, EdgeList *pathEdgeList, int *visited, int n,
            int vert, int start, int *count, int triangles) {

    // Mark the current vertex as visited
    visited[vert] = 1;

    // If the path of length (n-1) is found, check if it forms a cycle and save it if necessary
    if (n == 0) {

        // Mark vert as unvisited to make it usable again.
        visited[vert] = 0;

        // Check if the current path is a cycle by verifying that the last edge connects back to the starting vertex
        if (verifyEdgeExistence(vert, start, g->edges)) {

            EdgeListNode *edgeNode = newEdgeListNode(vert, start);
            edgeNode->next = pathEdgeList->head;
            pathEdgeList->head = edgeNode;
            *count = *count + 1;

            // Add the cycle (path) to the respective lists (trianglelist or tetraList)
            if (triangles) {

                EdgeListNode *edgeNode1 = pathEdgeList->head;
                TriangleListNode *node = newTriangleListNode(edgeNode1->edge, edgeNode1->next->edge,
                                                             edgeNode1->next->next->edge);

                /*node->next = trianglelist->head;
                trianglelist->head = node;*/
                insertNodeInTriangleList(node, trianglelist);
                popFromEdgeList(pathEdgeList);
                /*EdgeListNode *  ptr = pathEdgeList->head;
                pathEdgeList->head =pathEdgeList->head->next;
                free(ptr);*/
            } else {
                EdgeListNode *edgeNode1 = pathEdgeList->head;
                TetraListNode *node = newTetraListNode(edgeNode1->edge, edgeNode1->next->edge,
                                                       edgeNode1->next->next->edge, edgeNode1->next->next->next->edge);
                /*node->next = tetraList->head;
                tetraList->head = node;*/
                insertNodeInTetraList(tetraList, node);
/*
                EdgeListNode *  ptr = pathEdgeList->head;
                pathEdgeList->head =pathEdgeList->head->next;
                free(ptr);*/
                popFromEdgeList(pathEdgeList);
                /*EdgeListNode *  ptr = pathEdgeList->head;
                ptr = pathEdgeList->head;
                pathEdgeList->head =pathEdgeList->head->next;
                free(ptr);*/
            }
            return;
        } else
            return;
    }

    // For searching every possible path of length (n-1)
    for (int i = 0; i < g->current_nodes; i++)
        if (!visited[i] && verifyEdgeExistence(vert, i, g->edges)) {
            /* EdgeListNode* edgeNode = newEdgeListNode(i, vert);
             edgeNode->next = pathEdgeList->head;
             pathEdgeList->head= edgeNode;*/

            // Insert the current edge into the pathEdgeList and continue the DFS
            insertInEdgeList(pathEdgeList, vert, i);
            // Recursive call decreasing length by 1
            DFSforCycle(g, trianglelist, tetraList, pathEdgeList, visited, n - 1, i, start, count, triangles);
            // Backtrack by removing the last inserted edge from the pathEdgeList
            popFromEdgeList(pathEdgeList);
        }

    // Mark vert as unvisited to make it usable again.
    visited[vert] = 0;

}

/**
 * Find and store triangles in the graph 'g'.
 *
 * @param g The Graph in which triangles are to be found.
 * @param update A flag (1 or 0) indicating whether to update the 'trianglelist' or not.
 * @param trianglelist The TriangleList to store the found triangles.
 */
void findTriangles(Graph *g, int update, TriangleList *trianglelist) {
    if (!update) {
        EdgeList *edgeList = (EdgeList *) malloc(sizeof(EdgeList));
        int visited[g->V];
        int n = 3, count = 0, triangles = 1;

        for (int i = 0; i < g->V; i++)
            visited[i] = 0;

        // Start DFS from each vertex in the graph to find triangles
        for (int i = 0; i < g->current_nodes - (n - 1); i++) {
            DFSforCycle(g, trianglelist, NULL, edgeList, visited, n - 1, i, i, &count, triangles);
            visited[i] = 1;
        }
        free(edgeList);
    } else {
        // If 'update' is 1, it means that 'trianglelist' is already updated, so we do nothing and return.
        return;
    }
}

/**
 * Find and enumerate tetragons in the graph 'g'.
 * The tetragons found are added to the 'tetraList'.
 *
 * @param g The Graph in which to find tetragons.
 * @param update A flag to indicate whether to update the 'tetraList' or not.
 *              If update=1, it means that the 'tetraList' already contains some tetragons,
 *              and we don't need to find new ones, so the function will return without doing anything.
 * @param tetraList The TetraList to which the found tetragons are added.
 */
void findTetragons(Graph *g, int update, TetraList *tetraList) {

    if (!update) {
        EdgeList *edgeList = (EdgeList *) malloc(sizeof(EdgeList));
        int visited[g->V];
        int n = 4, count = 0, tetra = 0;

        for (int i = 0; i < g->V; i++)
            visited[i] = 0;

        // Perform DFS for each vertex to find cycles of length (n-1)
        for (int i = 0; i < g->current_nodes - (n - 1); i++) {
            DFSforCycle(g, NULL, tetraList, edgeList, visited, n - 1, i, i, &count, tetra);
            visited[i] = 1;
        }
        free(edgeList);

    } else {
        // If the 'tetraList' is already updated, return without doing anything
        return;
    }

}

/**
 * Compares two triangles for equality.
 *
 * This function compares two triangles represented by 't1' and 't2' for equality. It checks if any of the four possible
 * combinations of edges between the triangles are equal. If any combination is found to be equal, the function returns 1,
 * indicating that the triangles are equal. Otherwise, it returns 0, indicating inequality.
 *
 * @param t1 A pointer to the first Triangle structure to be compared.
 * @param t2 A pointer to the second Triangle structure to be compared.
 * @return 1 if the triangles are equal (regardless of the order of edges), 0 otherwise.
 */
int compareTriangle(Triangle *t1, Triangle *t2) {
    return compareEdge(t1->e1, t2->e1) && compareEdge(t1->e2, t2->e2) && compareEdge(t1->e3, t2->e3) ||
           compareEdge(t1->e1, t2->e2) && compareEdge(t1->e2, t2->e1) && compareEdge(t1->e3, t2->e3) ||
           compareEdge(t1->e1, t2->e1) && compareEdge(t1->e2, t2->e3) && compareEdge(t1->e3, t2->e2) ||
           compareEdge(t1->e1, t2->e3) && compareEdge(t1->e2, t2->e2) && compareEdge(t1->e3, t2->e1);
}

/**
 * Compares two tetragons for equality.
 *
 * This function compares two tetragons represented by 't1' and 't2' for equality.
 * It checks if any of the seven possible combinations of edges between the tetragons are equal.
 * If any combination is found to be equal, the function returns 1, indicating that the tetragons are equal.
 * Otherwise, it returns 0, indicating inequality.
 *
 * @param t1 A pointer to the first Tetragon structure to be compared.
 * @param t2 A pointer to the second Tetragon structure to be compared.
 * @return 1 if the tetragons are equal (regardless of the order of edges), 0 otherwise.
 */
int compareTetragon(Tetragon *t1, Tetragon *t2) {
    return compareEdge(t1->e1, t2->e1) && compareEdge(t1->e2, t2->e2) && compareEdge(t1->e3, t2->e3) &&
           compareEdge(t1->e4, t2->e4) ||
           compareEdge(t1->e1, t2->e2) && compareEdge(t1->e2, t2->e1) && compareEdge(t1->e3, t2->e3) &&
           compareEdge(t1->e4, t2->e4) ||
           compareEdge(t1->e1, t2->e1) && compareEdge(t1->e2, t2->e3) && compareEdge(t1->e3, t2->e2) &&
           compareEdge(t1->e4, t2->e4) ||
           compareEdge(t1->e1, t2->e3) && compareEdge(t1->e2, t2->e2) && compareEdge(t1->e3, t2->e1) &&
           compareEdge(t1->e4, t2->e4) ||
           compareEdge(t1->e1, t2->e4) && compareEdge(t1->e2, t2->e2) && compareEdge(t1->e3, t2->e3) &&
           compareEdge(t1->e4, t2->e1) ||
           compareEdge(t1->e1, t2->e4) && compareEdge(t1->e2, t2->e3) && compareEdge(t1->e3, t2->e2) &&
           compareEdge(t1->e4, t2->e1) ||
           compareEdge(t1->e1, t2->e1) && compareEdge(t1->e2, t2->e4) && compareEdge(t1->e3, t2->e3) &&
           compareEdge(t1->e4, t2->e2);
}

/**
 * Checks if a Triangle is valid.
 *
 * This function checks if the given Triangle 't' is a valid triangle, i.e., it consists of three distinct vertices.
 * The function creates a boolean array 'v' to mark the occurrence of vertices in the triangle's edges.
 * If exactly three vertices are marked as present in the triangle, the function returns 1, indicating that it is a valid triangle.
 * Otherwise, it returns 0, indicating that it is not a valid triangle.
 *
 * @param t A pointer to the Triangle structure to be checked for validity.
 * @return 1 if 't' is a valid triangle, 0 otherwise.
 */
int isATriangle(Triangle *t) {
    int v[MAX_VERTICES];
    for (int i = 0; i < MAX_VERTICES; i++)
        v[i] = 0;
    v[t->e1->src] = 1;
    v[t->e1->dest] = 1;
    v[t->e2->src] = 1;
    v[t->e2->dest] = 1;
    v[t->e3->src] = 1;
    v[t->e3->dest] = 1;
    int count = 0;

    for (int i = 0; i < MAX_VERTICES; i++)
        count += v[i];
    if (count == 3)
        return 1;
    return 0;
}

/**
 * Checks if a Tetragon is valid.
 *
 * This function checks if the given Tetragon 't' is a valid tetragon, i.e., it consists of four distinct vertices.
 * The function creates a boolean array 'v' to mark the occurrence of vertices in the tetragon's edges.
 * If exactly four vertices are marked as present in the tetragon, the function returns 1, indicating that it is a valid tetragon.
 * Otherwise, it returns 0, indicating that it is not a valid tetragon.
 *
 * @param t A pointer to the Tetragon structure to be checked for validity.
 * @return 1 if 't' is a valid tetragon, 0 otherwise.
 */
int isATetragon(Tetragon *t) {
    int v[MAX_VERTICES];
    for (int i = 0; i < MAX_VERTICES; i++)
        v[i] = 0;
    v[t->e1->src] = 1;
    v[t->e1->dest] = 1;
    v[t->e2->src] = 1;
    v[t->e2->dest] = 1;
    v[t->e3->src] = 1;
    v[t->e3->dest] = 1;
    v[t->e4->src] = 1;
    v[t->e4->dest] = 1;
    int count = 0;
    for (int i = 0; i < MAX_VERTICES; i++)
        count += v[i];
    if (count == 4)
        return 1;
    return 0;
}

/**
 * Verifies if a Triangle exists in the TriangleList.
 *
 * This function checks if the given Triangle 'triangle' exists in the TriangleList represented by 'triangleList'.
 * It iterates through each Triangle in the list and compares 'triangle' for equality using 'isATriangle' and 'compareTriangle' functions.
 * If an equal Triangle is found, the function returns 1, indicating existence. Otherwise, it returns 0, indicating non-existence.
 *
 * @param triangle A pointer to the Triangle structure to be verified.
 * @param triangleList A pointer to the TriangleList structure representing the list of Triangles.
 * @return 1 if 'triangle' exists in the list, 0 otherwise.
 */
int verifyTriangleExistence(Triangle *triangle, TriangleList *triangleList) {
    TriangleListNode *triangleListNode = triangleList->head;
    while (triangleListNode != NULL) {
        if (isATriangle(triangle))
            if (compareTriangle(triangle, triangleListNode->triangle))
                return 1;
        triangleListNode = triangleListNode->next;
    }
    return 0;
}

/**
 * Verifies if a Tetragon exists in the TetraList.
 *
 * This function checks if the given Tetragon 'tetra' exists in the TetraList represented by 'tetraList'.
 * It iterates through each Tetragon in the list and compares 'tetra' for validity using 'isATetragon' function
 * and for equality using 'compareTetragon' function. If a valid and equal Tetragon is found, the function returns 1,
 * indicating existence. Otherwise, it returns 0, indicating non-existence.
 *
 * @param tetraList A pointer to the TetraList structure representing the list of Tetragons.
 * @param tetra A pointer to the Tetragon structure to be verified.
 * @return 1 if 'tetra' exists in the list, 0 otherwise.
 */
int verifyTetragonExistence(TetraList *tetraList, Tetragon *tetra) {
    TetraListNode *tetraListNode = tetraList->head;
    while (tetraListNode != NULL) {
        if (isATetragon(tetra))
            if (compareTetragon(tetra, tetraListNode->tetragon))
                return 1;
        tetraListNode = tetraListNode->next;
    }
    return 0;
}

/**
 * Insert triangles from 'tList' into the graph 'g'.
 * The edges in the graph that are part of a triangle are updated to point to the corresponding triangle.
 *
 * @param g The Graph where triangles are to be inserted.
 * @param tList The TriangleList containing the triangles to be inserted.
 */
void insertTriangles(Graph *g, TriangleList *tList) {

    g->triangles = tList;
    EdgeListNode *node = g->edges->head;
    TriangleListNode *triangleNode = tList->head;

    while (triangleNode) {

        Triangle *t = triangleNode->triangle;
        while (node) {
            // Check if the edge is part of the current triangle 't'
            if (compareEdge(node->edge, t->e1) || compareEdge(node->edge, t->e2) || compareEdge(node->edge, t->e3))
                node->triangle = t;

            node = node->next;
        }
        triangleNode = triangleNode->next;
    }
}

/**
 * Insert tetragons from 'tList' into the graph 'g'.
 * The edges in the graph that are part of a tetragon are updated to point to the corresponding tetragon.
 *
 * @param g The Graph where tetragons are to be inserted.
 * @param tList The TetraList containing the tetragons to be inserted.
 */
void insertTetragons(Graph *g, TetraList *tList) {
    g->tetragons = tList;
    EdgeListNode *node = g->edges->head;
    TetraListNode *tetraNode = tList->head;

    while (tetraNode) {

        while (node) {
            Tetragon *t = tetraNode->tetragon;
            // Check if the edge is part of the current tetragon 't'
            if (compareEdge(node->edge, t->e1) || compareEdge(node->edge, t->e2) || compareEdge(node->edge, t->e3) ||
                compareEdge(node->edge, t->e4))
                if (node->tetragon == NULL)
                    node->tetragon = t;

            node = node->next;
        }
        tetraNode = tetraNode->next;
    }
}

/**
 * Inserts a TriangleListNode into the TriangleList if the Triangle does not already exist in the list.
 *
 * This function inserts a TriangleListNode 'node' into the TriangleList represented by 'triangleList',
 * but only if the Triangle in the node does not already exist in the list. It first checks if 'triangleList'
 * is not NULL, and then uses the 'verifyTriangleExistence' function to check for the existence of the Triangle.
 * If the Triangle is not found in the list, the node is added to the front of the list.
 *
 * @param node A pointer to the TriangleListNode structure to be inserted into the list.
 * @param triangleList A pointer to the TriangleList structure representing the list of Triangles.
 */
void insertNodeInTriangleList(TriangleListNode *node, TriangleList *triangleList) {
    if (triangleList == NULL)
        return;
    if (!verifyTriangleExistence(node->triangle, triangleList)) {
        node->next = triangleList->head;
        triangleList->head = node;
    }
}

/**
 * Inserts a TetraListNode into the TetraList if the Tetragon does not already exist in the list.
 *
 * This function inserts a TetraListNode 'node' into the TetraList represented by 'tetraList',
 * but only if the Tetragon in the node does not already exist in the list. It first checks if 'tetraList'
 * is not NULL, and then uses the 'verifyTetragonExistence' function to check for the existence of the Tetragon.
 * If the Tetragon is not found in the list, the node is added to the front of the list.
 *
 * @param tetraList A pointer to the TetraList structure representing the list of Tetragons.
 * @param node A pointer to the TetraListNode structure to be inserted into the list.
 */
void insertNodeInTetraList(TetraList *tetraList, TetraListNode *node) {
    if (tetraList == NULL)
        return;
    if (!verifyTetragonExistence(tetraList, node->tetragon)) {
        node->next = tetraList->head;
        tetraList->head = node;
    }
}

/* Copy Graph */

/**
 * "Deep" memcpy for graphs.
 *
 * This function copies the graph 'graphSrc' to 'graphDest' by performing a deep copy. It ensures that all the graph
 * elements are correctly duplicated, including the adjacency lists and edge lists. It allocates memory for the destination
 * graph and copies each element individually to avoid shallow copying.
 *
 * @param graphDest A pointer to the destination Graph where the 'graphSrc' will be copied.
 * @param graphSrc A pointer to the source Graph that needs to be copied.
 */
void copyGraph(Graph *graphDest, Graph *graphSrc) {
    memcpy(graphDest, graphSrc, sizeof(Graph));

    //Copy adjlist

    graphDest->nodes = (struct AdjList *) malloc(graphSrc->V * sizeof(struct AdjList));
    memcpy(graphDest->nodes, graphSrc->nodes, sizeof(struct AdjList));
    for (int i = 0; i < graphSrc->current_nodes; i++) {
        graphDest->nodes[i].head = (struct AdjListNode *) malloc(sizeof(struct AdjListNode));
        memcpy(graphDest->nodes[i].head, graphSrc->nodes[i].head, sizeof(struct AdjListNode));

        struct AdjListNode *node = graphSrc->nodes[i].head;
        struct AdjListNode *copy = graphDest->nodes[i].head;

        while (node->next != NULL) {
            copy->next = (struct AdjListNode *) malloc(sizeof(struct AdjListNode));
            memcpy(copy->next, node->next, sizeof(struct AdjListNode));
            copy = copy->next;

            node = node->next;
        }
    }

    // Copy edges list
    graphDest->edges = (struct EdgeList *) malloc(sizeof(struct EdgeList));
    memcpy(graphDest->edges, graphSrc->edges, sizeof(struct AdjList));
    graphDest->edges->head = (struct EdgeListNode *) malloc(sizeof(struct EdgeListNode));
    memcpy(graphDest->edges->head, graphSrc->edges->head, sizeof(struct EdgeListNode));

    struct EdgeListNode *EdgeNode = graphSrc->edges->head;
    struct EdgeListNode *copyEdge = graphDest->edges->head;

    while (EdgeNode->next != NULL) {
        copyEdge->next = (struct EdgeListNode *) malloc(sizeof(struct EdgeListNode));
        memcpy(copyEdge->next, EdgeNode->next, sizeof(struct EdgeListNode));
        copyEdge = copyEdge->next;
        EdgeNode = EdgeNode->next;
    }

}

/* Graph Tree Operations */

/**
 * Check whether the graph 'g' exists in the tree of graphs represented by 'treeNode'.
 *
 * @param treeNode The root of the tree of graphs.
 * @param g The graph to check for existence in the tree.
 * @return 1 if the graph exists in the tree, 0 otherwise.
 */
int check_graph_existence(PrimeGraphTreeNode *treeNode, Graph *g) {
    int find = 0;
    sparsegraph sg1, sg2;
    options.getcanon = TRUE;

    if (treeNode->graph->current_nodes == g->current_nodes) {
        SG_INIT(sg1);
        SG_INIT(sg2);
        SG_ALLOC(sg1, MAX_VERTICES, MAX_DEGREE * MAX_VERTICES, "malloc");
        SG_ALLOC(sg2, MAX_VERTICES, MAX_DEGREE * MAX_VERTICES, "malloc");

        copy_sparse_graph(g);
        nauty((graph *) &sg, lab, ptn, NULL, orbits, &options, &stats, workspace, WORKSIZE, MAXM, g->current_nodes,
              (graph *) &sg1);
        copy_sparse_graph(treeNode->graph);
        nauty((graph *) &sg, lab, ptn, NULL, orbits, &options, &stats, workspace, WORKSIZE, MAXM, g->current_nodes,
              (graph *) &sg2);

        // Check if the graphs 'g' and the graph at 'treeNode' are isomorphic
        find = aresame_sg(&sg1, &sg2);

        SG_FREE(sg1);
        SG_FREE(sg2);
    } else
        // If the number of nodes does not match, the graph cannot exist in the tree
        find = 0;

    // Recursive call to check existence in the children of the current node
    for (int i = 0; i < treeNode->num_children; i++) {
        find = find || check_graph_existence(treeNode->children[i], g);
    }
    return find;
}

/**
 * Extend an irreducible graph 'graphSrc' by adding an adjacent diamond (K4-) to it.
 * The extended graph is stored in 'graphDest'.
 *
 * @param graphSrc The original irreducible graph to be extended.
 * @param graphDest The graph that will store the extended version of 'graphSrc'.
 * @param edge The edge to be used for the adjacent diamond insertion.
 * @return 1 if the extension is successful, 0 otherwise (when there are not enough spare vertices).
 */
int extendIrreducibleGraphWithAdjDiamond(Graph *graphSrc, Graph *graphDest, Edge *edge) {
    copyGraph(graphDest, graphSrc);
    return adjDiamondInsertion(graphDest, edge->src, edge->dest);
}

/**
 * Extend an irreducible graph 'graphSrc' by adding a lollipop structure to it.
 * The extended graph is stored in 'graphDest'.
 *
 * @param graphSrc The original irreducible graph to be extended.
 * @param graphDest The graph that will store the extended version of 'graphSrc'.
 * @param edge The edge to be used for the lollipop insertion.
 * @return 1 if the extension is successful, 0 otherwise (when there are not enough spare vertices).
 */
int extendIrreducibleGraphWithLollipop(Graph *graphSrc, Graph *graphDest, Edge *edge) {
    copyGraph(graphDest, graphSrc);
    return lollipopInsertion(graphDest, edge->src, edge->dest);
}

/**
 * Extend a graph 'graphSrc' by adding two new nodes and an edge connecting them.
 * The extended graph is stored in 'graphDest'.
 *
 * @param graphSrc The original graph to be extended.
 * @param graphDest The graph that will store the extended version of 'graphSrc'.
 * @param edge1 The first edge to be used for the extension.
 * @param edge2 The second edge to be used for the extension.
 * @return
 */
int extendGraphWithEdgeInsertion(Graph *graphSrc, Graph *graphDest, Edge *edge1, Edge *edge2) {
    copyGraph(graphDest, graphSrc);
    return edgeInsertion(graphDest, edge1->src, edge1->dest, edge2->src, edge2->dest);
}

/**
 * Extend an irreducible graph 'graphSrc' by adding a non-adjacent diamond structure to it.
 * The extended graph is stored in 'graphDest'.
 *
 * @param graphSrc The original irreducible graph to be extended.
 * @param graphDest The graph that will store the extended version of 'graphSrc'.
 * @param edge1 The first edge to be used for the non-adjacent diamond insertion.
 * @param edge2 The second edge to be used for the non-adjacent diamond insertion.
 * @return 1 if the extension is successful, 0 otherwise (when there are not enough spare vertices).
 */
int extendIrreducibleGraphWithNonAdjDiamond(Graph *graphSrc, Graph *graphDest, Edge *edge1, Edge *edge2) {
    copyGraph(graphDest, graphSrc);
    return nonAdjDiamondInsertion(graphDest, edge1->src, edge1->dest, edge2->src, edge2->dest);
}

/**
 * Generate children of the graph 'currentGraph' using the "adjacent diamond" insertion.
 * Each child graph will be added as a child node in the PrimeGraphTree 'tree'.
 *
 * @param currentGraph The current graph for which children are generated.
 * @param treeNode The corresponding PrimeGraphTreeNode of 'currentGraph' in the PrimeGraphTree.
 * @param tree The PrimeGraphTree where children are added.
 */
void generateChildrenWithAdjDiamond(Graph *currentGraph, PrimeGraphTreeNode *treeNode, PrimeGraphTree *tree) {
    printf("Generating children with adjacent diamond...\n");

    EdgeListNode *currentEdgeNode = currentGraph->edges->head;
    while (currentEdgeNode) {
        Graph *graphChild = createGraph(MAX_VERTICES);

        if (!extendIrreducibleGraphWithAdjDiamond(currentGraph, graphChild, currentEdgeNode->edge))
            continue;

        if (isAIrreducibleGraph(graphChild)) {
            int find = 0;
            find = check_graph_existence(tree->head, graphChild);

            if (!find) {
                addGraphToTree(treeNode, graphChild, tree, 1, "AdjDiamond");
                if (graphChild->current_nodes == MAX_VERTICES)
                    cubicGraphs[++cubicIndex] = graphChild;
            }
        }
        currentEdgeNode = currentEdgeNode->next;
    }
}

/**
 * Generate children of the graph 'currentGraph' using the "non-adjacent diamond" insertion.
 * Each child graph will be added as a child node in the PrimeGraphTree 'tree'.
 *
 * @param currentGraph The current graph for which children are generated.
 * @param treeNode The corresponding PrimeGraphTreeNode of 'currentGraph' in the PrimeGraphTree.
 * @param tree The PrimeGraphTree where children are added.
 */
void generateChildrenWithNonAdjDiamond(Graph *currentGraph, PrimeGraphTreeNode *treeNode, PrimeGraphTree *tree) {
    printf("Generating children with non-adjacent diamond...\n");
    EdgeListNode *currentEdgeNode1 = currentGraph->edges->head;
    EdgeListNode *currentEdgeNode2 = NULL;

    while (currentEdgeNode1) {
        currentEdgeNode2 = currentEdgeNode1->next; // Start from the next edge

        while (currentEdgeNode2) {
            Graph *graphChild = createGraph(MAX_VERTICES);
            if (!extendIrreducibleGraphWithNonAdjDiamond(currentGraph, graphChild, currentEdgeNode2->edge,
                                                         currentEdgeNode1->edge))
                continue;

            if (isAIrreducibleGraph(graphChild)) {
                int find = 0;
                find = check_graph_existence(tree->head, graphChild);

                if (!find) {
                    addGraphToTree(treeNode, graphChild, tree, 1, "NoAdjDiamond");
                    if (graphChild->current_nodes == MAX_VERTICES)
                        cubicGraphs[++cubicIndex] = graphChild;
                }
            }
            currentEdgeNode2 = currentEdgeNode2->next;
        }
        currentEdgeNode1 = currentEdgeNode1->next;
    }
}

/**
 * Generate children of the graph 'currentGraph' using the "lollipop" insertion.
 * Each child graph will be added as a child node in the PrimeGraphTree 'tree'.
 *
 * @param currentGraph The current graph for which children are generated.
 * @param treeNode The corresponding PrimeGraphTreeNode of 'currentGraph' in the PrimeGraphTree.
 * @param tree The PrimeGraphTree where children are added.
 */
void generateChildrenWithLollipop(Graph *currentGraph, PrimeGraphTreeNode *treeNode, PrimeGraphTree *tree) {
    printf("Generating children with lollipop...\n");

    EdgeListNode *currentEdgeNode = currentGraph->edges->head;

    while (currentEdgeNode) {
        Graph *graphChild = createGraph(MAX_VERTICES);
        if (!extendIrreducibleGraphWithLollipop(currentGraph, graphChild, currentEdgeNode->edge))
            continue;

        if (isAIrreducibleGraph(graphChild)) {
            int find = 0;
            find = check_graph_existence(tree->head, graphChild);

            if (!find) {
                addGraphToTree(treeNode, graphChild, tree, 1, "Lollipop");
                if (graphChild->current_nodes == MAX_VERTICES)
                    cubicGraphs[++cubicIndex] = graphChild;
            }
        }
        currentEdgeNode = currentEdgeNode->next;
    }
}

/**
 * Generate children of the PrimeGraphTreeNode 'treeNode' using the "edge" insertion.
 * Each child graph will be added as a child node in the PrimeGraphTree 'tree'.
 *
 * @param treeNode The PrimeGraphTreeNode for which children are generated.
 * @param tree The PrimeGraphTree where children are added.
 */
void recursiveGenFromPrime(PrimeGraphTreeNode *treeNode, PrimeGraphTree *tree) {
    printf("Generating children with edge insertion...\n");

    EdgeListNode *currentEdgeNode1 = treeNode->graph->edges->head;
    EdgeListNode *currentEdgeNode2 = NULL;

    while (currentEdgeNode1) {

        currentEdgeNode2 = currentEdgeNode1->next;
        while (currentEdgeNode2) {
            Graph *graphChild = createGraph(MAX_VERTICES);
            extendGraphWithEdgeInsertion(treeNode->graph, graphChild, currentEdgeNode1->edge, currentEdgeNode2->edge);

            int find = 0;
            find = check_graph_existence(tree->head, graphChild);
            if (!find) {
                addGraphToTree(treeNode, graphChild, tree, 0, "Edge");
                if (graphChild->current_nodes == MAX_VERTICES)
                    cubicGraphs[++cubicIndex] = graphChild;
            }
            currentEdgeNode2 = currentEdgeNode2->next;
        }
        currentEdgeNode1 = currentEdgeNode1->next;
    }

    for (int i = 0; i < treeNode->num_children; i++)
        if (treeNode->children[i]->graph->current_nodes <= MAX_VERTICES - 2 && !treeNode->children[i]->isPrime)
            recursiveGenFromPrime(treeNode->children[i], tree);

}

/**
 * Generate prime cubic graphs with random transformations (lollipop, adjacency diamond, and non-adjacency diamond).
 * The generated prime graphs are added as children of the PrimeGraphTreeNode 'treeNode' in the PrimeGraphTree 'tree'.
 *
 * @param treeNode The PrimeGraphTreeNode for which prime cubic graphs are generated.
 * @param tree The PrimeGraphTree where the prime graphs are added as children.
 */
void generatePrimeTreesFullRandom(PrimeGraphTreeNode *treeNode, PrimeGraphTree *tree) {
    Graph *currentGraph = treeNode->graph;

    int choice = 0; // 0: No possible insertion, 1: Only AdjDiamond, 2: All transformations possible

    if (currentGraph->current_nodes <= MAX_VERTICES - 4)
        choice = 1; // Only AdjDiamond possible
    if (currentGraph->current_nodes <= MAX_VERTICES - 6)
        choice = 2; // All transformations possible

    EdgeListNode *currentEdgeNode1 = currentGraph->edges->head;
    EdgeListNode *currentEdgeNode2 = NULL;

    if (choice == 0)
        return;

    else if (choice == 1) {
        generateChildrenWithAdjDiamond(currentGraph, treeNode, tree);
        return;
    }

    int randomChoice = 1;

    while (currentEdgeNode1) {

        currentEdgeNode2 = currentEdgeNode1->next;

        while (currentEdgeNode2) {

            Graph *graphChild = NULL;
            randomChoice = rand() % 3 + 1;

            switch (randomChoice) {
                case 1: // Try Lollipop Insertion
                    graphChild = createGraph(MAX_VERTICES);
                    if (!extendIrreducibleGraphWithLollipop(currentGraph, graphChild, currentEdgeNode1->edge))
                        continue;

                    // If the child graph is irreducible and not isomorphic to an existing graph, add it to the tree
                    if (isAIrreducibleGraph(graphChild)) {
                        int find = 0;
                        find = check_graph_existence(tree->head, graphChild);

                        if (!find) {
                            addGraphToTree(treeNode, graphChild, tree, 1, "Lollipop");
                            if (graphChild->current_nodes == MAX_VERTICES)
                                cubicGraphs[++cubicIndex] = graphChild;
                            //continue;
                        }
                    }
                    continue;
                case 2: // Try AdjDiamond Insertion
                    graphChild = createGraph(MAX_VERTICES);
                    if (!extendIrreducibleGraphWithAdjDiamond(currentGraph, graphChild, currentEdgeNode1->edge))
                        continue;

                    // If the child graph is irreducible and not isomorphic to an existing graph, add it to the tree
                    if (isAIrreducibleGraph(graphChild)) {
                        int find = 0;
                        find = check_graph_existence(tree->head, graphChild);

                        if (!find) {
                            addGraphToTree(treeNode, graphChild, tree, 1, "AdjDiamond");
                            if (graphChild->current_nodes == MAX_VERTICES)
                                cubicGraphs[++cubicIndex] = graphChild;
                            //continue;
                        }
                    }
                    continue;
                case 3: // Try NoAdjDiamond Insertion
                    graphChild = createGraph(MAX_VERTICES);
                    if (!extendIrreducibleGraphWithNonAdjDiamond(currentGraph, graphChild, currentEdgeNode2->edge,
                                                                 currentEdgeNode1->edge))
                        continue;

                    // If the child graph is irreducible and not isomorphic to an existing graph, add it to the tree
                    if (isAIrreducibleGraph(graphChild)) {
                        int find = 0;
                        find = check_graph_existence(tree->head, graphChild);

                        if (!find) {
                            addGraphToTree(treeNode, graphChild, tree, 1, "NoAdjDiamond");
                            if (graphChild->current_nodes == MAX_VERTICES)
                                cubicGraphs[++cubicIndex] = graphChild;
                            //continue;
                        }
                    }
                    break;
            }
            currentEdgeNode2 = currentEdgeNode2->next;
        }
        currentEdgeNode1 = currentEdgeNode1->next;
    }

    // Recursive call to generate prime graphs for each child node of the current node
    for (int i = 0; i < treeNode->num_children; i++)
        generatePrimeTreesFullRandom(treeNode->children[i], tree);

}

/**
 * Generate prime cubic graphs using adjacency diamond, non-adjacency diamond, and lollipop insertions.
 * The generated prime graphs are added as children of the PrimeGraphTreeNode 'treeNode' in the PrimeGraphTree 'tree'.
 *
 * @param treeNode The PrimeGraphTreeNode for which prime cubic graphs are generated.
 * @param tree The PrimeGraphTree where the prime graphs are added as children.
 */
void generatePrimeTrees(PrimeGraphTreeNode *treeNode, PrimeGraphTree *tree) {
    Graph *currentGraph = treeNode->graph;
    //EdgeListNode *currentEdgeNode = currentGraph->edges->head;

    int choice = 0;
    if (currentGraph->current_nodes <= MAX_VERTICES - 6)
        choice = 2;
    else if (currentGraph->current_nodes <= MAX_VERTICES - 4)
        choice = 1;

    if (choice == 0)
        return;

    else if (choice == 1)
        generateChildrenWithAdjDiamond(currentGraph, treeNode, tree);
    else {
        generateChildrenWithAdjDiamond(currentGraph, treeNode, tree);
        generateChildrenWithNonAdjDiamond(currentGraph, treeNode, tree);
        generateChildrenWithLollipop(currentGraph, treeNode, tree);
    }

    // Recursive call
    for (int i = 0; i < treeNode->num_children; i++)
        generatePrimeTrees(treeNode->children[i], tree);
}

/**
 * Generate all possible cubic graphs by recursively exploring the PrimeGraphTree rooted at 'treeNode'.
 * The generated cubic graphs are added as children of the PrimeGraphTreeNode 'treeNode' in the PrimeGraphTree 'tree'.
 * The function stops generating cubic graphs when the current graph's number of vertices exceeds 'MAX_VERTICES - 2'.
 *
 * @param treeNode The PrimeGraphTreeNode for which cubic graphs are generated.
 * @param tree The PrimeGraphTree where the cubic graphs are added as children.
 */
void generateCubicGraphTree(PrimeGraphTreeNode *treeNode, PrimeGraphTree *tree) {
    if (!treeNode->isPrime)
        return;
    if (treeNode->graph->current_nodes > MAX_VERTICES - 2)
        return;
    recursiveGenFromPrime(treeNode, tree);

    // Recursive call
    for (int i = 0; i < treeNode->num_children; i++)
        generateCubicGraphTree(treeNode->children[i], tree);
}

/**
 * Add a new graph as a child to the given parent node in the PrimeGraphTree.
 *
 * @param parent The parent node to which the new graph node will be added as a child.
 * @param g The graph to be added as a child node.
 * @param tree The PrimeGraphTree to which the graph will be added.
 * @param isPrime A flag indicating whether the graph is a prime graph (1) or not (0).
 * @param op The operation used to generate the graph (e.g., "AdjDiamond", "NoAdjDiamond", "Lollipop", "Edge").
 */
void addGraphToTree(PrimeGraphTreeNode *parent, Graph *g, PrimeGraphTree *tree, int isPrime, char *op) {
    PrimeGraphTreeNode *newNode = newPrimeGraphTreeNode(tree, g, g->current_nodes, isPrime, op);
    newNode->parent = parent;
    parent->num_children++;
    //parent->children = realloc( parent->children, (parent->num_children * sizeof(PrimeGraphTreeNode))); // la realloc funge???? controllare meglio!
    parent->children[parent->num_children - 1] = newNode;

}

/* Random Graph Generation */

//Function Not Used
EdgeListNode *chooseRandomEligible(Graph *graph) {
    int i = 0;
    int max = graph->current_edges - 1, min = 0;
    int r = rand() % (max + 1 - min) + min;
    printf("random! %d\n", r);

    EdgeListNode *randEdge = graph->edges->head;
    for (; i < r; i++) {
        if (randEdge == NULL) {
            r = rand() % (max + 1 - min) + min;
            i = 0;
        }
        if (!randEdge->eligible)
            i++;
        randEdge = randEdge->next;
    }
    return randEdge;
}

/**
 * Generates a random cubic graph from the array of cubic graphs (cubicGraphs).
 * The function randomly selects a graph from the array and returns it.
 *
 * @return Graph* Pointer to the randomly selected cubic graph.
 */
Graph *generateUniform() {
    int max = cubicIndex, min = 1;
    int r = rand() % (max + 1 - min) + min;

    return cubicGraphs[r];
}

/* Printing Functions */

/**
 * Prints the list of edges in the graph.
 *
 * This function prints the list of edges in the graph 'graph'.
 * It traverses the EdgeList and prints the source, destination, and 'eligible' value of each edge.
 *
 * @param graph A pointer to the Graph structure representing the graph.
 */
void printEdgeList(Graph *graph) {
    printf("Edge List: \n");
    struct EdgeListNode *edgeNode = graph->edges->head;
    while (edgeNode != NULL) {
        printf("%d %d %d\n", edgeNode->edge->src, edgeNode->edge->dest, edgeNode->eligible);
        edgeNode = edgeNode->next;
    }
}

/**
 * Prints the list of triangles in the TriangleList.
 *
 * This function prints the list of triangles in the TriangleList 'triangleList'.
 * It traverses the TriangleList and prints the three edges of each triangle.
 *
 * @param triangleList A pointer to the TriangleList structure representing the list of triangles.
 */
void printTriangleList(TriangleList *triangleList) {
    printf("Printing triangle list...\n");
    TriangleListNode *triangleListNode = triangleList->head;
    while (triangleListNode != NULL) {
        printf("(%d ,%d) ", triangleListNode->triangle->e1->src, triangleListNode->triangle->e1->dest);
        printf("(%d ,%d) ", triangleListNode->triangle->e2->src, triangleListNode->triangle->e2->dest);
        printf("(%d ,%d)\n ", triangleListNode->triangle->e3->src, triangleListNode->triangle->e3->dest);
        triangleListNode = triangleListNode->next;
    }
}

/**
 * Prints the list of triangles in the graph.
 *
 * This function prints the list of triangles that are stored in the edges of the graph 'g'.
 * It traverses the EdgeList and prints the three edges of each triangle found in the graph.
 *
 * @param g A pointer to the Graph structure representing the graph.
 */
void printTriangleListAlt(Graph *g) {
    EdgeListNode *node = g->edges->head;
    while (node != NULL) {
        if (node->triangle != NULL) {
            printf("(%d ,%d) ", node->triangle->e1->src, node->triangle->e1->dest);
            printf("(%d ,%d) ", node->triangle->e2->src, node->triangle->e2->dest);
            printf("(%d ,%d)\n ", node->triangle->e3->src, node->triangle->e3->dest);
        }
        node = node->next;
    }
}

/**
 * Prints the list of tetragons in the TetraList.
 *
 * This function prints the list of tetragons in the TetraList 'tetraList'.
 * It traverses the TetraList and prints the four edges of each tetragon.
 *
 * @param tetraList A pointer to the TetraList structure representing the list of tetragons.
 */
void printTetragonList(TetraList *tetraList) {
    printf("Printing tetragon list...\n");
    TetraListNode *tetraListNode = tetraList->head;
    while (tetraListNode != NULL) {
        printf("(%d ,%d) ", tetraListNode->tetragon->e1->src, tetraListNode->tetragon->e1->dest);
        printf("(%d ,%d) ", tetraListNode->tetragon->e2->src, tetraListNode->tetragon->e2->dest);
        printf("(%d ,%d) ", tetraListNode->tetragon->e3->src, tetraListNode->tetragon->e3->dest);
        printf("(%d ,%d)\n ", tetraListNode->tetragon->e4->src, tetraListNode->tetragon->e4->dest);
        tetraListNode = tetraListNode->next;
    }
}

/**
 * Prints the list of tetragons in the graph.
 *
 * This function prints the list of tetragons that are stored in the edges of the graph 'g'.
 * It traverses the EdgeList and prints the four edges of each tetragon found in the graph.
 *
 * @param g A pointer to the Graph structure representing the graph.
 */
void printTetragonListAlt(Graph *g) {
    EdgeListNode *node = g->edges->head;
    while (node != NULL) {
        if (node->tetragon != NULL) {
            printf("(%d ,%d) ", node->tetragon->e1->src, node->tetragon->e1->dest);
            printf("(%d ,%d) ", node->tetragon->e2->src, node->tetragon->e2->dest);
            printf("(%d ,%d) ", node->tetragon->e3->src, node->tetragon->e3->dest);
            printf("(%d ,%d)\n ", node->tetragon->e4->src, node->tetragon->e4->dest);
        }
        node = node->next;
    }
}

/**
 * Prints the adjacency list of each vertex in the graph.
 *
 * This function prints the adjacency list of each vertex in the graph 'graph'.
 * It traverses the graph's nodes and prints the adjacent vertices of each vertex.
 *
 * @param graph A pointer to the Graph structure representing the graph.
 */
void printGraph(struct Graph *graph) {
    int v;
    for (v = 0; v < graph->V; ++v) {
        struct AdjListNode *pCrawl = graph->nodes[v].head;
        printf("Adjacency list of vertex %d head: ", v);
        while (pCrawl) {
            printf("-> %d", pCrawl->label);
            pCrawl = pCrawl->next;
        }
        printf("\n");
    }
}

/**
 * Gets spare vertices in the graph.
 *
 * This function retrieves the unused (spare) vertices in the graph 'graph' and stores them in the array 'spareVertices'.
 * The array 'spareVertices' should have a size of 'size' to accommodate the spare vertices.
 * It iterates through the graph's vertices and adds the index of any vertex with no adjacent nodes to the 'spareVertices' array.
 * The function returns the number of spare vertices found and stored in the array.
 *
 * @param graph A pointer to the Graph structure representing the graph.
 * @param spareVertices An array to store the indices of spare vertices.
 * @param size The size of the 'spareVertices' array.
 * @return The number of spare vertices found and stored in the 'spareVertices' array.
 */
int getSpareVertices(Graph *graph, int *spareVertices, int size) {
    int j = 0;
    for (int i = 0; i < graph->V && j < size; i++)
        if (graph->nodes[i].head == NULL) {
            spareVertices[j] = i;
            j++;
        }
    if (j < size)
        spareVertices[j] = -1;
    return j;
}

/* Tree Saving */

/**
 * Recursively prints the nodes of the PrimeGraphTree to a file in CSV format.
 *
 * @param node Pointer to the current PrimeGraphTreeNode being processed.
 * @param file Pointer to the file where the nodes will be printed.
 */
void recursivePrintNodesTree(PrimeGraphTreeNode *node, FILE *file) {
    if (node == NULL) {
        return;
    }

    fprintf(file, "%d,%d\n", node->id, node->graph->current_nodes);

    for (int i = 0; i < node->num_children; i++) {
        recursivePrintNodesTree(node->children[i], file);
    }
}

/**
 * Recursively print the edges of the PrimeGraphTree to a CSV file.
 *
 * @param node Pointer to the current PrimeGraphTreeNode in the recursion.
 * @param file Pointer to the file where the edges will be saved in CSV format.
 */
void recursivePrintEdgeTree(PrimeGraphTreeNode *node, FILE *file) {
    if (node == NULL) {
        return;
    }

    if (node->parent == NULL)
        fprintf(file, "%d,%d,%s,%d\n", 0, node->id, node->op, 0);
    else {
        int add = node->graph->current_nodes - node->parent->graph->current_nodes;
        fprintf(file, "%d,%d,%s,%d\n", node->parent->id, node->id, node->op, 0);
    }


    for (int i = 0; i < node->num_children; i++) {
        recursivePrintEdgeTree(node->children[i], file);
    }
}

/**
 * Save the PrimeGraphTree to a CSV file.
 *
 * @param root Pointer to the root PrimeGraphTreeNode of the tree.
 * @param filename The name of the file where the tree will be saved in CSV format.
 */
void saveTreeToCSV(PrimeGraphTreeNode *root, const char *filename) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        printf("Impossibile aprire il file %s per la scrittura.\n", filename);
        return;
    }

    // Write CSV header
    fprintf(file, "ID,CURR_NODES\n");

    recursivePrintNodesTree(root, file);

    fclose(file);
}

/**
 * Save the edges of the PrimeGraphTree to a CSV file.
 *
 * @param root Pointer to the root PrimeGraphTreeNode of the PrimeGraphTree.
 * @param filename The name of the file where the edges will be saved in CSV format.
 */
void saveEdgeTreeToCSV(PrimeGraphTreeNode *root, const char *filename) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        printf("Impossibile aprire il file %s per la scrittura.\n", filename);
        return;
    }

    // Write CSV header
    fprintf(file, "SRC_ID,DST_ID,OP,ADD\n");

    recursivePrintEdgeTree(root, file);

    fclose(file);
}

/* Main Function */

int main(int argc, char *argv[]) {

    if (argc != 4) {
        printf("Usage: %s <node_file_path> <edge_file_path> <full_random: 0/1>\n", argv[0]);
        return 1;
    }
    // Convert the string argument to a boolean value
    bool full_random = atoi(argv[3]) != 0;

    // Memory allocation
    cubicGraphs = (Graph **) malloc(10000 * sizeof(Graph *));
    srand(time(NULL));

    //For the timing
    struct tms TMS;
    unsigned int oldtime = 0;

    // Initialization
    init_nauty_options();
    Graph *g = createGraph(2 * MAX_VERTICES);
    PrimeGraphTree *primeTree = createPrimeGraphTree();
    initIrreducibleGraph(g);

    // Create a new Prime Graph
    char *op = "K4";
    PrimeGraphTreeNode *node = newPrimeGraphTreeNode(primeTree, g, 4, 1, op);
    primeTree->head = node;

    if (full_random)
    generatePrimeTreesFullRandom(node, primeTree);
    else
        generatePrimeTrees(node, primeTree);

    generateCubicGraphTree(node, primeTree);

    printf("Total number of generated cubic graphs: %d\n", primeTree->current_graphs);
    printf("Number of generated cubic graphs with vertices %d: %d\n", MAX_VERTICES, cubicIndex);
    printf("Random graph:\n");
    printEdgeList(generateUniform());

    // Calculate CPU time used and print the result
    times(&TMS);
    unsigned long savetime = oldtime + (unsigned int) TMS.tms_utime;
    fprintf(stderr, "CPU time: %.1f seconds.\n", (double) savetime / time_factor);

    saveTreeToCSV(primeTree->head, argv[1]);
    saveEdgeTreeToCSV(primeTree->head, argv[2]);

    return 0;
}