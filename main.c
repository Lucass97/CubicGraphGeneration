#include <stdio.h>
#include <stdlib.h>
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
Graph** cubicGraphs = NULL;


void init_nauty_options() {
    //options.getcanon = TRUE;
    //options.userautomproc = save_generators;

    /* Init the nauty datastructures */
    SG_INIT(sg);
    SG_ALLOC(sg, MAX_VERTICES, 3 * MAX_VERTICES, "malloc");

    //sg.v and sg.d only have to be set once
    int i;
    for(i = 0; i < MAX_VERTICES; i++) {
        sg.v[i] = i * MAX_DEGREE;
        sg.d[i] = 3;
    }

    SG_INIT(sg_canon);
    SG_ALLOC(sg_canon, MAX_VERTICES, 3 * MAX_VERTICES, "malloc");
}


int countVerticeDegree(Graph* graph, int label){
    int count = 0;
    struct AdjListNode* nodeDest = graph->nodes[label].head;
    while(nodeDest != NULL){
        count++;
        nodeDest = nodeDest->next;
    }
    return count;
}

void addEdge(Graph* graph, int src, int dest)
{
    // Check
    if (countVerticeDegree(graph, src) >= MAX_DEGREE || countVerticeDegree(graph, dest) >= MAX_DEGREE)
        return;
    if(graph->nodes[src].head == NULL)
        graph->current_nodes++;
    if(graph->nodes[dest].head == NULL)
        graph->current_nodes++;

    //Add Edge to AdjList
    struct AdjListNode* nodeDest = graph->nodes[src].head;
    while (nodeDest){
        if(nodeDest->label == dest)
            return;
        nodeDest = nodeDest->next;
    }

    struct AdjListNode* newNode = newAdjListNode(dest);
    newNode->next = graph->nodes[src].head;
    graph->nodes[src].head = newNode;

    newNode = newAdjListNode(src);
    newNode->next = graph->nodes[dest].head;
    graph->nodes[dest].head = newNode;

    //Add Edge to EdgeList
    EdgeListNode*  edge = newEdgeListNode(src, dest);
    edge->next = graph->edges->head;
    graph->edges->head = edge;

    //Update total edges
    graph->current_edges++;
}

void removeNodeFromAdjList(Graph* graph, int src, int dest){
    struct AdjListNode* nodeDest = graph->nodes[src].head;
    struct  AdjListNode* placeHolder = NULL;
    if (nodeDest->label == dest)
    {
        placeHolder = nodeDest;
        graph->nodes[src].head = nodeDest->next;
        free(placeHolder);
        return; //prima non c'era
    }
    while(nodeDest->next != NULL) {
        if (nodeDest->next->label == dest){
            placeHolder = nodeDest->next;
            nodeDest->next = nodeDest->next->next;
            free(placeHolder);
            return; //prima c'era break
        }
        nodeDest = nodeDest->next;
    }
}

void removeEdge(Graph* graph, int src, int dest){

    //Remove Edge from AdjList
    removeNodeFromAdjList(graph, src, dest);
    removeNodeFromAdjList(graph, dest, src);

    //Remove Edge from edgeList
    removeEdgeFromEdgeList(graph->edges, src, dest);
    removeEdgeFromEdgeList(graph->edges, dest, src);

    //Update current edges
    graph->current_edges--;
}

void removeEdgeFromEdgeList(EdgeList* edgeList, int src, int dest){
    EdgeListNode* edgeNode = edgeList->head;
    EdgeListNode* placeHolder = NULL;
    if (edgeNode->edge->src == src && edgeNode->edge->dest == dest)
    {
        placeHolder = edgeNode;
        edgeList->head = edgeNode->next;
        //free(placeHolder->edge);
        free(placeHolder);
        return; //non c'era prima
    }
    while(edgeNode->next != NULL) {
        if (edgeNode->next->edge->src == src && edgeNode->next->edge->dest == dest){
            EdgeListNode* placeHolder2 = edgeNode->next;
            edgeNode->next = edgeNode->next->next;
            //free(placeHolder2->edge);
            free(placeHolder2);
            return; //c'era break prima
        }
        edgeNode = edgeNode->next;
    }
}

void popFromEdgeList(EdgeList* edgeList) {
    EdgeListNode* headEdgeListNode = edgeList->head;
    if (headEdgeListNode != NULL) {
        edgeList->head = headEdgeListNode->next;
        //free(headEdgeListNode->edge);
        free(headEdgeListNode);
        return;
    }
    return;
}

void insertInEdgeList(EdgeList* edgeList, int src, int dest) {
    if (!verifyEdgeExistence(src, dest, edgeList)) {
        EdgeListNode *edge = newEdgeListNode(src, dest);
        edge->next = edgeList->head;
        edgeList->head = edge;
    }
}


/*Funzione per le Triangle List */
int compareTriangle(Triangle* t1, Triangle* t2){
    return compareEdge(t1->e1, t2->e1) && compareEdge(t1->e2, t2->e2) && compareEdge(t1->e3, t2->e3) ||
            compareEdge(t1->e1, t2->e2) && compareEdge(t1->e2, t2->e1) && compareEdge(t1->e3, t2->e3) ||
            compareEdge(t1->e1, t2->e1) && compareEdge(t1->e2, t2->e3) && compareEdge(t1->e3, t2->e2) ||
            compareEdge(t1->e1, t2->e3) && compareEdge(t1->e2, t2->e2) && compareEdge(t1->e3, t2->e1);
}

int isATriangle(Triangle* t){
    int v[MAX_VERTICES];
    for(int i=0; i<MAX_VERTICES; i++)
        v[i] = 0;
    v[t->e1->src] = 1;
    v[t->e1->dest] = 1;
    v[t->e2->src] = 1;
    v[t->e2->dest] = 1;
    v[t->e3->src] = 1;
    v[t->e3->dest] = 1;
    int count = 0;

    for(int i=0; i<MAX_VERTICES; i++)
        count += v[i];
    if(count == 3)
        return 1;
    return 0;
}

int verifyTriangleExistence(Triangle* triangle, TriangleList* triangleList){
    TriangleListNode* triangleListNode = triangleList->head;
    while(triangleListNode != NULL){
        if(isATriangle(triangle))
            if (compareTriangle(triangle, triangleListNode->triangle))
                return 1;
        triangleListNode = triangleListNode->next;
    }
    return 0;
}

void insertNodeInTriangleList(TriangleListNode * node, TriangleList* triangleList){
    if(triangleList == NULL)
        return;
    if (!verifyTriangleExistence(node->triangle, triangleList)) {
        node->next = triangleList->head;
        triangleList->head = node;
    }
}

/*Funzione per le Tetra List */
int compareTetragon(Tetragon* t1, Tetragon* t2){
    return compareEdge(t1->e1, t2->e1) && compareEdge(t1->e2, t2->e2) && compareEdge(t1->e3, t2->e3) && compareEdge(t1->e4, t2->e4) ||
           compareEdge(t1->e1, t2->e2) && compareEdge(t1->e2, t2->e1) && compareEdge(t1->e3, t2->e3) && compareEdge(t1->e4, t2->e4) ||
           compareEdge(t1->e1, t2->e1) && compareEdge(t1->e2, t2->e3) && compareEdge(t1->e3, t2->e2) && compareEdge(t1->e4, t2->e4)||
           compareEdge(t1->e1, t2->e3) && compareEdge(t1->e2, t2->e2) && compareEdge(t1->e3, t2->e1) && compareEdge(t1->e4, t2->e4) ||
           compareEdge(t1->e1, t2->e4) && compareEdge(t1->e2, t2->e2) && compareEdge(t1->e3, t2->e3) && compareEdge(t1->e4, t2->e1) ||
           compareEdge(t1->e1, t2->e4) && compareEdge(t1->e2, t2->e3) && compareEdge(t1->e3, t2->e2) && compareEdge(t1->e4, t2->e1) ||
           compareEdge(t1->e1, t2->e1) && compareEdge(t1->e2, t2->e4) && compareEdge(t1->e3, t2->e3) && compareEdge(t1->e4, t2->e2);
}

int isATetragon(Tetragon* t){
    int v[MAX_VERTICES];
    for(int i=0; i<MAX_VERTICES; i++)
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
    for(int i=0; i<MAX_VERTICES; i++)
        count += v[i];
    if(count == 4)
        return 1;
    return 0;
}

int verifyTetragonExistence(TetraList* tetraList, Tetragon* tetra){
    TetraListNode* tetraListNode = tetraList->head;
    while(tetraListNode != NULL){
        if(isATetragon(tetra))
            if (compareTetragon(tetra, tetraListNode->tetragon))
                return 1;
        tetraListNode = tetraListNode->next;
    }
    return 0;
}

void insertNodeInTetraList(TetraList* tetraList, TetraListNode * node){
    if(tetraList == NULL)
        return;
    if (!verifyTetragonExistence(tetraList, node->tetragon)) {
        node->next = tetraList->head;
        tetraList->head = node;
    }
}
void dfsForBridge(Graph* graph, int start, int* visited){
    visited[start] = 1;
    EdgeListNode* edgeListNode =  graph->edges->head;
    while(edgeListNode){
        if(edgeListNode->edge->src == start && !visited[edgeListNode->edge->dest])
            dfsForBridge(graph, edgeListNode->edge->dest, visited);
        if(edgeListNode->edge->dest == start && !visited[edgeListNode->edge->src])
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

int isABridge(Graph* graph, Edge* edge){
    removeEdge(graph, edge->src,edge->dest);
    int visited[graph->V];
    for(int i=0;  i < (graph->V); i++){
        visited[i] = 0;
    }
    dfsForBridge(graph, edge->src, visited);

    addEdge(graph,edge->src,edge->dest);
    return !visited[edge->dest];
}

//inutilizzata
int isABridgeList(Graph* graph){
    struct EdgeListNode* edgeNode = graph->edges->head;
    while(edgeNode != NULL){
        if (!isABridge(graph, edgeNode->edge))
            return 0;
        edgeNode = edgeNode->next;
    }
    return 1;
}

void printEdgeList(Graph* graph)
{
    struct EdgeListNode* edgeNode = graph->edges->head;
    while(edgeNode != NULL){
        printf("%d %d %d\n", edgeNode->edge->src, edgeNode->edge->dest, edgeNode->eligible);
        edgeNode = edgeNode->next;
    }
}

void printTriangle(TriangleList* triangleList){
    printf("Printing triangle list...\n");
    TriangleListNode* triangleListNode = triangleList->head;
    while(triangleListNode != NULL){
        printf("(%d ,%d) ", triangleListNode->triangle->e1->src, triangleListNode->triangle->e1->dest);
        printf("(%d ,%d) ", triangleListNode->triangle->e2->src, triangleListNode->triangle->e2->dest);
        printf("(%d ,%d)\n ", triangleListNode->triangle->e3->src, triangleListNode->triangle->e3->dest);
        triangleListNode = triangleListNode->next;
    }
}
void printTriangleNew(Graph* g){
    EdgeListNode * node = g->edges->head;
    while(node != NULL){
        if(node->triangle !=NULL){
        printf("(%d ,%d) ", node->triangle->e1->src, node->triangle->e1->dest);
        printf("(%d ,%d) ", node->triangle->e2->src, node->triangle->e2->dest);
        printf("(%d ,%d)\n ", node->triangle->e3->src, node->triangle->e3->dest);}
        node = node->next;
    }
}

void printTetragon(TetraList* tetraList){
    printf("Printing tetragon list...\n");
    TetraListNode* tetraListNode = tetraList->head;
    while(tetraListNode != NULL){
        printf("(%d ,%d) ", tetraListNode->tetragon->e1->src, tetraListNode->tetragon->e1->dest);
        printf("(%d ,%d) ", tetraListNode->tetragon->e2->src, tetraListNode->tetragon->e2->dest);
        printf("(%d ,%d) ", tetraListNode->tetragon->e3->src, tetraListNode->tetragon->e3->dest);
        printf("(%d ,%d)\n ", tetraListNode->tetragon->e4->src, tetraListNode->tetragon->e4->dest);
        tetraListNode = tetraListNode->next;
    }
}
void printTetragonNew(Graph* g){
    EdgeListNode * node = g->edges->head;
    while(node != NULL){
        if(node->tetragon !=NULL){
        printf("(%d ,%d) ", node->tetragon->e1->src, node->tetragon->e1->dest);
        printf("(%d ,%d) ", node->tetragon->e2->src, node->tetragon->e2->dest);
        printf("(%d ,%d) ", node->tetragon->e3->src, node->tetragon->e3->dest);
        printf("(%d ,%d)\n ", node->tetragon->e4->src, node->tetragon->e4->dest);}
        node = node->next;
    }
}

void printGraph(struct Graph* graph)
{
    int v;
    for (v = 0; v < graph->V; ++v)
    {
        struct AdjListNode* pCrawl = graph->nodes[v].head;
        printf("Adjacency list of vertex %d head: ", v);
        while (pCrawl)
        {
            printf("-> %d", pCrawl->label);
            pCrawl = pCrawl->next;
        }
        printf("\n");
    }
}

int getSpareVertices(Graph* graph, int* spareVertices, int size){
    int j = 0;
    for(int i=0;  i < graph->V && j < size; i++)
        if(graph->nodes[i].head == NULL) {
            spareVertices[j] = i;
            j++;
        }
    if (j < size)
        spareVertices[j] = -1;
    return j;
}

int edgeInsertion(Graph* graph, int src1, int dest1, int src2, int dest2)
{
    int size = 2;
    int spareVertices[size];
    int retrVertices = getSpareVertices(graph, spareVertices, size);
    if (retrVertices < 2){
        return 0; // non ci sono abbastanza vertici spare
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
int nonAdjDiamondInsertion(Graph* graph, int src1, int dest1, int src2, int dest2)
{
    int size = 6;
    int spareVertices[size];
    int retrVertices = getSpareVertices(graph, spareVertices, size);
    if (retrVertices < size){
        return 0; // non ci sono abbastanza vertici spare
    }
    //create K4-
    addEdge(graph, spareVertices[0], spareVertices[1]);
    addEdge(graph, spareVertices[1], spareVertices[2]);
    addEdge(graph, spareVertices[2], spareVertices[3]);
    addEdge(graph, spareVertices[3], spareVertices[0]);
    addEdge(graph, spareVertices[0], spareVertices[2]);

    //
    addEdge(graph, spareVertices[1], spareVertices[4]);
    addEdge(graph, spareVertices[3], spareVertices[5]);
    //
    removeEdge(graph, src1, dest1);
    removeEdge(graph, src2, dest2);
    addEdge(graph, src1, spareVertices[4]);
    addEdge(graph, spareVertices[4], dest1);
    addEdge(graph, src2,spareVertices[5] );
    addEdge(graph, spareVertices[5], dest2);

    return 1;
}
int adjDiamondInsertion(Graph* graph, int src, int dest)
{
    int size = 4;
    int spareVertices[size];
    int retrVertices = getSpareVertices(graph, spareVertices, size);
    if (retrVertices < size){
        return 0; // non ci sono abbastanza vertici spare
    }
    //create K4-
    addEdge(graph, spareVertices[0], spareVertices[1]);
    addEdge(graph, spareVertices[1], spareVertices[2]);
    addEdge(graph, spareVertices[2], spareVertices[3]);
    addEdge(graph, spareVertices[3], spareVertices[0]);
    addEdge(graph, spareVertices[0], spareVertices[2]);

    //
    removeEdge(graph, src, dest);
    addEdge(graph, src, spareVertices[3]);
    addEdge(graph, spareVertices[1], dest);

    return 1;
}
int lollipopInsertion(Graph* graph, int src, int dest)
{
    int size = 6;
    int spareVertices[size];
    int retrVertices = getSpareVertices(graph, spareVertices, size);
    if (retrVertices < size){
        printf(" ma cos\n");
        return 0; // non ci sono abbastanza vertici spare
    }

    //create K4+
    addEdge(graph, spareVertices[0], spareVertices[1]);
    addEdge(graph, spareVertices[1], spareVertices[2]);
    addEdge(graph, spareVertices[3], spareVertices[0]);
    addEdge(graph, spareVertices[0], spareVertices[2]);
    addEdge(graph, spareVertices[3], spareVertices[1]);
    addEdge(graph, spareVertices[3], spareVertices[4]);
    addEdge(graph, spareVertices[2], spareVertices[4]);

    //
    addEdge(graph, spareVertices[4], spareVertices[5]);
    removeEdge(graph, src, dest);
    addEdge(graph, src, spareVertices[5]);
    addEdge(graph, spareVertices[5], dest);

    return 1;
}
//inutilizzata
EdgeListNode* chooseRandomEligible(Graph* graph){
    int i = 0;
    int max = graph->current_edges-1, min = 0;
    int r = rand() % (max + 1 - min) + min;
    printf("random! %d\n",r);

    EdgeListNode* randEdge = graph->edges->head;
    for(;i<r;i++)
    {
      if (randEdge == NULL) {
          r = rand() % (max + 1 - min) + min;
          i=0;
      }
        if(!randEdge->eligible)
            i++;
        randEdge = randEdge->next;
    }
    return randEdge;
}


int compareEdge(Edge* e1, Edge* e2){
    return ((e1->src == e2->src && e1->dest == e2->dest) || (e1->src == e2->dest && e1->dest == e2->src));
}

// one of the edge endpoints is part of a triangle and the edge isn't
int irreducibilityCondition2(Graph* graph, Edge* edge){
    TriangleList* triangleList = graph->triangles;
    int src = edge->src;
    int dest = edge->dest;
    TriangleListNode* triangleNode = triangleList->head;
    while(triangleNode){
        if(!compareEdge(edge,triangleNode->triangle->e1)&& !compareEdge(edge,triangleNode->triangle->e2) && !compareEdge(edge,triangleNode->triangle->e3))
            if(src == triangleNode->triangle->e1->src || src == triangleNode->triangle->e1->dest || src == triangleNode->triangle->e2->src || src == triangleNode->triangle->e2->dest|| src == triangleNode->triangle->e3->src || src == triangleNode->triangle->e3->dest ||
                    dest == triangleNode->triangle->e1->src ||  dest == triangleNode->triangle->e1->dest ||  dest == triangleNode->triangle->e2->src ||  dest == triangleNode->triangle->e2->dest||  dest == triangleNode->triangle->e3->src ||  dest == triangleNode->triangle->e3->dest)
                return 1;

        triangleNode = triangleNode->next;

    }
    return 0;

}

// both of the edge endpoints are part of the same 4-gon and the edge isn't
int irreducibilityCondition3(Graph* graph, Edge* edge){
    TetraList* tList = graph->tetragons;
    int src = edge->src;
    int dest = edge->dest;
    TetraListNode* tNode = tList->head;
    while(tNode){
        if(!compareEdge(edge,tNode->tetragon->e1)&& !compareEdge(edge,tNode->tetragon->e2) && !compareEdge(edge,tNode->tetragon->e3)  && !compareEdge(edge,tNode->tetragon->e4))
            if((src == tNode->tetragon->e1->src || src == tNode->tetragon->e1->dest || src == tNode->tetragon->e2->src || src == tNode->tetragon->e2->dest|| src == tNode->tetragon->e3->src || src == tNode->tetragon->e3->dest  || src == tNode->tetragon->e4->src || src == tNode->tetragon->e4->dest) &&
                    (dest == tNode->tetragon->e1->src ||  dest == tNode->tetragon->e1->dest ||  dest == tNode->tetragon->e2->src ||  dest == tNode->tetragon->e2->dest||  dest == tNode->tetragon->e3->src ||  dest == tNode->tetragon->e3->dest || dest == tNode->tetragon->e4->src || dest == tNode->tetragon->e4->dest))
                return 1;

            tNode = tNode->next;
    }
    return 0;

}
int isAIrreducibleGraph(Graph* g){

    TriangleList* trianglelist = (TriangleList *) malloc(sizeof(TriangleList));
    trianglelist->head = NULL;
    findTriangles(g, 0,trianglelist);
    insertTriangles(g, trianglelist);
    TetraList* tetralist = (TetraList *) malloc(sizeof(TetraList));
    tetralist->head = NULL;
    findTetragons(g,0, tetralist);
    insertTetragons(g, tetralist);



    EdgeListNode* edgeListNode = g->edges->head;
    EdgeListNode* edgeListNodeSucc = NULL;
    while(edgeListNode){
        edgeListNodeSucc = edgeListNode->next;
        Edge* edge = edgeListNode->edge;

            int b = isABridge(g, edge);
            //edge = edgeListNode->edge;
            int irr2 = irreducibilityCondition2(g, edge);
            int irr3 = irreducibilityCondition3(g, edge);
            //printf("Edge: %d %d ", edge->src, edge->dest);
            //printf("c1:%d, c2:%d, c3:%d\n", b, irr2, irr3);

            if(!(b || irr2||  irr3))
                return 0;

        edgeListNode = edgeListNodeSucc;
    }
    return 1;

}

void initIrreducibleGraph(Graph* graph){

    int size = 4;
    int spareVertices[size];
    int retrVertices = getSpareVertices(graph, spareVertices, size);
    if (retrVertices < size){
        return; // non ci sono abbastanza vertici spare
    }
    addEdge(graph, spareVertices[0], spareVertices[1]);
    addEdge(graph, spareVertices[1], spareVertices[2]);
    addEdge(graph, spareVertices[2], spareVertices[3]);
    addEdge(graph, spareVertices[3], spareVertices[0]);
    addEdge(graph, spareVertices[0], spareVertices[2]);
    addEdge(graph, spareVertices[1], spareVertices[3]);

}

// "deep" memcpy per grafi
void copyGraph(Graph* graphDest, Graph* graphSrc)
{
    memcpy(graphDest, graphSrc, sizeof(Graph));

    //copy adjlist
    graphDest->nodes = (struct AdjList*) malloc(graphSrc->V * sizeof(struct AdjList));
    memcpy(graphDest->nodes, graphSrc->nodes, sizeof(struct AdjList));
    for(int i=0; i<graphSrc->current_nodes; i++) {
        graphDest->nodes[i].head = (struct AdjListNode *) malloc(sizeof(struct AdjListNode));
        memcpy(graphDest->nodes[i].head, graphSrc->nodes[i].head, sizeof(struct AdjListNode));

        struct AdjListNode* node = graphSrc->nodes[i].head;
        struct AdjListNode* copy = graphDest->nodes[i].head;

        while (node->next != NULL) {
            copy->next = (struct AdjListNode *) malloc(sizeof(struct AdjListNode));
            memcpy(copy->next, node->next, sizeof(struct AdjListNode));
            copy = copy->next;

            node = node->next;
        }
    }

    // copy edges list
    graphDest->edges = (struct EdgeList*) malloc(sizeof(struct EdgeList));
    memcpy(graphDest->edges, graphSrc->edges, sizeof(struct AdjList));
    graphDest->edges->head = (struct EdgeListNode*) malloc(sizeof(struct EdgeListNode));
    memcpy(graphDest->edges->head, graphSrc->edges->head, sizeof(struct EdgeListNode));

    struct EdgeListNode* EdgeNode = graphSrc->edges->head;
    struct EdgeListNode* copyEdge= graphDest->edges->head;

    while(EdgeNode->next !=NULL)
    {
        copyEdge->next = (struct EdgeListNode*) malloc(sizeof(struct EdgeListNode));
        memcpy(copyEdge->next, EdgeNode->next, sizeof(struct EdgeListNode));
        copyEdge = copyEdge ->next;
        EdgeNode = EdgeNode->next;
    }


}

int verifyEdgeExistence(int src, int dest, EdgeList* edgeList){
    struct EdgeListNode* edgeNode = edgeList->head;
    while(edgeNode != NULL){
        if((edgeNode->edge->src == src && edgeNode->edge->dest == dest)
            || (edgeNode->edge->src == dest && edgeNode->edge->dest == src))
            return 1;
        edgeNode = edgeNode->next;
    }
    return 0;
}


void DFSforCycle(Graph* g, TriangleList* trianglelist, TetraList* tetraList, EdgeList* pathEdgeList, int* visited, int n, int vert, int start, int* count, int triangles){
    visited[vert] = 1;
    // if the path of length (n-1) is found
    if (n == 0) {
        // mark vert as un-visited to make it usable again.
        visited[vert] = 0;

        // salvare le triplette/quadruple nelle liste
        if (verifyEdgeExistence(vert, start, g->edges))
        {

            EdgeListNode* edgeNode = newEdgeListNode(vert, start);
            edgeNode->next = pathEdgeList->head;
            pathEdgeList->head= edgeNode;
            *count = *count + 1;
            if(triangles){

                EdgeListNode * edgeNode1 = pathEdgeList->head;
                TriangleListNode* node = newTriangleListNode(edgeNode1->edge, edgeNode1->next->edge, edgeNode1->next->next->edge);

                /*node->next = trianglelist->head;
                trianglelist->head = node;*/
                insertNodeInTriangleList(node, trianglelist);
                popFromEdgeList(pathEdgeList);
                /*EdgeListNode *  ptr = pathEdgeList->head;
                pathEdgeList->head =pathEdgeList->head->next;
                free(ptr);*/
            }
            else {
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
     if (!visited[i] && verifyEdgeExistence(vert,i,g->edges)) {
        /* EdgeListNode* edgeNode = newEdgeListNode(i, vert);
         edgeNode->next = pathEdgeList->head;
         pathEdgeList->head= edgeNode;*/

            //ci perdiamo le info degli altri ma tanto per ora vediamo solo il count
            insertInEdgeList(pathEdgeList, vert, i);
            //y decreasing length by 1
            DFSforCycle(g, trianglelist, tetraList, pathEdgeList, visited, n - 1, i, start, count, triangles);
            //backtracking -> si rimuove l'ultimo arco inserito
            popFromEdgeList(pathEdgeList);
         }
     // marking vert as unvisited to make it usable again.
     visited[vert] = 0;

 }

 void findTriangles(Graph* g, int update, TriangleList* trianglelist)
 {
     if(!update) {
         EdgeList *edgeList = (EdgeList* ) malloc(sizeof(EdgeList));
         int visited[g->V];
         int n = 3, count = 0, triangles = 1;

         for (int i = 0; i < g->V; i++)
             visited[i] = 0;

         for (int i = 0; i < g->current_nodes - (n - 1); i++) {
             DFSforCycle(g, trianglelist, NULL, edgeList, visited, n - 1 , i, i, &count,triangles);
             visited[i] = 1;
         }
         free(edgeList);
         //printf("%d\n", count/2 );
     }
     else{
        return;
     }
 }

 void insertTriangles(Graph* g, TriangleList* tList){
    g->triangles = tList;
    EdgeListNode* node = g->edges->head;
    TriangleListNode* triangleNode = tList->head;
    while(triangleNode){
        Triangle* t = triangleNode->triangle;
        while(node)
        {
            if(compareEdge(node->edge,t->e1) || compareEdge(node->edge,t->e2) || compareEdge(node->edge,t->e3))
                node->triangle = t;

            node= node->next;

        }
        triangleNode = triangleNode->next;

    }


}

void insertTetragons(Graph* g, TetraList* tList){
    g->tetragons = tList;
    EdgeListNode* node = g->edges->head;
    TetraListNode* tetraNode = tList->head;
    while(tetraNode) {
        while (node) {
            Tetragon *t = tetraNode->tetragon;
            if (compareEdge(node->edge, t->e1) || compareEdge(node->edge, t->e2) || compareEdge(node->edge, t->e3) ||
                compareEdge(node->edge, t->e4))
                if(node->tetragon == NULL)
                    node->tetragon = t;

            node = node->next;
        }
        tetraNode = tetraNode->next;
    }
}

void findTetragons(Graph* g, int update, TetraList* tetraList){

    if(!update)
    {
        EdgeList *edgeList = (EdgeList *) malloc(sizeof(EdgeList));
        int visited[g->V];
        int n = 4, count = 0, tetra = 0;

        for (int i = 0; i < g->V; i++)
            visited[i] = 0;

        for (int i = 0; i < g->current_nodes - (n - 1); i++) {
            DFSforCycle(g, NULL, tetraList, edgeList, visited, n - 1 , i, i, &count,tetra);
            visited[i] = 1;
        }
        free(edgeList);
        //printf("%d\n", count/2 );

    }
    else return;

}

 int extendIrreducibleGraphWithAdjDiamond(Graph* graphSrc, Graph* graphDest, Edge* edge){
     copyGraph(graphDest, graphSrc);
     return adjDiamondInsertion(graphDest, edge->src, edge->dest);
 }

 int extendIrreducibleGraphWithLollipop(Graph* graphSrc, Graph* graphDest, Edge* edge){
    copyGraph(graphDest, graphSrc);
    printf("src:%d\n", edge->src);
    return lollipopInsertion(graphDest, edge->src, edge->dest);
    //return graphDest == graphSrc;
}

int extendGraphWithEdgeInsertion(Graph* graphSrc, Graph* graphDest, Edge* edge1, Edge* edge2){
    copyGraph(graphDest, graphSrc);
    return edgeInsertion(graphDest,edge1->src,edge1->dest,edge2->src,edge2->dest);
}

int extendIrreducibleGraphWithNonAdjDiamond(Graph* graphSrc, Graph* graphDest, Edge* edge1, Edge* edge2){
    copyGraph(graphDest, graphSrc);
    return nonAdjDiamondInsertion(graphDest, edge1->src, edge1->dest, edge2->src, edge2->dest);
    //return graphDest == graphSrc;
}

int check_graph_existence(PrimeGraphTreeNode* treeNode, Graph* g){
    int find = 0;
    sparsegraph sg1, sg2;
    options.getcanon = TRUE;

    if (treeNode->graph->current_nodes == g->current_nodes) {
            SG_INIT(sg1);
            SG_INIT(sg2);
            SG_ALLOC(sg1, MAX_VERTICES, MAX_DEGREE * MAX_VERTICES, "malloc");
            SG_ALLOC(sg2, MAX_VERTICES, MAX_DEGREE * MAX_VERTICES, "malloc");
            copy_sparse_graph(g);
            nauty((graph *) &sg, lab, ptn, NULL, orbits, &options, &stats, workspace, WORKSIZE, MAXM, g->current_nodes, (graph *) &sg1);
            copy_sparse_graph(treeNode->graph);
            nauty((graph *) &sg, lab, ptn, NULL, orbits, &options, &stats, workspace, WORKSIZE, MAXM, g->current_nodes, (graph *) &sg2);
            find = aresame_sg(&sg1, &sg2);
            SG_FREE(sg1);
            SG_FREE(sg2);
        }
    else
        find = 0;

    // chiamata ricorsiva
    for(int i=0; i<treeNode->num_children; i++) {
        find = find || check_graph_existence(treeNode->children[i], g);
    }
    return find;
}

void generateChildrenWithAdjDiamond(Graph* currentGraph, PrimeGraphTreeNode* treeNode, PrimeGraphTree* tree) {
    printf("Genero figli adj diamond...\n");
    EdgeListNode *currentEdgeNode = currentGraph->edges->head;
    while (currentEdgeNode) {
        Graph *graphChild = createGraph(MAX_VERTICES); // 50 c'era prima!
        if (!extendIrreducibleGraphWithAdjDiamond(currentGraph, graphChild, currentEdgeNode->edge))
            continue;
        printf("Possibile figlio del grafo con %d nodi -> grafoChild con %d nodi \n", currentGraph->current_nodes,
               graphChild->current_nodes);
        if (isAIrreducibleGraph(graphChild)) {
            int find = 0;
            find = check_graph_existence(tree->head, graphChild);
            printf("find = %d\n", find);
            if (!find) {
                addGraphToTree(treeNode, graphChild, tree, 1, "AdjDiamond");
                if(graphChild->current_nodes == MAX_VERTICES)
                    cubicGraphs[++cubicIndex] = graphChild;
            }
        }
        currentEdgeNode = currentEdgeNode->next;
    }
}

void generateChildrenWithNonAdjDiamond(Graph* currentGraph, PrimeGraphTreeNode* treeNode, PrimeGraphTree* tree) {
    printf("Genero figli non adj diamond...\n");
    EdgeListNode *currentEdgeNode1 = currentGraph->edges->head;
    EdgeListNode *currentEdgeNode2 = currentEdgeNode1;

    while (currentEdgeNode1) {
        currentEdgeNode2 = currentEdgeNode1->next; //parti dal successivo
        while (currentEdgeNode2) {
            Graph *graphChild = createGraph(MAX_VERTICES);
            if (!extendIrreducibleGraphWithNonAdjDiamond(currentGraph, graphChild, currentEdgeNode2->edge, currentEdgeNode1->edge))
                continue;
            printf("Possibile figlio del grafo con %d nodi -> grafoChild con %d nodi \n", currentGraph->current_nodes,
                   graphChild->current_nodes);
            if (isAIrreducibleGraph(graphChild)) {
                int find = 0;
                find = check_graph_existence(tree->head, graphChild);
                printf("find = %d\n", find);
                if (!find) {
                    addGraphToTree(treeNode, graphChild, tree, 1, "NoAdjDiamond");
                    if(graphChild->current_nodes == MAX_VERTICES)
                        cubicGraphs[++cubicIndex] = graphChild;
                }
            }
            currentEdgeNode2 = currentEdgeNode2->next;
        }
        currentEdgeNode1 = currentEdgeNode1->next;
    }
}

void generateChildrenWithLollipop(Graph* currentGraph, PrimeGraphTreeNode* treeNode, PrimeGraphTree* tree) {
    printf("Genero figli con lollipop...\n");
    EdgeListNode* currentEdgeNode = currentGraph->edges->head;
    while (currentEdgeNode) {
            Graph *graphChild = createGraph(MAX_VERTICES); //50 c'era prima!
            if (!extendIrreducibleGraphWithLollipop(currentGraph, graphChild, currentEdgeNode->edge))
                continue;
            printf("Possibile figlio del grafo con %d nodi -> grafoChild con %d nodi \n", currentGraph->current_nodes,
                   graphChild->current_nodes);
            if (isAIrreducibleGraph(graphChild)) {
                int find = 0;
                find = check_graph_existence(tree->head, graphChild);
                printf("find = %d\n", find);
                if (!find) {
                    //printf("è entrato");
                    addGraphToTree(treeNode, graphChild, tree, 1, "Lollipop");
                    if(graphChild->current_nodes == MAX_VERTICES)
                        cubicGraphs[++cubicIndex] = graphChild;
                }
            }
            currentEdgeNode = currentEdgeNode->next;
        }
    }

void recursiveGenFromPrime(PrimeGraphTreeNode* treeNode, PrimeGraphTree* tree){
    printf("Genero figli con edge...\n");
    EdgeListNode* currentEdgeNode1 = treeNode->graph->edges->head;
    EdgeListNode* currentEdgeNode2 = currentEdgeNode1;

    while(currentEdgeNode1) {
        currentEdgeNode2 = currentEdgeNode1->next;
        while (currentEdgeNode2) {
            Graph *graphChild = createGraph(MAX_VERTICES);
            extendGraphWithEdgeInsertion(treeNode->graph, graphChild, currentEdgeNode1->edge, currentEdgeNode2->edge);

            int find = 0;
            find = check_graph_existence(tree->head, graphChild);
            if (!find) {
                addGraphToTree(treeNode, graphChild, tree, 0, "Edge");
                if(graphChild->current_nodes == MAX_VERTICES)
                    cubicGraphs[++cubicIndex] = graphChild;
            }
            currentEdgeNode2 = currentEdgeNode2->next;
        }
        currentEdgeNode1 = currentEdgeNode1->next;
    }

    for(int i=0; i<treeNode->num_children; i++)
        if(treeNode->children[i]->graph->current_nodes <= MAX_VERTICES -2 && !treeNode->children[i]->isPrime)
            recursiveGenFromPrime(treeNode->children[i], tree);

}

void generatePrimeTrees(PrimeGraphTreeNode* treeNode, PrimeGraphTree* tree){
    Graph* currentGraph = treeNode->graph;
    EdgeListNode* currentEdgeNode = currentGraph->edges->head;
    int scelta = 0;
    if (currentGraph->current_nodes <= MAX_VERTICES - 6 ) //prima c'era l'or
        scelta = 2;
    else if(currentGraph->current_nodes <= MAX_VERTICES - 4)
        scelta = 1;

    printf("Grafo con %d nodi, La mia scelta: %d\n",currentGraph->current_nodes, scelta);

    if(scelta == 0)
        return;

    else if (scelta == 1 )
        generateChildrenWithAdjDiamond(currentGraph, treeNode, tree);
    else{
        generateChildrenWithAdjDiamond(currentGraph, treeNode, tree);
        generateChildrenWithNonAdjDiamond(currentGraph, treeNode, tree);
        generateChildrenWithLollipop(currentGraph, treeNode, tree);
    }


    /*if (currentGraph->current_nodes == MAX_VERTICES)
        return;*/

  /* while(currentEdgeNode){
        Graph* graphChild = createGraph(MAX_VERTICES);
        extendIrreducibleGraphWithAdjDiamond(currentGraph, graphChild, currentEdgeNode->edge);
        //extendIrreducibleGraphWithLollipop(currentGraph,graphChild,currentEdgeNode->edge);
        if(isAIrreducibleGraph(graphChild))
            addGraphToTree(treeNode, graphChild, tree);
        currentEdgeNode = currentEdgeNode->next;
   }*/
    // chiamata ricorsiva
    for(int i=0; i<treeNode->num_children; i++)
        generatePrimeTrees(treeNode->children[i], tree);
}

void generateCubicGraphTree(PrimeGraphTreeNode* treeNode, PrimeGraphTree* tree){
    if(!treeNode->isPrime)
        return;
    if(treeNode->graph->current_nodes > MAX_VERTICES - 2)
        return;
    recursiveGenFromPrime(treeNode, tree);

    for(int i=0; i<treeNode->num_children; i++)
        generateCubicGraphTree(treeNode->children[i], tree);
}

void addGraphToTree(PrimeGraphTreeNode* parent, Graph* g, PrimeGraphTree* tree, int isPrime, char* op)
 {
     PrimeGraphTreeNode* newNode = newPrimeGraphTreeNode(tree, g, g->current_nodes, isPrime, op);
     newNode->parent = parent;
     parent->num_children++;
     //parent->children = realloc( parent->children, (parent->num_children * sizeof(PrimeGraphTreeNode))); // la realloc funge???? controllare meglio!
     parent->children[parent->num_children-1] = newNode;

 }

void print_sparse_graph_nauty(sparsegraph sparse_graph) {
    int i, j;
    fprintf(stderr, "Printing sparse graph nauty:\n");
    for(i = 0; i < sparse_graph.nv; i++) {
        fprintf(stderr, "%d :", i);
        for(j = 0; j < sparse_graph.d[i];j++) {
            //fprintf(stderr, " %d", sparse_graph.e[i * REG + j]);
            fprintf(stderr, " %d", sparse_graph.e[sparse_graph.v[i] + j]);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "Number of directed edges: %lu\n", (unsigned long) sparse_graph.nde);
}
void copy_sparse_graph(Graph* graph) {

    sg.nv = graph->current_nodes;
    sg.nde = 3 * graph->current_nodes;

    int i, j;
    for(i = 0; i < graph->V; i++) {
        //These values were already set in init_nauty_options()
        //sg.v[i] = i * REG;
        //sg.d[i] = 3;
        struct AdjListNode* adjListNode = graph->nodes[i].head;
        for(j = 0; j < 3 && adjListNode != NULL ; j++) {
            sg.e[i * 3 + j] = adjListNode->label;
            adjListNode = adjListNode->next;
        }
    }
}

Graph* generateUniform(){
    int i = 0;
    int max = cubicIndex, min = 1;
    int r = rand() % (max + 1 - min) + min;

    return cubicGraphs[r];
}
void recursivePrintEdgeTreeForVis(PrimeGraphTreeNode* node) {
    if (node == NULL) {
        return;
    }

    /* Caso radice che non ha parent */
    if (node->parent == NULL)
        printf("%d%s#:%d,\n", node->id, node->op, node->graph->current_nodes);
    else{
        printf("%d%s#:%d,", node->parent->id, node->parent->op, node->parent->graph->current_nodes);
        printf("%d%s#:%d\n", node->id, node->op, node->graph->current_nodes);
    }

    for(int i = 0; i < node->num_children; i++) {
        recursivePrintEdgeTreeForVis(node->children[i]);
    }
}

void recursivePrintNodesTreeForVis(PrimeGraphTreeNode* node) {
    if (node == NULL) {
        return;
    }

    printf("%d%s#:%d,%d%s#:%d\n", node->id, node->op, node->graph->current_nodes,node->id, node->op, node->graph->current_nodes);


    for(int i = 0; i < node->num_children; i++) {
        recursivePrintNodesTreeForVis(node->children[i]);
    }
}
int main() {
    cubicGraphs = (Graph**) malloc(10000 * sizeof(Graph*));
    srand(time(NULL));
    //For the timing
    struct tms TMS;
    unsigned int oldtime = 0;
    init_nauty_options();
    Graph* g = createGraph(50);
    PrimeGraphTree* primeTree = createPrimeGraphTree();
    initIrreducibleGraph(g);
    char* op = "K4";
    PrimeGraphTreeNode* node = newPrimeGraphTreeNode(primeTree, g, 4, 1, op);
    primeTree->head = node;
    generatePrimeTrees(node, primeTree);
    generateCubicGraphTree(node, primeTree);
    int parent = 0;
    int id = 0;
    recursivePrintEdgeTreeForVis(primeTree->head);
    printf("--------------------------------------------------------------------\n");
    recursivePrintNodesTreeForVis(primeTree->head);
    printf("--------------------------------------------------------------------\n");
    printf("Numero totale di grafi cubici generati: %d\n", primeTree->current_graphs);
    printf("--------------------------------------------------------------------\n");
    printf("Numero di grafi cubici generati con vertici %d : %d\n", MAX_VERTICES, cubicIndex);
    printf("--------------------------------------------------------------------\n");
    printf("Grafo random:\n" );
    printEdgeList(generateUniform());
    printf("--------------------------------------------------------------------\n");


    times(&TMS);
    unsigned long savetime = oldtime + (unsigned int) TMS.tms_utime;
    fprintf(stderr, "CPU time: %.1f seconds.\n", (double) savetime / time_factor);

    return 0;
}