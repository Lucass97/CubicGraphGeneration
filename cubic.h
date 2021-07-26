//
// Created by luca and the great leader on 14/07/21.
//

#ifndef CUBICGRAPHGENERATION_CUBIC_H
#define CUBICGRAPHGENERATION_CUBIC_H

#include "nausparse.h"
#include <string.h>

#define MAX_DEGREE 3
#define MAX_VERTICES 14

/* ==================================== */
/* Strutture */


typedef struct Edge{
    int src;
    int dest;
}Edge;

typedef struct Triangle{
    Edge* e1;
    Edge* e2;
    Edge* e3;
}Triangle;

typedef struct Tetragon{
    Edge* e1;
    Edge* e2;
    Edge* e3;
    Edge* e4;
}Tetragon;

typedef struct TetraListNode{
    Tetragon* tetragon;
    struct TetraListNode* next;
}TetraListNode;

typedef struct TetraList
{
    struct TetraListNode *head;
}TetraList;

typedef struct TriangleListNode{
    Triangle* triangle;
    struct TriangleListNode* next;

}TriangleListNode;

typedef struct TriangleList
{
    struct TriangleListNode *head;
}TriangleList;


struct TetraListNode* newTetraListNode(Edge* e1, Edge* e2, Edge* e3,  Edge* e4)
{
    TetraListNode* newNode = ( TetraListNode*) malloc(sizeof(TetraListNode));
    newNode->tetragon = (Tetragon *) malloc(sizeof(Tetragon));
    newNode->tetragon->e1 = e1;
    newNode->tetragon->e2 = e2;
    newNode->tetragon->e3 = e3;
    newNode->tetragon->e4 = e4;

    newNode->next =NULL;
    return newNode;
};
struct TriangleListNode* newTriangleListNode(Edge* e1, Edge* e2, Edge* e3)
{
    TriangleListNode* newNode = (TriangleListNode*) malloc(sizeof(TriangleListNode));
    newNode->triangle = (Triangle*) malloc(sizeof(Triangle));
    newNode->triangle->e1 = e1;
    newNode->triangle->e2 = e2;
    newNode->triangle->e3 = e3;

    newNode->next =NULL;

    return newNode;
};
typedef struct AdjListNode
{
    int label;
    struct AdjListNode* next;
}AdjListNode;

typedef struct AdjList
{
    struct AdjListNode *head;
}AdjList;

struct AdjListNode* newAdjListNode(int dest)
{
    struct AdjListNode* newNode = (struct AdjListNode*) malloc(sizeof(struct AdjListNode));
    newNode->label = dest;
    newNode->next = NULL;
    return newNode;
}

/* ==================================== */
/* Strutture per edgeList */

typedef struct EdgeListNode
{
    Edge* edge;
    int eligible;
    struct EdgeListNode* next;
    Triangle* triangle;
    Tetragon* tetragon;
}EdgeListNode;

typedef struct EdgeList
{
    struct EdgeListNode* head;
}EdgeList;


struct EdgeListNode* newEdgeListNode(int src, int dest)
{
    struct EdgeListNode* newNode = (struct EdgeListNode*) malloc(sizeof(struct EdgeListNode));
    Edge* edge = (Edge*) malloc(sizeof(Edge));
    edge->src = src;
    edge->dest = dest;
    newNode->edge = edge;
    newNode->eligible = 1;
    newNode->next = NULL;
    newNode->triangle = NULL;
    newNode->tetragon = NULL;
    return newNode;
}

/* ==================================== */
/* Strutture per grafi */

typedef struct Graph
{
    int current_nodes;
    int current_edges;
    int V;
    struct AdjList* nodes;
    struct EdgeList* edges;
    struct TriangleList* triangles;
    struct TetraList* tetragons;
}Graph;


struct Graph* createGraph(int V)
{
    struct Graph* graph = (struct Graph*) malloc(sizeof(struct Graph));
    graph->V = V;
    graph->nodes = (struct AdjList*) malloc(V * sizeof(struct AdjList));
    graph->edges = (struct EdgeList*) malloc(sizeof(struct EdgeList));

    int i;
    for (i = 0; i < V; ++i)
        graph->nodes[i].head = NULL;
    graph->edges->head = NULL;
    graph->current_edges = 0;
    graph->current_edges = 0;
    graph->triangles = NULL;
    graph->tetragons = NULL;

    return graph;
}

typedef struct PrimeGraphTreeNode
{
    Graph* graph;
    int num_children;
    int isPrime;
    int id;
    char op[15];
    struct PrimeGraphTreeNode** children;
    struct PrimeGraphTreeNode* parent;
}PrimeGraphTreeNode;

typedef struct PrimeGraphTree
{
    int current_graphs;
    int current_max_n;
    struct PrimeGraphTreeNode* head;
}PrimeGraphTree;

struct PrimeGraphTree* createPrimeGraphTree()
{
    struct PrimeGraphTree* graphList = (PrimeGraphTree*) malloc(sizeof(PrimeGraphTree));
    graphList->head = NULL;
    graphList->current_graphs =0;
    graphList->current_max_n = 0;
    return graphList;
}
struct PrimeGraphTreeNode* newPrimeGraphTreeNode(PrimeGraphTree* tree, Graph* g, int n, int isPrime, char* op)
{
    struct PrimeGraphTreeNode* newNode = (PrimeGraphTreeNode*) malloc(sizeof(PrimeGraphTreeNode));
    newNode->graph = g;
    tree->current_graphs++;
    newNode->id = tree->current_graphs;
    if (n > tree->current_max_n)
        tree->current_max_n = n;
    newNode->children = (PrimeGraphTreeNode**) malloc(50 * sizeof(PrimeGraphTreeNode));// era NULL ma avevamo una realloc ogni volta.
    newNode->num_children = 0;
    newNode->parent = NULL;
    newNode->isPrime = isPrime;
    strcpy(newNode->op, op);
    return newNode;
}

/* ==================================== */
/* Funzioni */
int countVerticeDegree(Graph* graph, int label);
void removeNodeFromAdjList(Graph* graph, int src, int dest);
void addEdge(Graph* graph, int src, int dest);
void removeEdge(Graph* graph, int src, int dest);
void removeEdgeFromEdgeList(EdgeList* edgeList, int src, int dest);
void popFromEdgeList(EdgeList* edgeList);
void insertInEdgeList(EdgeList* edgeList, int src, int dest);
int compareEdge(Edge* e1, Edge* e2);
int verifyEdgeExistence(int src, int dest, EdgeList* edgeList);

int compareTriangle(Triangle* t1, Triangle* t2);
int verifyTriangleExistence(Triangle* triangle, TriangleList* triangleList);

void dfsForBridge(Graph* graph, int start, int* visited);
int isABridge(Graph* graph, Edge* edge);
int isABridgeList(Graph* graph);
void printEdgeList(Graph* graph);
void printTriangle(TriangleList* triangleList);
void printTriangleNew(Graph* g);
void printTetragon(TetraList* tetraList);
void printTetragonNew(Graph* g);
void printGraph(Graph* graph);
int getUnlabeledVertices(Graph* graph, int* unlabeledVertices, int size);
int edgeInsertion(Graph* graph, int src1, int dest1, int src2, int dest2);
int nonAdjDiamondInsertion(Graph* graph, int src1, int dest1, int src2, int dest2);
int adjDiamondInsertion(Graph* graph, int src, int dest);
int lollipopInsertion(Graph* graph, int src, int dest);
EdgeListNode* chooseRandomEligible(Graph* graph);
int irreducibilityCondition2(Graph* graph, Edge* edge);
int irreducibilityCondition3(Graph* graph, Edge* edge);
int isAIrreducibleGraph(Graph* g);
void initIrreducibleGraph(Graph* graph);
void copyGraph(Graph* graphDest, Graph* graphSrc);
void DFSforCycle(Graph* g, TriangleList* trianglelist, TetraList* tetraList, EdgeList* pathEdgeList, int* visited, int n, int vert, int start, int* count, int triangles);
void findTriangles(Graph* g, int update, TriangleList* trianglelist);
void insertTriangles(Graph* g, TriangleList* tList);
void insertTetragons(Graph* g, TetraList* tList);
void findTetragons(Graph* g, int update, TetraList* tetraList);
void insertNodeInTriangleList(TriangleListNode * node, TriangleList* triangleList);
int compareTetragon(Tetragon* t1, Tetragon* t2);
int isATetragon(Tetragon* t);
int verifyTetragonExistence(TetraList* tetraList, Tetragon* tetra);
void insertNodeInTetraList(TetraList* tetraList, TetraListNode * node);

void addGraphToTree(PrimeGraphTreeNode* parent, Graph* g, PrimeGraphTree* tree, int isPrime, char* op);
int check_graph_existence(PrimeGraphTreeNode* treeNode, Graph* g);
int extendIrreducibleGraphWithAdjDiamond(Graph* graphSrc, Graph* graphDest, Edge* edge);
int extendIrreducibleGraphWithNonAdjDiamond(Graph* graphSrc, Graph* graphDest, Edge* edge1, Edge* edge2);
int extendIrreducibleGraphWithLollipop(Graph* graphSrc, Graph* graphDest, Edge* edge);

void generateChildrenWithAdjDiamond(Graph* currentGraph, PrimeGraphTreeNode* treeNode, PrimeGraphTree* tree);
void generateChildrenWithNonAdjDiamond(Graph* currentGraph, PrimeGraphTreeNode* treeNode, PrimeGraphTree* tree);
void generateChildrenWithLollipop(Graph* currentGraph, PrimeGraphTreeNode* treeNode, PrimeGraphTree* tree);
void generatePrimeTrees(PrimeGraphTreeNode* treeNode, PrimeGraphTree* tree);

/* Funzione varie per nauty e grafi sparsi*/
void copy_sparse_graph(Graph* graph);
void print_sparse_graph_nauty(sparsegraph sparse_graph);

void init_nauty_options();

#endif //CUBICGRAPHGENERATION_CUBIC_H
