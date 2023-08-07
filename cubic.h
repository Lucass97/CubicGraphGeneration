// File: cubic.h
// Author: Luca Gregori & Alessandro Wood
// Date: 14th July 2021
// Description: Contains functions definitions for graph operations and manipulation and related data structures.

#ifndef CUBICGRAPHGENERATION_CUBIC_H
#define CUBICGRAPHGENERATION_CUBIC_H

#include "nausparse.h"
#include <string.h>

#define MAX_DEGREE 3
#define MAX_VERTICES 16

/* ==================================== */
/*|          Data structures           |*/
/* ==================================== */

typedef struct Edge {
    int src;
    int dest;
} Edge;

typedef struct Triangle {
    Edge* e1;
    Edge* e2;
    Edge* e3;
} Triangle;

typedef struct Tetragon {
    Edge* e1;
    Edge* e2;
    Edge* e3;
    Edge* e4;
} Tetragon;

typedef struct TetraListNode {
    Tetragon* tetragon;
    struct TetraListNode* next;
} TetraListNode;

typedef struct TetraList {
    struct TetraListNode* head;
} TetraList;

typedef struct TriangleListNode {
    Triangle* triangle;
    struct TriangleListNode* next;
} TriangleListNode;

typedef struct TriangleList {
    struct TriangleListNode* head;
} TriangleList;

typedef struct AdjListNode {
    int label;
    struct AdjListNode* next;
} AdjListNode;

typedef struct AdjList {
    struct AdjListNode* head;
} AdjList;

typedef struct EdgeListNode {
    Edge* edge;
    int eligible;
    struct EdgeListNode* next;
    Triangle* triangle;
    Tetragon* tetragon;
} EdgeListNode;

typedef struct EdgeList
{
    struct EdgeListNode* head;
}EdgeList;

// Graph Structure

typedef struct Graph {
    int current_nodes;
    int current_edges;
    int V;
    struct AdjList* nodes;
    struct EdgeList* edges;
    struct TriangleList* triangles;
    struct TetraList* tetragons;
} Graph;

// Prime Graph Tree Node Structure

typedef struct PrimeGraphTreeNode{
    Graph* graph;
    int num_children;
    int isPrime;
    int id;
    char op[15];
    struct PrimeGraphTreeNode** children;
    struct PrimeGraphTreeNode* parent;
} PrimeGraphTreeNode;

typedef struct PrimeGraphTree{
    int current_graphs;
    int current_max_n;
    struct PrimeGraphTreeNode* head;
} PrimeGraphTree;

/* ==================================== */
/*| Functions for Edge,                |*/
/*| Triangle, and Tetragon Nodes       |*/
/* ==================================== */

struct TetraListNode* newTetraListNode(Edge* e1, Edge* e2, Edge* e3,  Edge* e4){
    TetraListNode* newNode = ( TetraListNode*) malloc(sizeof(TetraListNode));
    newNode->tetragon = (Tetragon *) malloc(sizeof(Tetragon));
    newNode->tetragon->e1 = e1;
    newNode->tetragon->e2 = e2;
    newNode->tetragon->e3 = e3;
    newNode->tetragon->e4 = e4;

    newNode->next =NULL;
    return newNode;
};

struct TriangleListNode* newTriangleListNode(Edge* e1, Edge* e2, Edge* e3){
    TriangleListNode* newNode = (TriangleListNode*) malloc(sizeof(TriangleListNode));
    newNode->triangle = (Triangle*) malloc(sizeof(Triangle));
    newNode->triangle->e1 = e1;
    newNode->triangle->e2 = e2;
    newNode->triangle->e3 = e3;

    newNode->next =NULL;

    return newNode;
};

struct AdjListNode* newAdjListNode(int dest){
    struct AdjListNode* newNode = (struct AdjListNode*) malloc(sizeof(struct AdjListNode));
    newNode->label = dest;
    newNode->next = NULL;
    return newNode;
}

struct EdgeListNode* newEdgeListNode(int src, int dest){
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
/* |   Graph Creation and Operations   |*/
/* ==================================== */

struct Graph* createGraph(int V){
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

struct PrimeGraphTree* createPrimeGraphTree(){
    struct PrimeGraphTree* graphList = (PrimeGraphTree*) malloc(sizeof(PrimeGraphTree));
    graphList->head = NULL;
    graphList->current_graphs =0;
    graphList->current_max_n = 0;
    return graphList;
}

struct PrimeGraphTreeNode* newPrimeGraphTreeNode(PrimeGraphTree* tree, Graph* g, int n, int isPrime, char* op){
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
/* |             Functions             |*/
/* ==================================== */

// Nauty and Sparse Graph Functions
void init_nauty_options();
void copy_sparse_graph(Graph* graph);
void print_sparse_graph_nauty(sparsegraph sparse_graph);

// Basic Graph Operations
int countVertexDegree(Graph* graph, int label);
void removeNodeFromAdjList(Graph* graph, int src, int dest);
void addEdge(Graph* graph, int src, int dest);
void removeEdge(Graph* graph, int src, int dest);
void removeEdgeFromEdgeList(EdgeList* edgeList, int src, int dest);
void insertInEdgeList(EdgeList* edgeList, int src, int dest);
void popFromEdgeList(EdgeList* edgeList);
int compareEdge(Edge* e1, Edge* e2);
int verifyEdgeExistence(int src, int dest, EdgeList* edgeList);

// Graph Manipulation
int edgeInsertion(Graph* graph, int src1, int dest1, int src2, int dest2);
int nonAdjDiamondInsertion(Graph* graph, int src1, int dest1, int src2, int dest2);
int adjDiamondInsertion(Graph* graph, int src, int dest);
int lollipopInsertion(Graph* graph, int src, int dest);

// Irreducibility Checks
void dfsForBridge(Graph* graph, int start, int* visited);
int isABridgeList(Graph* graph);
int isABridge(Graph* graph, Edge* edge);
int irreducibilityCondition2(Graph* graph, Edge* edge);
int irreducibilityCondition3(Graph* graph, Edge* edge);
int isAIrreducibleGraph(Graph* g);
void initIrreducibleGraph(Graph* graph);

// Triangle and Tetragon Operations
void DFSforCycle(Graph* g, TriangleList* trianglelist, TetraList* tetraList, EdgeList* pathEdgeList, int* visited, int n, int vert, int start, int* count, int triangles);
void findTriangles(Graph* g, int update, TriangleList* trianglelist);
void findTetragons(Graph* g, int update, TetraList* tetraList);
int compareTriangle(Triangle* t1, Triangle* t2);
int compareTetragon(Tetragon* t1, Tetragon* t2);
int isATriangle(Triangle* t);
int isATetragon(Tetragon* t);
int verifyTriangleExistence(Triangle* triangle, TriangleList* triangleList);
int verifyTetragonExistence(TetraList* tetraList, Tetragon* tetra);
void insertTriangles(Graph* g, TriangleList* tList);
void insertTetragons(Graph* g, TetraList* tList);
void insertNodeInTriangleList(TriangleListNode* node, TriangleList* triangleList);
void insertNodeInTetraList(TetraList* tetraList, TetraListNode* node);

// Copy Graph
void copyGraph(Graph* graphDest, Graph* graphSrc);

// Graph Tree Operations
void addGraphToTree(PrimeGraphTreeNode* parent, Graph* g, PrimeGraphTree* tree, int isPrime, char* op);
int check_graph_existence(PrimeGraphTreeNode* treeNode, Graph* g);
int extendIrreducibleGraphWithAdjDiamond(Graph* graphSrc, Graph* graphDest, Edge* edge);
int extendIrreducibleGraphWithNonAdjDiamond(Graph* graphSrc, Graph* graphDest, Edge* edge1, Edge* edge2);
int extendIrreducibleGraphWithLollipop(Graph* graphSrc, Graph* graphDest, Edge* edge);
void generateChildrenWithAdjDiamond(Graph* currentGraph, PrimeGraphTreeNode* treeNode, PrimeGraphTree* tree);
void generateChildrenWithNonAdjDiamond(Graph* currentGraph, PrimeGraphTreeNode* treeNode, PrimeGraphTree* tree);
void generateChildrenWithLollipop(Graph* currentGraph, PrimeGraphTreeNode* treeNode, PrimeGraphTree* tree);
void recursiveGenFromPrime(PrimeGraphTreeNode* treeNode, PrimeGraphTree* tree);
void generatePrimeTreesFullRandom(PrimeGraphTreeNode* treeNode, PrimeGraphTree* tree);
void generateCubicGraphTree(PrimeGraphTreeNode* treeNode, PrimeGraphTree* tree);
int getSpareVertices(Graph* graph, int* spareVertices, int size);

// Random Graph Generation
Graph* generateUniform();

// Printing Functions
void printEdgeList(Graph* graph);
void printTriangleList(TriangleList* triangleList);
void printTriangleListAlt(Graph* g);
void printTetragonList(TetraList* tetraList);
void printTetragonListAlt(Graph* g);
void printGraph(Graph* graph);

// Tree Saving
void recursivePrintNodesTree(PrimeGraphTreeNode* node, FILE* file);
void saveTreeToCSV(PrimeGraphTreeNode* root, const char* filename);
void recursivePrintEdgeTree(PrimeGraphTreeNode* node, FILE* file);
void saveEdgeTreeToCSV(PrimeGraphTreeNode* root, const char* filename);


#endif //CUBICGRAPHGENERATION_CUBIC_H
