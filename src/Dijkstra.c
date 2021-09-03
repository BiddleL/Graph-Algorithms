// Dijkstra API implementation
// COMP2521 Assignment 2
// z5311885

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "Dijkstra.h"
#include "Graph.h"
#include "PQ.h"

// function 
static ShortestPaths setupSPS(Graph g, Vertex src);
static PredNode* insertPred(PredNode *head, PredNode *insert);
static void freePredList(PredNode *head);
static PredNode *newPred(Vertex v);
static void printPred(PredNode *node);

// this function is an implementation of dijkstra's algo, using a priority queue
// it returns a struct ShortestPaths which contains everything you could need
// for finding paths, including the shortest path in a graph from a source vertex
ShortestPaths dijkstra(Graph g, Vertex src) {
	if(g == NULL) {
		fprintf(stderr, "Invalid graph\n");
	} else if (src < 0 || src >= GraphNumVertices(g)) {
		fprintf(stderr, "Invalid source vertex\n");
	}
	
	ShortestPaths sps = setupSPS(g, src);
	PQ q = PQNew();

	// O(numNodes)
	// Adding all outNodes to PQ
	for(int i = 0; i < sps.numNodes; i++) {
		PQInsert(q, i, sps.dist[i]);
	}
	
	// O(nV^2)
	while(!PQIsEmpty(q)) {
		Vertex u = PQDequeue(q);
		AdjList curr = GraphOutIncident(g, u);
		
		// u is the current node
		// curr is the viewing node
		while(curr != NULL) {
			// calculating the alt path length
			int alt = sps.dist[u] + curr->weight;
			// checking that there isn't an integer overflow
			if(sps.dist[u] != INFINITY) {
				if(alt < sps.dist[curr->v]) {
					sps.dist[curr->v] = alt;
					// replacing the path as there is a shorter one
					freePredList(sps.pred[curr->v]);
					sps.pred[curr->v] = insertPred(NULL , newPred(u));
					PQUpdate(q, curr->v, alt);
				} else if (alt == sps.dist[curr->v]) {
					sps.pred[curr->v] = insertPred(sps.pred[curr->v], newPred(u));
					PQUpdate(q, curr->v, alt);
				}
			}
			curr = curr->next;
		}
		
	}
	PQFree(q);

	return sps;
}

void showShortestPaths(ShortestPaths sps) {
	printf("Number of Node: %d\n", sps.numNodes);
	printf("Source Vertex: %d\n", sps.src);
	printf("=== Distance Array ===\n");
	for(int i = 0; i < sps.numNodes; i++) {
		printf("%d -> %d:  ", sps.src, i);
		if(sps.dist[i] == INFINITY) {
			printf("INF\n");
		} else if (sps.dist[i] >= 0) {
			printf("%d\n", sps.dist[i]);
		}
	}
	for(int i = 0; i < sps.numNodes; i++) {
		printf("%d : ", i);
		printPred(sps.pred[i]);
	}
}

void freeShortestPaths(ShortestPaths sps) {
	for(int i = 0; i < sps.numNodes; i++) {
		freePredList(sps.pred[i]);
	}
	free(sps.pred);
	free(sps.dist);
}




// ========================
// === HELPER FUNCTIONS ===
// ========================

// This function sets up the ShortestPaths struct
// O(numNodes)
static ShortestPaths setupSPS(Graph g, Vertex src) {
	ShortestPaths sps;
	sps.numNodes = GraphNumVertices(g);
	sps.src = src;
	sps.dist = malloc(sps.numNodes * sizeof(int));
	sps.pred = malloc(sps.numNodes * sizeof(PredNode *));
	if(sps.dist != NULL && sps.pred != NULL) {
		for(int i = 0; i < sps.numNodes; i++) {
				sps.dist[i] = INFINITY;
				sps.pred[i] = NULL;
		}
		sps.dist[src] = 0;
	}
	return sps;
}

// this function just inserts a Prednode into a linked list
// a bit pointless but it works
// O(length of linked list)
static PredNode* insertPred(PredNode *head, PredNode *insert) {
	if (head == NULL) {
		return insert;
	} else { 
		insert->next = head;
		return insert;
	}
}

// recursive free function
// O(1)
static void freePredList(PredNode *head) {
    if (head == NULL) {
        return;
    }
    freePredList(head->next);
    free(head);
}

// sets up and return a new PredNode 
// O(1)
static PredNode *newPred(Vertex v) {
	PredNode *newNode = malloc(sizeof(struct PredNode));
    newNode->v = v;
    newNode->next = NULL;
    return newNode;
}
// prints PredNodes
static void printPred(PredNode *node) {
	if (node == NULL) {
		printf("NULL\n");
	} else {
		printf("[%d]->", node->v);
	}
	printPred(node->next);
}