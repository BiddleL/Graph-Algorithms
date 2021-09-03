// Lance-Williams Algorithm for Hierarchical Agglomerative Clustering
// COMP2521 Assignment 2

#include <assert.h>
#include <float.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#include "Graph.h"
#include "LanceWilliamsHAC.h"

#define INFINITY DBL_MAX


// === DEV DEFINES ===
#define v(node) node->v
#define weight(node) node->weight
#define next(node) node->next

#define left(node) node->left
#define right(node) node->right
#define vertex(node) node->vertex

// === HELPER FUNCTIONS ===
static Dendrogram newDend(int v);
static double **calDistance(Graph g, double **dist);
//static double **removeCluster(Graph g, double **dist, int v1, int v2);

/**
 * Generates  a Dendrogram using the Lance-Williams algorithm (discussed
 * in the spec) for the given graph  g  and  the  specified  method  for
 * agglomerative  clustering. The method can be either SINGLE_LINKAGE or
 * COMPLETE_LINKAGE (you only need to implement these two methods).
 * 
 * The function returns a 'Dendrogram' structure.
 */
Dendrogram LanceWilliamsHAC(Graph g, int method) {
	if(method != SINGLE_LINKAGE && method != COMPLETE_LINKAGE) {
		fprintf(stderr, "Wrong Method, must be Single or Complete Linkage\n");
		return NULL;
	}

	int nV = GraphNumVertices(g);
	double **dist = malloc(nV * sizeof(double *));
	Dendrogram *dend = malloc(nV * sizeof(Dendrogram));
	for(int i = 0; i =< nV; i++) {
		dist[i] = calloc(nV, sizeof(double));

		// initisation of Dendrogram array
		Dendrogram new = newDend(i);
        dend[i] = new;
	}

	dist = calDistance(g, dist);

	// finding shorter distances 
	for (Vertex v = 1; v < nV; v++) {
		double min = INFINITY;
		int v1, v2 = 0;

		for(Vertex i = 0; i < nV + 1 - v; i++) {
			for(Vertex j = 0; j < nV + 1 - v; j++) {
				if (dist[i][j] <= min && i != j) {
					min = dist[i][j];
					v1 = i;
					v2 = j;
				}
			} 
		}
		Dendrogram d1 = dend[v1];
		Dendrogram d2 = dend[v2];
		// shifting dendrogram array

		for(Vertex i = v1; i < nV - 1; i++) {
			dend[i] = dend[i + 1];
		}
		int second = (v2 > v1 ? v2 - 1: v2);
		for(Vertex i = second; i < nV - 1; i++) {
			dend[i] = dend[i + 1];
		}
		Dendrogram new = newDend(1);
		int index = nV - 1 - v;
		free(dend[index]);
		dend[index] = new;
		vertex(dend[index]) = -1;
		left(dend[index]) = d1;
		right(dend[index]) = d2;
	
		// calulating the distance to the new cluster formed above
		int dex = 0;
		double *newDist = calloc(nV - v, sizeof(double));
		for(Vertex i = 0; i < nV - v + 1; i++) {
			if(i != v1 && i !=v2) {
				if (method == SINGLE_LINKAGE) {
					if (dist[v1][i] < dist[v2][i]) {
						newDist[dex] = dist[v1][i];
					} else {
						newDist[dex] = dist[v2][i];
					}
				} else if (method == COMPLETE_LINKAGE) {
					if (dist[v1][i] < dist[v2][i]) {
						newDist[dex] = dist[v1][i];
					} else {
						newDist[dex] = dist[v2][i];
					}
				}
				dex++;
			}
		}
		newDist[nV - v - 1] = INFINITY;
		
		// removing cluster
		for(int i = v1; i <= nV; i++) {
			for(int j = 0; j < nV; j++) {
				dist[i][j] = dist[i + 1][j];
				dist[j][i] = dist[j][i + 1];
			}
		}
		int offset = (v2 > v1 ? v2 - 1: v2);
		for(int i = offset; i <= nV ; i++) {
			for(int j = 0; j < nV; j++) {
				dist[i][j] = dist[i + 1][j];
				dist[j][i] = dist[j][i + 1];
			}
		}
	
		for(int i = 0; i < nV - i; i++) {
			int index = nV - i - 1;
			dist[i][index] = newDist[i];
			dist[index][i] = newDist[i];	
		}
	}
	
	return dend[0];
}

/**
 * Frees all memory associated with the given Dendrogram structure.
 */
void freeDendrogram(Dendrogram d) {
	if (d == NULL) return;
	freeDendrogram(d->left);
	freeDendrogram(d->right);
	free(d);
}
// ======================================
// ======== HELPER FUNCTIONS ============
// ======================================

static Dendrogram newDend(int v) {
	Dendrogram new = malloc(sizeof(Dendrogram));
	left(new) = NULL;
	right(new) = NULL;
	vertex(new) = v;
	return new;
}

static double **calDistance(Graph g, double **dist) {
	int nV = GraphNumVertices(g);
	for(int i = 0; i < nV; i++) {
		for(int j = 0; j < nV; j++) {
			// checking for a direct link between nodes of the graph
			if(GraphIsAdjacent(g, i, j) || GraphIsAdjacent(g, i ,j)) {
				// calulating the distances
				int weightIn = 0;
				int weightOut = 0;
				AdjList in = GraphInIncident(g, i);
				AdjList out = GraphOutIncident(g, i);

				// for in links
				while(in != NULL) {
					if (v(in) == j) {
						weightIn = weight(in);
					}
					in = next(in);
				}
				// for out links
				while(out != NULL) {
					if (v(out) == j) {
						weightOut = weight(out);
					}
					out = next(out);
				}

				double value;
				
				// setting distance between vertices
				if(weightOut - weightIn > 0) {
					value = 1 / (double)weightOut;
				} else {
					value = 1/ (double)weightIn;
				}

				dist[i][j] = value;
				dist[j][i] = value;
			} else {
				// if vertices i and j aren't connected
				dist[i][j] = INFINITY;
				dist[j][i] = INFINITY;
			} 
		}	
	}
	return dist;
}

/*static double **removeCluster(Graph g, double **dist, int v1, int v2) {
	int nV = GraphNumVertices(g) - 1;
	

	return dist;
}*/