// Lance-Williams Algorithm for Hierarchical Agglomerative Clustering
// COMP2521 Assignment 2
// z5311885

#include <assert.h>
#include <float.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>


#include "Graph.h"
#include "LanceWilliamsHAC.h"

#define INFINITY DBL_MAX

#define MERGED -1

#define vertex(node) node->vertex
#define left(node) node->left
#define right(node) node->right

// struct that makes everything look nicer
// :)
typedef struct shape {
    double **distance;
    Dendrogram *dendro;
    int method;
    Graph graph;
} Shape;


// === HELPER FUNCTIONS ===
static double max(double x, double y);
static double min(double x, double y);
static Dendrogram newDNode (int v);
static Dendrogram newLinkNode (Shape s, int a, int b);
static double calDistance(double x, double y);
static int cluster(Shape s);
static void merge(Shape s, int v1, int v2);
static void singleLinkage(Shape s, int v1, int v2, int count);
static void completeLinkage(Shape s, int v1, int v2, int count);

/**
 * Generates  a Dendrogram using the Lance-Williams algorithm (discussed
 * in the spec) for the given graph  g  and  the  specified  method  for
 * agglomerative  clustering. The method can be either SINGLE_LINKAGE or
 * COMPLETE_LINKAGE (you only need to implement these two methods).
 * 
 * The function returns a 'Dendrogram' structure.
 */

// O( nV * nV * (nV - 1))
Dendrogram LanceWilliamsHAC(Graph g, int method){
    if(g == NULL) {
		fprintf(stderr, "Invalid graph\n");
	}
    if (method != SINGLE_LINKAGE && method != COMPLETE_LINKAGE) {
        fprintf(stderr, "Invalid method: must be SINGLE or COMPLETE linkage\n");
    } 
    
    
    int nV = GraphNumVertices(g);
    // distance calulated array
    double** dist = malloc(nV * sizeof(double*));
    // dendrogram array
    Dendrogram *dendro = malloc(nV * sizeof(Dendrogram));
    // direct distance array, directly from graph
    int** directDistance = malloc(nV * sizeof(int*));
    // looping to initalise all the memory for the stuff
    // O(nV)
    for (int i = 0; i < nV; i++){
        dist[i] = calloc(nV, sizeof(double));
        dendro[i] = newDNode(i);
        directDistance[i] = calloc(nV, sizeof(int));
    }

    // populating the directDistance array
    // O(nV^2)
    for (int i = 0; i < nV; i++){
        AdjList out = GraphOutIncident(g,i);
        AdjList head;
        for(head = out; head != NULL; head = head->next) {
            directDistance[i][head->v] = head->weight;
        }
    }
   
   
    // O(nV^2)
    for (int i = 0; i < nV; i++){
    // populating dist matrix with the calulated distance
        for (int j = i + 1; j < nV; j++){
            if(directDistance[i][j] == 0 && directDistance[j][i] == 0){
            // if no edges exist
                dist[i][j] = INFINITY;
            } else{
                double a = (double)directDistance[j][i];
                double b = (double)directDistance[i][j];
                dist[i][j] = calDistance(a, b);
            }
        }
    }
    Shape s;
    s.distance = dist;
    s.dendro = dendro;
    s.graph = g;
    s.method = method;


    int index = 0;
    
    int check = cluster(s);
    
    // finding all remaining clusters
    // O(nV^2 * nV - 1)
    while(check != MERGED){

        index = check;
        check = cluster(s);
    }
    Dendrogram final = dendro[index];
    // O(nV)
    for(int i = 0; i < nV; i++) {
        free(dist[i]);
        free(directDistance[i]);
    }
    free(dendro);
    free(dist);
    free(directDistance);

    

    return final;
    // return the last head node 
}


/**
 * Frees all memory associated with the given Dendrogram structure.
 */
// O(1)
void freeDendrogram(Dendrogram d) {
    if(d == NULL) return;
    freeDendrogram(left(d));
    freeDendrogram(right(d));
    free(d);
}

// ========================
// === HELPER FUNCTIONS ===
// ========================

// O(1)
// returns min of two numbers
static double min(double x,double y) {
    if(x == INFINITY && y == INFINITY) {
        return INFINITY;
    } else if(x == INFINITY)  {
        return y;
    } else if (y == INFINITY) {
        return x;
    }
    return (x > y ? y : x);
}

// O(1)
// returns max of two nums
static double max(double x, double y) {
    if(x == INFINITY && y == INFINITY) {
        return INFINITY;
    } else if(x == INFINITY)  {
        return y;
    } else if (y == INFINITY) {
        return x;
    }
    return (x > y ? x : y);
}

// O(1)
// Calculates distance according to spec
static double calDistance(double x, double y) {
    double ret = max(x, y);
    return 1 / ret;
}

// O(1)
// new linking node
static Dendrogram newLinkNode(Shape s, int a, int b) {
    Dendrogram new = malloc(sizeof(DNode));
    
    vertex(new) = -1;
    left(new) = s.dendro[a];
    right(new) = s.dendro[b];
    
    return new;
}

// O(1)
// new normal node
// used in initalisation
static Dendrogram newDNode(int v) {
    Dendrogram new = malloc(sizeof(DNode));
    vertex(new) = v;
    left(new) = NULL;
    right(new) = NULL;
    
    return new;
}


// O(nV^2)
// function that clusters nodes
static int cluster(Shape s){
    int nV = GraphNumVertices(s.graph);
    double minimum  = INFINITY;
    int a = -1;
    int b = -1;

    // finding a shortest path
    for (int i = 0; i < nV; i++){
        for (int j = i + 1; j < nV; j++){
            if(s.distance[i][j] != 0){
                if(s.distance[i][j] <= minimum){
                minimum = s.distance[i][j];
                a = i;
                b = j;
            }
            } 
        }
    }

    // merging is complete
    if(a == MERGED && b == MERGED){
        return MERGED;
    } else{
        merge(s, a, b);
    }
    return a;
}

// O(nV)
// this function merges the nodes and node clusters into bigger clusters
// updating the distance between new nodes
static void merge(Shape s, int v1, int v2) {
    int nV = GraphNumVertices(s.graph);
    
    // new link node
    Dendrogram new = newLinkNode(s, v1, v2);
    //update the array and inserting new cluster node
    s.dendro[v1] = new;
    s.dendro[v2] = NULL;

    //update the matrix
    for(int i = 0; i < nV; i++) {
        // since the array is a ladder so there are several conditions for q
        if(i == v1){
            // renew dist[A][B] to be set as found
            s.distance[i][v2] = 0;
        }
        if (s.method == SINGLE_LINKAGE) {
            singleLinkage(s, v1, v2, i);
        } else if(s.method == COMPLETE_LINKAGE) {
            completeLinkage(s, v1, v2, i);
        }
        
        if(i < v1) {
            s.distance[i][v2] = 0;
        } else if( i > v1 && i < v2) {
            s.distance[i][v2] = 0;
        } else if (i > v2) {
            s.distance[v2][i] = 0;
        }
    }
}

// O(1)
// updating distance for the new cluster
// according to using a single linkage method
static void singleLinkage(Shape s, int v1, int v2, int count) {
        if(count < v1) {
            s.distance[count][v1] = min(s.distance[count][v1],s.distance[count][v2]);
        } else if(count > v1 && count < v2) {
            s.distance[v1][count] = min(s.distance[v1][count],s.distance[count][v2]);
        } else if (count > v2) {
            s.distance[v1][count] = min(s.distance[v1][count],s.distance[v2][count]);
        }
}

// O(1)
// updating distance for the new cluster
// according to using a single linkage method
static void completeLinkage(Shape s, int v1, int v2, int count) {
    if(count < v1) {
        s.distance[count][v1] = max(s.distance[count][v1],s.distance[count][v2]);
    } else if(count > v1 && count < v2) {
        s.distance[v1][count] = max(s.distance[v1][count],s.distance[count][v2]);
    } else if (count > v2) {
        s.distance[v1][count] = max(s.distance[v1][count],s.distance[v2][count]);
    }
}