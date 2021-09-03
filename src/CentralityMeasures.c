// Centrality Measures API implementation
// COMP2521 Assignment 2

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "CentralityMeasures.h"
#include "Dijkstra.h"
#include "PQ.h"

// struct that makes everything look nicer
// :)
typedef struct path {
	Vertex src;
	Vertex middle;
	Vertex dest;
} Path;


static double calCloseness(double distanceSum, int numVisted, int nV);
static NodeValues createNodeValues(int numNodes);
static bool pathCheck(Path p);
static double calBetweeness(ShortestPaths sps, Path p);
static int numPaths(ShortestPaths sps, Path p);
static int numPathFinder(ShortestPaths sps, Path p);


// O(n^3)
// This function calculates the closeness of centrality of every node in a graph
// While using a Wasserman and Faust formula on a directed graph
// Closeness is a measure of how long it would take for information to spread from the 
// src node to all other nodes
NodeValues closenessCentrality(Graph g) {
	if(g == NULL) {
		fprintf(stderr, "Invalid graph\n");
	}
	int nV = GraphNumVertices(g);
	NodeValues nvs = createNodeValues(nV);
	// looping through all nodes of the graph
	for(int i = 0; i < nV; i++) {
		// running the shortest path algo on the g from src vertex: i
		ShortestPaths sps = dijkstra(g, i);
		// finding num of nodes reached and sum of distances
		int numVisted = 1; // as src node is reached
		double sumDistance = 0;
		
		for(int j = 0; j < nV; j++) {
			// checking that the node isn't src or unreachable
			if(sps.dist[j] != 0 && sps.dist[j] != INFINITY) {
				numVisted++;
				sumDistance = sumDistance + sps.dist[j];
			}
			
		}
		// checking if the node has no edges to it
		AdjList in = GraphInIncident(g,i);
		AdjList out = GraphOutIncident(g, i);
		
		// the specific cases don't really need to be checked
		// as the nvs.value is calloced, and closeness will return 0
		// just to be safe
		if (in == NULL && out == NULL) {
			nvs.values[i] = 0;
		} else if (sumDistance == 0) { // if the node is the src
			nvs.values[i] = 0;
		} else { // calculating the closeness
			nvs.values[i] = calCloseness(sumDistance, numVisted, nV);
		}
		freeShortestPaths(sps);
	}
	return nvs;
}

// O(n^4)
// This function calulates the betweenes of centrality
// it is similar to the closeness but uses a method
// of a path through the graph with a source, middle and destination node

// Betweeness is the measurement of how many times the middle is a bridge along a 
// path between different src, and dest nodes
NodeValues betweennessCentrality(Graph g) {
	if(g == NULL) {
		fprintf(stderr, "Invalid graph\n");
	}

	int nV = GraphNumVertices(g);
	NodeValues nvs = createNodeValues(nV);
	
	// Looping through g
	// middle node
	for(int middleN = 0; middleN < nV; middleN++) {
		// source node
		for(int srcN = 0; srcN < nV; srcN++) {
			ShortestPaths sps = dijkstra(g, srcN);
			
			// destination node
			for(int destN = 0; destN < nV; destN++) {
				Path p;
				p.src = srcN;
				p.middle = middleN;
				p.dest = destN;
				
				if(pathCheck(p)) {
					nvs.values[middleN] = nvs.values[middleN] + calBetweeness(sps, p);
				}
				
			}
			freeShortestPaths(sps);
		}
	}
	return nvs;
}

// Same as betweennessCentrality, but normalises the value of the betweeness
// The betweeness normalise is just to rescale the value for nodes not including v
NodeValues betweennessCentralityNormalised(Graph g) {
	double nV = GraphNumVertices(g);
	NodeValues nvs = betweennessCentrality(g);
	double normalise = 1 / ((nV - 1) * (nV - 2));
	for(Vertex v = 0; v < nV; v++) {
		// nV can't be less than 2, as function becomes undefined
		if (nV > 2) {
			double value = normalise * nvs.values[v];
			nvs.values[v] = value;
		} else {
			nvs.values[v] = 0;
		} 
	}

	return nvs;
}

void showNodeValues(NodeValues nvs) {
	printf("nV : %d\n", nvs.numNodes);
	printf("=== Node Values ===\n");
	for(int i = 0; i < nvs.numNodes; i++) {
		printf("%d : %lf\n", i, nvs.values[i]);
	}
}

void freeNodeValues(NodeValues nvs) {	
	free(nvs.values);
}

// ========================
// === HELPER FUNCTIONS ===
// ========================

// helper that creates a NodeValues array
static NodeValues createNodeValues(int numNodes) {
	NodeValues new;
	new.numNodes = numNodes;
	new.values = calloc(numNodes, sizeof(double));
	return new;
}


// calulates closeness using the Wasserman and Faust formula
// this function assumes that the node is connected and not isolated
static double calCloseness(double distanceSum, int numVisted, int nV) {
	double numVist = (numVisted - 1) * (numVisted - 1);
	double nnV = nV - 1;
	double sum = 1 / distanceSum;
	
	double value = (numVist / nnV) * (sum);
	
	return value;
}

// checks to see if the path could be a valid path
// i.e src != dest, src != middle, middle != dest
static bool pathCheck(Path p) {
	bool first = (p.src != p.middle);
	bool second = (p.src != p.dest);
	bool third = (p.middle != p.dest);
	return (first && second && third);
}

// helper function that calulates the betweeness of the graph
static double calBetweeness(ShortestPaths sps, Path p) {
	double numberPaths = numPaths(sps, p);
	// isolated node
	if (numberPaths == 0) {
		return 0;
	}

	double numThrough = numPathFinder(sps, p);

	double calulated = numThrough / numberPaths;
	return calulated;
}

// calulcates number of paths p.src to p.dest 
static int numPaths(ShortestPaths sps, Path p) {
	int numP = 0;
	PredNode *head = sps.pred[p.dest];
	if(p.src == p.dest) {
		return 1;
	} else {
		// calculating
		while(head != NULL) {
			Path s = p;
			s.dest = head->v;
			numP = numP + numPaths(sps, s);
			head = head->next;
		}
	}
	return numP;
}

// counting how many times p.middle appears
static int numPathFinder(ShortestPaths sps, Path p) {
	int count = 0;

	if(p.middle != p.dest) {
		PredNode *head = sps.pred[p.dest];
		while(head != NULL) {
			Path s = p;
			s.dest = head->v;
			count = count + numPathFinder(sps, s);
			head = head->next;
		}
	} else {
		count = numPaths(sps, p);
	}

	return count;
}