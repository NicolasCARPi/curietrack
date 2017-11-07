/*******************************************************************************************************************
									Node Class
								  Node representing one structure
									Institut Curie
									UMR - 144
									by Victor Racine
									12 july 2002
*******************************************************************************************************************/

#ifndef __nodes__
#define __nodes__


#include "general.h"
#include "../generalParam/NoError.h"
#include "info_node.h"

#include "chain_neighbor.h"
#include <string.h>


#define DIST_SIGNIFICANT	10.
#define DIST_INFINITY		1.e10



class Node{
protected:
public:
	InfoNode* info;
	//int n;
	Node(int _time);
	~Node();

	float distance(Node *node);
	float distance(Node *n1, Node *n2);
	float distancedxdy(Node *node);
	float distancedxdy(Node *n1, Node *n2);
	//InfoData2D	info;

	int time;
	Node*	connectionFuture;
	Node*	connectionBlinkFuture;
	Node*	connectionPast;
	Node*	connectionBlinkPast;
	Node*	connectionFission1Future;
	Node*	connectionFission2Future;
	Node*	connectionFissionPast;
	Node*	connectionFusion1Past;
	Node*	connectionFusion2Past;
	Node*	connectionFusionFuture;
	unsigned long trackedStructureIndex;
	unsigned long trackedStructureAmiraIndex;
	unsigned short colorIndex[2];

	ChainNeighbor* neighborsPresent;
	long neighborsPresentCount;
	long neighborsPresentCountCluster;
	ChainNeighbor* neighborsPast;
	long neighborsPastCount;
	ChainNeighbor* neighborsFuture;
	long neighborsFutureCount;

	ChainNeighbor* neighborsBlinkPast;
	long neighborsBlinkPastCount;
	ChainNeighbor* neighborsBlinkFuture;
	long neighborsBlinkFutureCount;


	int sortByIntersityNeighbors(ChainNeighbor** neighborsHead);
	virtual double likeness(Node* n){exit(T_ERR_OOP);return(0);}
	virtual double likenessMeanIntensity(Node* n){exit(T_ERR_OOP);return(0);}
	virtual double likenessIntensity(Node* n){exit(T_ERR_OOP);return(0);}
	virtual double likenessDistance(Node* n){exit(T_ERR_OOP);return(0);}
	virtual double likenessFusion(Node* n1, Node* n2){exit(T_ERR_OOP);return(0);}
	virtual double likenessIntensityFusion(Node* n1, Node* n2){exit(T_ERR_OOP);return(0);}
	virtual double likenessMeanIntensityFusion(Node* n1, Node* n2){exit(T_ERR_OOP);return(0);}
	virtual double likenessDistanceFusion(Node* n1, Node* n2){exit(T_ERR_OOP);return(0);}
	virtual double likenessFission(Node* n1, Node* n2){exit(T_ERR_OOP);return(0);}
	virtual double likenessBlink(Node* n){exit(T_ERR_OOP);return(0);}


	double distanceToBorder(void);
	double stericEffect(void);
};



#endif//	__nodes__
