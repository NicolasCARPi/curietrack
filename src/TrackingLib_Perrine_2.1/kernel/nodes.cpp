/*******************************************************************************************************************
									Node Class
								  Node representing one structure
									Institut Curie
									UMR - 144
									by Victor Racine
									12 july 2002
*******************************************************************************************************************/

/********************************************************************************************************************
Class Node contains the caracteristic off the node, and the chain of the other nodes (time after and time before)
links to the current objet.
********************************************************************************************************************/

#include "nodes.h"
#include "../generalParam/generalParam.h"
#include <algorithm>
float globalDistMax;

extern FILE *outputStream;
extern PublicParameters *publicParameters;


using namespace std;

/*********************************************************************************************
double distanceToBorder(void)
  get the distance from this to the closer border
*********************************************************************************************/
double Node::distanceToBorder(void){
	double mini;
	mini=info->getX(); // distance with y= 0
	mini=min(mini, publicParameters->maxX-info->getX());
	mini=min(mini, info->getY()); //distance with x=0 border
	mini=min(mini, publicParameters->maxY-info->getY());

	return(mini);

}

/*********************************************************************************************
double stericEffect(void)
  Compute the mean of radius of the node this
*********************************************************************************************/
double Node::stericEffect(void){

	double ret;
	double val;

    val=info->getSurface();
    //surface of a sphere is PI * r^2
    ret=sqrt(val/PI);

	return(ret);
}

/********************************************************************************************************************
Constructuor(int _time) of Node class.
The current objets is at the time "time".
********************************************************************************************************************/
Node::Node(int _time){
	time=_time;
	neighborsPresent=NULL;
	neighborsPresentCount=0;
	neighborsPresentCountCluster=0;
	neighborsPast=NULL;
	neighborsPastCount=0;
	neighborsFuture=NULL;
	neighborsFutureCount=0;
	neighborsBlinkPast=NULL;
	neighborsBlinkPastCount=0;
	neighborsBlinkFuture=NULL;
	neighborsBlinkFutureCount=0;

	connectionFuture=NULL;
	connectionPast=NULL;
	connectionBlinkPast=NULL;
	connectionBlinkFuture=NULL;
	connectionFission1Future=NULL;
	connectionFission2Future=NULL;
	connectionFissionPast=NULL;
	connectionFusion1Past=NULL;
	connectionFusion2Past=NULL;
	connectionFusionFuture=NULL;
	if(publicParameters->nbDim==2)
		info=new InfoNode2D();

	publicParameters->nbNodeAllocated++;
}

/********************************************************************************************************************
Destructor() of Node class.
Delete all the linkage objets.
********************************************************************************************************************/
Node::~Node(){

	ChainNeighbor*	l;
	while(neighborsFuture!=NULL){
		l=neighborsFuture->next;
		delete neighborsFuture;
		neighborsFuture=l;
	}

	while(neighborsPast!=NULL){
		l=neighborsPast->next;
		delete neighborsPast;
		neighborsPast=l;
	}

	while(neighborsPresent!=NULL){
		l=neighborsPresent->next;
		delete neighborsPresent;
		neighborsPresent=l;
	}
	delete info;
	publicParameters->nbNodeDeallocated++;
	trackedStructureIndex=0;
}

/********************************************************************************************************************
distance(Node *node) gives the distance between "node" and the current object.
********************************************************************************************************************/
float Node::distance(Node *node){	//Distance 1
	return(info->distance(node->info));
}

/********************************************************************************************************************
distance(Node *n1, Node *n2) gives the distance between the mass center of ( n1 and n2), to the current object.
********************************************************************************************************************/
float Node::distance(Node *n1, Node *n2){
	return(info->distance(n1->info,n2->info));
}

float Node::distancedxdy(Node *node){	//Distance 1
	return(info->distancedxdy(node->info));
}

/********************************************************************************************************************
distance(Node *n1, Node *n2) gives the distance between the mass center of ( n1 and n2), to the current object.
********************************************************************************************************************/
float Node::distancedxdy(Node *n1, Node *n2){
	return(info->distancedxdy(n1->info,n2->info));
}

/********************************************************************************************************************
sortByIntersityNeighbors(ChainNeighbor* neighborsHead) sort the chained list begining at neighborsHead by the intensity.
********************************************************************************************************************/
int Node::sortByIntersityNeighbors(ChainNeighbor** neighborsHead){
	ChainNeighbor* l;
	ChainNeighbor* l2;
	ChainNeighbor* lOld;
	bool change;

	if(*neighborsHead!=NULL)
		if((*neighborsHead)->next!=NULL)
			do{
				change=false;
				l=*neighborsHead;
				lOld=NULL;
				while(l->next!=NULL){

					if(l->next->neighbor->info->getIntensity()>l->neighbor->info->getIntensity()){
						l2=l->next;
						if(lOld==NULL)
							*neighborsHead=l2;
						else
							lOld->next=l2;
						l->next=l2->next;
						l2->next=l;
						change=true;
						lOld=l2;
					}else{
						l=l->next;
						if(lOld==NULL)
							lOld=*neighborsHead;
						else
							lOld=lOld;
					}
				}

			}while(change);

	return(RET_OK);

}

