/*******************************************************************************************************************
    TimeNodes Class
    Array containing all the Node s of a particular instant
    Institut Curie
    UMR - 144
    by Victor Racine
    12 july 2002
*******************************************************************************************************************/

/********************************************************************************************************************
Class TimeNodes contains the Node objets for one time.
********************************************************************************************************************/

#include "time_nodes.h"

extern FILE *outputStream;
extern PublicParameters *publicParameters;
/********************************************************************************************************************
Constructuor(int _time, int _nbNodes) of TimeNodes class.
The current objets is at the time "time".
Create nbNodes objets Node.
********************************************************************************************************************/
TimeNodes::TimeNodes(int _time, int _nbNodes){
	nbNodes=_nbNodes;
	time=_time;
	tabNodes= new Node*[nbNodes];
}


/********************************************************************************************************************
Destructor() of TimeNodes class.
Delete all objets Node.
********************************************************************************************************************/
TimeNodes::~TimeNodes(){
	for(int i=0;i<nbNodes;i++)
		delete tabNodes[i];
	delete [] tabNodes;
}

/********************************************************************************************************************
operator()(int i) to access to all Node objets.
********************************************************************************************************************/
Node *TimeNodes::operator()(int i){
	return(tabNodes[i]);
}

/********************************************************************************************************************
linkNeighbors() build the chained list of all the neighbours of each nodes presents at this time.
********************************************************************************************************************/
int TimeNodes::linkNeighborsPresent(){
	long i,j;
	for(i=0;i<nbNodes;i++)
		for(j=0;j<nbNodes;j++)
			if(i!=j){
				if(tabNodes[i]->distancedxdy(tabNodes[j])<=publicParameters->distMaxNeighbour){
					tabNodes[i]->neighborsPresent=new ChainNeighbor(tabNodes[j], tabNodes[i]->neighborsPresent);// check the neigbors of this node
					tabNodes[i]->neighborsPresentCount++;
				}
				if(tabNodes[i]->distancedxdy(tabNodes[j])<=publicParameters->distMaxNeighbour/publicParameters->ratioDistCluster){
					tabNodes[i]->neighborsPresentCountCluster++;
					}
				}

	return(RET_OK);
}


/********************************************************************************************************************
linkBeforeAndAfter() for each structure at the time t, it builds the chained list of all the structures of the next
time (t+1) close to it.
********************************************************************************************************************/
int TimeNodes::linkNeighborsBeforeAndAfter(TimeNodes & futureTimeNodes){
	long presentNode,futureNode;
	for(presentNode=0;presentNode<nbNodes;presentNode++)
		for(futureNode=0;futureNode<futureTimeNodes.nbNodes;futureNode++)
			if(tabNodes[presentNode]->distancedxdy(futureTimeNodes.tabNodes[futureNode])<=publicParameters->distMaxPastAndFuture){
				tabNodes[presentNode]->neighborsFuture=new ChainNeighbor(futureTimeNodes.tabNodes[futureNode],
					tabNodes[presentNode]->neighborsFuture);
				tabNodes[presentNode]->neighborsFutureCount++;
				futureTimeNodes.tabNodes[futureNode]->neighborsPast=new ChainNeighbor(tabNodes[presentNode],
					futureTimeNodes.tabNodes[futureNode]->neighborsPast);
				futureTimeNodes.tabNodes[futureNode]->neighborsPastCount++;
			}

	return(RET_OK);
}

/********************************************************************************************************************
linkBlinkBeforeAndAfter() for each structure at the time t, it builds the chained list of all the structures of the next time (t+dt) close to it.
********************************************************************************************************************/
int TimeNodes::linkNeighborsBlinkBeforeAndAfter(TimeNodes & futureTimeNodes){
	long presentNode,futureNode;

	for(presentNode=0;presentNode<nbNodes;presentNode++)
		for(futureNode=0;futureNode<futureTimeNodes.nbNodes;futureNode++)
			if(tabNodes[presentNode]->distancedxdy(futureTimeNodes.tabNodes[futureNode])<=publicParameters->distMaxPastAndFuture){ // dist max past and future or speed*dt?
				tabNodes[presentNode]->neighborsBlinkFuture=new ChainNeighbor(futureTimeNodes.tabNodes[futureNode],
					tabNodes[presentNode]->neighborsBlinkFuture);
				tabNodes[presentNode]->neighborsBlinkFutureCount++;
				futureTimeNodes.tabNodes[futureNode]->neighborsBlinkPast=new ChainNeighbor(tabNodes[presentNode],
					futureTimeNodes.tabNodes[futureNode]->neighborsBlinkPast);
				futureTimeNodes.tabNodes[futureNode]->neighborsBlinkPastCount++;
			}

	return(RET_OK);
}

/********************************************************************************************************************
sortByIntensityAllNeighbors() Sort by intersity all the neighbors of each nodes.
********************************************************************************************************************/
int TimeNodes::sortByIntensityAllNeighbors(){
	int ret;
	for(int i=0;i<nbNodes;i++){
		if((ret=tabNodes[i]->sortByIntersityNeighbors(&(tabNodes[i]->neighborsPast)))<0) return(ret);
		if((ret=tabNodes[i]->sortByIntersityNeighbors(&(tabNodes[i]->neighborsPresent)))<0) return(ret);
		if((ret=tabNodes[i]->sortByIntersityNeighbors(&(tabNodes[i]->neighborsFuture)))<0) return(ret);
	}
	return(RET_OK);
}
