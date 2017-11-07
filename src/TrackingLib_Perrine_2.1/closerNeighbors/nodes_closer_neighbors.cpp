/*******************************************************************************************************************
									Nodes closer neighbors Class
									It is derived from the class Nodes
								This class find the closer neighbors for each structure.
								Institut Curie
									UMR - 144
									by Victor Racine
									2004 06 21
*******************************************************************************************************************/


#include "nodes_closer_neighbors.h"


/********************************************************************************************************************
Constructuor(int _time) of NodeCloserNeighbors class.
********************************************************************************************************************/
NodeCloserNeighbors::NodeCloserNeighbors(int _time):Node(_time){

}





void NodeCloserNeighbors::makeConnection(void){
	ChainNeighbor*c;
	Node*nMin=NULL;
	double distMin=-1.;
	for(c=neighborsFuture;c!=NULL;c=c->next){
		if(c->neighbor->connectionPast==NULL){
			if(distMin<0.){
				distMin=distance(c->neighbor);
				nMin=c->neighbor;
			}else
				if(distance(c->neighbor)<distMin){
					distMin=distance(c->neighbor);
					nMin=c->neighbor;
				}
		}
	}
	if(nMin!=NULL){
		nMin->connectionPast=this;
		connectionFuture=nMin;
	}
}
