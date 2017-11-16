/*******************************************************************************************************************
    ChainNeighbor Class
    This class allows to make a chained list to chain all the present neighbor nodes
    Institut Curie
    UMR - 144
    by Victor Racine
    2004 06 17
*******************************************************************************************************************/

#include "chain_neighbor.h"

extern PublicParameters *publicParameters;

/********************************************************************************************************************
ChainNeighbor(Node*	_neighbor, ChainNeighbor* _next) of ChainNeighbor class.
********************************************************************************************************************/
ChainNeighbor::ChainNeighbor(Node*	_neighbor, ChainNeighbor* _next){
	neighbor=_neighbor;
	next=_next;

	publicParameters->nbChainNeighborAllocated++;
}

/********************************************************************************************************************
Destructuor() of ChainLink class.
********************************************************************************************************************/
ChainNeighbor::~ChainNeighbor(){

	publicParameters->nbChainNeighborDeallocated++;
}
