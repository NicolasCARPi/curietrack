/*******************************************************************************************************************
									ChainNeighbor Class
								This class allows to make a chained list to chain all the presnet neighbor nodes
								Institut Curie
									UMR - 144
									by Victor Racine
									2004 06 17
*******************************************************************************************************************/

#ifndef __chain_neighbor__
		#define __chain_neighbor__
		#include "general.h"


		class Node;
		class ChainNeighbor{
		public:
			ChainNeighbor(Node*	_neighbor, ChainNeighbor* _next);
			~ChainNeighbor();
			Node*	neighbor;
			ChainNeighbor*	next;
		};

#endif	//__chain_neighbor__
