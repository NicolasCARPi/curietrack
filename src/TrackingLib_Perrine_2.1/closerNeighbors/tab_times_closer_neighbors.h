/*******************************************************************************************************************
    TabTimesCloserNeighbors Class
    It is derived from the class TabTimes
    This class find the closer neighbors for each structure.
    Institut Curie
    UMR - 144
    by Victor Racine
    2004 06 21
*******************************************************************************************************************/
#ifndef __tab_times_closer_neighbors__
		#define __tab_times_closer_neighbors__
		#
		#include "../kernel/tab_times.h"
		#include "time_nodes_closer_neighbors.h"

	class TabTimesCloserNeighbors:public TabTimes{
	protected:
		TimeNodes* allocTimeNodes(long time, long nbNodes);
		void allocNode();
	public:
		TabTimesCloserNeighbors();
		~TabTimesCloserNeighbors();

		int runModel(std::vector<std::vector<double> > *outputVector);
	};


#endif	//__tab_times_closer_neighbors__
