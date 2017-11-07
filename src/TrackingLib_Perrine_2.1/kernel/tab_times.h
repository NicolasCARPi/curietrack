/*******************************************************************************************************************
									TabTimes Class
							Array containing all the TimeNodes for all instants
									Institut Curie
									UMR - 144
									by Victor Racine
									12 july 2002
									modified by Perrine August 2011
*******************************************************************************************************************/


#ifndef __tab_times__
		#define __tab_times__
		#include <vector>

		#include <stdlib.h>
		#include <time.h>  //for rand()

		#include "general.h"
		#include "time_nodes.h"

		using namespace std;

		class TabTimes{
		protected:
			TimeNodes**	timeNodes;
			int linkNeighborsPresent();
			int linkNeighborsBeforeAndAfter();
			int linkNeighborsBlinkBeforeAndAfter();
			int sortByIntensityAllNeighbors();
			int	readStatParameters(double** objectstotrack, int allobjectsinframes, int nbframesinmovie, int sizeimageX, int sizeimageY);
			//int printConnections();
			//int printConnectionsASCII();
			int printConnectionsASCIIwithoutput( std::vector<std::vector<double> > *outputVector);
			void setIndexToTrackedStructures();
			virtual TimeNodes* allocTimeNodes(long time, long nbNodes){return(NULL);};
			void allocTabTimes();
		public:

			int giveResult(std::vector<std::vector<double> > *output);
			TabTimes();
			~TabTimes();
			TimeNodes *operator()(int i);
			TimeNodes **operator[](int i);
			int nbTimes;

			int initialisation(double** objectstotrack, int allobjectsinframes, int nbframesinmovie, int sizeimageX, int sizeimageY);
			virtual int runModel(){return(T_ERR_OOP);};

		};





#endif	//__tab_times__
