/*******************************************************************************************************************
									TimeNodesDistanceIntensity Class
									It is derived from the class TimeNodes
								This class find the neighbor with the best likeness.
								Likeness depend on the distance and the intensity.
								Institut Curie
									UMR - 144
									by Victor Racine
									2004 06 24
*******************************************************************************************************************/

#ifndef __time_nodes_dist_inten__
		#define __time_nodes_dist_inten__

		#include "../kernel/time_nodes.h"
		#include "nodes_dist_inten.h"


class TimeNodesDistanceIntensity: public TimeNodes{
protected:
	Node* allocNode(void);
public:
	TimeNodesDistanceIntensity(int _time, int _nbNodes);
	void computeConstantSigmas(double &stdIntensityComputing, long &countStdInt,
		double &stdMeanComputing, long &countStdMean, double &stdDistanceComputing, long &countStdDist, bool sigmaClipping, long*hist=NULL);
	double makeAllInitializationConnections(void);
	double tryRandomChange(double limit);

	double getLikenessTime(long lastTime);

};

#endif //__time_nodes_dist_inten__
