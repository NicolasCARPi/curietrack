// TrackingLib.h
#pragma once
#include <stdio.h>

#include "general.h"
#include "tab_times.h"
#include "../closerNeighbors/tab_times_closer_neighbors.h"
#include "../likenessdistInten/tab_times_dist_inten.h"


namespace TrackingFunc{

    class CRTracking
    {

	public :
		static int TrackingCellRace(int temperaturedecrease, double normalizedIntensityWeight, double maxDistance,double dydX, int disapearancetime, double costBirthDeath ,double** objectstotrack, int allobjectsinframes, int nbframesinmovie, std::vector<std::vector<double> > *output, int sizeimageX, int sizeimageY, double ratioCluster) ;

    };

}
