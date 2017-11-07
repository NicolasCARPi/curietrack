/*******************************************************************************************************************
								General Structures used for the tracking
									Institut Curie
									UMR - 144
									by Victor Racine
									12 july 2002
									modified by Perrine 2011
*******************************************************************************************************************/

#pragma once
#include <fcntl.h>
#include <sys/stat.h>
#include <stdio.h>
#include <math.h>
#include <wchar.h>
#include <stdlib.h>



class PublicParameters{
public:
    PublicParameters(){};
    ~PublicParameters(){};
    double distMaxNeighbour;
    double ratioDistCluster;
    double distMaxPastAndFuture;
    double dydX;
    int nbNodeAllocated;
    int nbNodeDeallocated;
    int nbChainNeighborAllocated;
    int nbChainNeighborDeallocated;
    char statFile[1024];
    bool inputAscii;
    char fileRes[1024];
    char fileResASCII[1024];

    char fileResASCIIdivision[1024];

    bool verbose;
    long model; //the number of the used model
    double likenessIntensityMin;
    double temperatureDecrease;
    long nbDim;
    long nbMinSuccessiveStructures;
    bool allowFusionAndFission;
    long nbTemporalWindows;
    double costBirthDeath;
    double stdIntensity;
    double stdMean;
    double stdDistance;
    double covMeanIntensity;
    double maxX;
    double maxY;

    double weightDistance;
    double weightIntensity;
};
