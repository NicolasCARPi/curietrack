// MathFuncsLib.cpp
// compile with: /c /EHsc
// post-build command: lib MathFuncsLib.obj

#include "TrackingLib.h"

using namespace std;





FILE *outputStream;

PublicParameters *publicParameters;

namespace TrackingFunc
{

int getOpt(int temperaturedecrease, double normalizedIntensityWeight, double maxDistance,double dydX, int maxDisapearanceTime, double costBirthDeath, double ratiocluster){


	publicParameters->likenessIntensityMin=0.2;
	publicParameters->verbose=false;
	publicParameters->distMaxNeighbour=maxDistance;
	publicParameters->ratioDistCluster=ratiocluster;
	publicParameters->distMaxPastAndFuture=maxDistance;
	publicParameters->dydX=dydX;

	publicParameters->temperatureDecrease=temperaturedecrease;
	publicParameters->nbDim=2;
	publicParameters->nbMinSuccessiveStructures=1;
	publicParameters->allowFusionAndFission=true;
	publicParameters->nbTemporalWindows=maxDisapearanceTime;
	publicParameters->costBirthDeath=costBirthDeath;
    publicParameters->model=2;

	publicParameters->weightIntensity =normalizedIntensityWeight;
	publicParameters->weightDistance =1-normalizedIntensityWeight;



    return(RET_OK);
	// todo : add a watch event for blinking (missed segmenattion, fusion: there is a problem, and false division
}

#define RETURN(param) {if(publicParameters->fileRes[0]!=L'\0'){fclose(outputStream);}return(param);}

int CRTracking::TrackingCellRace(int temperaturedecrease, double normalizedIntensityWeight, double maxDistance,double dydX, int maxdisapear, double costBirthDeath ,double** objectstotrack, int allobjectsinframes, int nbframesinmovie,std::vector<std::vector<double> > *outputVectorResults, int sizeimageX, int sizeimageY, double ratioCluster){
	TabTimes*	tabTimes;
	int ret;

	publicParameters=new PublicParameters();


	if((ret=getOpt( temperaturedecrease, normalizedIntensityWeight, maxDistance, dydX, maxdisapear,  costBirthDeath, ratioCluster ))<0) return(ret);
	//printf("error code: %d \n",ret);

	outputStream=stdout;

	//printf("\n result copy est :%25s \n",publicParameters->statFile);




	tabTimes=new TabTimesDistanceIntensity();



	if((ret=tabTimes->initialisation( objectstotrack,allobjectsinframes, nbframesinmovie, sizeimageX,sizeimageY))<0)
		RETURN(ret);
	if((ret=tabTimes->runModel())<0)
		RETURN(ret);
	if((ret=tabTimes->giveResult(outputVectorResults))<0)
		RETURN(ret);

	delete tabTimes;

	if(publicParameters->fileRes[0]!=L'\0'){
		fclose(outputStream);
	}

	delete publicParameters;



	return(RET_OK);
}
}
#undef RETURN
