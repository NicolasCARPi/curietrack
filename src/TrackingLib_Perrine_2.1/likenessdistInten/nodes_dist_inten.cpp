/*******************************************************************************************************************
									NodeDistanceIntensity  Class
									It is derived from the class Nodes
								This class find the neighbor with the best likeness.
								Likeness depend on the distance and the intensity.
								Institut Curie
									UMR - 144
									by Victor Racine
									2004 06 24
									modified by Perrine August-September 2011 +Jannuary2012
*******************************************************************************************************************/


#include "nodes_dist_inten.h"
#include "../generalParam/generalParam.h"
#include <algorithm>

extern PublicParameters *publicParameters;
extern FILE *outputStream;
#define _weightDistancediv_ 0.6  //0.6
#define _weightIntensitydiv_ 0.2  //0.2
#define _weightMeanIntensitydiv_ 0.2  //0.2
//#define _weightDistance_ 0.8  //0.6
//#define _weightIntensity_ 0.1  //0.2
#define _weightMeanIntensity_ 0.1  //0.2
#define _minEnergyFusion_ 0
#define _minEnergyBlink_ 0


using namespace std;

/********************************************************************************************************************
bool metropolis(double deltaE, double temperature)
Metropolis algorithm. metropolis returns a boolean variable that issues a verdict on whether
to accept a reconfiguration that leads to a change deltaE in the objective function. If deltaE>0,
metropolis = true, while if deltaE<0, metropolis is only true with probability exp(-deltaE/temperature), where
temperature is a temperature determined by the annealing schedule.
********************************************************************************************************************/
bool NodeDistanceIntensity::metropolis(double deltaE, double temperature)
{

    if(deltaE>=0.)
        return(true);
    double test=(double)(rand()+1)/(double)RAND_MAX;
    if( test> exp(deltaE*100/temperature))
        return(false);
    else
        return(true);



    /*if(deltaE>=-temperature)  //same convergence!!!
    	return(true);
    else
    	return(false);*/

}


/********************************************************************************************************************
Constructuor(int _time) of NodeDistanceIntensity class.
********************************************************************************************************************/
NodeDistanceIntensity::NodeDistanceIntensity(int _time):Node(_time)
{

}

/********************************************************************************************************************
double getLikeness(double difference, double sigma)
compute the likeness for a difference of value compared to the sigma
likeness=1 if difference=0
likeness=0 if difference>=3*sigma
********************************************************************************************************************/
double NodeDistanceIntensity::getLikeness(double difference, double sigma)
{
    if (difference==0) // if difference is null, then very likely
        return (1.);
    if (sigma==0) // if standard deviation is null, then if difference is not null, very unlikely
    {
        if (difference==0)//
            return (1.);
        else
            return(0.);
    }
// likeliness should be between 0 (unlikely) and 1 (very likely)- The more difference there is, the less likelu it is.
//function erf approximé par sqrt(1-exp(-4*X^2/pi)) (approximation Chu55) when x>0
    //return 1-erf
    return( 1-(sqrt(1-exp((-4*((difference)/(3.*sigma))*(difference)/(3.*sigma))/PI))));

    //return(-difference*difference/sigma/sigma);

    //return(max(1.- fabs(difference)/(3.*sigma), 3.- fabs(difference)/sigma));

    //return(1.- sqrt(difference/(3.*sigma)) );


    //return(1.-log(fabs(difference)+1)/log(3.*sigma));
}





/*********************************************************************************************
double NodeDistanceIntensity::likenessDistance(Node* n)
  Compute the likeness from the distance between this and n.
  New version from 2005/01/12
*********************************************************************************************/
double NodeDistanceIntensity::likenessDistance(Node* n)
{
    if(n==NULL)
        return(getLikeness(min(distanceToBorder(), publicParameters->stdDistance), publicParameters->stdDistance)); // attention je rajoute le *3 pour les dooonées d'aqi(autoquant)
    return( getLikeness( distance(n), publicParameters->stdDistance));
}




/*********************************************************************************************
double NodeDistanceIntensity::stericEffect(Node* n1, Node* n2)
  Compute the sum of the minimum radius of the nodes n1 and n2
  New version from 2005/01/12
*********************************************************************************************/
double NodeDistanceIntensity::stericEffect(Node* n1, Node* n2)
{

    return(n1->stericEffect()+n2->stericEffect());
}

/*********************************************************************************************
double NodeDistanceIntensity::likenessDistanceFusion(Node* n1, Node* n2)
  Compute the likeness from the distance of a fusion between this and n1 and n2.
  New version from 2005/01/12
*********************************************************************************************/
double NodeDistanceIntensity::likenessDistanceFusion(Node* n1, Node* n2)
{
    if(n1==NULL)
        return(0.);
    if(n2==NULL)
        return(0.);
    //Likeness =likeness in distance of n and (n&,n2) gravity center+ likeness of n1,n2 beeing near enough (less stddistance)
    return( getLikeness(distance(n1,n2),publicParameters->stdDistance)+ getLikeness( fabs(n1->distance(n2)-stericEffect(n1,n2)), publicParameters->stdDistance/2));
    //return( getLikeness( fabs(distance(n1,n2)+n1->distance(n2)-stericEffect(n1,n2)), publicParameters->stdDistance));
}





/*********************************************************************************************
double NodeDistanceIntensity::likenessMeanIntensity(Node* n)
  Compute the likeness from the mean of intensity between this and n.
  New version from 2005/01/12

  We use the hypothesis that
	sigma(mean intensity) proportional to surface^-1/6 simga_noise / Mean in 2D.
*********************************************************************************************/
double NodeDistanceIntensity::likenessMeanIntensity(Node* n)
{

    if(n==NULL)
        return( getLikeness( info->getMean(), publicParameters->stdMean*
                             max(1., 3.*publicParameters->stdDistance/this->distanceToBorder())) );

    else
        return( getLikeness( info->getMean()-n->info->getMean(), publicParameters->stdMean));

}

/*********************************************************************************************
double NodeDistanceIntensity::likenessMeanIntensityFusion(Node* n1, Node* n2)
  Compute the likeness from the mean of intensity in fusion between this and n1 and n2 .
  New version from 2005/01/12

*********************************************************************************************/
double NodeDistanceIntensity::likenessMeanIntensityFusion(Node* n1, Node* n2)
{

    if(n1==NULL)
        return(0.);	//Error ???
    if(n2==NULL)
        return(0.);	//Error ???

    return( getLikeness( info->getMean()-
                         (n1->info->getIntensity()+n2->info->getIntensity())
                         /(n1->info->getSurface()+n2->info->getSurface()), publicParameters->stdMean));

}



/*********************************************************************************************
double NodeDistanceIntensity::likenessIntensity(Node* n)
  Compute the likeness from the total intensity between this and n.
  New version from 2005/01/12

  We use the hypothesis that sigma(intensity) proportional to simga_noise/Mean.
*********************************************************************************************/
double NodeDistanceIntensity::likenessIntensity(Node* n)
{

    if(n==NULL)
        return( getLikeness( info->getIntensity(), publicParameters->stdIntensity*
                             max(1., 3.*publicParameters->stdDistance/this->distanceToBorder()))); //max(1., 3.*publicParameters->stdDistance/this->distanceToBorder())));
    else
        return( getLikeness( info->getIntensity()-n->info->getIntensity(), publicParameters->stdIntensity));

}




/*********************************************************************************************
double NodeDistanceIntensity::likenessIntensityFusion(Node* n1, Node* n2
  Compute the likeness for the intensity in a fusion between this and n1 and n2.
  New version from 2005/01/12

  We use the hypothesis that sigma(intensity) proportional to simga_noise/Mean.
*********************************************************************************************/
double NodeDistanceIntensity::likenessIntensityFusion(Node* n1, Node* n2)
{
    if(n1==NULL)
        return(0.);	//inutil
    if(n2==NULL)
        return(0.);	//inutil


    return( getLikeness( info->getIntensity()-n1->info->getIntensity()-n2->info->getIntensity(), publicParameters->stdIntensity));
}



/*********************************************************************************************
double NodeDistanceIntensity::likeness(Node* n)
  Compute the likeness between this and n.
  New version from 2005/01/12
*********************************************************************************************/
double NodeDistanceIntensity::likeness(Node* n)
{


    if(n==NULL)
    {
        return(max((publicParameters->weightDistance*likenessDistance(n)+
                    publicParameters->weightIntensity*likenessIntensity(n)+
                    _weightMeanIntensity_*likenessMeanIntensity(n))/2. - publicParameters->costBirthDeath,0.0));// the /2 is because the sum of
        // the likeness of 2 unconnected nodes has to be compared to the likeness of the
        // connection between those two nodes.
    }

    else
    {
        return(
                  publicParameters->weightDistance*likenessDistance(n)+
                  publicParameters->weightIntensity*likenessIntensity(n)+
                  _weightMeanIntensity_*likenessMeanIntensity(n));
    }
}


/*********************************************************************************************
double NodeDistanceIntensity::likeness(Node* n)
  Compute the likeness between this and n.
  New version from 2005/01/12
*********************************************************************************************/
double NodeDistanceIntensity::likenessBlink(Node* n)
{

    double alpha=0.7;
    long dt;
    dt=abs(time-n->time);
    //return(likeness(n)* getLikeness(dt-1, 1)-_minEnergyBlink_);
    return(likeness(n)*alpha + (likeness(NULL)+ n->likeness(NULL))*(1.-alpha/2.));
}




/*********************************************************************************************
double NodeDistanceIntensity::likenessFusion(Node* n1, Node* n2)
  Compute the likeness in fusion between this and n1 and n2.
  New version from 2005/01/12
*********************************************************************************************/

double NodeDistanceIntensity::likenessFusion(Node* n1, Node* n2)
{
    if(n1==NULL)
        return(likeness(n2));
    if(n2==NULL)
        return(likeness(n1));




    return((
               _weightDistancediv_*likenessDistanceFusion(n1,n2)+
               _weightIntensitydiv_*likenessIntensityFusion(n1,n2)+
               _weightMeanIntensitydiv_*likenessMeanIntensityFusion(n1,n2)+
               _minEnergyFusion_

           )*3./2.); // the *3/2 is because the sum of
    // the likeness of 1 unconnected node and 2 linked nodes has to be compared to the likeness of the
    // fusion connection between those 3 nodes.
}



/*********************************************************************************************
double NodeDistanceIntensity::likenessFission(Node* n1, Node* n2)
  Compute the likeness in fission between this and n1 and n2.
  New version from 2005/01/12
*********************************************************************************************/
double NodeDistanceIntensity::likenessFission(Node* n1, Node* n2)
{
    return(likenessFusion(n1, n2));
}



/*********************************************************************************************
void NodeDistanceIntensity::computeConstantSigmas(double &stdIntensityComputing,
	double &stdMeanComputing, double &stdDistanceComputing, long &countStd)
  Find in the vicinity of each node the one with the biggest likeness to compute the constant
  kInt and kMean and kDist
*********************************************************************************************/
void NodeDistanceIntensity::computeConstantSigmas(double &stdIntensityComputing, long &countStdInt,
        double &stdMeanComputing, long &countStdMean, double &stdDistanceComputing, long &countStdDist, bool sigmaClipping, long *hist)
{

    ChainNeighbor*c;
    Node*nMax=NULL;
    double likenessMax=-1e10;
    double likenessVal;

    for(c=neighborsFuture; c!=NULL; c=c->next)
    {
        likenessVal=likeness(c->neighbor);
        if(likenessVal>likenessMax)
        {
            likenessMax=likenessVal;
            nMax=c->neighbor;
        }
    }


    if(nMax!=NULL)
    {

        if(hist!=NULL)
        {
            hist[(long)(distance(nMax))]++;
        }
        else
        {
            if(sigmaClipping)
            {
                if(fabs(info->getIntensity()-nMax->info->getIntensity())<3.*publicParameters->stdIntensity)
                {
                    stdIntensityComputing+=(info->getIntensity()-nMax->info->getIntensity())*(info->getIntensity()-nMax->info->getIntensity());
                    countStdInt++;
                }
                if(fabs(info->getMean()-nMax->info->getMean())<3.*publicParameters->stdMean)
                {
                    stdMeanComputing+=(info->getMean()-nMax->info->getMean())*(info->getMean()-nMax->info->getMean());
                    countStdMean++;
                }

                if(distance(nMax)<3.*publicParameters->stdDistance)
                {
                    stdDistanceComputing+=(distance(nMax))*(distance(nMax));
                    countStdDist++;
                }
            }
            else
            {
                stdIntensityComputing+=(info->getIntensity()-nMax->info->getIntensity())*(info->getIntensity()-nMax->info->getIntensity());
                stdMeanComputing+=(info->getMean()-nMax->info->getMean())*(info->getMean()-nMax->info->getMean());
                stdDistanceComputing+=(distance(nMax))*(distance(nMax));
                countStdInt++;
                countStdMean++;
                countStdDist++;
            }
        }
    }

}

/*********************************************************************************************
double NodeDistanceIntensity::makeInitializationConnection(void)
  Make the initial connections. For each strucutres, chose the one with the highest likeness
  if it is avaliable. If it is not possible no connection.

  It returns the likeness.

  New version from 2005/01/12
*********************************************************************************************/
double NodeDistanceIntensity::makeInitializationConnection(void)
{
    ChainNeighbor*c;
    Node*nMax=NULL;
    double likenessMax=-1e10;
    double likenessVal;

    for(c=neighborsFuture; c!=NULL; c=c->next)
    {
        if(c->neighbor->connectionPast==NULL)
        {
            likenessVal=likeness(c->neighbor);

            if(likenessVal>likenessMax)
            {
                likenessMax=likenessVal;
                nMax=c->neighbor;
            }
        }
    }
    if(nMax!=NULL)
    {
        nMax->connectionPast=this;
        connectionFuture=nMax;
        return(likenessMax);
    }


    return(likeness(NULL));
    //old return(0.);
}


/*********************************************************************************************
double NodeDistanceIntensity::getLikenessOfLocalNode(int lastTime)
  Compute the likeness according to the connections.

  It returns the likeness.

  New version from 2005/01/13
*********************************************************************************************/
double NodeDistanceIntensity::getLikenessOfLocalNode(int lastTime)
{

    double _likeness=0.;

    if(time>0)
    {
        //to the past
        if((connectionPast==NULL) && (connectionFissionPast==NULL) && (connectionFusion1Past==NULL) )
            _likeness=likeness(NULL);
        if(connectionFusion1Past!=NULL)
            _likeness=+likenessFusion(connectionFusion1Past,connectionFusion2Past);
    }

    if(time<lastTime-1)
    {
        //to the future
        if((connectionFuture==NULL) && (connectionFission1Future==NULL) && (connectionFusionFuture==NULL) )
            _likeness+=likeness(NULL);
        if(connectionFuture!=NULL)
            _likeness+=likeness(connectionFuture);
        if(connectionFission1Future!=NULL)
            _likeness+=likenessFission(connectionFission1Future,connectionFission2Future);
    }
    return(_likeness);
}

double NodeDistanceIntensity::makeChange(double limitLoss)
{
    //choose an operation between:
    //simple linkage
    //strucutre fusion
    //strucutre dissociation



    if(publicParameters->allowFusionAndFission)
    {
        if(connectionFuture!=NULL)
        {
            switch(rand() % 4)
            {
            //case 0:
            //return(makeChangeFusion(limitLoss));
            case 1:
                return(makeChangeFission(limitLoss));
            default:
                return(makeChangeSimpleLinkage(limitLoss));
            }
        }

        if((publicParameters->nbTemporalWindows>1)&&(connectionBlinkFuture!=NULL))
        {
            return(makeChangeBlinkLinkage(limitLoss));
        }

        if(connectionFission1Future!=NULL)
            return(makeChangeFission(limitLoss));

        if(connectionFusionFuture!=NULL)
            return(makeChangeFusion(limitLoss));

        if((publicParameters->nbTemporalWindows>1)&&(rand()%2==0))
            return(makeChangeBlinkLinkage(limitLoss));

        return(makeChangeSimpleLinkage(limitLoss));

    }

    if(connectionFuture!=NULL)
        return(makeChangeSimpleLinkage(limitLoss));

    if(connectionBlinkFuture!=NULL)
        return(makeChangeBlinkLinkage(limitLoss));


    if((publicParameters->nbTemporalWindows>2)&&(rand()%2==0))
        return(makeChangeBlinkLinkage(limitLoss));


    return(makeChangeSimpleLinkage(limitLoss));
}

double NodeDistanceIntensity::makeChangeFission0(double limitLoss)
{
    //no fission in old configuration
    //try to find 2 structures to make fission:
    //From
    //1(this)-------2
    //4-------------3(cRand)

    //To
    //1(this)-------2
    //        \-----3(cRand)
    //4-------------NULL


    double deltaLikeness;
    long idRandNeighbor;
    ChainNeighbor*cRand;
    long i;

    if(connectionBlinkFuture!=NULL)
        return(-1); //ERROR


    if(neighborsFutureCount<2)
        return(0.); //no neighbors like 3
    idRandNeighbor=rand() % neighborsFutureCount;
    //set the node cRand
    for(i=0, cRand=neighborsFuture; (i<idRandNeighbor)&&(cRand!=NULL); cRand=cRand->next, i++);
    //test if the node is good
    if(cRand==NULL)
        return(-1.); //ERROR
    if(cRand->neighbor==connectionFuture) // nothing to do, the new random link is the same than connectionFuture.
        return(0.);
    if(cRand->neighbor->connectionFusion1Past!=NULL)// crand is already implicated in a fusion.
        return(0.);
    if(cRand->neighbor->connectionFissionPast!=NULL)// crand is already implicated in a fission.
        return(0.);
    if(cRand->neighbor->connectionBlinkPast!=NULL) // crand is already implicated in a past blink.
        return(0.);

    if(distancedxdy(cRand->neighbor)>publicParameters->distMaxPastAndFuture) //distance is not respected
        return(0.);

    deltaLikeness=	+ likenessFission(cRand->neighbor, connectionFuture)
                    + (cRand->neighbor->connectionPast==NULL ? 0. : cRand->neighbor->connectionPast->likeness(NULL) )
                    - likeness(connectionFuture)
                    - cRand->neighbor->likeness(cRand->neighbor->connectionPast);
    if(!metropolis(deltaLikeness, limitLoss)) //deltaLikeness<limitLoss) // if the new solution is to bad
        return(0.);
    else  // change the configuration
    {


        connectionFission1Future=connectionFuture;
        connectionFission2Future=cRand->neighbor;
        connectionFuture->connectionFissionPast=this;
        cRand->neighbor->connectionFissionPast=this;
        connectionFuture->connectionPast=NULL;
        connectionFuture=NULL;

        if(cRand->neighbor->connectionPast!=NULL)
            cRand->neighbor->connectionPast->connectionFuture=NULL;
        cRand->neighbor->connectionPast=NULL;



        return(deltaLikeness);
    }


}


double NodeDistanceIntensity::makeChangeFission1(double limitLoss, ChainNeighbor*cRand)
{
    //From
    //1(this)---------2(connectionFission1Future)
    //          \-----3(connectionFission2Future)
    //4(cRand)--------5
    //possiblilty 1
    //To
    //1(this)---------2(connectionFission1Future)
    //4(cRand)--------3(connectionFission2Future)
    //NULL------------5

    if(cRand->neighbor->distancedxdy(connectionFission2Future)>publicParameters->distMaxPastAndFuture) //distance is not respected
        return(0.);

    if(cRand->neighbor->connectionFission1Future!=NULL)// cRand is already implicated in a fission.
        return(0);
    if(cRand->neighbor->connectionFusionFuture!=NULL)// cRand is already implicated in a fusion.
        return(0);
    if(cRand->neighbor->connectionBlinkFuture!=NULL) // crand is already implicated in a future blink.
        return(0.);

    double deltaLikeness=  + likeness(connectionFission1Future)
                           + cRand->neighbor->likeness(connectionFission2Future)
                           + ( cRand->neighbor->connectionFuture==NULL ? 0 : cRand->neighbor->connectionFuture->likeness(NULL) )
                           - likenessFusion(connectionFission1Future, connectionFission2Future)
                           - cRand->neighbor->likeness(cRand->neighbor->connectionFuture);
    if(!metropolis(deltaLikeness, limitLoss))//deltaLikeness<limitLoss) // if the new solution is to bad
        return(0.);
    else  // change the configuration
    {





        connectionFuture=connectionFission1Future;
        connectionFuture->connectionPast=this;
        if(cRand->neighbor->connectionFuture!=NULL)
            cRand->neighbor->connectionFuture->connectionPast=NULL;
        cRand->neighbor->connectionFuture=connectionFission2Future;
        connectionFission2Future->connectionPast=cRand->neighbor;
        connectionFission1Future->connectionFissionPast=NULL;
        connectionFission2Future->connectionFissionPast=NULL;
        connectionFission1Future=NULL;
        connectionFission2Future=NULL;


        return(deltaLikeness);
    }
}

double NodeDistanceIntensity::makeChangeFission2(double limitLoss, ChainNeighbor*cRand)
{
    //From
    //1(this)---------2(connectionFission1Future)
    //          \-----3(connectionFission2Future)
    //4(cRand)--------5
    //possiblilty 2
    //To
    //1---------------3(connectionFission2Future)
    //4(cRand)--------2(connectionFission1Future)
    //NULL------------5

    if(cRand->neighbor->distancedxdy(connectionFission1Future)>publicParameters->distMaxPastAndFuture) //distance is not respected
        return(0.);

    if(cRand->neighbor->connectionFission2Future!=NULL)// cRand is already implicated in a fission.
        return(0);
    if(cRand->neighbor->connectionFusionFuture!=NULL)// cRand is already implicated in a fusion.
        return(0);

    if(cRand->neighbor->connectionBlinkFuture!=NULL) // crand is already implicated in a future blink.
        return(0.);

    double deltaLikeness=	+ likeness(connectionFission2Future)
                            + cRand->neighbor->likeness(connectionFission1Future)
                            + ( cRand->neighbor->connectionFuture==NULL ? 0 : cRand->neighbor->connectionFuture->likeness(NULL) )
                            - likenessFusion(connectionFission2Future, connectionFission1Future)
                            - cRand->neighbor->likeness(cRand->neighbor->connectionFuture);
    if(!metropolis(deltaLikeness, limitLoss))//deltaLikeness<limitLoss) // if the new solution is to bad
        return(0.);
    else  // change the configuration
    {




        connectionFuture=connectionFission2Future;
        connectionFuture->connectionPast=this;
        if(cRand->neighbor->connectionFuture!=NULL)
            cRand->neighbor->connectionFuture->connectionPast=NULL;
        cRand->neighbor->connectionFuture=connectionFission1Future;
        connectionFission1Future->connectionPast=cRand->neighbor;
        connectionFission2Future->connectionFissionPast=NULL;
        connectionFission1Future->connectionFissionPast=NULL;
        connectionFission2Future=NULL;
        connectionFission1Future=NULL;

        return(deltaLikeness);
    }


}

double NodeDistanceIntensity::makeChangeFission3(double limitLoss)
{
    //From
    //1(this)---------2(connectionFission1Future)
    //          \-----3(connectionFission2Future)

    //possiblilty 3
    //To
    //1(this)---------2(connectionFission1Future)
    //NULL------------3(connectionFission2Future)

    double deltaLikeness=	+ likeness(connectionFission1Future)
                            + connectionFission2Future->likeness(NULL)
                            - likenessFission(connectionFission1Future, connectionFission2Future);
    if(!metropolis(deltaLikeness, limitLoss))//deltaLikeness<limitLoss) // if the new solution is to bad
        return(0.);
    else  // change the configuration
    {






        connectionFuture=connectionFission1Future;
        connectionFuture->connectionPast=this;
        connectionFuture->connectionFissionPast=NULL;
        connectionFission2Future->connectionFissionPast=NULL;
        connectionFission1Future=NULL;
        connectionFission2Future=NULL;


        return(deltaLikeness);
    }
}

double NodeDistanceIntensity::makeChangeFission4(double limitLoss)
{
    //possiblilty 4
    //To
    //1(this)---------3(connectionFission2Future)
    //NULL------------2(connectionFission1Future)
    double deltaLikeness=	+ likeness(connectionFission2Future)
                            + connectionFission1Future->likeness(NULL)
                            - likenessFission(connectionFission1Future, connectionFission2Future);
    if(!metropolis(deltaLikeness, limitLoss)) //deltaLikeness<limitLoss) // if the new solution is to bad
        return(0.);
    else  // change the configuration
    {




        connectionFuture=connectionFission2Future;
        connectionFuture->connectionPast=this;
        connectionFuture->connectionFissionPast=NULL;
        connectionFission1Future->connectionFissionPast=NULL;
        connectionFission1Future=NULL;
        connectionFission2Future=NULL;



        return(deltaLikeness);
    }
}

double NodeDistanceIntensity::makeChangeFission5(double limitLoss, ChainNeighbor*cRand)
{
    //From
    //1(this)---------2(connectionFission1Future)
    //          \-----3(connectionFission2Future)
    //4(cRand)--------5

    //possibility 5
    //1(this)---------2(connectionFission1Future)
    //           /----3(connectionFission2Future)
    //4(cRand)--/-----5

    if(cRand->neighbor->distancedxdy(connectionFission2Future)>publicParameters->distMaxPastAndFuture) //distance is not respected
        return(0.);
    if(cRand->neighbor->connectionBlinkFuture!=NULL) // crand is already implicated in a future blink.
        return(0.);

    double deltaLikeness=  + likeness(connectionFission1Future)
                           + cRand->neighbor->likenessFusion(connectionFission2Future, cRand->neighbor->connectionFuture)
                           - likenessFusion(connectionFission1Future, connectionFission2Future)
                           - cRand->neighbor->likeness(cRand->neighbor->connectionFuture);
    if(!metropolis(deltaLikeness, limitLoss))//deltaLikeness<limitLoss) // if the new solution is to bad
        return(0.);
    else  // change the configuration
    {



        connectionFuture=connectionFission1Future;
        connectionFuture->connectionPast=this;
        cRand->neighbor->connectionFission1Future=connectionFission2Future;
        cRand->neighbor->connectionFission1Future->connectionFissionPast=cRand->neighbor;
        cRand->neighbor->connectionFission2Future=cRand->neighbor->connectionFuture;
        cRand->neighbor->connectionFission2Future->connectionFissionPast=cRand->neighbor;

        connectionFission1Future->connectionFissionPast=NULL;
        connectionFission1Future=NULL;
        connectionFission2Future=NULL;
        cRand->neighbor->connectionFuture->connectionPast=NULL;
        cRand->neighbor->connectionFuture=NULL;




        return(deltaLikeness);
    }

}


double NodeDistanceIntensity::makeChangeFission8(double limitLoss, ChainNeighbor*cRand)
{
    //From
    //1(this)---------2(connectionFission1Future)
    //          \-----3(connectionFission2Future)
    //4(cRand)--------5

    //possibility 8
    //1(this)---------3(connectionFission2Future)
    //           /----2(connectionFission1Future)
    //4(cRand)--/-----5

    if(cRand->neighbor->distancedxdy(connectionFission1Future)>publicParameters->distMaxPastAndFuture) //distance is not respected
        return(0.);
    if(cRand->neighbor->connectionBlinkFuture!=NULL) // crand is already implicated in a future blink.
        return(0.);

    double deltaLikeness=  + likeness(connectionFission2Future)
                           + cRand->neighbor->likenessFusion(connectionFission1Future, cRand->neighbor->connectionFuture)
                           - likenessFusion(connectionFission1Future, connectionFission2Future)
                           - cRand->neighbor->likeness(cRand->neighbor->connectionFuture);
    if(!metropolis(deltaLikeness, limitLoss))//deltaLikeness<limitLoss) // if the new solution is to bad
        return(0.);
    else  // change the configuration
    {



        connectionFuture=connectionFission2Future;
        connectionFuture->connectionPast=this;
        cRand->neighbor->connectionFission2Future=connectionFission1Future;
        cRand->neighbor->connectionFission2Future->connectionFissionPast=cRand->neighbor;
        cRand->neighbor->connectionFission1Future=cRand->neighbor->connectionFuture;
        cRand->neighbor->connectionFission1Future->connectionFissionPast=cRand->neighbor;

        connectionFission2Future->connectionFissionPast=NULL;
        connectionFission2Future=NULL;
        connectionFission1Future=NULL;
        cRand->neighbor->connectionFuture->connectionPast=NULL;
        cRand->neighbor->connectionFuture=NULL;



        return(deltaLikeness);
    }

}

double NodeDistanceIntensity::makeChangeFission9(double limitLoss, ChainNeighbor*cRand)
{
    //From
    //1(this)---------2(connectionFission1Future)
    //          \-----3(connectionFission2Future)
    //4(cRand)--------5

    //possibility 8
    //1(this)---------5
    //           /----2(connectionFission1Future)
    //4(cRand)--/-----3(connectionFission2Future)
    if(cRand->neighbor->connectionBlinkFuture!=NULL) // crand is already implicated in a future blink.
        return(0.);

    if(cRand->neighbor->connectionFuture!=NULL)
        if(distancedxdy(cRand->neighbor->connectionFuture)>publicParameters->distMaxPastAndFuture) //distance is not respected
            return(0.);

    if(cRand->neighbor->distancedxdy(connectionFission1Future)>publicParameters->distMaxPastAndFuture) //distance is not respected
        return(0.);
    if(cRand->neighbor->distancedxdy(connectionFission2Future)>publicParameters->distMaxPastAndFuture) //distance is not respected
        return(0.);

    if(cRand->neighbor->connectionFission1Future!=NULL) return(0.);
    if(cRand->neighbor->connectionFusionFuture!=NULL) return(0.);

    double deltaLikeness=	+ likeness(cRand->neighbor->connectionFuture)
                            + cRand->neighbor->likenessFusion(connectionFission1Future, connectionFission2Future)
                            - likenessFusion(connectionFission1Future, connectionFission2Future)
                            - cRand->neighbor->likeness(cRand->neighbor->connectionFuture);
    if(!metropolis(deltaLikeness, limitLoss))//deltaLikeness<limitLoss) // if the new solution is to bad
        return(0.);
    else  // change the configuration
    {



        cRand->neighbor->connectionFission1Future=connectionFission1Future;
        connectionFission1Future->connectionFissionPast=cRand->neighbor;
        cRand->neighbor->connectionFission2Future=connectionFission2Future;
        connectionFission2Future->connectionFissionPast=cRand->neighbor;
        connectionFuture=cRand->neighbor->connectionFuture;
        if(connectionFuture!=NULL) connectionFuture->connectionPast=this;
        connectionFission1Future=NULL;
        connectionFission2Future=NULL;
        cRand->neighbor->connectionFuture=NULL;


        return(deltaLikeness);
    }

}

double NodeDistanceIntensity::makeChangeFission6(double limitLoss, ChainNeighbor*cRand)
{
    //From
    //1(this)---------2(connectionFission1Future)
    //          \-----3(connectionFission2Future)
    //4(cRand)--------5(cRand->neighbor->connectionFuture)

    //To
    //possibility 6
    //1(this)--\------2(connectionFission1Future)
    //          \-----5
    //4(cRand)--------3(connectionFission2Future)

    if(cRand->neighbor->connectionBlinkFuture!=NULL) // crand is already implicated in a future blink.
        return(0.);

    if(cRand->neighbor->connectionFuture!=NULL)
        if(distancedxdy(cRand->neighbor->connectionFuture)>publicParameters->distMaxPastAndFuture) //distance is not respected
            return(0.);
    if(cRand->neighbor->distancedxdy(connectionFission2Future)>publicParameters->distMaxPastAndFuture) //distance is not respected
        return(0.);

    double deltaLikeness=  + likenessFusion(connectionFission1Future, cRand->neighbor->connectionFuture)
                           + cRand->neighbor->likeness(connectionFission2Future)
                           - likenessFusion(connectionFission1Future, connectionFission2Future)
                           - cRand->neighbor->likeness(cRand->neighbor->connectionFuture);
    if(!metropolis(deltaLikeness, limitLoss))//deltaLikeness<limitLoss) // if the new solution is to bad
        return(0.);
    else  // change the configuration
    {



        connectionFission2Future->connectionPast=cRand->neighbor;
        cRand->neighbor->connectionFuture->connectionFissionPast=this;
        Node*localN=cRand->neighbor->connectionFuture;
        cRand->neighbor->connectionFuture=connectionFission2Future;
        connectionFission2Future->connectionFissionPast=NULL;
        connectionFission2Future=localN;
        connectionFission2Future->connectionPast=NULL;


        return(deltaLikeness);
    }
}

double NodeDistanceIntensity::makeChangeFission7(double limitLoss, ChainNeighbor*cRand)
{
    //From
    //1(this)---------2(connectionFission1Future)
    //          \-----3(connectionFission2Future)
    //4(cRand)--------5(cRand->neighbor->connectionFuture)

    //To
    //possibility 7
    //1(this)-\-------3(connectionFission2Future)
    //          \-----5
    //4(cRand)--------2(connectionFission1Future)

    if(cRand->neighbor->connectionBlinkFuture!=NULL) // crand is already implicated in a future blink.
        return(0.);

    if(cRand->neighbor->connectionFuture!=NULL)
        if(distancedxdy(cRand->neighbor->connectionFuture)>publicParameters->distMaxPastAndFuture) //distance is not respected
            return(0.);
    if(cRand->neighbor->distancedxdy(connectionFission1Future)>publicParameters->distMaxPastAndFuture) //distance is not respected
        return(0.);


    double deltaLikeness=  + likenessFusion(connectionFission2Future, cRand->neighbor->connectionFuture)
                           + cRand->neighbor->likeness(connectionFission1Future)
                           - likenessFusion(connectionFission2Future, connectionFission1Future)
                           - cRand->neighbor->likeness(cRand->neighbor->connectionFuture);
    if(!metropolis(deltaLikeness, limitLoss))//deltaLikeness<limitLoss) // if the new solution is to bad
        return(0.);
    else  // change the configuration
    {


        connectionFission1Future->connectionPast=cRand->neighbor;
        cRand->neighbor->connectionFuture->connectionFissionPast=this;
        Node*localN=cRand->neighbor->connectionFuture;
        cRand->neighbor->connectionFuture=connectionFission1Future;
        connectionFission1Future->connectionFissionPast=NULL;
        connectionFission1Future=localN;
        connectionFission1Future->connectionPast=NULL;


        return(deltaLikeness);
    }
}


double NodeDistanceIntensity::makeChangeFission(double limitLoss)
{


    //is the node already makes a Fission with an other structure
    if((connectionFuture==NULL)&&(connectionFission1Future==NULL)) //nothing connected in the future!
        return(0.);

    if(connectionBlinkFuture!=NULL) // crand is already implicated in a future blink.
        return(0.);

    if(connectionFusionFuture!=NULL)
        return(0.); //already implicated in a Fusion

    if(connectionFission1Future==NULL)
    {
        //no fission in old configuration
        //try to find 2 structures to make fission:
        //From
        //1(this)-------2
        //4-------------3(cRand)

        //To
        //1(this)-------2
        //        \-----3(cRand)
        //4-------------NULL

        return(makeChangeFission0(limitLoss));

    }
    else
    {
        //there is a fission in old configuration
        //Try to break or change the fission:
        //From
        //1(this)---------2(connectionFission1Future)
        //          \-----3(connectionFission2Future)
        //4(cRand)--------5

        //To
        //possiblilty 1
        //1---------------2(connectionFission1Future)
        //4(cRand)--------3(connectionFission2Future)
        //NULL------------5
        //or
        //possiblilty 2
        //1---------------3(connectionFission2Future)
        //4(cRand)--------2(connectionFission1Future)
        //NULL------------5
        //or
        //possibility 3
        //1(this)---------2(connectionFission1Future)
        //NULL------------3(connectionFission2Future)
        //or
        //possibility 4
        //1(this)---------3(connectionFission2Future)
        //NULL------------2(connectionFission1Future)

        //possibility 5
        //1(this)---------2(connectionFission1Future)
        //           /----3(connectionFission2Future)
        //4(cRand)--/-----5

        //possibility 6
        //1(this)--\------2(connectionFission1Future)
        //          \-----5
        //4(cRand)--------3(connectionFission2Future)

        //possibility 7
        //1(this)-\-------3(connectionFission2Future)
        //          \-----5
        //4(cRand)--------2(connectionFission1Future)

        //possibility 8
        //1(this)---------3(connectionFission2Future)
        //           /----2(connectionFission1Future)
        //4(cRand)--/-----5

        //possibility 9
        //1(this)---------5
        //           /----2(connectionFission1Future)
        //4(cRand)--/-----3(connectionFission2Future)


        switch(rand()%8)
        {
        case 0:
            return(makeChangeFission3(limitLoss));
        case 1:
            return(makeChangeFission4(limitLoss));
        }

        long idRandNeighbor;
        ChainNeighbor*cRand;
        long i;


        if(connectionFission2Future->neighborsPastCount<2) //no neighbors like "4"
            switch(rand()%2)
            {
            case 0:
                return(makeChangeFission3(limitLoss));
            case 1:
                return(makeChangeFission4(limitLoss));
            }

        //try to find a neighbor "4" in the neighborhood of "3"
        idRandNeighbor=rand() % connectionFission2Future->neighborsPastCount;
        //set the node cRand
        for(i=0, cRand=connectionFission2Future->neighborsPast; (i<idRandNeighbor)&&(cRand!=NULL); cRand=cRand->next, i++)
        {
            //test if the node is good

            if(cRand->neighbor==this) // nothing to do, the new random link is the same than this.
                return(0);
            if(cRand->neighbor->connectionBlinkFuture!=NULL) // crand is already implicated in a future blink.
                return(0.);

            if(cRand->neighbor->connectionFuture!=NULL)
                switch(rand()%7)
                {
                case 0:
                    return(makeChangeFission5(limitLoss, cRand));
                case 1:
                    return(makeChangeFission6(limitLoss, cRand));
                case 2:
                    return(makeChangeFission7(limitLoss, cRand));
                case 3:
                    return(makeChangeFission8(limitLoss, cRand));
                }

            switch(rand()%3)
            {
            case 0:
                return(makeChangeFission1(limitLoss, cRand));
            case 1:
                return(makeChangeFission2(limitLoss, cRand));
            default:
                return(makeChangeFission9(limitLoss, cRand));

            }
        }
    }
}

double NodeDistanceIntensity::makeChangeFusion0(double limitLoss)
{
    double deltaLikeness;
    long idRandNeighbor;
    ChainNeighbor*cRand;
    long i;

    //no fusion in old configuration
    //try to find 2 structures to make fusion:
    //From
    //1(this)-------2(connectionFuture)
    //3(cRand)------4(can be NULL)

    //To
    //1-------------2
    //3(cRand)--/
    //NULL----------4

    if(connectionFuture->neighborsPastCount<2) return(0.); //no neighbors like 3
    idRandNeighbor=rand() % connectionFuture->neighborsPastCount;
    //set the node cRand
    for(i=0, cRand=connectionFuture->neighborsPast; (i<idRandNeighbor)&&(cRand!=NULL); cRand=cRand->next, i++);
    //test if the node is good
    if(cRand==NULL)
        return(-1.); //ERROR
    if(cRand->neighbor==this) // nothing to do, the new random link is the same than this.
        return(0);
    if(cRand->neighbor->connectionFission1Future!=NULL)// cRand is already implicated in a fission.
        return(0);
    if(cRand->neighbor->connectionFusionFuture!=NULL)// cRand is already implicated in a fusion.
        return(0);
    if(cRand->neighbor->connectionBlinkFuture!=NULL)// cRand is already implicated in a blink connection
        return(0);

    if(connectionFuture->distancedxdy(cRand->neighbor)>publicParameters->distMaxPastAndFuture) //distance is not respected
        return(0.);

    deltaLikeness=	+ connectionFuture->likenessFusion(cRand->neighbor, this)
                    + ((cRand->neighbor->connectionFuture!=NULL) ? cRand->neighbor->connectionFuture->likeness(NULL) : 0.)
                    - likeness(connectionFuture)
                    - cRand->neighbor->likeness(cRand->neighbor->connectionFuture);

    //if(cRand->neighbor->connectionFuture!=NULL)
    //	deltaLikeness+=	cRand->neighbor->connectionFuture->likeness(NULL);

    if(!metropolis(deltaLikeness, limitLoss)) //deltaLikeness<limitLoss) // if the new solution is to bad
        return(0.);
    else  // change the configuration
    {
        //unlink connectionFuture



        connectionFusionFuture=connectionFuture;
        connectionFuture->connectionFusion1Past=this;
        cRand->neighbor->connectionFusionFuture=connectionFuture;
        connectionFuture->connectionFusion2Past=cRand->neighbor;

        if(cRand->neighbor->connectionFuture!=NULL)
            cRand->neighbor->connectionFuture->connectionPast=NULL;

        cRand->neighbor->connectionFuture=NULL;
        connectionFuture=NULL;
        connectionFusionFuture->connectionPast=NULL;


        return(deltaLikeness);
    }
}





double NodeDistanceIntensity::makeChangeFusion1(double limitLoss, ChainNeighbor*cRand)
{
    Node*node;

    //find "node"
    if(connectionFusionFuture->connectionFusion1Past==this)
        node=connectionFusionFuture->connectionFusion2Past;
    else
        node=connectionFusionFuture->connectionFusion1Past;
    //there is a fusion in old configuration
    //Try to break the fusion:


    //From
    //1(this)---------2(connectionFusionFuture)
    //3(node)----/
    //5---------------4(cRand)

    //To
    //possiblilty 1
    //1---------------2
    //3---------------4(cRand)
    //5---------------NULL


    if(cRand->neighbor->connectionFissionPast!=NULL)// cRand is already implicated in a fission.
        return(0);
    if(cRand->neighbor->connectionFusion1Past!=NULL)// cRand is already implicated in a fusion.
        return(0);

    if(cRand->neighbor->distancedxdy(node)>publicParameters->distMaxPastAndFuture) //distance is not respected
        return(0.);


    double deltaLikeness=	+ likeness(connectionFusionFuture)
                            + node->likeness(cRand->neighbor)
                            + (cRand->neighbor->connectionPast==NULL ? 0 : cRand->neighbor->connectionPast->likeness(NULL) )
                            - connectionFusionFuture->likenessFusion(node, this)
                            - cRand->neighbor->likeness(cRand->neighbor->connectionPast);
    if(!metropolis(deltaLikeness, limitLoss))//deltaLikeness<limitLoss) // if the new solution is to bad
        return(0.);
    else  // change the configuration
    {

        connectionFusionFuture->connectionPast=this;
        connectionFusionFuture->connectionFusion1Past=NULL;
        connectionFusionFuture->connectionFusion2Past=NULL;
        node->connectionFuture=cRand->neighbor;
        node->connectionFusionFuture=NULL;
        connectionFuture=connectionFusionFuture;
        connectionFusionFuture=NULL;

        if(cRand->neighbor->connectionPast!=NULL)
            cRand->neighbor->connectionPast->connectionFuture=NULL;
        cRand->neighbor->connectionPast=node;


        return(deltaLikeness);
    }
}

double NodeDistanceIntensity::makeChangeFusion2(double limitLoss, ChainNeighbor*cRand)
{

    Node*node;

    //find "node"
    if(connectionFusionFuture->connectionFusion1Past==this)
        node=connectionFusionFuture->connectionFusion2Past;
    else
        node=connectionFusionFuture->connectionFusion1Past;

    //From
    //1(this)---------2(connectionFusionFuture)
    //3(node)----/
    //5---------------4(cRand)

    //possiblilty 2
    //To
    //1---------------4(cRand)
    //3(node)---------2(connectionFusionFuture)
    //5---------------NULL
    //try to find a neighbor "4" in the neighborhood of "1"

    if(cRand->neighbor->connectionFissionPast!=NULL)// cRand is already implicated in a fission.
        return(0);
    if(cRand->neighbor->connectionFusion1Past!=NULL)// cRand is already implicated in a fusion.
        return(0);

    if(distancedxdy(cRand->neighbor)>publicParameters->distMaxPastAndFuture) //distance is not respected
        return(0.);



    double deltaLikeness=	+ likeness(cRand->neighbor)
                            + node->likeness(connectionFusionFuture)
                            + (cRand->neighbor->connectionPast==NULL ? 0 : cRand->neighbor->connectionPast->likeness(NULL))
                            - connectionFusionFuture->likenessFusion(node, this)
                            - cRand->neighbor->likeness(cRand->neighbor->connectionPast);
    if(!metropolis(deltaLikeness, limitLoss))//deltaLikeness<limitLoss) // if the new solution is to bad
        return(0.);
    else  // change the configuration
    {
        //unlink connectionFuture


        connectionFusionFuture->connectionPast=node;
        connectionFusionFuture->connectionFusion1Past=NULL;
        connectionFusionFuture->connectionFusion2Past=NULL;
        node->connectionFusionFuture=NULL;
        node->connectionFuture=connectionFusionFuture;
        connectionFusionFuture=NULL;
        connectionFuture=cRand->neighbor;
        if(cRand->neighbor->connectionPast!=NULL)
            cRand->neighbor->connectionPast->connectionFuture=NULL;
        cRand->neighbor->connectionPast=this;

        return(deltaLikeness);
    }
}

double NodeDistanceIntensity::makeChangeFusion3(double limitLoss)
{
    double deltaLikeness;
    Node*node;

    //find "node"
    if(connectionFusionFuture->connectionFusion1Past==this)
        node=connectionFusionFuture->connectionFusion2Past;
    else
        node=connectionFusionFuture->connectionFusion1Past;
    //From
    //1(this)---------2(connectionFusionFuture)
    //3(node)----/

    //possiblilty 3
    //To
    //1(this)---------2(connectionFusionFuture)
    //3(node)---------NULL

    deltaLikeness=	+ likeness(connectionFusionFuture)
                    + node->likeness(NULL)
                    - connectionFusionFuture->likenessFusion(node, this);
    if(!metropolis(deltaLikeness, limitLoss))//deltaLikeness<limitLoss) // if the new solution is to bad
        return(0.);
    else  // change the configuration
    {


        connectionFusionFuture->connectionPast=this;
        connectionFusionFuture->connectionFusion1Past=NULL;
        connectionFusionFuture->connectionFusion2Past=NULL;
        connectionFuture=connectionFusionFuture;
        connectionFusionFuture=NULL;
        node->connectionFusionFuture=NULL;
        if(node->connectionFuture!=NULL)
            return(-1.);//impossible inutile


        return(deltaLikeness);
    }
}

double NodeDistanceIntensity::makeChangeFusion4(double limitLoss)
{
    double deltaLikeness;
    Node*node;

    //find "node"
    if(connectionFusionFuture->connectionFusion1Past==this)
        node=connectionFusionFuture->connectionFusion2Past;
    else
        node=connectionFusionFuture->connectionFusion1Past;
    //From
    //1(this)---------2(connectionFusionFuture)
    //3(node)----/

    //possiblilty 4
    //1(this)---------NULL
    //3(node)---------2(connectionFusionFuture)
    deltaLikeness=	+ node->likeness(connectionFusionFuture)
                    + likeness(NULL)
                    - connectionFusionFuture->likenessFusion(node, this);
    if(!metropolis(deltaLikeness, limitLoss))//deltaLikeness<limitLoss) // if the new solution is to bad
        return(0.);
    else  // change the configuration
    {


        connectionFusionFuture->connectionPast=node;
        connectionFusionFuture->connectionFusion1Past=NULL;
        connectionFusionFuture->connectionFusion2Past=NULL;
        node->connectionFuture=connectionFusionFuture;
        connectionFusionFuture=NULL;
        node->connectionFusionFuture=NULL;

        if(connectionFuture!=NULL)
            return(-1.);//impossible inutile
        return(deltaLikeness);
    }
}



double NodeDistanceIntensity::makeChangeFusion5(double limitLoss, ChainNeighbor*cRand)
{
    Node*node;

    //find "node"
    if(connectionFusionFuture->connectionFusion1Past==this)
        node=connectionFusionFuture->connectionFusion2Past;
    else
        node=connectionFusionFuture->connectionFusion1Past;

    //From
    //1(this)-----/---2(connectionFusionFuture)
    //3(node)----/
    //5---------------4(cRand)


    //possibility 5
    //1(this)---------2(connectionFusionFuture)
    //3(node)----\
    //5-----------\---4(cRand)

    if(node->distancedxdy(cRand->neighbor)>publicParameters->distMaxPastAndFuture) //distance is not respected
        return(0.);

    if(cRand->neighbor->connectionFissionPast!=NULL)// cRand is already implicated in a fission.
        return(0);
    if(cRand->neighbor->connectionFusion1Past!=NULL)// cRand is already implicated in a fusion.
        return(0);

    if(cRand->neighbor->connectionPast==NULL)// not possible to reconnection the fusion.
        return(0);

    double deltaLikeness=	+ likeness(connectionFusionFuture)
                            + cRand->neighbor->likenessFusion(node, cRand->neighbor->connectionPast)
                            - connectionFusionFuture->likenessFusion(node, this)
                            - cRand->neighbor->likeness(cRand->neighbor->connectionPast);
    if(!metropolis(deltaLikeness, limitLoss))//deltaLikeness<limitLoss) // if the new solution is to bad
        return(0.);
    else  // change the configuration
    {


        connectionFuture=connectionFusionFuture;
        connectionFuture->connectionPast=this;

        node->connectionFusionFuture=cRand->neighbor;
        cRand->neighbor->connectionFusion1Past=node;

        cRand->neighbor->connectionPast->connectionFusionFuture=cRand->neighbor;
        cRand->neighbor->connectionFusion2Past=cRand->neighbor->connectionPast;

        connectionFuture->connectionFusion1Past=NULL;
        connectionFusionFuture=NULL;

        connectionFuture->connectionFusion2Past=NULL;

        cRand->neighbor->connectionPast->connectionFuture=NULL;
        cRand->neighbor->connectionPast=NULL;

        return(deltaLikeness);
    }
}


double NodeDistanceIntensity::makeChangeFusion7(double limitLoss, ChainNeighbor*cRand)
{
    Node*node;

    //find "node"
    if(connectionFusionFuture->connectionFusion1Past==this)
        node=connectionFusionFuture->connectionFusion2Past;
    else
        node=connectionFusionFuture->connectionFusion1Past;

    //From
    //1(this)-----/---2(connectionFusionFuture)
    //3(node)----/
    //5---------------4(cRand)


    //possibility 7
    //3(node)---------2(connectionFusionFuture)
    //1(this)----\
    //5-----------\---4(cRand)

    if(cRand->neighbor->connectionFissionPast!=NULL)// cRand is already implicated in a fission.
        return(0);
    if(cRand->neighbor->connectionFusion1Past!=NULL)// cRand is already implicated in a fusion.
        return(0);

    if(cRand->neighbor->connectionPast==NULL)// not possible to reconnection the fusion.
        return(0);

    if(distancedxdy(cRand->neighbor)>publicParameters->distMaxPastAndFuture) //distance is not respected
        return(0.);
    double deltaLikeness=	+ node->likeness(connectionFusionFuture)
                            + cRand->neighbor->likenessFusion(this, cRand->neighbor->connectionPast)
                            - connectionFusionFuture->likenessFusion(node, this)
                            - cRand->neighbor->likeness(cRand->neighbor->connectionPast);
    if(!metropolis(deltaLikeness, limitLoss))//deltaLikeness<limitLoss) // if the new solution is to bad
        return(0.);
    else  // change the configuration
    {



        node->connectionFuture=node->connectionFusionFuture;
        node->connectionFuture->connectionPast=node;

        connectionFusionFuture=cRand->neighbor;
        cRand->neighbor->connectionFusion1Past=this;

        cRand->neighbor->connectionPast->connectionFusionFuture=cRand->neighbor;
        cRand->neighbor->connectionFusion2Past=cRand->neighbor->connectionPast;

        node->connectionFuture->connectionFusion1Past=NULL;
        node->connectionFusionFuture=NULL;

        node->connectionFuture->connectionFusion2Past=NULL;

        cRand->neighbor->connectionPast->connectionFuture=NULL;
        cRand->neighbor->connectionPast=NULL;


        return(deltaLikeness);
    }
}


double NodeDistanceIntensity::makeChangeFusion8(double limitLoss, ChainNeighbor*cRand)
{
    Node*node;

    //find "node"
    if(connectionFusionFuture->connectionFusion1Past==this)
        node=connectionFusionFuture->connectionFusion2Past;
    else
        node=connectionFusionFuture->connectionFusion1Past;
    //From
    //1(this)---------2(connectionFusionFuture)
    //3(node)----/
    //5---------------4(cRand)

    //possibility 8
    //5---------------2(connectionFusionFuture)
    //1(this)----\
    //3(node)-----\---4(cRand)

    if(cRand->neighbor->connectionFissionPast!=NULL)// cRand is already implicated in a fission.
        return(0);
    if(cRand->neighbor->connectionFusion1Past!=NULL)// cRand is already implicated in a fusion.
        return(0);

    if(cRand->neighbor->connectionPast==NULL)// not possible to reconnection the fusion.
        return(0);

    if(distancedxdy(cRand->neighbor)>publicParameters->distMaxPastAndFuture) //distance is not respected
        return(0.);
    if(node->distancedxdy(cRand->neighbor)>publicParameters->distMaxPastAndFuture) //distance is not respected
        return(0.);
    if(cRand->neighbor->connectionPast->distancedxdy(connectionFusionFuture)>publicParameters->distMaxPastAndFuture) //distance is not respected
        return(0.);

    double deltaLikeness=	+ cRand->neighbor->connectionPast->likeness(connectionFusionFuture)
                            + cRand->neighbor->likenessFusion(this, node)
                            - connectionFusionFuture->likenessFusion(node, this)
                            - cRand->neighbor->likeness(cRand->neighbor->connectionPast);
    if(!metropolis(deltaLikeness, limitLoss))//deltaLikeness<limitLoss) // if the new solution is to bad
        return(0.);
    else  // change the configuration
    {


        connectionFusionFuture->connectionPast=cRand->neighbor->connectionPast;
        connectionFusionFuture->connectionPast->connectionFuture=connectionFusionFuture;

        cRand->neighbor->connectionFusion1Past=connectionFusionFuture->connectionFusion1Past;
        cRand->neighbor->connectionFusion2Past=connectionFusionFuture->connectionFusion2Past;
        cRand->neighbor->connectionPast=NULL;
        //connectionFusionFuture->connectionPast->connectionFuture=NULL;
        node->connectionFusionFuture=cRand->neighbor;
        connectionFusionFuture->connectionFusion1Past=NULL;
        connectionFusionFuture->connectionFusion2Past=NULL;
        connectionFusionFuture=cRand->neighbor;


        return(deltaLikeness);
    }
}

double NodeDistanceIntensity::makeChangeFusion6(double limitLoss, ChainNeighbor*cRand)
{
    Node*node;

    //find "node"
    if(connectionFusionFuture->connectionFusion1Past==this)
        node=connectionFusionFuture->connectionFusion2Past;
    else
        node=connectionFusionFuture->connectionFusion1Past;

    //From
    //1(this)-----/---2(connectionFusionFuture)
    //3(node)----/
    //5---------------4(cRand)


    //possibility 6
    //1(this)-----/---2(connectionFusionFuture)
    //5----------/
    //3(node)---------4(cRand)

    if(node->distancedxdy(cRand->neighbor)>publicParameters->distMaxPastAndFuture) //distance is not respected
        return(0.);
    if(cRand->neighbor->connectionPast!=NULL)
        if(cRand->neighbor->connectionPast->distancedxdy(connectionFusionFuture)>publicParameters->distMaxPastAndFuture) //distance is not respected
            return(0.);


    double deltaLikeness=	+ connectionFusionFuture->likenessFusion(this, cRand->neighbor->connectionPast)
                            + cRand->neighbor->likeness(node)
                            - connectionFusionFuture->likenessFusion(node, this)
                            - cRand->neighbor->likeness(cRand->neighbor->connectionPast);
    if(!metropolis(deltaLikeness, limitLoss))//deltaLikeness<limitLoss) // if the new solution is to bad
        return(0.);
    else  // change the configuration
    {


        if(connectionFusionFuture->connectionFusion1Past==this)
        {
            cRand->neighbor->connectionPast->connectionFusionFuture=connectionFusionFuture;
            cRand->neighbor->connectionPast->connectionFuture=NULL;
            connectionFusionFuture->connectionFusion2Past=cRand->neighbor->connectionPast;
            node->connectionFusionFuture=NULL;
            node->connectionFuture=cRand->neighbor;
            cRand->neighbor->connectionPast=node;
        }
        else
        {
            cRand->neighbor->connectionPast->connectionFusionFuture=connectionFusionFuture;
            cRand->neighbor->connectionPast->connectionFuture=NULL;
            connectionFusionFuture->connectionFusion1Past=cRand->neighbor->connectionPast;
            node->connectionFusionFuture=NULL;
            node->connectionFuture=cRand->neighbor;
            cRand->neighbor->connectionPast=node;

        }


        return(deltaLikeness);
    }
}

double NodeDistanceIntensity::makeChangeFusion9(double limitLoss, ChainNeighbor*cRand)
{
    Node*node;

    //find "node"
    if(connectionFusionFuture->connectionFusion1Past==this)
        node=connectionFusionFuture->connectionFusion2Past;
    else
        node=connectionFusionFuture->connectionFusion1Past;
    Node*node2;
    //find "node2"
    if(cRand->neighbor->connectionFissionPast->connectionFission1Future==cRand->neighbor)
        node2=cRand->neighbor->connectionFissionPast->connectionFission2Future;
    else
        node2=cRand->neighbor->connectionFissionPast->connectionFission1Future;

    //possibility 9
    //From
    //1(this)-----/---2(connectionFusionFuture)
    //3(node)----/
    //5---------\-----4(cRand)
    //           \----6(node2)
    //To
    //1(this)---------2(connectionFusionFuture)
    //3(node)---------4(cRand)
    //5---------------6(node2)

    if(node->distancedxdy(cRand->neighbor)>publicParameters->distMaxPastAndFuture) //distance is not respected
        return(0.);


    double deltaLikeness=	+ likeness(connectionFusionFuture)
                            + node->likeness(cRand->neighbor)
                            + cRand->neighbor->connectionFissionPast->likeness(node2)
                            - connectionFusionFuture->likenessFusion(node, this)
                            - cRand->neighbor->connectionFissionPast->likenessFission(cRand->neighbor, node2);
    if(!metropolis(deltaLikeness, limitLoss))//deltaLikeness<limitLoss) // if the new solution is to bad
        return(0.);
    else  // change the configuration
    {


        connectionFusionFuture->connectionPast=this;
        connectionFuture=connectionFusionFuture;
        connectionFusionFuture->connectionFusion1Past=NULL;
        connectionFusionFuture->connectionFusion2Past=NULL;
        node->connectionFusionFuture=NULL;
        connectionFusionFuture=NULL;
        node->connectionFuture=cRand->neighbor;
        cRand->neighbor->connectionPast=node;
        cRand->neighbor->connectionFissionPast=NULL;
        node2->connectionFissionPast->connectionFission1Future=NULL;
        node2->connectionFissionPast->connectionFission2Future=NULL;
        node2->connectionFissionPast->connectionFuture=node2;
        node2->connectionPast=node2->connectionFissionPast;
        node2->connectionFissionPast=NULL;




        return(deltaLikeness);
    }
}

double NodeDistanceIntensity::makeChangeFusion(double limitLoss)
{


    //is the node already making a Fusion with an other structure
    if((connectionFuture==NULL)&&(connectionFusionFuture==NULL)) //nothing connected in the future!
        return(0.);

    if(connectionFission1Future!=NULL)
        return(0.); //already implicated in a Fission

    if(connectionFusionFuture==NULL)
    {
        //no fusion in old configuration
        //try to find 2 structures to make fusion:
        //From
        //1(this)-------2
        //3-------------4(cRand)

        //To
        //1-------------2
        //3-------/
        //NULL----------4(cRand)

        return(makeChangeFusion0(limitLoss));
    }
    else
    {
        //there is a fusion in old configuration
        //Try to break the fusion:
        //From
        //1(this)---------2(connectionFusionFuture)
        //3(node)----/
        //5---------------4(cRand)

        //To
        //possiblilty 1
        //1---------------2
        //3---------------4(cRand)
        //5---------------NULL
        //or
        //possiblilty 2
        //1---------------4(cRand)
        //3(node)---------2(connectionFusionFuture)
        //5---------------NULL
        //or
        //possibility 3
        //1(this)---------2(connectionFusionFuture)
        //3(node)---------NULL
        //or
        //possibility 4
        //1(this)---------NULL
        //3(node)---------2(connectionFusionFuture)

        //possibility 5
        //1(this)---------2(connectionFusionFuture)
        //3(node)----\
        //5-----------\---4(cRand)

        //possibility 6
        //1(this)----/----2(connectionFusionFuture)
        //5---------/
        //3(node)---------4(cRand)

        //possibility 7
        //3(node)---------2(connectionFusionFuture)
        //1(this)----\
        //5-----------\---4(cRand)

        //possibility 8
        //5---------------2(connectionFusionFuture)
        //1(this)----\
        //3(node)-----\---4(cRand)

        //possibility 9
        //From
        //1(this)---------2(connectionFusionFuture)
        //3(node)----/
        //5---------------4(cRand)
        //           \----6
        //To
        //1(this)---------2(connectionFusionFuture)
        //3(node)---------4(cRand)
        //5---------------6


        //just unlink the fusion

        switch(rand()%8)
        {
        case 0:
            return(makeChangeFusion3(limitLoss));
        case 1:
            return(makeChangeFusion4(limitLoss));
        }


        long idRandNeighbor;
        ChainNeighbor*cRand;
        long i;

        //get the cRand
        if(neighborsFutureCount<2)  //no neighbors like "4", no cRand
            switch(rand()%2)
            {
            case 0:
                return(makeChangeFusion3(limitLoss));
            case 1:
                return(makeChangeFusion4(limitLoss));
            }

        idRandNeighbor=rand() % neighborsFutureCount;
        //set the node cRand
        for(i=0, cRand=neighborsFuture; (i<idRandNeighbor)&&(cRand!=NULL); cRand=cRand->next, i++);
        //test if the node is good
        if(cRand==NULL)
            return(-1.); //ERROR

        if(cRand->neighbor->connectionFissionPast!=NULL)
        {
            return(makeChangeFusion9(limitLoss, cRand));
        }
        if(cRand->neighbor->connectionBlinkPast!=NULL)// cRand is already implicated in a blink connection
            return(0);

        if(cRand->neighbor->connectionPast!=NULL)
            switch(rand()%7)
            {
            case 0:
                return(makeChangeFusion5(limitLoss, cRand));
            case 1:
            case 2:
                return(makeChangeFusion6(limitLoss, cRand));
            case 3:
                return(makeChangeFusion7(limitLoss, cRand));
            case 4:
                return(makeChangeFusion8(limitLoss, cRand));
            }

        switch(rand()%2)
        {
        case 0:
            return(makeChangeFusion1(limitLoss, cRand));
        default:
            return(makeChangeFusion2(limitLoss, cRand));
        }


    }

}

double NodeDistanceIntensity::makeChangeSimpleLinkage0(double limitLoss)  //simple deconnection
{

    //From
    //Present	Future
    //1(this)-----2(connectionFuture)

    //To
    //1(this)-----NULL
    //NULL--------2(connectionFuture)

    double deltaLikeness=	+ likeness(NULL)
                            + connectionFuture->likeness(NULL)
                            - likeness(connectionFuture);
    if(!metropolis(deltaLikeness, limitLoss))//deltaLikeness<limitLoss) // if the new solution is to bad
        return(0.);
    else  // change the configuration
    {



        connectionFuture->connectionPast=NULL;
        connectionFuture=NULL;

        return(deltaLikeness);
    }

}


double NodeDistanceIntensity::makeChangeSimpleLinkage1(double limitLoss, ChainNeighbor*cRand)
{
    //	From
    //	Present		Future
    //	1(this)-----2
    //	NULL--------3(cRand)

    //	To
    //	Present		Future
    //	1(this)-----3(cRand)
    //	NULL--------2

    if(cRand->neighbor->connectionBlinkPast!=NULL)// cRand is already implicated in a blink connection
        return(0);

    if(distancedxdy(cRand->neighbor)>publicParameters->distMaxPastAndFuture) //distance is not respected
        return(0.);



    double deltaLikeness=	+ likeness(cRand->neighbor)
                            + (connectionFuture!=NULL ? connectionFuture->likeness(NULL) : 0)
                            - cRand->neighbor->likeness(NULL)
                            - likeness(connectionFuture);
    if(!metropolis(deltaLikeness, limitLoss))//deltaLikeness<limitLoss) // if the new solution is to bad
        return(0.);
    else  // change the configuration
    {


        //unlink connectionFuture
        if(connectionFuture!=NULL) connectionFuture->connectionPast=NULL;
        connectionFuture=cRand->neighbor;
        connectionFuture->connectionPast=this;

        return(deltaLikeness);
    }


}

double NodeDistanceIntensity::makeChangeSimpleLinkage2(double limitLoss, ChainNeighbor*cRand)
{

    // try to break the link
    //	From
    //	Present		Future
    //	1(this)-----2
    //	4-----------3(cRand)
    //	To
    //	Present		Future
    //	1-----------3(cRand)
    //	NULL--------2
    //	4-----------NULL

    if(distancedxdy(cRand->neighbor)>publicParameters->distMaxPastAndFuture) //distance is not respected
        return(0.);

    double deltaLikeness=	+ likeness(cRand->neighbor)
                            + (connectionFuture==NULL ? 0:connectionFuture->likeness(NULL) )
                            + (cRand->neighbor->connectionPast==NULL ? 0 : cRand->neighbor->connectionPast->likeness(NULL) )
                            - likeness(connectionFuture)
                            - cRand->neighbor->likeness(cRand->neighbor->connectionPast);
    if(!metropolis(deltaLikeness, limitLoss))//deltaLikeness<limitLoss) // if the new solution is to bad
        return(0.);
    else  // change the configuration
    {
        //unlink connectionFuture


        if(connectionFuture!=NULL) connectionFuture->connectionPast=NULL;
        connectionFuture=cRand->neighbor;
        connectionFuture->connectionPast->connectionFuture=NULL;
        connectionFuture->connectionPast=this;

        return(deltaLikeness);
    }
}

double NodeDistanceIntensity::makeChangeSimpleLinkage3(double limitLoss, ChainNeighbor*cRand, long count)
{
    //	From
    //	Present		Future
    //	1-----------2
    //	4-----------3(cRand)
    //	NULL--------5(c)

    long idRandNeighbor;
    ChainNeighbor*c;
    //find a random one
    idRandNeighbor=rand() % count;
    int iCount=0;

    c=cRand->neighbor->connectionPast->neighborsFuture;
    do
    {
        if((c->neighbor->connectionFusion1Past==NULL)&&(c->neighbor->connectionFissionPast==NULL))
        {
            if(c->neighbor->connectionPast==NULL)
            {
                iCount++;
            }
            else
            {
                if(c->neighbor==connectionFuture)
                {
                    iCount++;
                }
            }
        }
        if(iCount<=idRandNeighbor)
            c=c->next;

    }
    while(iCount<=idRandNeighbor);


    //	To
    //	Present		Future
    //	1-----------3(cRand)
    //	4-----------5(c)
    //	NULL--------2
    if(c==NULL)
        return(-1.); //ERROR
    if(c->neighbor->connectionFusion1Past!=NULL)
        return(-1); //error
    if(c->neighbor->connectionFissionPast!=NULL)
        return(-1); //error
    if(c->neighbor!=connectionFuture)
        if(c->neighbor->connectionPast!=NULL)
            return(-1); //error

    if(distancedxdy(cRand->neighbor)>publicParameters->distMaxPastAndFuture) //distance is not respected
        return(0.);
    if(cRand->neighbor->connectionBlinkPast!=NULL)// cRand is already implicated in a blink connection
        return(0);
    if(c->neighbor->connectionBlinkPast!=NULL)// c is already implicated in a blink connection
        return(0);
    if(cRand->neighbor->connectionPast!=NULL)
        if(cRand->neighbor->connectionPast->distancedxdy(c->neighbor)>publicParameters->distMaxPastAndFuture) //distance is not respected
            return(0.);

    double deltaLikeness=	+ likeness(cRand->neighbor)
                            + cRand->neighbor->connectionPast->likeness(c->neighbor)
                            + (connectionFuture==NULL ? 0 : connectionFuture->likeness(NULL))
                            - likeness(connectionFuture)
                            - cRand->neighbor->likeness(cRand->neighbor->connectionPast)
                            - c->neighbor->likeness(NULL);
    if(!metropolis(deltaLikeness, limitLoss)) //deltaLikeness<limitLoss) // if the new solution is to bad
        return(0.);
    else  // change the configuration
    {
        //unlink connectionFuture


        if(connectionFuture!=NULL) connectionFuture->connectionPast=NULL; //????
        connectionFuture=cRand->neighbor;
        connectionFuture->connectionPast->connectionFuture=c->neighbor;
        connectionFuture->connectionPast->connectionFuture->connectionPast=connectionFuture->connectionPast;
        connectionFuture->connectionPast=this;

        return(deltaLikeness);
    }

}

double NodeDistanceIntensity::makeChangeSimpleLinkage4(double limitLoss, ChainNeighbor*cRand)
{
    //	From
    //	Present		Future
    //	1(this)-----2
    //	4-----------3(cRand)
    //	To
    //	Present		Future
    //	1-----------3(cRand)
    //	4-----------2


    if(cRand->neighbor->connectionPast!=NULL)
        if(cRand->neighbor->connectionPast->distancedxdy(connectionFuture)>publicParameters->distMaxPastAndFuture) //distance is not respected
            return(0.);
    if(distancedxdy(cRand->neighbor)>publicParameters->distMaxPastAndFuture) //distance is not respected
        return(0.);

    double deltaLikeness=	+ likeness(cRand->neighbor)
                            + connectionFuture->likeness(cRand->neighbor->connectionPast)
                            - cRand->neighbor->likeness(cRand->neighbor->connectionPast)
                            - likeness(connectionFuture);
    if(!metropolis(deltaLikeness, limitLoss))//deltaLikeness<limitLoss) // if the new solution is to bad
        return(0.);
    else  // change the configuration
    {


        cRand->neighbor->connectionPast->connectionFuture=connectionFuture;
        connectionFuture->connectionPast=cRand->neighbor->connectionPast;
        connectionFuture=cRand->neighbor;
        cRand->neighbor->connectionPast=this;

        return(deltaLikeness);
    }
}


double NodeDistanceIntensity::makeChangeBlinkLinkage0(double limitLoss)
{
    //transform :
    //	From
    //	t			t+1		t+2
    //	1(this)-------------NULL
    //	NULL----------------2


    //	To
    //	t			t+1		t+2
    //	1(this)-------------2(connectionBlinkFuture)
    long dt,i;
    ChainNeighbor*cRand;

    //if(connectionFuture!=NULL)
    //	return(0);

    dt=rand()%(publicParameters->nbTemporalWindows-1)+2;
    long nbNodeAtGoodTime=0;
    cRand=NULL; //inutile
    bool exitFor=false;
    long idRandNeighbor;

    for(i=0, cRand=neighborsBlinkFuture; cRand!=NULL; cRand=cRand->next)
        //	if(cRand->neighbor->time-time>dt)
        //		;//exitFor=true;
        //	else
        if((cRand->neighbor->time-time==dt)
                &&(cRand->neighbor->connectionPast==NULL)
                &&(cRand->neighbor->connectionBlinkPast==NULL)
                &&(cRand->neighbor->connectionFissionPast==NULL)
                &&(cRand->neighbor->connectionFusion1Past==NULL))
            nbNodeAtGoodTime++;
    if(nbNodeAtGoodTime==0)
        return(0);

    idRandNeighbor=rand() % nbNodeAtGoodTime;
    i=0;
    exitFor=false;
    for(cRand=neighborsBlinkFuture; (!exitFor)&&(cRand!=NULL);)
    {
        if((cRand->neighbor->time-time==dt)
                &&(cRand->neighbor->connectionPast==NULL)
                &&(cRand->neighbor->connectionBlinkPast==NULL)
                &&(cRand->neighbor->connectionFissionPast==NULL)
                &&(cRand->neighbor->connectionFusion1Past==NULL))
        {

            if(i==idRandNeighbor)
                exitFor=true;
            i++;
        }
        if(!exitFor)
            cRand=cRand->next;
    }


    if(cRand==NULL)
        return(-1.); //ERROR

    double deltaLikeness=+ likenessBlink(cRand->neighbor) - cRand->neighbor->likeness(NULL) - likeness(NULL);
    if(!metropolis(deltaLikeness, limitLoss))//deltaLikeness<limitLoss) // if the new solution is to bad
        return(0.);
    else  // change the configuration
    {


        cRand->neighbor->connectionBlinkPast=this;
        connectionBlinkFuture=cRand->neighbor;

        return(deltaLikeness);
    }
}

double NodeDistanceIntensity::makeChangeBlinkLinkage1(double limitLoss)
{
    //transform :
    //	From
    //	t			t+1		t+2
    //	1(this)-------------2(connectionBlinkFuture)

    //	To
    //	t			t+1		t+2
    //	1(this)-------------NULL
    //	NULL----------------2
    double deltaLikeness= + connectionBlinkFuture->likeness(NULL) + likeness(NULL) - likenessBlink(connectionBlinkFuture);
    if(!metropolis(deltaLikeness, limitLoss)) // deltaLikeness<limitLoss) // if the new solution is to bad
        return(0.);
    else  // change the configuration
    {



        connectionBlinkFuture->connectionBlinkPast=NULL;
        connectionBlinkFuture=NULL;

        return(deltaLikeness);
    }
}

double NodeDistanceIntensity::makeChangeBlinkLinkage2(double limitLoss, ChainNeighbor*cRand)
{
    //transform :
    //	From
    //	t			t+1					t+n
    //	1(this)-------------------------2(connectionBlinkFuture)
    //	4---------3(cRand)----...
    //  1=this
    //  2=connectionBlinkFuture
    //	3=cRand->neighbor
    //	4=cRand->neighbor->connectionPast


    //	To
    //	t			t+1					t+n
    //	4-------------------------------2(connectionBlinkFuture)
    //	1(this)---3(cRand)----...


    double deltaLikeness= + likeness(cRand->neighbor)+ cRand->neighbor->connectionPast->likenessBlink(connectionBlinkFuture)
                          - likenessBlink(connectionBlinkFuture) - cRand->neighbor->connectionPast->likeness(cRand->neighbor);
    if(!metropolis(deltaLikeness, limitLoss)) // deltaLikeness<limitLoss) // if the new solution is to bad
        return(0.);
    else  // change the configuration
    {


        connectionBlinkFuture->connectionBlinkPast=cRand->neighbor->connectionPast;
        cRand->neighbor->connectionPast->connectionBlinkFuture=connectionBlinkFuture;
        cRand->neighbor->connectionPast->connectionFuture=NULL;
        cRand->neighbor->connectionPast=this;
        connectionFuture=cRand->neighbor;
        connectionBlinkFuture=NULL;


        return(deltaLikeness);
    }
}

double NodeDistanceIntensity::makeChangeBlinkLinkage3(double limitLoss, ChainNeighbor*cRand)
{
    //transform :
    //	From
    //	t					 t+n-1		t+n
    //	1(this)-------------------------2(connectionBlinkFuture)
    //				...------3(cRand)---4
    //  1=this
    //  2=connectionBlinkFuture
    //	3=cRand->neighbor
    //	4=cRand->neighbor->connectionFuture


    //	To
    //	t					 t+n-1		t+n
    //	1(this)-------------------------4
    //				...------3(cRand)---2(connectionBlinkFuture)
    //  makeChangeBlinkLinkage3

    double deltaLikeness= + likenessBlink(cRand->neighbor->connectionFuture)+ cRand->neighbor->likeness(connectionBlinkFuture)
                          - likenessBlink(connectionBlinkFuture)-cRand->neighbor->connectionFuture->likeness(cRand->neighbor);
    if(!metropolis(deltaLikeness, limitLoss)) // deltaLikeness<limitLoss) // if the new solution is to bad
        return(0.);
    else  // change the configuration
    {



        cRand->neighbor->connectionFuture->connectionBlinkPast=this;
        cRand->neighbor->connectionFuture->connectionPast=NULL;
        connectionBlinkFuture->connectionBlinkPast=NULL;
        connectionBlinkFuture->connectionPast=cRand->neighbor;
        Node * tmp=cRand->neighbor->connectionFuture;
        cRand->neighbor->connectionFuture=connectionBlinkFuture;
        connectionBlinkFuture=tmp;



        return(deltaLikeness);
    }
}

double NodeDistanceIntensity::makeChangeBlinkLinkage4(double limitLoss, ChainNeighbor*cRand)
{

    //transform :
    //	From
    //	t			t+1					t+n
    //	1(this)-------------------------2(connectionBlinkFuture)
    //	NULL--------3(cRand2)----...
    //  1=this
    //  2=connectionBlinkFuture
    //	3=cRand2->neighbor


    //	To
    //	t			t+1					t+n
    //	NULL----------------------------2(connectionBlinkFuture)
    //	1(this)---3(cRand2)----...
    //  makeChangeBlinkLinkage4

    double deltaLikeness= + likeness(cRand->neighbor)+ connectionBlinkFuture->likeness(NULL)
                          - likenessBlink(connectionBlinkFuture) - cRand->neighbor->likeness(NULL);
    if(!metropolis(deltaLikeness, limitLoss)) // deltaLikeness<limitLoss) // if the new solution is to bad
        return(0.);
    else  // change the configuration
    {

        connectionBlinkFuture->connectionBlinkPast=NULL;
        connectionBlinkFuture=NULL;
        cRand->neighbor->connectionPast=this;
        connectionFuture=cRand->neighbor;


        return(deltaLikeness);
    }
}

double NodeDistanceIntensity::makeChangeBlinkLinkage5(double limitLoss, ChainNeighbor*cRand)
{

    //transform :
    //	From
    //	t					 t+n-1		t+n
    //	1(this)-------------------------2(connectionBlinkFuture)
    //				...------3(cRand2)---NULL
    //  1=this
    //  2=connectionBlinkFuture
    //	3=cRand2->neighbor


    //	To
    //	t					 t+n-1		t+n
    //	1(this)-------------------------NULL
    //				...------3(cRand2)---2(connectionBlinkFuture)
    //  makeChangeBlinkLinkage5


    double deltaLikeness= + likeness(NULL)+ cRand->neighbor->likeness(connectionBlinkFuture)
                          - likenessBlink(connectionBlinkFuture)-cRand->neighbor->likeness(NULL);
    if(!metropolis(deltaLikeness, limitLoss)) // deltaLikeness<limitLoss) // if the new solution is to bad
        return(0.);
    else  // change the configuration
    {

        connectionBlinkFuture->connectionBlinkPast=NULL;
        connectionBlinkFuture->connectionPast=cRand->neighbor;
        cRand->neighbor->connectionFuture=connectionBlinkFuture;
        connectionBlinkFuture=NULL;


        return(deltaLikeness);
    }
}

double NodeDistanceIntensity::makeChangeBlinkLinkage(double limitLoss)
{

    int i;
    if(connectionBlinkFuture==NULL)   //nothing connected in the future
    {
        //transform :
        //	From
        //	t			t+1		t+n
        //	1(this)-------------NULL
        //	NULL----------------2


        //	To
        //	t			t+1		t+n
        //	1(this)-------------2(connectionBlinkFuture)
        return(makeChangeBlinkLinkage0(limitLoss));
    }
    else   //connectionBlinkFuture!=NULL
    {
        //transform :
        //	From
        //	t			t+1		t+n
        //	1(this)-------------2(connectionBlinkFuture)

        //	To
        //	t			t+1		t+n
        //	1(this)-------------NULL
        //	NULL----------------2
        //  makeChangeBlinkLinkage1

        //transform :
        //	From
        //	t			t+1					t+n
        //	1(this)-------------------------2(connectionBlinkFuture)
        //	4---------3(cRand)----...
        //  1=this
        //  2=connectionBlinkFuture
        //	3=cRand->neighbor
        //	4=cRand->neighbor->connectionPast


        //	To
        //	t			t+1					t+n
        //	4-------------------------------2(connectionBlinkFuture)
        //	1(this)---3(cRand)----...
        //  makeChangeBlinkLinkage2

        //OR

        //transform :
        //	From
        //	t					 t+n-1		t+n
        //	1(this)-------------------------2(connectionBlinkFuture)
        //				...------3(cRand)---4
        //  1=this
        //  2=connectionBlinkFuture
        //	3=cRand->neighbor
        //	4=cRand->neighbor->connectionFuture


        //	To
        //	t					 t+n-1		t+n
        //	1(this)-------------------------4
        //				...------3(cRand)---2(connectionBlinkFuture)
        //  makeChangeBlinkLinkage3



        //OR

        //transform :
        //	From
        //	t			t+1					t+n
        //	1(this)-------------------------2(connectionBlinkFuture)
        //	NULL--------3(cRand2)----...
        //  1=this
        //  2=connectionBlinkFuture
        //	3=cRand2->neighbor


        //	To
        //	t			t+1					t+n
        //	NULL----------------------------2(connectionBlinkFuture)
        //	1(this)---3(cRand2)----...
        //  makeChangeBlinkLinkage4

        //OR

        //transform :
        //	From
        //	t					 t+n-1		t+n
        //	1(this)-------------------------2(connectionBlinkFuture)
        //				...------3(cRand2)---NULL
        //  1=this
        //  2=connectionBlinkFuture
        //	3=cRand2->neighbor


        //	To
        //	t					 t+n-1		t+n
        //	1(this)-------------------------NULL
        //				...------3(cRand2)---2(connectionBlinkFuture)
        //  makeChangeBlinkLinkage5

        if(rand()%2)
            return(makeChangeBlinkLinkage1(limitLoss));

        //find a potential cRand
        //cRand has to be either just after this or just before connectionBlinkFuture

        if(rand()%2)  //cRand just after this
        {
            ChainNeighbor*cRand;
            ChainNeighbor*cRand2;
            long nbNodeAtGoodTime;
            long nbNodeAtGoodTime2;
            nbNodeAtGoodTime=0;
            nbNodeAtGoodTime2=0;
            cRand=NULL;
            cRand2=NULL;
            bool exitFor=false;

            long idRandNeighbor;

            for(i=0, cRand=neighborsFuture; cRand!=NULL; cRand=cRand->next)
                if(cRand->neighbor->time-time==1)
                {
                    if(cRand->neighbor->connectionPast!=NULL)
                    {
                        nbNodeAtGoodTime++;
                    }
                    else
                    {
                        if((cRand->neighbor->connectionFusion1Past==NULL)&&
                                (cRand->neighbor->connectionFissionPast==NULL)&&
                                (cRand->neighbor->connectionBlinkPast==NULL))
                            nbNodeAtGoodTime2++;

                    }
                }
            if((nbNodeAtGoodTime==0)&&(nbNodeAtGoodTime2==0))
                return(makeChangeBlinkLinkage1(limitLoss));

            if(((rand()%2)&&(nbNodeAtGoodTime!=0))||(nbNodeAtGoodTime2==0))
            {
                idRandNeighbor=rand() % nbNodeAtGoodTime;
                i=0;
                exitFor=false;
                for(cRand=neighborsFuture; (!exitFor)&&(cRand!=NULL);)
                {
                    if((cRand->neighbor->time-time==1)&&(cRand->neighbor->connectionPast!=NULL))
                    {
                        if(i==idRandNeighbor)
                            exitFor=true;
                        i++;
                    }
                    if(!exitFor)
                        cRand=cRand->next;
                }
                if(cRand==NULL)
                    return(-1.); //ERROR

                if(distancedxdy(cRand->neighbor)>publicParameters->distMaxPastAndFuture) //cRand not linked to this!
                    return(0);
                if(connectionBlinkFuture->distancedxdy(cRand->neighbor->connectionPast)>publicParameters->distMaxPastAndFuture) //connectionBlinkFuture not linked to cRand->neighbor->connectionPast
                    return(0);
                return(makeChangeBlinkLinkage2(limitLoss, cRand));
            }
            else
            {
                idRandNeighbor=rand() % nbNodeAtGoodTime2;
                i=0;
                exitFor=false;
                for(cRand2=neighborsFuture; (!exitFor)&&(cRand2!=NULL);)
                {
                    if((cRand2->neighbor->time-time==1)
                            &&(cRand2->neighbor->connectionPast==NULL)
                            &&(cRand2->neighbor->connectionFusion1Past==NULL)
                            &&(cRand2->neighbor->connectionFissionPast==NULL)
                            &&(cRand2->neighbor->connectionBlinkPast==NULL))
                    {
                        if(i==idRandNeighbor)
                            exitFor=true;
                        i++;
                    }
                    if(!exitFor)
                        cRand2=cRand2->next;
                }
                if(cRand2==NULL)
                    return(-1.); //ERROR

                if(distancedxdy(cRand2->neighbor)>publicParameters->distMaxPastAndFuture) //cRand2 not linked to this!
                    return(0);

                return(makeChangeBlinkLinkage4(limitLoss, cRand2));

            }
        }
        else   // cRand is just before connectionBlinkFuture
        {
            ChainNeighbor*cRand;
            ChainNeighbor*cRand2;
            long nbNodeAtGoodTime;
            nbNodeAtGoodTime=0;
            long nbNodeAtGoodTime2;
            nbNodeAtGoodTime2=0;
            cRand=NULL;
            cRand2=NULL;
            bool exitFor=false;

            long idRandNeighbor;

            for(i=0, cRand=connectionBlinkFuture->neighborsPast; cRand!=NULL; cRand=cRand->next)
                if(connectionBlinkFuture->time-cRand->neighbor->time==1)
                {
                    if(cRand->neighbor->connectionFuture!=NULL)
                    {
                        nbNodeAtGoodTime++;
                    }
                    else
                    {
                        if((cRand->neighbor->connectionFission1Future==NULL)
                                &&(cRand->neighbor->connectionFusionFuture==NULL)
                                &&(cRand->neighbor->connectionBlinkFuture==NULL))
                            nbNodeAtGoodTime2++;
                    }
                }
            if((nbNodeAtGoodTime==0)&&(nbNodeAtGoodTime2==0))
                return(makeChangeBlinkLinkage1(limitLoss));

            if(((rand()%2)&&(nbNodeAtGoodTime!=0))||(nbNodeAtGoodTime2==0))
            {

                idRandNeighbor=rand() % nbNodeAtGoodTime;
                i=0;
                exitFor=false;
                for(cRand=connectionBlinkFuture->neighborsPast; (!exitFor)&&(cRand!=NULL);)
                {
                    if((connectionBlinkFuture->time-cRand->neighbor->time==1)
                            &&(cRand->neighbor->connectionFuture!=NULL))
                    {
                        if(i==idRandNeighbor)
                            exitFor=true;
                        i++;
                    }
                    if(!exitFor)
                        cRand=cRand->next;
                }
                if(cRand==NULL)
                    return(-1.); //ERROR
                if(connectionBlinkFuture->distancedxdy(cRand->neighbor)>publicParameters->distMaxPastAndFuture) //cRand not linked to connectionBlinkFuture!
                    return(0);
                if(distancedxdy(cRand->neighbor->connectionFuture)>publicParameters->distMaxPastAndFuture) //this  not linked to cRand->neighbor->connectionFuture!
                    return(0);


                return(makeChangeBlinkLinkage3(limitLoss, cRand));
            }
            else
            {
                idRandNeighbor=rand() % nbNodeAtGoodTime2;
                i=0;
                exitFor=false;
                for(cRand2=connectionBlinkFuture->neighborsPast; (!exitFor)&&(cRand2!=NULL);)
                {
                    if((connectionBlinkFuture->time-cRand2->neighbor->time==1)
                            &&(cRand2->neighbor->connectionFuture==NULL)
                            &&(cRand2->neighbor->connectionFission1Future==NULL)
                            &&(cRand2->neighbor->connectionFusionFuture==NULL)
                            &&(cRand2->neighbor->connectionBlinkFuture==NULL))
                    {
                        if(i==idRandNeighbor)
                            exitFor=true;
                        i++;
                    }
                    if(!exitFor)
                        cRand2=cRand2->next;
                }
                if(cRand2==NULL)
                    return(-1.); //ERROR
                if(connectionBlinkFuture->distancedxdy(cRand2->neighbor)>publicParameters->distMaxPastAndFuture) //cRand2 not linked to connectionBlinkFuture!
                    return(0);
                return(makeChangeBlinkLinkage5(limitLoss, cRand2));
            }
        }
    }
}
double NodeDistanceIntensity::makeChangeSimpleLinkage(double limitLoss)
{
    int i;
    ChainNeighbor*cRand;
    ChainNeighbor*c;
    // is the difference between the new configuration likeness and the old.
    //if deltaLikeness>0 then the new solution is better

    /*//test if this is implicated in a fusion or a fission
    if(connectionFission1Future!=NULL) //! Fission
    	return(0.);
    if(connectionFusionFuture!=NULL) //! Fusion
    	return(0.);*/



    if((connectionFuture==NULL)||(rand()%2==0))
    {

        //cRand is a random node in the vicinity of "this" in the future.

        //compute index of cRand
        long idRandNeighbor;
        if(neighborsFutureCount==0)
            return(0.); //no neighbors
        idRandNeighbor=rand() % neighborsFutureCount;

        //set the node cRand
        cRand=NULL; //inutile
        for(i=0, cRand=neighborsFuture; (i<idRandNeighbor)&&(cRand!=NULL); cRand=cRand->next, i++);
        //test if the node is good
        if(cRand==NULL)
            return(-1.); //ERROR
        if(cRand->neighbor==connectionFuture) // the new random link is the same than the connectionFuture.
            return(makeChangeSimpleLinkage0(limitLoss)); //try simple deconnection


        if(cRand->neighbor->connectionFissionPast!=NULL) // cRand already implicated in a fission
            return(0);
        if(cRand->neighbor->connectionFusion1Past!=NULL) // cRand already implicated in a fusion
            return(0);
        if(cRand->neighbor->connectionBlinkPast!=NULL)// cRand is already implicated in a blink connection
            return(0);



        //is there a conflict between this and the pastConnextion of cRand->neighbor?
        if(cRand->neighbor->connectionPast==NULL)  //no conflict cRand->neighbor can be linked with this.
        {
            //	From
            //	Present		Future
            //	1(this)-----2
            //	NULL--------3(cRand)

            //	To
            //	Present		Future
            //	1(this)-----3(cRand)
            //	NULL--------2

            return(makeChangeSimpleLinkage1(limitLoss, cRand));

        }
        else   //there is a conflict, the node cRand->neighbor is already connected to a node at the same time than this
        {
            if((rand()%2)&&(connectionFuture!=NULL))  //just flip the 2 nodes
            {
                //	From
                //	Present		Future
                //	1(this)-----2
                //	4-----------3(cRand)
                //	To
                //	Present		Future
                //	1-----------3(cRand)
                //	4-----------2
                return(makeChangeSimpleLinkage4(limitLoss, cRand));


            }
            else   //find a new node to be connected
            {
                // try to find a free node to connect the connectionPast of the conflict node.
                //count the number of free node
                //	From
                //	Present		Future
                //	1(this)-----2
                //	4-----------3(cRand)
                int count;
                count=0;
                for(c=cRand->neighbor->connectionPast->neighborsFuture; c!=NULL; c=c->next)
                {
                    if((c->neighbor->connectionFusion1Past==NULL)&&(c->neighbor->connectionFissionPast==NULL))
                    {
                        if(c->neighbor->connectionPast==NULL)
                        {
                            count++;
                        }
                        else
                        {
                            if(c->neighbor==connectionFuture)
                            {
                                count++;
                            }
                        }
                    }
                }
                if(count==0)    //no node available to relink the conflicted node
                {
                    // try to break the link
                    //	From
                    //	Present		Future
                    //	1(this)-----2
                    //	4-----------3(cRand)
                    //	To
                    //	Present		Future
                    //	1-----------3(cRand)
                    //	NULL--------2
                    //	4-----------NULL
                    return(makeChangeSimpleLinkage2(limitLoss, cRand));

                }
                else    //there is several available nodes to relink the conflicted node
                {
                    //	From
                    //	Present		Future
                    //	1-----------2
                    //	4-----------3(cRand)
                    //	NULL--------5(c)

                    //	To
                    //	Present		Future
                    //	1-----------3(cRand)
                    //	4-----------5(c)
                    //	NULL--------2
                    return(makeChangeSimpleLinkage3(limitLoss, cRand, count));
                }
            }
        }
    }
    else
    {
        return(makeChangeSimpleLinkage0(limitLoss));
    }
}

