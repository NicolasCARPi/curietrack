/*******************************************************************************************************************
    TabTimes Class
    Array containing all the TimeNodes for all instants
    Institut Curie
    UMR - 144
    by Victor Racine
    12 july 2002
    modified in 2011 Perrine Gilloteaux
*******************************************************************************************************************/

/********************************************************************************************************************
///Class TabTimes contains several TabTimes objets for each time (nbTimes).
********************************************************************************************************************/
#include "tab_times.h"
#include <stdio.h>
using namespace std;

extern FILE *outputStream;

extern PublicParameters *publicParameters;

#ifndef _countof
#define _countof(X) (sizeof(X) / sizeof(X[0]))
#endif

/********************************************************************************************************************
///Constructuor(int _nbTimes) of TabTimes class.
Create nbTimes objets TimeNodes.
********************************************************************************************************************/
TabTimes::TabTimes(){
	publicParameters->nbNodeAllocated=0;
	publicParameters->nbNodeDeallocated=0;
	publicParameters->nbChainNeighborDeallocated=0;
	publicParameters->nbChainNeighborAllocated=0;
	timeNodes=NULL;
	nbTimes=0;
}


/********************************************************************************************************************
Destructuor() of TabTimes class.
Delete all objets TimeNodes.
********************************************************************************************************************/
TabTimes::~TabTimes(){
	for(int i=0;i<nbTimes;i++)
		delete timeNodes[i];
	delete [] timeNodes;
	if((publicParameters->verbose) || (publicParameters->nbNodeAllocated!=publicParameters->nbNodeDeallocated)
			|| (publicParameters->nbChainNeighborAllocated!=publicParameters->nbChainNeighborDeallocated)){
	}
}

/********************************************************************************************************************
allocTabTimes() alloc the timeNodes array of nbTimes long
********************************************************************************************************************/
void TabTimes::allocTabTimes(){
	timeNodes= new TimeNodes*[nbTimes];
}



/********************************************************************************************************************
initialisation() initialisation of the data in memory
********************************************************************************************************************/
int TabTimes::initialisation(double** objectstotrack, int allobjectsinframes,int nbframesinmovie,int sizeimageX, int sizeimageY){
	int ret;
	if((ret=readStatParameters( objectstotrack, allobjectsinframes,nbframesinmovie,sizeimageX,sizeimageY))<0) return(ret);
	if((ret=linkNeighborsPresent())<0) return(ret);
	if((ret=linkNeighborsBeforeAndAfter())<0) return(ret);
	if((ret=sortByIntensityAllNeighbors())<0) return(ret);

	if(publicParameters->nbTemporalWindows>1)
		if((ret=linkNeighborsBlinkBeforeAndAfter())<0)
			return(ret);


	srand( (unsigned)time( NULL ) );
	return(RET_OK);
}


int TabTimes::printConnectionsASCIIwithoutput( std::vector<std::vector<double> > *outputVector){
	//ccFILE* stream;

	std::vector<double> structurecharacteristics;

	for(int t=0;t<nbTimes;t++){
		TimeNodes* _localTimeNodes=timeNodes[t];
		for(int n=0;n<_localTimeNodes->nbNodes;n++){
			Node* node=_localTimeNodes->tabNodes[n];
			if((node->connectionPast==NULL)&&(node->trackedStructureIndex>0)){
				for(Node* nodeConnection=node;nodeConnection!=NULL;nodeConnection=nodeConnection->connectionFuture){
					structurecharacteristics.push_back(node->trackedStructureIndex); //color tag
					structurecharacteristics.push_back(nodeConnection->time); //frame
					structurecharacteristics.push_back(nodeConnection->info->getIndex());

					if  (nodeConnection->connectionFissionPast!=NULL){
					    structurecharacteristics.push_back( nodeConnection->connectionFissionPast->trackedStructureIndex);
						}
					else
						{
						structurecharacteristics.push_back( 0);
						}
					structurecharacteristics.push_back(node->neighborsPresentCountCluster);

			    outputVector->push_back(structurecharacteristics);
				structurecharacteristics.clear();
					}

				}
			}
		}
	return(RET_OK);
}

/********************************************************************************************************************
readStatParameters() reads the file specified in publicParameters and allocates the corresponding TimeNodes
********************************************************************************************************************/
int TabTimes::readStatParameters(double** objectstotrack, int allobjectsinframes, int nbframesinmovie, int sizeimageX, int sizeimageY){
  // some frames can have no objects segmented inside.
    //unsigned long nbTimes;
    unsigned long nbNodes;
    int time;
    int count=0;
    int i,j;
    TimeNodes * localTimeNodes;

    int node=0;

    nbTimes=nbframesinmovie;
    allocTabTimes(); //Allocation of tabTimes

    long incre=0;

    publicParameters->maxX=sizeimageX;
    publicParameters->maxY=sizeimageY;

    for(time=0;time<(long)nbTimes;time++){ //for time #1
        nbNodes=0;
        for (i=0;i<allobjectsinframes;i++){ //for i#1
            // count the nodes in this frame
            // initialization: we do not take the max number of nodes since taggs can be arbitrarily set
            if (objectstotrack[i][0]==time)
                nbNodes++;
        }// end for i#1
		node=0;

        //ccprintf( "\nnbnode: %i\n",nbNodes);
        timeNodes[time]=allocTimeNodes(time,nbNodes);  //Allocation of timeNodes
                    //new TimeNodes(time,(int)nbNodes);
        localTimeNodes=timeNodes[time];

        if(localTimeNodes==NULL){
           //fprintf(outputStream, "Error when allocation the nodes of the time %u\n", time);
            return(T_ERR_ALLOC);
        }//end if



        int test1,test2;
        float test3,test4,test7, test5,test6;
        for (i=0;i<allobjectsinframes;i++){//for i#2

            if (objectstotrack[i][0]==time){//if objtotrack
                    test1=int(objectstotrack[i][0]); // time
                    test2=int(objectstotrack[i][1]); //index
                    test3=float(objectstotrack[i][2]); //x
                    test4=float(objectstotrack[i][3]); //y
                    test5=float(objectstotrack[i][4]); //area
                    test6=float(objectstotrack[i][4]); //area again, not used
                    test7=float(objectstotrack[i][6]); //question INtegrated intensity or not???



                    (*localTimeNodes)(node)->info->convertinfoAsciiinfo(test1,test2,test3,test4,test5,test6,test7);

                    if((*localTimeNodes)(node)->info->getIntensity()<0){
                        printf("%ld : <0\n", ++incre);
                    }

// max X and max Y should be the image size for the border computation

               /*     if(publicParameters->maxX<(*localTimeNodes)(node)->info->getX())
                        publicParameters->maxX=(*localTimeNodes)(node)->info->getX();
                    if(publicParameters->maxY<(*localTimeNodes)(node)->info->getY())
                        publicParameters->maxY=(*localTimeNodes)(node)->info->getY();*/
					node++;
            }//end if objtotrack


                ////if(publicParameters->verbose)
                    //**Global Likeness:
                   //// fprintf(outputStream, "Time %u: %u nodes loaded\n", time, nbNodes);
			}// end for i#2

            //Sorting the nodes by the intensity for one time
            Node* nodeMax;
            int indexMax;
            double maxIntensity;
             for(i=0;i<timeNodes[time]->nbNodes-1;i++){ //for i#3
                    maxIntensity=-1.;
                    for(j=i;j<localTimeNodes->nbNodes;j++)
                        if(maxIntensity< (*localTimeNodes)(j)->info->getIntensity()){
                            maxIntensity=(*localTimeNodes)(j)->info->getIntensity();
                            indexMax=j;
                        }
                    nodeMax=localTimeNodes->tabNodes[indexMax];
                    localTimeNodes->tabNodes[indexMax]=localTimeNodes->tabNodes[i];
                    localTimeNodes->tabNodes[i]=nodeMax;
               } //end for i#3

   } //end for time#1
return(RET_OK);

}


/********************************************************************************************************************
operator()(int i) to access to all TimeNodes objets.
********************************************************************************************************************/
TimeNodes *TabTimes::operator()(int i){
	return(timeNodes[i]);
}

/********************************************************************************************************************
operator[](int i) to access to pointors to all TimeNodes objets.
********************************************************************************************************************/
TimeNodes **TabTimes::operator[](int i){
	return(timeNodes+i);
}

/********************************************************************************************************************
int linkNeighborsPresent() to link of the neighbors in the presents of all nodes of all times
********************************************************************************************************************/
int TabTimes::linkNeighborsPresent(){
	int ret;
	for(long time=0;time<nbTimes;time++)
		if((ret=timeNodes[time]->linkNeighborsPresent())<0) return(ret);
	return(RET_OK);
}


/********************************************************************************************************************
int linkNeighborsBeforeAndAfter() to link of the neighbors in the past and in the future of all nodes of all times
********************************************************************************************************************/
int TabTimes::linkNeighborsBeforeAndAfter(){
	int ret;
	for(long time=0;time<nbTimes-1;time++)
		if((ret=timeNodes[time]->linkNeighborsBeforeAndAfter(*(timeNodes[time+1])))<0) return(ret);
	return(RET_OK);
}

/********************************************************************************************************************
int linkNeighborsBlinkBeforeAndAfter() to link of the neighbors in the past and in the future of all nodes of
all times for blinking purpose
********************************************************************************************************************/
int TabTimes::linkNeighborsBlinkBeforeAndAfter(){
	int ret;
	for(long time=0;time<nbTimes;time++)
		for(long dt=2;dt<=publicParameters->nbTemporalWindows;dt++)
			if(time+dt<nbTimes)
				if((ret=timeNodes[time]->linkNeighborsBlinkBeforeAndAfter(*(timeNodes[time+dt])))<0) return(ret);
	return(RET_OK);
}

/********************************************************************************************************************
int sortByIntensityAllNeighbors() sort by intensity of the node neighbors of each time
********************************************************************************************************************/
int TabTimes::sortByIntensityAllNeighbors(){
	int ret;
	for(long time=0;time<nbTimes;time++)
		if((ret=timeNodes[time]->sortByIntensityAllNeighbors())<0) return(ret);
	return(RET_OK);

}

/********************************************************************************************************************
void setIndexToTrackedStructures() set a unique index to each tracked structures
the fusion and fission split the strucutres.
********************************************************************************************************************/
void TabTimes::setIndexToTrackedStructures(){
	unsigned long newIndex=1;
	long nbSuccessiveNodes;
	Node* node;
	Node* oldNode;
	bool isNewIndex;
	for(long time=0;time<nbTimes;time++){
		for(long indexNode=0;indexNode<timeNodes[time]->nbNodes;indexNode++){
			if (timeNodes[time]->tabNodes[indexNode]->connectionBlinkFuture!=NULL) {
				timeNodes[time]->tabNodes[indexNode]->connectionFuture=timeNodes[time]->tabNodes[indexNode]->connectionBlinkFuture;
				timeNodes[time]->tabNodes[indexNode]->connectionFuture->connectionPast=timeNodes[time]->tabNodes[indexNode];
				// blinking taken into account
				}
			if(timeNodes[time]->tabNodes[indexNode]->connectionPast==NULL){
				for(node=timeNodes[time]->tabNodes[indexNode], oldNode=NULL,nbSuccessiveNodes=0;
					(node!=NULL)&&(nbSuccessiveNodes<publicParameters->nbMinSuccessiveStructures);
					node=node->connectionFuture,nbSuccessiveNodes++){//count number of successive nodes
						oldNode=node;
						node->trackedStructureIndex=0;
					}

				isNewIndex=false;
				if(nbSuccessiveNodes>=publicParameters->nbMinSuccessiveStructures)//it was ==
					isNewIndex=true;
				if(timeNodes[time]->tabNodes[indexNode]->connectionFissionPast!=NULL)
					isNewIndex=true;
				if(timeNodes[time]->tabNodes[indexNode]->connectionFusion1Past!=NULL)
					isNewIndex=true;
				//if(timeNodes[time]->tabNodes[indexNode]->connectionBlinkPast!=NULL)
				//	isNewIndex=true; //si la cellule n'avait pas de passé , alors elle a un nouveau numero uniquement si elle vient d'une division
				if(oldNode!=NULL){
					if(oldNode->connectionFission1Future!=NULL)
						isNewIndex=true;
					if(oldNode->connectionFusionFuture!=NULL)
						isNewIndex=true;
					if(oldNode->connectionBlinkFuture!=NULL)
						isNewIndex=true;
				}

				if(isNewIndex){
					timeNodes[time]->tabNodes[indexNode]->trackedStructureIndex=newIndex;
					newIndex++;
				}
            }else{
				timeNodes[time]->tabNodes[indexNode]->trackedStructureIndex=
					timeNodes[time]->tabNodes[indexNode]->connectionPast->trackedStructureIndex;
            }
		}
    }
}

/********************************************************************************************************************
int giveResult() just call the function to display the result
********************************************************************************************************************/
int TabTimes::giveResult(std::vector<std::vector<double> > *output){
	int ret;

	setIndexToTrackedStructures();

	if((ret=printConnectionsASCIIwithoutput(output))<0) return(ret);

	return(RET_OK);
}
