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
		////fprintf(outputStream,"node alloc:%d,\tdealloc:%d\n", publicParameters->nbNodeAllocated,
		////	publicParameters->nbNodeDeallocated);
		////fprintf(outputStream,"neighbors alloc:%d,\tdealloc:%d\n", publicParameters->nbChainNeighborAllocated,
		////	publicParameters->nbChainNeighborDeallocated);
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


/********************************************************************************************************************
printConnectionsASCII() print the result of the connections in a ASCII file
Colonne 1         numero de la trace
Colonne 2         numero de l’image, ou temps
Colonne 3		   Former tag
Colonne 4         X
Colonne 5         Y
Colonne 6         Intensity
Colonne 7         Surface
Colonne 8			Mother
Colonne 9			 Daughter 1
Colonne 10			Daughter 2
+ print the result in division ascii:
Frame_Mother Tag_ori_M Tag_new_M Frame_D1 Tag_ori_D1 Tag_New_D2 Tag_ori_D2
********************************************************************************************************************/
/*int TabTimes::printConnectionsASCII(){
	FILE* stream;
	printf("\n check the end: result file is :%25s \n",publicParameters->fileResASCII);
	if( (stream=fopen(publicParameters->fileResASCII,"w"))==NULL){
		printf("The file %ls can not be opened\n", publicParameters->fileResASCII );
		return(T_ERR_OPENNING_FILE);
	}
	FILE* streamdivision;
	if( (streamdivision=fopen(publicParameters->fileResASCIIdivision,"w"))==NULL){
		printf("The file %ls can not be opened\n", publicParameters->fileResASCIIdivision );
		return(T_ERR_OPENNING_FILE);
	}

	fprintf(stream,"StructureTag\tTime\tFormerTag\tX\tY\tIntensity\tArea\tMother Tag\tDaugther1 Tag\tDaughter2 tag \n");
	fprintf(streamdivision,"Time\tMother Tag\t Mother Former Taf \t Daughter 1 tag\t Daughter formertag \t Daughter2 tag \t Daughter 2 formertag \n");
	for(int t=0;t<nbTimes;t++){
		TimeNodes* _localTimeNodes=timeNodes[t];
		for(int n=0;n<_localTimeNodes->nbNodes;n++){
			Node* node=_localTimeNodes->tabNodes[n];
			if((node->connectionPast==NULL)&&(node->trackedStructureIndex>0)){
				for(Node* nodeConnection=node;nodeConnection!=NULL;nodeConnection=nodeConnection->connectionFuture){
					fprintf(stream, "%d\t", node->trackedStructureIndex);
					fprintf(stream, "%d\t", nodeConnection->time);
					fprintf(stream, "%d\t", nodeConnection->info->getIndex()); // former index for reference
					fprintf(stream, "%f\t", nodeConnection->info->getY());
					fprintf(stream, "%f\t", nodeConnection->info->getX());
					fprintf(stream, "%f\t", nodeConnection->info->getIntensity());

					fprintf(stream, "%d\t", nodeConnection->info->getSurface());
					if  (nodeConnection->connectionFissionPast!=NULL)
						fprintf(stream, "%d\t", nodeConnection->connectionFissionPast->trackedStructureIndex); // add the mother as well
					else
						fprintf(stream,"0\t");

					if (nodeConnection->connectionFission1Future!=NULL) {//there is a fission division here Daughter 1 anbd Daughter 2

								fprintf(stream,"%d\t",nodeConnection->connectionFission1Future->trackedStructureIndex);
								fprintf(stream,"%d\n",nodeConnection->connectionFission2Future->trackedStructureIndex);
								// Other file specific to division Frame_Mother Tag_ori_M Tag_new_M Frame_D1 Tag_ori_D1 Tag_New_D2 Tag_ori_D2
								fprintf(streamdivision, "%d\t", nodeConnection->time);
								fprintf(streamdivision, "%d\t", node->trackedStructureIndex);
								fprintf(streamdivision, "%d\t", nodeConnection->info->getIndex()); // former index for reference
								fprintf(streamdivision,"%d\t",nodeConnection->connectionFission1Future->trackedStructureIndex);
								fprintf(streamdivision,"%d\t",nodeConnection->connectionFission1Future->info->getIndex());
								fprintf(streamdivision,"%d\t",nodeConnection->connectionFission2Future->trackedStructureIndex);
								fprintf(streamdivision,"%d\n",nodeConnection->connectionFission2Future->info->getIndex());

						}
						else{
								fprintf(stream,"0\t");
								fprintf(stream,"0\n");
						}


					}

				}
			}
		}
	fclose(stream);
	fclose(streamdivision);
	return(RET_OK);
}
*/

/********************************************************************************************************************
printConnectionsASCIIwithoutput() print the result of the connections in a ASCII file AND return the result in a vector
 In an 2D array vector of double[10]or 1D vector (to see)vector.h of vector from stdlib

 Structure For this: TrajectorytagColor|frame| FormerTag| TrajectoryTagMother(0 if no mother) |Numberofneighborsinthisframe

Colonne 1         numero de la trace
Colonne 2         numero de l’image, ou temps
Colonne 3		   Former tag
Colonne 4         X
Colonne 5         Y
Colonne 6         Intensity
Colonne 7         Surface
Colonne 8			Mother
Colonne 9			 Daughter 1
Colonne 10			Daughter 2
+ print the result in division ascii:
Frame_Mother Tag_ori_M Tag_new_M Frame_D1 Tag_ori_D1 Tag_New_D2 Tag_ori_D2
********************************************************************************************************************/
int TabTimes::printConnectionsASCIIwithoutput( std::vector<std::vector<double> > *outputVector){
	//ccFILE* stream;

	std::vector<double> structurecharacteristics;

	//ccprintf("\n check the end: result file is :%25s \n",publicParameters->fileResASCII);
	//ccif( (stream=fopen(publicParameters->fileResASCII,"w"))==NULL){
	//cc	printf("The file %ls can not be opened\n", publicParameters->fileResASCII );
	//cc	return(T_ERR_OPENNING_FILE);
	//cc}
	//ccFILE* streamdivision;
	//ccif( (streamdivision=fopen(publicParameters->fileResASCIIdivision,"w"))==NULL){
	//cc	printf("The file %ls can not be opened\n", publicParameters->fileResASCIIdivision );
	//cc	return(T_ERR_OPENNING_FILE);
	//cc}

	//ccfprintf(stream,"StructureTag\tTime\tFormerTag\tX\tY\tIntensity\tArea\tMother Tag\tDaugther1 Tag\tDaughter2 tag \n");
	//ccfprintf(streamdivision,"Time\tMother Tag\t Mother Former Taf \t Daughter 1 tag\t Daughter formertag \t Daughter2 tag \t Daughter 2 formertag \n");
	for(int t=0;t<nbTimes;t++){
		TimeNodes* _localTimeNodes=timeNodes[t];
		for(int n=0;n<_localTimeNodes->nbNodes;n++){
			Node* node=_localTimeNodes->tabNodes[n];
			if((node->connectionPast==NULL)&&(node->trackedStructureIndex>0)){
				for(Node* nodeConnection=node;nodeConnection!=NULL;nodeConnection=nodeConnection->connectionFuture){
					//ccfprintf(stream, "%d\t", node->trackedStructureIndex);
					structurecharacteristics.push_back(node->trackedStructureIndex); //color tag
					//ccfprintf(stream, "%d\t", nodeConnection->time);
					structurecharacteristics.push_back(nodeConnection->time); //frame
					//ccfprintf(stream, "%d\t", nodeConnection->info->getIndex()); // former index for reference
					structurecharacteristics.push_back(nodeConnection->info->getIndex());
					//ccfprintf(stream, "%f\t", nodeConnection->info->getY());
					//ccfprintf(stream, "%f\t", nodeConnection->info->getX());
					//ccfprintf(stream, "%f\t", nodeConnection->info->getIntensity());

					//ccfprintf(stream, "%d\t", nodeConnection->info->getSurface());
					if  (nodeConnection->connectionFissionPast!=NULL){
						//ccfprintf(stream, "%d\t", nodeConnection->connectionFissionPast->trackedStructureIndex); // add the mother as well
					    structurecharacteristics.push_back( nodeConnection->connectionFissionPast->trackedStructureIndex);
						}
					else
						{
						//ccfprintf(stream,"0\t");
						structurecharacteristics.push_back( 0);
						}
					structurecharacteristics.push_back(node->neighborsPresentCountCluster);
					if (nodeConnection->connectionFission1Future!=NULL) {//there is a fission division here Daughter 1 anbd Daughter 2

							//cc	fprintf(stream,"%d\t",nodeConnection->connectionFission1Future->trackedStructureIndex);
								//ccfprintf(stream,"%d\n",nodeConnection->connectionFission2Future->trackedStructureIndex);
								// Other file specific to division Frame_Mother Tag_ori_M Tag_new_M Frame_D1 Tag_ori_D1 Tag_New_D2 Tag_ori_D2
								//ccfprintf(streamdivision, "%d\t", nodeConnection->time);
								//ccfprintf(streamdivision, "%d\t", node->trackedStructureIndex);
								//ccfprintf(streamdivision, "%d\t", nodeConnection->info->getIndex()); // former index for reference
								//ccfprintf(streamdivision,"%d\t",nodeConnection->connectionFission1Future->trackedStructureIndex);
								//ccfprintf(streamdivision,"%d\t",nodeConnection->connectionFission1Future->info->getIndex());
								//ccfprintf(streamdivision,"%d\t",nodeConnection->connectionFission2Future->trackedStructureIndex);
								//ccfprintf(streamdivision,"%d\n",nodeConnection->connectionFission2Future->info->getIndex());

						}
					else{
								//ccfprintf(stream,"0\t");
								//ccfprintf(stream,"0\n");
						}

			    outputVector->push_back(structurecharacteristics);
				//printf("\n nb vector is:%d \n",outputVector->size());
				//printf("\n nb vector char is:%d \n",structurecharacteristics.size());
				structurecharacteristics.clear();
				//printf("\n nb vector is after:%d \n",outputVector->size());
				//printf("\n nb vector char is:%d \n",structurecharacteristics.size());
					}

				}
			}
		}
	//ccfclose(stream);
	//ccfclose(streamdivision);
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
    ////printf( "\nNumber of frames: %i\n", nbTimes);
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
        /**
                    void InfoNode2D::convertinfoAsciiinfo(int i1,int i2,float i3,float i4, float i5,float i6,float i7){
                    info.moment.x=i3;
                    info.moment.y=i4;
                    info.surf=i5;
                    info.intensity=i7*i5;
                    info.mean=i7;
                    info.index=i2;
                    info.time=i1;
        */
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
	for(long time=0;time<nbTimes;time++)
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
				//if (timeNodes[time]->tabNodes[indexNode]->connectionPast!=NULL){
				timeNodes[time]->tabNodes[indexNode]->trackedStructureIndex=
					timeNodes[time]->tabNodes[indexNode]->connectionPast->trackedStructureIndex;//}
				//else{
				//timeNodes[time]->tabNodes[indexNode]->trackedStructureIndex=
				//	timeNodes[time]->tabNodes[indexNode]->connectionBlinkPast->trackedStructureIndex;}
			}
		}

	//ccfprintf(outputStream, "%d tracked structures found with a minimum of %d times long\n\n", newIndex-1, publicParameters->nbMinSuccessiveStructures);
	//ccfprintf(outputStream, "The standard variation of the intensity sum is %f\n", publicParameters->stdIntensity);
	//ccfprintf(outputStream, "The standard variation of the intensity average is %f\n", publicParameters->stdMean);
	//ccfprintf(outputStream, "The standard variation of the distance is %f\n", publicParameters->stdDistance);

	}






/********************************************************************************************************************
int giveResult() just call the function to display the result
********************************************************************************************************************/
int TabTimes::giveResult(std::vector<std::vector<double> > *output){
	int ret;

	setIndexToTrackedStructures();

	//if(publicParameters->fileResASCII[0]!=L'\0')
	//ccif(publicParameters->fileResASCII != L'\0')
	if((ret=printConnectionsASCIIwithoutput(output))<0) return(ret);

	return(RET_OK);
}

