/*******************************************************************************************************************
    NodeDistanceIntensity  Class
    It is derived from the class Nodes
    This class find the neighbor with the best likeness.
    Likeness depend on the distance and the intensity.
    Institut Curie
    UMR - 144
    by Victor Racine
    2004 06 24
*******************************************************************************************************************/

#ifndef __node_distance_intensity__
#define __node_distance_intensity__

#include "../kernel/nodes.h"
#include "../kernel/general.h"


#define STD_NOISE	40

class NodeDistanceIntensity:public Node{
private:
	double getLikeness(double difference, double sigma);
	double likenessDistance(Node* n);
	double likenessDistanceFusion(Node* n1, Node* n2);
	double likenessFusion(Node* n1, Node* n2);
	double likenessFission(Node* n1, Node* n2);
	double likenessIntensity(Node* n);
	double likenessIntensityFusion(Node* n1, Node* n2);
	double likenessMeanIntensity(Node* n);
	double likenessMeanIntensityFusion(Node* n1, Node* n2);
	double makeChangeSimpleLinkage(double limit);
	double makeChangeSimpleLinkage0(double limit);
	double makeChangeSimpleLinkage1(double limit, ChainNeighbor*cRand);
	double makeChangeSimpleLinkage2(double limit, ChainNeighbor*cRand);
	double makeChangeSimpleLinkage3(double limit, ChainNeighbor*cRand, long count);
	double makeChangeSimpleLinkage4(double limit, ChainNeighbor*cRand);
	double makeChangeFusion(double limit);
	double makeChangeFusion0(double limitLoss);
	double makeChangeFusion1(double limitLoss, ChainNeighbor*cRand);
	double makeChangeFusion2(double limitLoss, ChainNeighbor*cRand);
	double makeChangeFusion3(double limitLoss);
	double makeChangeFusion4(double limitLoss);
	double makeChangeFusion5(double limitLoss, ChainNeighbor*cRand);
	double makeChangeFusion6(double limitLoss, ChainNeighbor*cRand);
	double makeChangeFusion7(double limitLoss, ChainNeighbor*cRand);
	double makeChangeFusion8(double limitLoss, ChainNeighbor*cRand);
	double makeChangeFusion9(double limitLoss, ChainNeighbor*cRand);
	double makeChangeFission(double limitLoss);
	double makeChangeFission0(double limitLoss);
	double makeChangeFission1(double limitLoss, ChainNeighbor*cRand);
	double makeChangeFission2(double limitLoss, ChainNeighbor*cRand);
	double makeChangeFission3(double limitLoss);
	double makeChangeFission4(double limitLoss);
	double makeChangeFission5(double limitLoss, ChainNeighbor*cRand);
	double makeChangeFission6(double limitLoss, ChainNeighbor*cRand);
	double makeChangeFission7(double limitLoss, ChainNeighbor*cRand);
	double makeChangeFission8(double limitLoss, ChainNeighbor*cRand);
	double makeChangeFission9(double limitLoss, ChainNeighbor*cRand);
	double makeChangeBlinkLinkage(double limitLoss);
	double makeChangeBlinkLinkage0(double limitLoss);
	double makeChangeBlinkLinkage1(double limitLoss);
	double makeChangeBlinkLinkage2(double limitLoss, ChainNeighbor*cRand);
	double makeChangeBlinkLinkage3(double limitLoss, ChainNeighbor*cRand);
	double makeChangeBlinkLinkage4(double limitLoss, ChainNeighbor*cRand);
	double makeChangeBlinkLinkage5(double limitLoss, ChainNeighbor*cRand);
	double stericEffect(Node* n1, Node* n2);
	bool metropolis(double deltaE, double temperature);
public:
	void computeConstantSigmas(double &stdIntensityComputing, long &countStdInt,
		double &stdMeanComputing, long &countStdMean, double &stdDistanceComputing, long &countStdDist, bool sigmaClipping, long*hist);
	double getLikenessOfLocalNode(int lastTime);
	double makeChange(double limit);
	double likeness(Node* n);
	double likenessBlink(Node* n);
	NodeDistanceIntensity(int _time);
	double makeInitializationConnection(void);


};

#endif //__node_distance_intensity__
