/*******************************************************************************************************************
    InfoNode InfoNode2D Classes
    InfoDataXX contains the information for a given Node in 2D or in 3D
    Institut Curie
    UMR - 144
    by Victor Racine
    2004/06/24
    modified Perrine 2011
*******************************************************************************************************************/

#ifndef __info_node__
#define __info_node__

#include "general.h"
#include "info_2D.h"
#include "../generalParam/NoError.h"

class InfoNode{
public:
	InfoNode();
	~InfoNode(){};
	virtual float distance(InfoNode*infoData){return(T_ERR_OOP);};
	virtual float distance(InfoNode* i1, InfoNode*i2){return(T_ERR_OOP);};
	virtual float distancedxdy(InfoNode*infoData){return(T_ERR_OOP);};
	virtual float distancedxdy(InfoNode* i1, InfoNode*i2){return(T_ERR_OOP);};

	virtual double getMean(){return(T_ERR_OOP);};
	virtual long getVolume(){return(T_ERR_OOP);};
	virtual long getSurface(){return(T_ERR_OOP);};
	virtual double getX(){return(T_ERR_OOP);};
	virtual double getY(){return(T_ERR_OOP);};
	virtual double getIntensity(){return(T_ERR_OOP);};
	virtual unsigned long getIndex(){return(T_ERR_OOP);};
	virtual unsigned long getTime(){return(T_ERR_OOP);};
	virtual unsigned long sizeofInfoData(){return(T_ERR_OOP);};
	virtual unsigned long sizeofInfoDataAscii(){return(T_ERR_OOP);};
	virtual void* getPtrInfoData(){return(NULL);};
	virtual void* getPtrInfoDataAscii(){return(NULL);};
	virtual void convertinfoAsciiinfo(int i1,int i2,float i3,float i4, float i5,float i6,float i7){return;};
};

class InfoNode2D: public InfoNode{
public:
	InfoNode2D():InfoNode(){};
	float distance(InfoNode *infoData);
	float distance(InfoNode* i1, InfoNode*i2);
	float distancedxdy(InfoNode *infoData);
	float distancedxdy(InfoNode* i1, InfoNode*i2);
	InfoData2D info;
	InfoData2DAscii infoAscii;
	double getMean();
	double getX();
	double getY();
	double getIntensity();
	long   getSurface();
	unsigned long getIndex();
	unsigned long getTime();
	unsigned long sizeofInfoData();
	unsigned long sizeofInfoDataAscii();
	void* getPtrInfoData();
	void* getPtrInfoDataAscii();
	void convertinfoAsciiinfo(int i1,int i2,float i3,float i4, float i5,float i6,float i7);
};

#endif //__info_node__
