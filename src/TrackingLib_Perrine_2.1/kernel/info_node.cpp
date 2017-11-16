/*******************************************************************************************************************
    InfoNode Class
    InfoNode contains the information for a given Node in 2D or in 3D
    Institut Curie
    UMR - 144
    by Victor Racine
    2004/06/24
    modified by Perrine 2011
*******************************************************************************************************************/

#include "info_node.h"

extern PublicParameters *publicParameters;

InfoNode::InfoNode(){

}

float InfoNode2D::distancedxdy(InfoNode*infoNode){
	float res;
	double dX=getX()-infoNode->getX();
	double dY=(getY()-infoNode->getY())*publicParameters->dydX;// Perrine: Modified to favorize displacement in y
	res=(float)sqrt(dX*dX+dY*dY);
	if(res<0.) return(0.);
	return(res);
}
//compute the distance to the center of gravity
float InfoNode2D::distancedxdy(InfoNode* i1, InfoNode*i2){
	float res;
	double n1n2X=(i1->getX()+i2->getX())/2;
	double n1n2Y=(i1->getY()+i2->getY())/2;
	double dX=getX()-n1n2X;
	double dY=(getY()-n1n2Y)*publicParameters->dydX;
	res=(float)sqrt(dX*dX+dY*dY);
	if(res<0.)
		return(0.);
	if(res<1e-9)
		res+=0;
	return(res);
}

float InfoNode2D::distance(InfoNode*infoNode){
	float res;
	double dX=getX()-infoNode->getX();
	double dY=(getY()-infoNode->getY());// Perrine: Modified to favorize displacement in y
	res=(float)sqrt(dX*dX+dY*dY);
	if(res<0.) return(0.);
	return(res);
}
//compute the distance to the center of gravity
float InfoNode2D::distance(InfoNode* i1, InfoNode*i2){
	float res;
	double n1n2X=(i1->getX()+i2->getX())/2;
	double n1n2Y=(i1->getY()+i2->getY())/2;
	double dX=getX()-n1n2X;
	double dY=(getY()-n1n2Y);
	res=(float)sqrt(dX*dX+dY*dY);
	if(res<0.)
		return(0.);
	if(res<1e-9)
		res+=0;
	return(res);
}
double InfoNode2D::getMean(){
	return(info.mean);
}
double InfoNode2D::getIntensity(){
	return(info.intensity);
}
unsigned long InfoNode2D::getIndex(){
	return(info.index);
}
unsigned long InfoNode2D::getTime(){
	return(info.time);
}

long InfoNode2D::getSurface(){
	return(info.surf);
}

double InfoNode2D::getX(){
	return(info.moment.x);
}

double InfoNode2D::getY(){
	return(info.moment.y);
}

unsigned long InfoNode2D::sizeofInfoData(){
	return(sizeof(InfoData2D));
}

unsigned long InfoNode2D::sizeofInfoDataAscii(){
	return(sizeof(InfoData2DAscii));
}

void* InfoNode2D::getPtrInfoData(){
	return(&info);
}
void * InfoNode2D::getPtrInfoDataAscii(){
	return(&infoAscii);
}
void InfoNode2D::convertinfoAsciiinfo(int i1,int i2,float i3,float i4, float i5,float i6,float i7){
	info.moment.x=i3;
	info.moment.y=i4;
	info.surf=i5;
	info.intensity=i7*i5;
	info.mean=i7;
	info.index=i2;
	info.time=i1;

	return;
}
