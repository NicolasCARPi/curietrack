/*******************************************************************************************************************
    Link class to treat data from segmentation of images 2D (16bit only)
    Institut Curie
    UMR - 144
    by Victor Racine
    april 25th 2003
*******************************************************************************************************************/
#include "info_2D.h"

InfoStruct2D::InfoStruct2D(int index){
	data.index=index;
	next=NULL;

}

InfoStruct2D::~InfoStruct2D(){
}
