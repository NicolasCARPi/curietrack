/*******************************************************************************************************************
								Link class to treat data from segmentation of images 2D (16bit only)
									Institut Curie
									UMR - 144
									by Victor Racine
									april 25th 2003
*******************************************************************************************************************/

#ifndef __link_2D__
#define __link_2D__

#include <stdio.h>


typedef struct {
	double	x;
	double	y;
}	val2D;


typedef struct {
	long	x;
	long	y;
}	coor2D;

typedef struct {
	float	maxPotentialWithNeighbor;
	long	structureInContact;
	float	minPotentialOnTheStructure;
	float	maxPotentialOnTheStructure;
	long x;
	long y;
	bool kept;
	bool locked;
	long structureToBecaume;
	long newIndex;
}	wsFusion2D;


typedef struct {
	double			val;
	unsigned long	x;
	unsigned long	y;
}	maxStruct2D;

typedef struct {
	double	xx;
	double	yy;
	double	xy;
}	varianceStruct2D;


typedef struct {
	unsigned long		peri;
	unsigned long		surf;
	maxStruct2D			max;
	double				morpho;
	val2D				moment;
	val2D				momentSurf;
	double				mean;
	double				intensity;
	varianceStruct2D	var;
	val2D				sigma;
	double				angle;
	unsigned long		index;
	unsigned long		time;
	long		maxLength;
} InfoData2D;

typedef struct {
	int		time;
	int		index;
	val2D				moment;
	unsigned long		surf;
	unsigned long		peri;
	double				mean;
} InfoData2DAscii;


class InfoStruct2D{
public:
	InfoStruct2D(int index);
	~InfoStruct2D();
	InfoStruct2D* next;

	InfoData2D data;
};



#endif
