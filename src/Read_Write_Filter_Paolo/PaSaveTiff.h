/*#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/video/background_segm.hpp>
#include <opencv2/imgproc/imgproc_c.h>
#include <opencv2/opencv.hpp>*/
#include <gsl/gsl_statistics.h>
#include <fstream>
#include <iostream>
#include <list>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <time.h>
#include <vector>
#include <tiffio.h>
#include <tiff.h>

using namespace cv;
using namespace std;

/////////////////// Save th output 16bit name string!!
void SaveTiffOuputNameStr(string nameout, vector<Mat>& outframes, int& ntot,  int& imagewidth, int& imagelength)
{

	TIFF* OutputTif = TIFFOpen(nameout.c_str(), "w");
	for (int i=0; i<ntot; i++)
	{
		TIFFSetField(OutputTif, TIFFTAG_IMAGELENGTH, imagelength);
		TIFFSetField(OutputTif, TIFFTAG_IMAGEWIDTH, imagewidth);
		TIFFSetField(OutputTif, TIFFTAG_BITSPERSAMPLE, 16);
		TIFFSetField(OutputTif, TIFFTAG_SAMPLESPERPIXEL, 1);
		TIFFSetField(OutputTif, TIFFTAG_PLANARCONFIG, 1);
		TIFFSetField(OutputTif, TIFFTAG_COMPRESSION, 1);
		TIFFSetField(OutputTif, TIFFTAG_SAMPLEFORMAT, 1);
		for (int l=0; l<imagelength; l++)
			{
			unsigned short LineBuf[imagewidth];
			//LineBuf = (unsigned short *)_TIFFmalloc(imagewidth);
			for (int w=0; w<imagewidth; w++)
				LineBuf[w] = outframes[i].at<unsigned short>(l,w);
			TIFFWriteScanline(OutputTif, LineBuf, l, 0);
			}
		if (i < ntot-1)
			TIFFWriteDirectory(OutputTif);
		TIFFFlush(OutputTif);
	}
	TIFFClose(OutputTif);

}
/////////////////////// 8bit version
void SaveTiffOuputNameStr8bit(string nameout, vector<Mat>& outframes, int& ntot,  int& imagewidth, int& imagelength)
{

	TIFF* OutputTif = TIFFOpen(nameout.c_str(), "w");
	for (int i=0; i<ntot; i++)
	{
		TIFFSetField(OutputTif, TIFFTAG_IMAGELENGTH, imagelength);
		TIFFSetField(OutputTif, TIFFTAG_IMAGEWIDTH, imagewidth);
		TIFFSetField(OutputTif, TIFFTAG_BITSPERSAMPLE, 8);
		TIFFSetField(OutputTif, TIFFTAG_SAMPLESPERPIXEL, 1);
		TIFFSetField(OutputTif, TIFFTAG_PLANARCONFIG, 1);
		TIFFSetField(OutputTif, TIFFTAG_COMPRESSION, 1);
		TIFFSetField(OutputTif, TIFFTAG_SAMPLEFORMAT, 1);
		for (int l=0; l<imagelength; l++)
			{
			unsigned char LineBuf[imagewidth];
			//LineBuf = (unsigned short *)_TIFFmalloc(imagewidth);
			for (int w=0; w<imagewidth; w++)
				LineBuf[w] = outframes[i].at<unsigned char>(l,w);
			TIFFWriteScanline(OutputTif, LineBuf, l, 0);
			}
		if (i < ntot-1)
			TIFFWriteDirectory(OutputTif);
		TIFFFlush(OutputTif);
	}
	TIFFClose(OutputTif);

}
