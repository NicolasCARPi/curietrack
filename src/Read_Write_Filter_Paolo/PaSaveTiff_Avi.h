/*#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
//#include <opencv2/video/background_segm.hpp>
//#include <opencv2/imgproc/imgproc_c.h>
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
void SaveTiffOutputNameStr(string nameout, vector<Mat>& outframes, int& ntot,  int& imagewidth, int& imagelength)
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
void SaveTiffOutputNameStr8bit(string nameout, vector<Mat>& outframes, int& ntot,  int& imagewidth, int& imagelength)
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
///////////////////////////
void Save8bitcolor(string nameout, vector<Mat>& outframes)
{
	int ntot = outframes.size();
	int imagewidth, imagelength;
	imagewidth = outframes[0].cols;
	imagelength = outframes[0].rows;
	Size WriOut_size=Size(imagewidth, imagelength);
	double WriOut_FPS=5;
	if (ntot < 10)
		WriOut_FPS = 2;
	if (ntot < 5)
		WriOut_FPS = 1;

	//int WriOut_fourCC = CV_FOURCC('I','4','2','0');
	//int WriOut_fourCC = CV_FOURCC('I','Y','U','V');
	//int WriOut_fourCC = CV_FOURCC('D','I','B',' ');
	//int WriOut_fourCC = CV_FOURCC('P', 'I', 'M', '1'); // = MPEG-1 codec
	int WriOut_fourCC = CV_FOURCC('M', 'J', 'P', 'G'); // = motion-jpeg codec (does not work well)
	//int WriOut_fourCC = CV_FOURCC('M', 'P', '4', '2'); // = MPEG-4.2 codec
	//int WriOut_fourCC = CV_FOURCC('D', 'I', 'V', '3'); // = MPEG-4.3 codec
	//int WriOut_fourCC = CV_FOURCC('D', 'I', 'V', 'X'); // = MPEG-4 codec
	//int WriOut_fourCC = CV_FOURCC('U', '2', '6', '3'); // = H263 codec
	//int WriOut_fourCC = CV_FOURCC('I', '2', '6', '3'); // = H263I codec
	//int WriOut_fourCC = CV_FOURCC('F', 'L', 'V', '1'); // = FLV1 codec

	VideoWriter WriOut(nameout, WriOut_fourCC, WriOut_FPS, WriOut_size, true);
	int pos=0;
	for (pos=0; pos<ntot; pos++)
		{
		Mat temp_outframe = outframes[pos]; //, temp_outframe_color;
		WriOut << temp_outframe;
		}
}

void Save8bitcolorTiff(string nameout, vector<Mat>& outframes)
{
	int ntot = outframes.size();
	int imagewidth, imagelength;
	imagewidth = outframes[0].cols;
	imagelength = outframes[0].rows;

	TIFF* OutputTif = TIFFOpen(nameout.c_str(), "w");
	for (int i=0; i<ntot; i++)
	{
		TIFFSetField(OutputTif, TIFFTAG_IMAGELENGTH, imagelength);
		TIFFSetField(OutputTif, TIFFTAG_IMAGEWIDTH, imagewidth);
		TIFFSetField(OutputTif, TIFFTAG_BITSPERSAMPLE, 8);
		TIFFSetField(OutputTif, TIFFTAG_SAMPLESPERPIXEL, 3);
		TIFFSetField(OutputTif, TIFFTAG_PLANARCONFIG, 1);
		TIFFSetField(OutputTif, TIFFTAG_COMPRESSION, 1);
		TIFFSetField(OutputTif, TIFFTAG_SAMPLEFORMAT, 1);
		for (int l=0; l<imagelength; l++)
			{
			Vec3b LineBuf[imagewidth];
			//LineBuf = (unsigned short *)_TIFFmalloc(imagewidth);
			for (int w=0; w<imagewidth; w++)
				LineBuf[w] = outframes[i].at<Vec3b>(l,w);
			TIFFWriteScanline(OutputTif, LineBuf, l, 0);
			}
		if (i < ntot-1)
			TIFFWriteDirectory(OutputTif);
		TIFFFlush(OutputTif);
	}
	TIFFClose(OutputTif);
}
