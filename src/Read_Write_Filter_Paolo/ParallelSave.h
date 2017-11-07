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

/////////////////// Save th output
void ParallelSaveOuput(char* argv[], vector<Mat>& outframes, int& ntot,  int& imagewidth, int& imagelength, string& nameout)
{
	nameout=string(argv[1]);
	string stravi=".tif";
	nameout.replace(nameout.find(stravi), stravi.length(), "_out.avi");
	//nameout="/home/paolo/Documents/OpenCV-soft/ActinWave&TractionField/out.avi";
	Size WriOut_size=Size(imagewidth, imagelength);
	double WriOut_FPS=5;
	//int WriOut_fourCC = -1;
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
	#pragma omp parallel for ordered schedule(dynamic)
	for (pos=0; pos<ntot; pos++)
		{
		#pragma omp ordered
			{
			#pragma omp critical (write)
				{
				Mat temp_outframe, temp_outframe_color;
				normalize(outframes[pos], temp_outframe, 0, 255, NORM_MINMAX, -1);
				temp_outframe.convertTo(temp_outframe, CV_8UC1);
				//imshow("pre saved", temp_outframe);
				cvtColor(temp_outframe, temp_outframe_color, CV_GRAY2RGB);
				WriOut << temp_outframe_color;
				//imshow("saved", temp_outframe_color);
				//int keycode = waitKey(50);
				}
			}
		}
}
/////////////////// Save th output color !!
void ParallelSaveOuputColor(char* argv[], vector<Mat>& outframes, int& ntot,  int& imagewidth, int& imagelength, string& nameout)
{
	nameout=string(argv[1]);
	string stravi=".tif";
	nameout.replace(nameout.find(stravi), stravi.length(), "_out.avi");
	//nameout="/home/paolo/Documents/OpenCV-soft/ActinWave&TractionField/out.avi";
	Size WriOut_size=Size(imagewidth, imagelength);
	double WriOut_FPS=5;
	//int WriOut_fourCC = -1;
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
	#pragma omp parallel for ordered schedule(dynamic)
	for (pos=0; pos<ntot; pos++)
		{
		#pragma omp ordered
			{
			#pragma omp critical (write)
				{
				Mat temp_outframe = outframes[pos]; //, temp_outframe_color;
				//normalize(outframes[pos], temp_outframe, 0, 255, NORM_MINMAX, -1);
				//temp_outframe.convertTo(temp_outframe, CV_8UC1);
				//imshow("pre saved", temp_outframe);
				//cvtColor(temp_outframe, temp_outframe, CV_GRAY2RGB);
				//WriOut << temp_outframe_color;
				WriOut << temp_outframe;
				//imshow("saved", temp_outframe_color);
				//int keycode = waitKey(50);
				}
			}
		}
}

/////////////////// Save th output color name string!!
void ParallelSaveOuputColorNameStr(string nameout, vector<Mat>& outframes, int& ntot,  int& imagewidth, int& imagelength)
{
	/*nameout=string(argv[1]);
	string stravi=".tif";
	nameout.replace(nameout.find(stravi), stravi.length(), "_out.avi");*/
	//nameout="/home/paolo/Documents/OpenCV-soft/ActinWave&TractionField/out.avi";
	Size WriOut_size=Size(imagewidth, imagelength);
	double WriOut_FPS=5;
	//int WriOut_fourCC = -1;
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
	#pragma omp parallel for ordered schedule(dynamic)
	for (pos=0; pos<ntot; pos++)
		{
		#pragma omp ordered
			{
			#pragma omp critical (write)
				{
				Mat temp_outframe = outframes[pos]; //, temp_outframe_color;
				//normalize(outframes[pos], temp_outframe, 0, 255, NORM_MINMAX, -1);
				//temp_outframe.convertTo(temp_outframe, CV_8UC1);
				//imshow("pre saved", temp_outframe);
				//cvtColor(temp_outframe, temp_outframe, CV_GRAY2RGB);
				//WriOut << temp_outframe_color;
				WriOut << temp_outframe;
				//imshow("saved", temp_outframe_color);
				//int keycode = waitKey(50);
				}
			}
		}
}


/////////////////// Save th output color name string!!
void ParallelSaveOuputNameStr(string nameout, vector<Mat>& outframes, int& ntot,  int& imagewidth, int& imagelength)
{
	/*nameout=string(argv[1]);
	string stravi=".tif";
	nameout.replace(nameout.find(stravi), stravi.length(), "_out.avi");*/
	//nameout="/home/paolo/Documents/OpenCV-soft/ActinWave&TractionField/out.avi";
	Size WriOut_size=Size(imagewidth, imagelength);
	double WriOut_FPS=5;
	//int WriOut_fourCC = -1;
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
	#pragma omp parallel for ordered schedule(dynamic)
	for (pos=0; pos<ntot; pos++)
		{
		#pragma omp ordered
			{
			#pragma omp critical (write)
				{
				Mat temp_outframe = outframes[pos]; //, temp_outframe_color;
				normalize(outframes[pos], temp_outframe, 0, 255, NORM_MINMAX, -1);
				temp_outframe.convertTo(temp_outframe, CV_8UC1);
				//imshow("pre saved", temp_outframe);
				cvtColor(temp_outframe, temp_outframe, CV_GRAY2RGB);
				//WriOut << temp_outframe_color;
				WriOut << temp_outframe;
				//imshow("saved", temp_outframe_color);
				//int keycode = waitKey(50);
				}
			}
		}
}

/////////////////// Save th output color name string!!
void SaveOuputNameStr(string nameout, vector<Mat>& outframes, int& ntot,  int& imagewidth, int& imagelength)
{
	/*nameout=string(argv[1]);
	string stravi=".tif";
	nameout.replace(nameout.find(stravi), stravi.length(), "_out.avi");*/
	//nameout="/home/paolo/Documents/OpenCV-soft/ActinWave&TractionField/out.avi";
	Size WriOut_size=Size(imagewidth, imagelength);
	double WriOut_FPS=5;
	//int WriOut_fourCC = -1;
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
		Mat temp_outframe = outframes[pos];
		Mat temp_outframe_color;
		normalize(outframes[pos], temp_outframe, 0, 255, NORM_MINMAX, -1);
		temp_outframe.convertTo(temp_outframe, CV_8UC1);
		//imshow("pre saved", temp_outframe);
		cvtColor(temp_outframe, temp_outframe_color, CV_GRAY2RGB);
		WriOut << temp_outframe_color;
		//imshow("saved", temp_outframe_color);
		//int keycode = waitKey(50);
		}
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
	//int WriOut_fourCC = -1;
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
