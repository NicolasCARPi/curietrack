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
#include <dirent.h>
//#include </usr/include/tiffio.h>
//#include </usr/include/tiff.h>
#include <tiffio.h>
#include <tiff.h>
#include "../Read_Write_Filter_Paolo/FileListAlphaNumeric.h"

using namespace cv;
using namespace std;

//////////////////////////////////
typedef Vec<unsigned short, 3> Vec3us;
//////////////////// read tif (8-16bit color or grayscale), convert it to 16 bit grayscale and put in the vector
void ReadTIFF(char* argv, int& ntot, int& imagewidth, int& imagelength, vector<Mat>& input_frames)
{
    //cout << TIFFLIB_VERSION << endl;
    TIFF* tif = TIFFOpen(argv, "r");
    if (tif)
    {
        int dircount = 0;
        do
        {
            dircount++;
            //uint32 imagelength, imagewidth;
            uint16 config, spp, bps;
            unsigned int sf, comp;

            TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &imagelength);
            TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &imagewidth);
            TIFFGetField(tif, TIFFTAG_PLANARCONFIG, &config);
            TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &spp);
            TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bps);
            TIFFGetField(tif, TIFFTAG_SAMPLEFORMAT, &sf);
            TIFFGetField(tif, TIFFTAG_COMPRESSION, &comp);

            //printf("config:%u\t spp:%u\t bps:%u\t sf:%u\t comp:%u\n", config, spp, bps, sf, comp);
            if (config == PLANARCONFIG_SEPARATE)
            {
                printf("Not supported format: PLANARCONFIG_SEPARATE");
                exit(0);
            }
            Mat frame;
            if(spp == 1 && bps == 16)  //// 16 bit gray
            {
                //printf("16 bit gray\n");
                unsigned short* LineBuf;
                LineBuf = (unsigned short*)_TIFFmalloc(TIFFScanlineSize(tif));
                frame = Mat(imagelength, imagewidth, CV_16UC1);
                for (int row = 0; row < imagelength; row++)
                {
                    TIFFReadScanline(tif, LineBuf, row);
                    for (int col = 0; col < imagewidth; col++)
                    {
                        frame.at<unsigned short>(row,col)=LineBuf[col];
                    }
                }
                //normalize(frame, frame, 0, 65535, NORM_MINMAX, -1);
                _TIFFfree(LineBuf);
            }
            else if(spp ==1 && bps == 8)  ///// 8 bit gray
            {
                //printf("8 bit gray\n");
                unsigned char* LineBuf;
                LineBuf = (unsigned char*)_TIFFmalloc(TIFFScanlineSize(tif));
                frame = Mat(imagelength, imagewidth, CV_8UC1);
                for (int row = 0; row < imagelength; row++)
                {
                    TIFFReadScanline(tif, LineBuf, row);
                    for (int col = 0; col < imagewidth; col++)
                    {
                        frame.at<unsigned char>(row,col)=LineBuf[col];
                    }
                }
                //normalize(frame, frame, 0, 255, NORM_MINMAX, -1);
                frame.convertTo(frame, CV_16UC1);
                //normalize(frame, frame, 0, 65535, NORM_MINMAX, -1);
                _TIFFfree(LineBuf);
            }
            else if(spp ==3 && bps == 8)  ///// 8 bit color
            {
                //printf("8 bit color\n");
                unsigned char* LineBuf;
                LineBuf = (unsigned char *)_TIFFmalloc(TIFFScanlineSize(tif));
                Mat frameColor = Mat(imagelength, imagewidth, CV_8UC3);
                frame = Mat(imagelength, imagewidth, CV_8UC3);
                for (int row = 0; row < imagelength; row++)
                {
                    TIFFReadScanline(tif, LineBuf, row, 0);
                    int col, tmpcolorcol=0;
                    for (col = 0; col < imagewidth; col++)
                    {
                        tmpcolorcol=3*col;
                        frameColor.at<Vec3b>(row,col)[2]=LineBuf[tmpcolorcol];
                        frameColor.at<Vec3b>(row,col)[1]=LineBuf[tmpcolorcol+1];
                        frameColor.at<Vec3b>(row,col)[0]=LineBuf[tmpcolorcol+2];
                    }
                }
                cvtColor(frameColor, frame, CV_RGB2GRAY);
                frame.convertTo(frame, CV_16UC1);
                //normalize(frame, frame, 0, 65535, NORM_MINMAX, -1);
                _TIFFfree(LineBuf);
                frameColor.release();
            }
            else if(spp ==3 && bps == 16)  ///// 16 bit color
            {
                //printf("16 bit color\n");
                unsigned short* LineBuf;
                LineBuf = (unsigned short *)_TIFFmalloc(TIFFScanlineSize(tif));
                Mat frameColor = Mat(imagelength, imagewidth, CV_16UC3);
                frame = Mat(imagelength, imagewidth, CV_16UC1);
                for (int row = 0; row < imagelength; row++)
                {
                    TIFFReadScanline(tif, LineBuf, row, 0);
                    int col, tmpcolorcol=0;
                    for (col = 0; col < imagewidth; col++)
                    {
                        tmpcolorcol=3*col;
                        frameColor.at<Vec3us>(row,col)[2]=LineBuf[tmpcolorcol];
                        frameColor.at<Vec3us>(row,col)[1]=LineBuf[tmpcolorcol+1];
                        frameColor.at<Vec3us>(row,col)[0]=LineBuf[tmpcolorcol+2];
                    }
                }
                cvtColor(frameColor, frame, CV_RGB2GRAY);
                //normalize(frame, frame, 0, 65535, NORM_MINMAX, -1);
                _TIFFfree(LineBuf);
                frameColor.release();
            }
            else
            {
                printf("Not supported format\n");
                exit (0);
            }
            //imshow("open", frame);
            //waitKey();
            input_frames.push_back(frame);
            frame.release();
        }
        while (TIFFReadDirectory(tif));
        ntot=dircount;
    }
    TIFFClose(tif);
}

void ImportTIFFdir(char* argv, int& ntot, int& imagewidth, int& imagelength, vector<Mat>& input_frames)
{
    string dirname=argv;
    if ( dirname.compare(dirname.size()-1, 1, "/") != 0 ) dirname.append("/");
    //printf("dirname: %s\n", dirname.c_str());
    vector<string> filelist;
    listFileAphaNumeric(dirname, filelist);
    ntot = filelist.size();
    printf("ntot: %i\n", ntot);

    for (int i=0; i<ntot; i++)
    {
        int temp_n;
        int temp_imagewidth, temp_imagelength;
        Mat temp_input_frame;
        vector< Mat > temp_input_frames;
        string filenamestr;
        filenamestr.append(dirname);
        filenamestr.append(filelist[i].c_str());
        char *namefile = new char[strlen(filenamestr.c_str())];
        strcpy(namefile, filenamestr.c_str());
        //printf("#%i: %s\n",i,namefile);
        //cout << namefile << endl;


        ReadTIFF(namefile, temp_n, temp_imagewidth, temp_imagelength, temp_input_frames);
        temp_input_frames[0].copyTo(temp_input_frame);
        /*temp_input_frame= imread(namefile, 0);
        temp_input_frame.convertTo(temp_input_frame, CV_16UC1);*/
        if (i==0)
        {
            /*imagewidth=temp_imagewidth;
            imagelength=temp_imagelength;*/
            imagewidth=temp_input_frame.cols;
            imagelength=temp_input_frame.rows;
        }
        else
        {
            //if (imagewidth != temp_imagewidth || imagelength != temp_imagelength)
            if (imagewidth != temp_input_frame.cols || imagelength != temp_input_frame.rows)
            {
                printf("image width or length different!\n aborted\n");
                break;
            }
        }
        //cout << temp_n << endl;
        //for (int j=0; j<temp_input_frames.size(); j++)
        //Mat frame = temp_input_frames[0];
        //input_frames.push_back(frame);
        input_frames.push_back(temp_input_frame);
        printf("#%i: %s\n",i,namefile);
    }

}

