#include <opencv2/opencv.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/core/mat.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc_c.h>
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

void ThrGausDecay8bit2(const Mat& hist, int& max, int& sigma3)
{
    double minVal=0, maxVal=0;
    Point minLoc, maxLoc;
    minMaxLoc(hist, &minVal, &maxVal, &minLoc, &maxLoc);
    int di=0;
    while (hist.at<unsigned char>(0, maxLoc.x+di) > 0.5*maxVal)
    {
        di += 1;
    }
    if (di != 0)
        sigma3 = cvRound( 3*di/sqrt(2*log(2)) );
    else
        sigma3 =0;
    max=round(maxLoc.x);
    printf("\t\t8bit\t Max: %i\t 3sigma_gaussian: %i\n", max, sigma3);
}
void ThrGausDecay16bit2(const Mat& hist, int& max, int& sigma3)
{
    double minVal=0, maxVal=0;
    Point minLoc, maxLoc;
    minMaxLoc(hist, &minVal, &maxVal, &minLoc, &maxLoc);
    int di=0;
    while (hist.at<unsigned short>(0, maxLoc.x+di) > 0.5*maxVal)
    {
        di += 1;
    }
    if (di != 0)
        sigma3 = cvRound( 3*di/sqrt(2*log(2)) );
    else
        sigma3 =0;
    max=round(maxLoc.x);
    printf("\t\t16bit\t Max: %i\t 3sigma_gaussian: %i\n", max, sigma3);
}
/// build cumulative histogram + Integral in parallel and Identify max and sigma GAUSSIAN in a frame from hist e bin
void CumHistogramInt_ThrGausDecay8bit(vector<Mat>& inframes, int& ntot, Mat& CumHist, Mat& CumHistInt, int& max, int& sigma3)
{
    int bin=256;
    int pos;
    CumHist = Mat::zeros(1, bin, CV_32F);
    CumHistInt = Mat::zeros(1, bin, CV_32F);
#pragma omp parallel for ordered schedule(dynamic)
    for (pos=0; pos<ntot; pos++)
    {
        Mat frame;
#pragma omp critical (readcapHist)
        frame = inframes[pos];
        Mat hist=Mat::zeros(1, bin, CV_32F);
        for (int h=0; h<frame.rows; h++)
        {
            for (int k=0; k<frame.cols; k++)
            {
                int pixvalue = frame.at<unsigned char>(h,k);
                hist.at<float>(0, pixvalue) += 1;
            }
        }
#pragma omp critical (updateCumHist)
        {
            CumHist += hist;
            blur(CumHist, CumHist, Size(3,1));
        }
        frame.release();
    }

    /*for (int j=0; j<bin; j++)
        for (int z=0; z<j+1; z++)
            CumHistInt.at<float>(0,j) += CumHist.at<float>(0,z);*/
    #pragma omp parallel for schedule(dynamic)
        for (int j=0; j<bin; j++)
            for (int z=0; z<j+1; z++)
                #pragma omp atomic
                    CumHistInt.at<float>(0,j) += CumHist.at<float>(0,z);

/*#pragma omp parallel for num_threads(4) schedule(dynamic)
    for (int j=0; j<bin; j++)
    {
        Mat CumHistCrop;
        resize(CumHist, CumHistCrop, Size(j+1,1), 1, 1);
        CumHistInt.at<float>(0,j) = sum(CumHistCrop)[0];
    }*/

    normalize(CumHist, CumHist, 0, 255, NORM_MINMAX);
    CumHist.convertTo(CumHist, CV_8UC1);
    normalize(CumHistInt, CumHistInt, 0, 255, NORM_MINMAX);
    CumHistInt.convertTo(CumHistInt, CV_8UC1);

    /*double minVal=0, maxVal=0;
    Point minLoc, maxLoc;*/
    // To plot Histogram and integrated histogram
    /*
    Mat CumHistImg = Mat::zeros(600, bin, CV_8UC3);
    minMaxLoc(CumHist, &minVal, &maxVal, &minLoc, &maxLoc);
    for (int posbin=0; posbin<bin; posbin++)
    	rectangle(CumHistImg, Point(posbin+1, floor(600*CumHist.at<unsigned char>(0,posbin)/maxVal)), Point(posbin, 0), Scalar(100, 100, 100));
    	flip(CumHistImg, CumHistImg, 0);
    imshow("CumHist8bit", CumHistImg);
    waitKey(500);
    minMaxLoc(CumHistInt, &minVal, &maxVal, &minLoc, &maxLoc);
    CumHistImg = Mat::zeros(600, bin, CV_8UC3);
    for (int posbin=0; posbin<bin; posbin++)
    	rectangle(CumHistImg, Point(posbin+1, floor(600*CumHistInt.at<unsigned char>(0,posbin)/maxVal)), Point(posbin, 0), Scalar(100, 100, 100));
    	flip(CumHistImg, CumHistImg, 0);
    imshow("CumHist8bit", CumHistImg);
    waitKey(500);
    destroyAllWindows();
    */
    ThrGausDecay8bit2 (CumHist, max, sigma3);
}
/// bild cumulative histogram + Integral in parallel and Identify max and sigma GAUSSIAN in a frame from hist e bin 16bit
void CumHistogramInt_ThrGausDecay16bit(vector<Mat>& inframes, int& ntot, Mat& CumHist, Mat& CumHistInt, int& max, int& sigma3)
{
    int bin=65536;
    int pos;
    CumHist = Mat::zeros(1, bin, CV_32F);
    CumHistInt = Mat::zeros(1, bin, CV_32F);
#pragma omp parallel for ordered schedule(dynamic) num_threads(4)
    for (pos=0; pos<ntot; pos++)
    {
        Mat frame;
#pragma omp critical
        inframes[pos].copyTo(frame); //frame = inframes[pos];
        Mat hist=Mat::zeros(1, bin, CV_32F);
        for (int h=0; h<frame.rows; h++)
        {
            for (int k=0; k<frame.cols; k++)
            {
                int pixvalue = frame.at<unsigned short>(h,k);
                hist.at<float>(0, pixvalue) += 1;
            }
        }
#pragma omp critical
        {
            CumHist += hist;
            blur(CumHist, CumHist, Size(3,1));
        }
        frame.release();
    }

    printf("\tsingle element access start\n");

    #pragma omp parallel for schedule(dynamic)
        for (int j=0; j<bin; j++)
            for (int z=0; z<j+1; z++)
                #pragma omp atomic
                    CumHistInt.at<float>(0,j) += CumHist.at<float>(0,z);

// very fast!
/*
#pragma omp parallel for num_threads(5) schedule(dynamic)
    for (int j=0; j<bin; j++)
    {
        Mat CumHistCrop;
        resize(CumHist, CumHistCrop, Size(j+1,1), 1, 1);
        CumHistInt.at<float>(0,j) = sum(CumHistCrop)[0];
    }
*/
// on the GPU
    /*#pragma omp parallel for num_threads(4) schedule(dynamic)
    for (int j=0; j<bin; j++)
        {
            Mat CumHistCrop;
            GpuMat CumHistCrop_GPU;
            //printf("Resizing frame: %i\n",j);
            resize(CumHist, CumHistCrop, Size(j+1,1), 1, 1);
            CumHistCrop_GPU.upload(CumHistCrop);
            //printf("Frame: %i resized\n",j);
            CumHistInt.at<float>(0,j) = gpu::sum(CumHistCrop_GPU)[0];
        }*/

    printf("\tsingle element access end\n");

    normalize(CumHist, CumHist, 0, 65535, NORM_MINMAX);
    CumHist.convertTo(CumHist, CV_16UC1);
    normalize(CumHistInt, CumHistInt, 0, 65535, NORM_MINMAX);
    CumHistInt.convertTo(CumHistInt, CV_16UC1);
    ThrGausDecay16bit2 (CumHist, max, sigma3);
}

//// Parallel Filter ant convert to 8bit but not cumulative hist and equalize hist
void Scale_Convert8bit(vector<Mat> inframes, int ntot, vector<Mat>& outframes, int& MinStack, int& MaxStack)
{
    printf("\tScale and convert to 8bit\n");
    //printf("\t Start\n");

    double Max[ntot], Min[ntot];
    double MaxValue, MinValue;
    int pos=0;
#pragma omp parallel for ordered schedule(dynamic)
    for (pos=0; pos<ntot; pos++)
    {
        Mat frame0;
        double minval, maxval;
#pragma omp critical (readcap)
        //frame0=inframes[pos].clone();
        inframes[pos].copyTo(frame0);
        minMaxLoc(frame0, &minval, &maxval);
#pragma omp critical (minmax)
        {
            Min[pos] =  minval;
            Max[pos] =  maxval;
            frame0.release();
        }
    }
    /*MaxValue = 0.5*(gsl_stats_mean(Max, 1, ntot)+gsl_stats_max(Max, 1, ntot));
    MinValue = 0.5*(gsl_stats_mean (Min, 1, ntot)+gsl_stats_min (Min, 1, ntot));*/
    MaxValue = gsl_stats_mean(Max, 1, ntot);
    MinValue = gsl_stats_mean (Min, 1, ntot);
    MinStack = int(round(MinValue));
    MaxStack = int(round(MaxValue));

#pragma omp parallel for ordered schedule(dynamic)
    for (pos=0; pos<ntot; pos++)
    {
        Mat frame0;
#pragma omp critical (readtemp)
        {
            //frame0=inframes[pos].clone();
            inframes[pos].copyTo(frame0);
        }
        Mat frame1 = Mat::zeros(frame0.rows, frame0.cols, CV_8UC1);
        for (int i=0; i<frame0.rows; i++)
            for (int j=0; j<frame0.cols; j++)
            {
                int value;
                value=round((255/(MaxValue-MinValue))*frame0.at<unsigned short>(i,j)-255*MinValue/(MaxValue-MinValue));
                /*if ((255/(MaxValue-MinValue))*(frame0.at<unsigned short>(i,j)-MinValue) > 255)
                    frame1.at<unsigned char>(i,j) = 255;
                else
                    frame1.at<unsigned char>(i,j) = (255/(MaxValue-MinValue))*(frame0.at<unsigned short>(i,j)-MinValue);*/
                if (value > 255)
                    frame1.at<unsigned char>(i,j) = 255;
                else
                    frame1.at<unsigned char>(i,j) = value;
            }
#pragma omp critical (tempoutframespos)
        {
            //outframes[pos] = frame1.clone();
            frame1.copyTo(outframes[pos]);
            /*imshow("Original", frame0);
            imshow("8bit Scaled", frame1);
            waitKey(1555);*/
            frame1.release();
            frame0.release();
        }
    }
    //printf("\t End\n");
}
void Enhance16bit(vector<Mat> inframes, int ntot, vector<Mat>& outframes, int max, int sigma3)
{
    int pos=0;
#pragma omp parallel for ordered schedule(dynamic)
    for (pos=0; pos<ntot; pos++)
    {
        Mat frame0, frame1;
#pragma omp critical (readtemp)
        {
            //frame0=inframes[pos];
            //frame0=inframes[pos].clone();
            inframes[pos].copyTo(frame0);
        }
        /*Mat frame1 = Mat::zeros(frame0.rows, frame0.cols, CV_16UC1);
        for (int i=0; i<frame0.cols; i++)
            for (int j=0; j<frame0.rows; j++)
                frame1.at<unsigned short>(i,j) = frame0.at<unsigned short>(i,j)-max;*/
        frame1 = frame0 - max;
#pragma omp critical (tempoutframespos)
        {
            //outframes[pos] = frame1.clone();
            frame1.copyTo(outframes[pos]);
            frame0.release();
            frame1.release();
        }
    }
}
//// only filter
//// Filter to segment
void FilterSegmSingleFrame(Mat image, Mat& dest16, Mat& filter, int tamagno, Mat& MorphEllipseTamagno)
{
    //The last!!
    GaussianBlur(image, filter, Size(tamagno, tamagno), 0, 0, BORDER_DEFAULT);
    Mat image16, filter16;
    filter.convertTo(filter16, CV_16UC1);
    morphologyEx(filter, filter, MORPH_OPEN, MorphEllipseTamagno, Point(-1, -1), 1);
    medianBlur(filter, filter, tamagno);
    image.convertTo(image16, CV_16UC1);
    multiply(image16, filter16, dest16);
    image16.release();
    filter16.release();
}
void FilterSegmSingleFrame1(Mat image, Mat& filter, int tamagno, Mat& MorphEllipseTamagno)
{
    //The last ? ? !!
    GaussianBlur(image, filter, Size(tamagno, tamagno), 0, 0, BORDER_DEFAULT);
    morphologyEx(filter, filter, MORPH_OPEN, MorphEllipseTamagno, Point(-1, -1), 1);
    medianBlur(filter, filter, tamagno);
}
////////////////////////////////////////
void Filter_8bit(vector<Mat> inframes, int ntot, vector<Mat>& outframes, int tamagno)
{
    Mat MorphEllipseTamagno = getStructuringElement(MORPH_ELLIPSE, Size(tamagno,tamagno));
    int pos=0;
    vector<Mat> temp_outl(ntot);
    vector<Mat> temp_filter(ntot);
    vector<Mat> temp_outl_scaled(ntot);
#pragma omp parallel for ordered schedule(dynamic) shared(temp_outl, temp_filter)
    for (pos=0; pos<ntot; pos++)
    {
        Mat frame0, frame1, filter;
#pragma omp critical (readtemp)
        {
            //frame0=inframes[pos];
            inframes[pos].copyTo(frame0);
        }
        FilterSegmSingleFrame(frame0, frame1, filter, tamagno, MorphEllipseTamagno);
        //printf("8bit conversion\n");
#pragma omp critical (tempoutframespos)
        {
            //temp_outl[pos] = frame1;
            //temp_filter[pos] = filter;
            imshow("ter", frame1);
            waitKey();
            frame1.copyTo(temp_outl[pos]);
            filter.copyTo(temp_filter[pos]);
        }
        frame0.release();
        frame1.release();
        filter.release();
    }
    MorphEllipseTamagno.release();
    cout << "End I part" << endl;
    int tempMAx, tempMin;
    Scale_Convert8bit(temp_outl, ntot, temp_outl_scaled, tempMAx, tempMin);
    temp_outl.clear();
#pragma omp parallel for ordered schedule(dynamic) shared(inframes, temp_outl_scaled, temp_filter, outframes)
    for (pos=0; pos<ntot; pos++)
    {
        Mat image, frame0, frame1, filter;
#pragma omp critical (readtemp)
        {
            //image = inframes[pos];
            //frame0=temp_outl[pos];
            //filter=temp_filter[pos];
            inframes[pos].copyTo(image);
            temp_outl_scaled[pos].copyTo(frame0);
            temp_filter[pos].copyTo(filter);
        }
        subtract(frame0, filter, frame1);
        //frame1 = frame0 - filter;
#pragma omp critical (tempoutframespos)
        {
            //outframes[pos] = frame1;
            frame1.copyTo(outframes[pos]);
            /*imshow("image", image);
            imshow("filter", filter);
            imshow("FilterTamagno", frame1);
            waitKey(555);*/
        }
        image.release();
        frame1.release();
        frame0.release();
        filter.release();
    }
    temp_outl_scaled.clear();
    temp_filter.clear();
    printf("\tSubstract Background\n");
}
void Filter_8bit1(vector<Mat> inframes, int ntot, vector<Mat>& outframes, int tamagno)
{
    Mat MorphEllipseTamagno = getStructuringElement(MORPH_ELLIPSE, Size(tamagno,tamagno));
    int pos=0;
    vector<Mat> temp_filter(ntot);
#pragma omp parallel for ordered schedule(dynamic) shared(temp_filter)
    for (pos=0; pos<ntot; pos++)
    {
        Mat frame0, filter;
#pragma omp critical (readtemp)
        {
            //frame0=inframes[pos];
            inframes[pos].copyTo(frame0);
        }
        FilterSegmSingleFrame1(frame0, filter, tamagno, MorphEllipseTamagno);
        //printf("8bit conversion\n");
#pragma omp critical (tempoutframespos)
        {
            //temp_outl[pos] = frame1;
            //temp_filter[pos] = filter;
            filter.copyTo(temp_filter[pos]);
        }
        frame0.release();
        filter.release();
    }
    MorphEllipseTamagno.release();
    cout << "End I part" << endl;
#pragma omp parallel for ordered schedule(dynamic) shared(inframes, temp_filter, outframes)
    for (pos=0; pos<ntot; pos++)
    {
        Mat image, frame1, filter;
#pragma omp critical (readtemp)
        {
            //image = inframes[pos];
            //frame0=temp_outl[pos];
            //filter=temp_filter[pos];
            inframes[pos].copyTo(image);
            temp_filter[pos].copyTo(filter);
        }
        subtract(image, filter, frame1);
        //frame1 = frame0 - filter;
#pragma omp critical (tempoutframespos)
        {
            //outframes[pos] = frame1;
            frame1.copyTo(outframes[pos]);
            /*imshow("image", image);
            imshow("filter", filter);
            imshow("FilterTamagno", frame1);
            waitKey(555);*/
        }
        image.release();
        frame1.release();
        filter.release();
    }
    temp_filter.clear();
    printf("\tSubstract Background\n");
}
////////////////////////////////////////
////////////////////////////////////////
int FreeVectorMat( vector<Mat>& V)
{
    vector<Mat> tmp;
    V.swap( tmp );
    V.clear();
    return(0);
}
////////////////////////////////////////
/// only enhance
void Enhance8bit(vector<Mat>& inframes, int ntot, vector<Mat>& outframes, Mat& CumHistInt, int tamagno, double SubBackFactor)
{
    int pos=0;
    //vector <Mat> temp_output(ntot);
#pragma omp parallel for ordered schedule(dynamic)
    for (pos=0; pos<ntot; pos++)
    {
        Mat frame0, frame1;
#pragma omp critical (readtemp)
        {
            //frame0=inframes[pos];
            inframes[pos].copyTo(frame0);
            //imshow("8bit", frame0);
            //waitKey();
        }
        LUT(frame0, CumHistInt, frame1);
#pragma omp critical (tempoutframespos)
        {
            //temp_output[pos] = frame1;
            frame1.copyTo(outframes[pos]);
            frame1.release();
            frame0.release();
            /*imshow("8bit", frame0);
            imshow("Enhanced", frame1);
            waitKey(333);*/
        }
    }
    Mat CumHistInt_temp, CumHist_temp;
    int max_temp, sigma3_temp;
    CumHistogramInt_ThrGausDecay8bit(outframes, ntot, CumHist_temp, CumHistInt_temp, max_temp, sigma3_temp);
#pragma omp parallel for ordered schedule(dynamic)
    for (pos=0; pos<ntot; pos++)
    {
        Mat image, frame0, frame1;
#pragma omp critical (readtemp1)
        {
            //image = inframes[pos];
            //frame0=temp_output[pos];
            inframes[pos].copyTo(image);
            outframes[pos].copyTo(frame0);
        }
        frame1 = frame0 - SubBackFactor*(max_temp+sigma3_temp);
        //GaussianBlur(frame1, frame1, Size(3, 3), 0, 0, BORDER_DEFAULT);
        //GaussianBlur(frame1, frame1, Size(tamagno, tamagno), 0, 0, BORDER_DEFAULT);
        //medianBlur(frame1, frame1, tamagno);
#pragma omp critical (tempoutframespos1)
        {
            //temp_output[pos] = frame1;
            frame1.copyTo(outframes[pos]);
            /*imshow("Image", image);
            imshow("Enhanced-Background", frame1);
            waitKey(333);*/
            image.release();
            frame1.release();
            frame0.release();
        }
    }
    //Filter_8bit(temp_output, ntot, outframes, tamagno);
    //outframes = temp_output;

    printf("\tEnhanced\n");

}
////////////// background substraction
/*void BackGroundSubOpenCV (vector<Mat>& inframes, vector<Mat>& outframes)
{
    int ntot = inframes.size();
    BackgroundSubtractorMOG2 bg_model;
    int pos=0;
#pragma omp parallel for ordered schedule(dynamic)
    for (pos=0; pos<ntot; pos++)
    {
        Mat frame0, frame1, mask;
#pragma omp critical (readtemp)
        frame0=inframes[pos];
        bg_model(frame0, mask, -1);
        frame0.copyTo(frame1, mask);
#pragma omp critical (writetemp)
        {
            outframes[pos]=frame1;
            imshow("frame0", frame0);
            imshow("frame1", frame1);
            waitKey(555);
        }
    }
}*/
/////////////////
void GenerateColorOut(vector<Mat>& inframes, vector<Mat>& outcolor)
{
    int ntot = inframes.size();
    int pos=0;
#pragma omp parallel for ordered schedule(dynamic)
    for (pos=0; pos<ntot; pos++)
    {
        Mat frame0, frame1;
#pragma omp critical (readtemp)
        //frame0=inframes[pos];
        inframes[pos].copyTo(frame0);
        //normalize(frame0, frame1, 0, 255, NORM_MINMAX);
        //frame1.convertTo(frame1, CV_8UC1);
        cvtColor(frame0, frame1, CV_GRAY2RGB);
#pragma omp critical (writetemp)
        {
            //outcolor[pos]=frame1;
            frame1.copyTo(outcolor[pos]);
            /*imshow("frame0", frame0);
            imshow("frame1", frame1);
            waitKey(555); */
        }
        frame0.release();
        frame1.release();
    }
}

void GenerateColorOutMixChannels(vector<Mat>& inframes, vector<Mat>& outcolor)
{
    int from_to[] = {0,0};
    int ntot = inframes.size();
    int pos=0;
#pragma omp parallel for ordered schedule(dynamic)
    for (pos=0; pos<ntot; pos++)
    {
        Mat frame0, frame1;
#pragma omp critical (readtemp)
        inframes[pos].copyTo(frame0);
        //cvtColor(frame0, frame1, CV_GRAY2RGB);
        frame1 = Mat::zeros(frame0.rows, frame0.cols, CV_8UC3);
        mixChannels(&frame0, 1, &frame1, 1, from_to, 1);
#pragma omp critical (writetemp)
        {
            /*imshow("frame0", frame0);
            imshow("frame1", frame1);
            waitKey(555);*/
            frame1.copyTo(outcolor[pos]);
        }
        frame0.release();
        frame1.release();
    }
}
