#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_math.h>
#include <fstream>
#include <iostream>
#include <sstream>
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
#include <cstring>
#include <memory>
#include <new>
#include "../TrackingLib_Perrine_2.1/kernel/TrackingLib.h"

using namespace cv;
using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct Roi
{
    int origID; // originale identifier
    int frame; // frame
    Point2f center;  // Roi center
    double iI; // Integrated Intensity
    double peri; // Perimeter
    double area; // Area
    double convexity; // Area
    vector< Point> contour; // contour
    int trajtag;
    int mothertag;
    int neigbors;
};
double RoiDistance(Roi A, Roi B)
{
    double d=0;
    d = sqrt( pow(A.center.x-B.center.x,2) + pow(A.center.y-B.center.y,2) );
    return d;
}
double RoiDistanceX(Roi A, Roi B)
{
    double d=0;
    d = sqrt( pow(A.center.x-B.center.x,2));
    return d;
}
double RoiDistanceY(Roi A, Roi B)
{
    double d=0;
    d = sqrt(pow(A.center.y-B.center.y,2) );
    return d;
}
int RoiFrameDistance(Roi A, Roi B)
{
    int df=0;
    df = abs(A.frame-B.frame);
    return df;
}
int RoiiIDifference(Roi A, Roi B)
{
    int diI=0;
    diI = sqrt(pow(A.iI-B.iI,2));
    return diI;
}
int RoiPeriDifference(Roi A, Roi B)
{
    int d=0;
    d = sqrt(pow(A.peri-B.peri,2));
    return d;
}
int RoiAreaDifference(Roi A, Roi B)
{
    int d=0;
    d = sqrt(pow(A.area-B.area,2));
    return d;
}
int FreeVectorRoi( vector<Roi>& V)
{
    vector<Roi> tmp;
    V.swap( tmp );
    V.clear();
    return(0);
}
int FreeVectorVectorRoi( vector< vector<Roi> >& V)
{
    vector< vector<Roi> > tmp;
    V.swap( tmp );
    V.clear();
    return(0);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Define dimension of interest objects in a single frame
int DefiSizeObjsErodeSingleFrame2(const Mat& img, int& threshod)
{
    Mat tmp_img;
    threshold(img, tmp_img, threshod, 255, THRESH_BINARY);
    Scalar thr_mean;
    thr_mean = mean(tmp_img);
    //printf("thr_mean: %f\n", thr_mean[0]);
    int n_erode=0;
    while (thr_mean[0] != 0) // erode while the mean intensity is not 0
    {
        erode(tmp_img, tmp_img, getStructuringElement(MORPH_ELLIPSE, Size(3,3)), Point(-1,-1), 1);
        thr_mean = mean(tmp_img);
        //printf("\tthr_mean: %f\n", thr_mean[0]);
        n_erode += 1;
    }
    if (n_erode == 0)
        n_erode = 1;
    return n_erode;
}
// Define dimension of interest objects in all frames
int DefiSizeObjsErodeAll2(vector<Mat>& input_frames, int& ntot, int& max, int& sigma3)
{
    int Nerode=0, erode2=0, erodew=0;
    int expthreshold=0;
    expthreshold = max+0.33*sigma3;
    //printf("expthreshold: %i\n", expthreshold);
    if (expthreshold == 0)
    {
        Nerode=9;
        return Nerode;
    }
    int pos;
    double nerodecount[ntot];
    //printf("Dimension?\n");
#pragma omp parallel for ordered schedule(dynamic)
    for (pos=0; pos<ntot; pos++)
    {
        Mat frame;
#pragma omp critical (erode)
        {
            frame=input_frames[pos];
        }
        int tmp_erode=0;
        tmp_erode=DefiSizeObjsErodeSingleFrame2(frame, expthreshold);
        /*#pragma omp atomic
        		erode2 += pow(tmp_erode,2);
        #pragma omp atomic
        		erodew += tmp_erode;*/
#pragma omp critical (erodeout)
        nerodecount[pos]=tmp_erode;
    }
    //printf("erode2: %i\t erodew: %i\n", erode2, erodew);
    //Nerode = round(2*erode2/erodew);
    Nerode = int(round(2*gsl_stats_mean(nerodecount, 1, ntot)));
    if (GSL_IS_EVEN(Nerode)==1)
        Nerode -= 1;
    if (Nerode<9)
        Nerode = 9;
    else if (Nerode>45)
        Nerode = 45;
    printf("\tPutative dimension (minimum 9 - maximum 45): %i\n", Nerode);
    return Nerode;
}
int DefiSizeObjsErodeAll3(vector<Mat>& input_frames, int& ntot, int& max, int& sigma3)
{
    double nerodecount[ntot];
    int pos=0;
#pragma omp parallel for ordered schedule(dynamic)
    for (pos=0; pos<ntot; pos++)
    {
        Mat tmp_frame, bgmask;
        input_frames[pos].copyTo(tmp_frame);
        adaptiveThreshold(tmp_frame, bgmask, 255, ADAPTIVE_THRESH_GAUSSIAN_C, THRESH_BINARY, 5, -max);
        int temperode=0;
        Scalar thr_mean;
        thr_mean = mean(bgmask);
        while (thr_mean[0] != 0) // erode while the mean intensity is not 0
        {
            //printf("mean: %f\n", thr_mean[0]);
            erode(bgmask, bgmask, getStructuringElement(MORPH_ELLIPSE, Size(3,3)), Point(-1,-1), 1);
            thr_mean = mean(bgmask);
            temperode += 1;
        }
        //printf("temp erode: %i\n", temperode);
        nerodecount[pos]= temperode;
    }
    int Nerode = int(round(9*gsl_stats_mean(nerodecount, 1, ntot)));
    if (GSL_IS_EVEN(Nerode)==1)
        Nerode -= 1;
    if (Nerode<9)
        Nerode = 9;
    else if (Nerode>31)
        Nerode = 31;
    printf("\tPutative dimension (minimum 9 - maximum 31): %i\n", Nerode);
    return Nerode;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Refine Segments
void refineSegments_Cells(const Mat& img, Mat& mask, int tamagno, int nFrame, vector<Roi>& RoisFramei, char* argv[], double MinArea, double MaxArea, FILE* MyObjectsFile)
{
    vector<vector<Point> > contours;
    vector<Vec4i> hierarchy;

    findContours( mask, contours, hierarchy, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_NONE);

    if( contours.size() == 0 )
        return;

    int idx = 0, largestComp = 0;
    for( ; idx >= 0; idx = hierarchy[idx][0] )
    {
        const vector<Point>& c = contours[idx];
        double area = fabs(contourArea(Mat(c)));
        double peri = fabs(arcLength(Mat(c), true));
        double circularity = pow(peri, 2)/(2*M_PI*area);
        vector<Point> convexHullpoints;  // Convex hull points
        vector<Point> convexHullcontour;  // Convex hull contour points
        double convexity=1;
        double epsilon = 0.001; // Contour approximation accuracy
        // Calculate convex hull of original points (which points positioned on the boundary)
        convexHull(Mat(c), convexHullpoints, false);
        // Approximating polygonal curve to convex hull
        approxPolyDP(Mat(convexHullpoints), convexHullcontour, 0.001, true);
        convexity = (fabs(contourArea(Mat(convexHullcontour)))-area)/area;
        if
        (
            area > MinArea &&
            area < MaxArea /*&&
            circularity < 7*/
        )
        {
            Roi temp_roi;
            Mat object_mask = Mat::zeros(img.size(), CV_8UC1);
            drawContours(object_mask, contours, idx, 255, CV_FILLED);
            //double peri = arcLength(Mat(c), true);
            Moments ObjectMoments = moments(Mat(c));
            Scalar MeanIntensityObject = mean(img, object_mask);
            double MeanIntensityObject_double = MeanIntensityObject[0];
            temp_roi.origID = idx;
            temp_roi.frame = nFrame;
            //temp_roi.center = CvPoint2D64f(ObjectMoments.m10/ObjectMoments.m00, ObjectMoments.m01/ObjectMoments.m00);
            temp_roi.center.x = ObjectMoments.m10/ObjectMoments.m00;
            temp_roi.center.y = ObjectMoments.m01/ObjectMoments.m00;
            temp_roi.iI = area*MeanIntensityObject_double;
            temp_roi.peri = peri; //arcLength(Mat(c), true);
            temp_roi.area = area;
            temp_roi.convexity=convexity;
            temp_roi.contour = c;
            temp_roi.mothertag = 0;
            RoisFramei.push_back(temp_roi);

            fprintf(MyObjectsFile, "%i\t%i\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n", nFrame, idx+1, ObjectMoments.m10/ObjectMoments.m00, ObjectMoments.m01/ObjectMoments.m00, area, peri, convexity, MeanIntensityObject_double); /// remember remember idx+1
        }
        /*else
            printf("Objectet size < min or > max\n");*/
    }
    fflush (MyObjectsFile);
}
// Parallel Segmentation cells
void ParallelSegmentCells(vector<Mat>& inframes, vector<Mat>& inframesMask, vector<Mat>& output_segmented, int& ntot, int& tamagno, int& MaxAll, int& StdevAll, vector< vector<Roi> >& RoisFrames, char* argv[], int nameposargv, int nclose, int nopen, double threshold_Mask, double threshold_Visual, int GeomFilterLineSizeX, int GeomFilterLineSizeY, int GeomFilterChanSizeX, int GeomFilterChanSizeY, int GeomFilter2DSizeX, int GeomFilter2DSizeY, int VisOut)
{
    FILE* MyObjectsFile;
    string nameout;
    nameout=string(argv[nameposargv]);
    string strtif=".tif";
    size_t foundtif;
    foundtif=nameout.rfind(strtif);
    if (foundtif != string::npos)
        nameout.replace(nameout.find(strtif), strtif.length(), "_Objects.txt");
    else
    {
        string newnameout = nameout;
        if ( newnameout.compare(newnameout.size()-1, 1, "/") == 0 ) newnameout.resize(newnameout.size()-1);
        string delim="/";
        size_t founddelim;
        founddelim=newnameout.rfind(delim);
        if (founddelim != string::npos)
            newnameout = newnameout.substr(founddelim);
        newnameout.append("_Objects.txt");
        if ( nameout.compare(nameout.size()-1, 1, "/") != 0 ) nameout.append("/");
        nameout.append(newnameout);
    }
    //cout << "File Objects: " << nameout.c_str() << endl;
    MyObjectsFile = fopen (nameout.c_str(), "w");
    fprintf( MyObjectsFile, "Frame\t#\tx\ty\tA\tP\t\C\tI\n");
    int pos=0;
    double minArea=0, maxArea=0;
    if (string(argv[4]) == "LINE")
    {
        minArea = M_PI*pow(0.05*tamagno,2);
        maxArea = M_PI*pow(0.5*tamagno,2);
    }
    else if (string(argv[4]) == "CHANNEL")
    {
        minArea = M_PI*pow(0.1*tamagno,2);
        maxArea = M_PI*pow(0.5*tamagno,2);
    }
    else // 2D
    {
        minArea = M_PI*pow(0.1*tamagno,2);
        maxArea = M_PI*pow(0.7*tamagno,2);
    }
    printf("\t MinArea: %.0f \t MaxArea: %.0f\n", minArea, maxArea);

#pragma omp parallel for ordered schedule(dynamic)
    for (pos=0; pos<ntot; pos++)
    {
        Mat tmp_frame0, tmp_frame, bgmask0, bgmask1, bgmask2, bgmask3;
        vector<Roi> RoisFramei;
//#pragma omp critical (readcap)
        {
            //tmp_frame = inframesMask[pos];
            //tmp_frame0 = inframes[pos];
            inframesMask[pos].copyTo(tmp_frame);
            inframes[pos].copyTo(tmp_frame0);
        }
//#pragma omp ordered
        {
            if (string(argv[4]) == "LINE")
            {
                //GaussianBlur(tmp_frame, tmp_frame, Size(tamagno, tamagno), 0, 0, BORDER_DEFAULT); // ONLY Harvard
                threshold(tmp_frame, bgmask0, round(threshold_Mask*(MaxAll+StdevAll)), 255, THRESH_TOZERO);
                threshold(tmp_frame0, tmp_frame0, round(threshold_Visual*(MaxAll+StdevAll)), 255, THRESH_TOZERO);
                bgmask0.copyTo(bgmask1);
                /*int max_ec=nopen;
                //* Not OPEN for Havard Datas
                for (int iter_ec=0; iter_ec<max_ec; iter_ec++)
                {
                    erode(bgmask1, bgmask1, getStructuringElement(MORPH_ELLIPSE, Size(GeomFilterLineSizeX,GeomFilterLineSizeY)), Point(-1,-1), 1);
                    dilate(bgmask1, bgmask1, getStructuringElement(MORPH_ELLIPSE, Size(GeomFilterLineSizeX,GeomFilterLineSizeY)), Point(-1,-1), 1);
                }*/
                erode(bgmask1, bgmask1, getStructuringElement(MORPH_ELLIPSE, Size(GeomFilterLineSizeX,GeomFilterLineSizeY)), Point(-1,-1), nopen);
                dilate(bgmask1, bgmask1, getStructuringElement(MORPH_ELLIPSE, Size(GeomFilterLineSizeX,GeomFilterLineSizeY)), Point(-1,-1), nopen);
                bgmask1.copyTo(bgmask2);
                /*max_ec=nclose;   // ONLY Harvard CLOSE
                for (int iter_ec=0; iter_ec<max_ec; iter_ec++)
                {
                    dilate(bgmask2, bgmask2, getStructuringElement(MORPH_ELLIPSE, Size(GeomFilterLineSizeX,GeomFilterLineSizeY)), Point(-1,-1), 1);
                    erode(bgmask2, bgmask2, getStructuringElement(MORPH_ELLIPSE, Size(GeomFilterLineSizeX,GeomFilterLineSizeY)), Point(-1,-1), 1);
                }*/
                dilate(bgmask2, bgmask2, getStructuringElement(MORPH_ELLIPSE, Size(GeomFilterLineSizeX,GeomFilterLineSizeY)), Point(-1,-1), nclose);
                erode(bgmask2, bgmask2, getStructuringElement(MORPH_ELLIPSE, Size(GeomFilterLineSizeX,GeomFilterLineSizeY)), Point(-1,-1), nclose);
                bgmask3 = min(tmp_frame, bgmask2);
                refineSegments_Cells(tmp_frame0, bgmask3, tamagno, pos, RoisFramei, argv, minArea, maxArea, MyObjectsFile);
            }
            else if (string(argv[4]) == "CHANNEL")
            {
                threshold(tmp_frame, bgmask0, round(threshold_Mask*(MaxAll+StdevAll)), 255, THRESH_TOZERO);
                threshold(tmp_frame0, tmp_frame0, round(threshold_Visual*(MaxAll+StdevAll)), 255, THRESH_TOZERO);
                bgmask0.copyTo(bgmask1);
                /*int max_ec=nopen;
                // Not OPEN for Havard Datas
                for (int iter_ec=0; iter_ec<max_ec; iter_ec++)
                {
                    erode(bgmask1, bgmask1, getStructuringElement(MORPH_ELLIPSE, Size(GeomFilterChanSizeX,GeomFilterChanSizeY)), Point(-1,-1), 2);
                    dilate(bgmask1, bgmask1, getStructuringElement(MORPH_ELLIPSE, Size(GeomFilterChanSizeX,GeomFilterChanSizeY)), Point(-1,-1), 2);
                }*/
                erode(bgmask1, bgmask1, getStructuringElement(MORPH_ELLIPSE, Size(GeomFilterChanSizeX,GeomFilterChanSizeY)), Point(-1,-1), nopen);
                dilate(bgmask1, bgmask1, getStructuringElement(MORPH_ELLIPSE, Size(GeomFilterChanSizeX,GeomFilterChanSizeY)), Point(-1,-1), nopen);
                bgmask1.copyTo(bgmask2);
                /*max_ec=nclose;   // ONLY Harvard CLOSE
                for (int iter_ec=0; iter_ec<max_ec; iter_ec++)
                {
                    dilate(bgmask2, bgmask2, getStructuringElement(MORPH_ELLIPSE, Size(GeomFilterChanSizeX,GeomFilterChanSizeY)), Point(-1,-1), 2);
                    erode(bgmask2, bgmask2, getStructuringElement(MORPH_ELLIPSE, Size(GeomFilterChanSizeX,GeomFilterChanSizeY)), Point(-1,-1), 2);
                }*/
                dilate(bgmask2, bgmask2, getStructuringElement(MORPH_ELLIPSE, Size(GeomFilterChanSizeX,GeomFilterChanSizeY)), Point(-1,-1), nclose);
                erode(bgmask2, bgmask2, getStructuringElement(MORPH_ELLIPSE, Size(GeomFilterChanSizeX,GeomFilterChanSizeY)), Point(-1,-1), nclose);
                //bgmask3 = min(tmp_frame, bgmask2);
                bgmask2.copyTo(bgmask3);
                refineSegments_Cells(tmp_frame0, bgmask3, tamagno, pos, RoisFramei, argv, minArea, maxArea, MyObjectsFile);
            }
            else // 2D
            {
                threshold(tmp_frame, bgmask0, round(threshold_Mask*(MaxAll+StdevAll)), 255, THRESH_TOZERO);
                threshold(tmp_frame0, tmp_frame0, round(threshold_Visual*(MaxAll+StdevAll)), 255, THRESH_TOZERO);
                bgmask0.copyTo(bgmask1);
                /*int max_ec=nopen;
                // Not OPEN for Havard Datas
                for (int iter_ec=0; iter_ec<max_ec; iter_ec++)
                {
                    erode(bgmask1, bgmask1, getStructuringElement(MORPH_ELLIPSE, Size(GeomFilter2DSizeX,GeomFilter2DSizeY)), Point(-1,-1), 1);
                    dilate(bgmask1, bgmask1, getStructuringElement(MORPH_ELLIPSE, Size(GeomFilter2DSizeX,GeomFilter2DSizeY)), Point(-1,-1), 1);
                }*/
                erode(bgmask1, bgmask1, getStructuringElement(MORPH_ELLIPSE, Size(GeomFilter2DSizeX,GeomFilter2DSizeY)), Point(-1,-1), nopen);
                dilate(bgmask1, bgmask1, getStructuringElement(MORPH_ELLIPSE, Size(GeomFilter2DSizeX,GeomFilter2DSizeY)), Point(-1,-1), nopen);
                bgmask1.copyTo(bgmask2);
                /*max_ec=nclose;   // ONLY Harvard CLOSE
                for (int iter_ec=0; iter_ec<max_ec; iter_ec++)
                {
                    dilate(bgmask2, bgmask2, getStructuringElement(MORPH_ELLIPSE, Size(GeomFilter2DSizeX,GeomFilter2DSizeY)), Point(-1,-1), 1);
                    erode(bgmask2, bgmask2, getStructuringElement(MORPH_ELLIPSE, Size(GeomFilter2DSizeX,GeomFilter2DSizeY)), Point(-1,-1), 1);
                }*/
                dilate(bgmask2, bgmask2, getStructuringElement(MORPH_ELLIPSE, Size(GeomFilter2DSizeX,GeomFilter2DSizeY)), Point(-1,-1), nclose);
                erode(bgmask2, bgmask2, getStructuringElement(MORPH_ELLIPSE, Size(GeomFilter2DSizeX,GeomFilter2DSizeY)), Point(-1,-1), nclose);
                bgmask3 = min(tmp_frame, bgmask2);
                //bgmask2.copyTo(bgmask3);
                refineSegments_Cells(tmp_frame0, bgmask3, tamagno, pos, RoisFramei, argv, minArea, maxArea, MyObjectsFile);
            }
        }
//#pragma omp critical (outframespos)
#pragma omp ordered
        {
            //RoisFrames.push_back(RoisFramei);
            RoisFrames[pos]=RoisFramei;
            //output_segmented[pos]=bgmask2;
            bgmask2.copyTo(output_segmented[pos]);
            if(VisOut==1)
            {
                imshow("Frame+", tmp_frame);
                imshow("Frame0", tmp_frame0);
                imshow("Mask0", bgmask0);
                imshow("Mask1", bgmask1);
                imshow("Mask2", bgmask2);
                waitKey(100);
            }
            bgmask0.release();
            bgmask1.release();
            bgmask2.release();
            bgmask3.release();
            tmp_frame0.release();
            tmp_frame.release();
            FreeVectorRoi(RoisFramei);
        }
    }
    //destroyAllWindows();
    fclose(MyObjectsFile);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Build trajectories
void BuildTrajecto (char* argv[], vector< vector<Roi> >& RoisFrames, vector< vector<Roi> >& Trajecto, int nbframesinmovie, double maxdisplacement, double dydX, int disapearancetime, double costBirthDeath, int Sizex, int Sizey, int temperaturedecrease, double normalizedIntensityWeight)
{
    printf("Build trajectories\n");
    int allobjectsinframes=0;
    for (int i=0; i<RoisFrames.size(); i++)
    {
        //printf("frame i=%i\t objects n=%i\n", i, RoisFrames[i].size());
        allobjectsinframes += RoisFrames[i].size();
        //printf("allobjectsinframes = %i\n", allobjectsinframes);
        /*if (i==RoisFrames.size()-1)
        	printf("\tend\n");*/
    }
    printf("\tall objects in all frames: %i\n", allobjectsinframes);
    //fflush(stdout);
    //double** objectstotrak;// [allobjectsinframes][7];
    double **objectstotrak;
    objectstotrak = new double*[allobjectsinframes];
    for (int oi = 0; oi <allobjectsinframes; oi++)
        objectstotrak[oi] = new double[7];
    //double** objectstotrak=(double**) malloc(allobjectsinframes * sizeof(double)); // working one!
//cout << "Start1!" << endl;
    int k=0;
    //printf("\tBuild Matrix - Start\n");
    for (int i=0; i<RoisFrames.size(); i++)
    {
        //printf("\tRoisFrames.size(): %i\n", RoisFrames.size());
        for (int j=0; j<RoisFrames[i].size(); j++)
        {
            //objectstotrak[k] = (double*) malloc(7 * sizeof(double));
            //printf("\t\tRoisFrames[i].size(): %i\n", RoisFrames[i].size());
            //	printf("\t\t\tRoisFrames[i][j].frame: %i\n", RoisFrames[i][j].frame);
            objectstotrak[k][0]=RoisFrames[i][j].frame;
            objectstotrak[k][1]=RoisFrames[i][j].origID+1;  // Remember remember +1 on the tag!
            objectstotrak[k][2]=RoisFrames[i][j].center.x;
            objectstotrak[k][3]=RoisFrames[i][j].center.y;
            objectstotrak[k][4]=RoisFrames[i][j].area;
            objectstotrak[k][5]=RoisFrames[i][j].peri;
            objectstotrak[k][6]=RoisFrames[i][j].iI;

            //printf("%i:\t\t %.0f\t %.0f\t %.3f\t %.3f\t %.3f\t %.3f\t %.3f\n", k, objectstotrak[k][0],objectstotrak[k][1],objectstotrak[k][2],objectstotrak[k][3],objectstotrak[k][4],objectstotrak[k][5],objectstotrak[k][6]);
            k++;
        }
    }
//cout << "End1!" << endl;
    printf("\t start tracking\n");
    vector<vector<double> > resultstrajectories;
    //int temperaturedecrease = 5;
    //double normalizedIntensityWeight = 0.12;
    int DoTrack= TrackingFunc::CRTracking::TrackingCellRace(temperaturedecrease, normalizedIntensityWeight, maxdisplacement, dydX,  disapearancetime, costBirthDeath, objectstotrak, allobjectsinframes, nbframesinmovie, &resultstrajectories, Sizex, Sizey, 5.0);
    printf("DoTrack: %i\n", DoTrack);
    //fflush(stdout);
    int ntrajecto=1;
//cout << "Start2!" << endl;
    if (resultstrajectories.size()>0)
    {
        vector <Roi> temp_trajectoi;
        for (int t=0; t<resultstrajectories.size(); t++)
        {
            int nframe, origID, giotag, trajtag, neigbors;
            nframe = round(resultstrajectories[t][1]);
            origID = round(resultstrajectories[t][2])-1;  // Remember remember +1 on the tag! from line 370~
            giotag = round(resultstrajectories[t][3]);
            trajtag = round(resultstrajectories[t][0]);
            neigbors = round(resultstrajectories[t][4]);
            /*for (int idx=0; idx<4; idx++)
            	printf("%.2f\t ",resultstrajectories[t][idx]);
            printf( "\n ");
            fflush(stdout);*/
            //printf("#: %i/%i\t Trajectory: %.0f\t Frame: %i\t FormerTag: %i\n", t, resultstrajectories.size(), resultstrajectories[t][0], nframe, origID);
            //printf("Inside1\n");
            Roi temp_roi;
            for(int h=0; h<RoisFrames[nframe].size(); h++)
            {
                if (RoisFrames[nframe][h].origID == origID)
                {
                    temp_roi = RoisFrames[nframe][h];
                    //printf("Roi copied\n");
                    temp_roi.mothertag = giotag;
                    temp_roi.trajtag = trajtag;
                    temp_roi.neigbors = neigbors;
                    break;
                }
            }
            //printf("\t\tRoi updated#\n");
            //printf("#: %i/%i\t Trajectory: %i\t Frame: %i\t FormerTag: %i\n\n", t, resultstrajectories.size(), ntrajecto, temp_roi.frame, temp_roi.origID);
            //printf("Inside2\n");
            if (t < resultstrajectories.size()-1)
            {
                if (resultstrajectories[t][0] == resultstrajectories[t+1][0])
                    temp_trajectoi.push_back(temp_roi);
                else if (resultstrajectories[t+1][0] == 1+resultstrajectories[t][0])
                {
                    temp_trajectoi.push_back(temp_roi);
                    Trajecto.push_back(temp_trajectoi);
                    temp_trajectoi.clear();
                    ntrajecto += 1;
                }
            }
            else
            {
                temp_trajectoi.push_back(temp_roi);
                Trajecto.push_back(temp_trajectoi);
            }
        }
        //printf("Inside3\n");

        for (int w=0; w<Trajecto.size(); w++)
            for(int z=1; z<Trajecto[w].size(); z++)
                Trajecto[w][z].mothertag = Trajecto[w][0].mothertag;
        //printf("Inside4\n");
    }
    printf("\tTrajectories: %lu\n", Trajecto.size());
    printf("\t end tracking\n");
    fflush(stdout);
}
//////////////////// Parallel build Trajectories
void ParallelBuildTrajecto (char* argv[], vector< vector<Roi> >& RoisFrames, vector< vector<Roi> >& Trajecto, int ntot, double maxdisplacement, double dydX, int disapearancetime, double costBirthDeath, int Sizex, int Sizey)
{
    // Split vector of rois per Frame in 4 sub-vectors
    int n4=floor(ntot/4);
    vector< vector< vector<Roi> > > SubRoisFrames(4);
    vector< vector< vector<Roi> > > SubTrajecto(4);
    int Frame0[4];
#pragma omp parallel for ordered schedule(dynamic) num_threads(4)
    for (int i=0; i<4; i++)
    {
        vector< vector<Roi> > SubRoisFrames4;
        //int Frame0=0;
        for(int j=0; j<RoisFrames.size(); j++)
        {
            if(j<(i+1)*n4+disapearancetime && j>=i*n4)
            {
                SubRoisFrames4.push_back(RoisFrames[j]);
                if(SubRoisFrames4.size()==1)
                    Frame0[i] = SubRoisFrames4.front()[0].frame;
                //printf("Frame=%i\t Frame0=%i\t size=%i\n", j, Frame0, SubRoisFrames4.back().size());
                for(int k=0; k<SubRoisFrames4.back().size(); k++)
                {
                    //printf("Frame0=%i\n", SubRoisFrames4[0][0].frame);
                    SubRoisFrames4.back()[k].frame=SubRoisFrames4.back()[k].frame-Frame0[i];
                }
            }
        }
        SubRoisFrames[i]=SubRoisFrames4;
    }

#pragma omp parallel for ordered schedule(dynamic) num_threads(4)
    for (int i=0; i<4; i++)
    {
        printf("Build Sub-trajectories:%i of size:%lu\n", i, SubRoisFrames[i].size());
        int allobjectsinframes=0;
        for (int j=0; j<SubRoisFrames[i].size(); j++)
        {
            allobjectsinframes += SubRoisFrames[i][j].size();
        }
        printf("\tall objects in all frmaes: %i\n", allobjectsinframes);

        double **objectstotrak;
        objectstotrak = new double*[allobjectsinframes];
        for (int oi = 0; oi <allobjectsinframes; oi++)
            objectstotrak[oi] = new double[7];
        int k=0;
        for (int j=0; j<SubRoisFrames[i].size(); j++)
        {
            for (int h=0; h<SubRoisFrames[i][j].size(); h++)
            {
                objectstotrak[k][0]=SubRoisFrames[i][j][h].frame;
                objectstotrak[k][1]=SubRoisFrames[i][j][h].origID+1;  // Remember remember +1 on the tag!
                objectstotrak[k][2]=SubRoisFrames[i][j][h].center.x;
                objectstotrak[k][3]=SubRoisFrames[i][j][h].center.y;
                objectstotrak[k][4]=SubRoisFrames[i][j][h].area;
                objectstotrak[k][5]=SubRoisFrames[i][j][h].peri;
                objectstotrak[k][6]=SubRoisFrames[i][j][h].iI;
                k++;
            }
        }
        printf("\t start tracking\n");
        vector<vector<double> > resultstrajectories;
        int temperaturedecrease = 5;
        double normalizedIntensityWeight = 0.12;
        int DoTrack;
#pragma omp critical
        DoTrack= TrackingFunc::CRTracking::TrackingCellRace(temperaturedecrease, normalizedIntensityWeight, maxdisplacement, dydX,  disapearancetime, costBirthDeath, objectstotrak, allobjectsinframes, SubRoisFrames[i].size(), &resultstrajectories, Sizex, Sizey, 5.0);
        printf("DoTrack: %i\n", DoTrack);
        fflush(stdout);
//#pragma omp barrier
        int ntrajecto=1;
        if (resultstrajectories.size()>0)
        {
            vector <Roi> temp_trajectoi;
            for (int t=0; t<resultstrajectories.size(); t++)
            {
                int nframe, origID, giotag, trajtag, neigbors;
                nframe = round(resultstrajectories[t][1]);
                origID = round(resultstrajectories[t][2])-1;  // Remember remember +1 on the tag! from line 370~
                giotag = round(resultstrajectories[t][3]);
                trajtag = round(resultstrajectories[t][0]);
                neigbors = round(resultstrajectories[t][4]);
                Roi temp_roi;
                for(int h=0; h<SubRoisFrames[i][nframe].size(); h++)
                {
                    if (SubRoisFrames[i][nframe][h].origID == origID)
                    {
                        temp_roi = SubRoisFrames[i][nframe][h];
                        temp_roi.frame = SubRoisFrames[i][nframe][h].frame + Frame0[i];
                        temp_roi.mothertag = giotag;
                        temp_roi.trajtag = trajtag;
                        temp_roi.neigbors = neigbors;
                        break;
                    }
                }
                if (t < resultstrajectories.size()-1)
                {
                    if (resultstrajectories[t][0] == resultstrajectories[t+1][0])
                        temp_trajectoi.push_back(temp_roi);
                    else if (resultstrajectories[t+1][0] == 1+resultstrajectories[t][0])
                    {
                        temp_trajectoi.push_back(temp_roi);
//#pragma omp critical
                        SubTrajecto[i].push_back(temp_trajectoi);
                        temp_trajectoi.clear();
                        ntrajecto += 1;
                    }
                }
                else
                {
                    temp_trajectoi.push_back(temp_roi);
                    SubTrajecto[i].push_back(temp_trajectoi);
                }
            }

            for (int w=0; w<Trajecto.size(); w++)
                for(int z=1; z<Trajecto[w].size(); z++)
                    Trajecto[w][z].mothertag = Trajecto[w][0].mothertag;
        }
        printf("\tTrajectories: %lu\n", Trajecto.size());
        printf("\t end tracking\n");
        fflush(stdout);
    }
    Trajecto = SubTrajecto[0];

}

// Serial build Trajectories
void SerialBuildTrajecto (char* argv[], vector< vector<Roi> >& RoisFrames, vector< vector<Roi> >& Trajecto, int ntot, double maxdisplacement, double dydX, int disapearancetime, double costBirthDeath, int Sizex, int Sizey, int temperaturedecrease, double normalizedIntensityWeight)
{
    //ParallelBuildTrajecto(argv, RoisFrames, Trajecto, ntot, maxdisplacement, dydX, disapearancetime, costBirthDeath, Sizex, Sizey);
    //printf("End Parallel!\n");
    BuildTrajecto(argv, RoisFrames, Trajecto, ntot, maxdisplacement, dydX, disapearancetime, costBirthDeath, Sizex, Sizey, temperaturedecrease, normalizedIntensityWeight);
}

// Draw trajectories
void DrawTrajecto (vector<Mat>& output_color, vector<vector<Roi> >& Trajectos)
{
    Mat tmp_output_color;
    vector <Scalar> ColorsVector;
    Scalar color;
    vector <int> TAG_trajecto;
    for (int i=0; i<Trajectos.size(); i++)
    {

        if (Trajectos[i][0].mothertag == 0 )
        {
            color = Scalar( rand()&255, rand()&255, rand()&255 );
            //color = Scalar( rand()&255, 0, 0);
            ColorsVector.push_back(color);
            TAG_trajecto.push_back(Trajectos[i][0].trajtag);
        }
        else
        {
            int giotag = Trajectos[i][0].mothertag;
            if (ColorsVector.size() >0 )
            {
                for (int cc=0; cc<TAG_trajecto.size(); cc++)
                {
                    if (TAG_trajecto[cc] == giotag)
                    {
                        color=ColorsVector[cc];
                        break;
                    }
                }
            }
            else
            {
                color = Scalar( rand()&255, rand()&255, rand()&255 );
                //color = Scalar( rand()&255, 0, 0);
                ColorsVector.push_back(color);
                TAG_trajecto.push_back(Trajectos[i][0].trajtag);
            }
        }

        for (int t=0; t<Trajectos[i].size(); t++)
        {
            tmp_output_color = output_color[Trajectos[i][t].frame];
            vector<vector<Point> > tmp_vector1_contour;
            tmp_vector1_contour.push_back(Trajectos[i][t].contour);
            if (Trajectos[i].size() > 5)
                drawContours(tmp_output_color, tmp_vector1_contour, -1, color, 1);
            output_color[Trajectos[i][t].frame] = tmp_output_color;
            int countpoint =0;
            if(t>0)
            {
                while(countpoint < 25 && t-countpoint-1 > 0)
                {
                    line(output_color[Trajectos[i][t].frame], Trajectos[i][t-countpoint].center, Trajectos[i][t-countpoint-1].center, color);
                    /*printf("countpoint: %i\n", countpoint);
                    fflush(stdout);*/
                    countpoint += 1;
                }
            }
        }
    }
}
void DrawTrajectoGreen (vector<Mat>& output_color, vector<vector<Roi> >& Trajectos)
{
    Mat tmp_output_color;
    Scalar color;
    for (int i=0; i<Trajectos.size(); i++)
    {
        color = Scalar( 0, rand()&255, 0);
        for (int t=0; t<Trajectos[i].size(); t++)
        {
            tmp_output_color = output_color[Trajectos[i][t].frame];
            vector<vector<Point> > tmp_vector1_contour;
            tmp_vector1_contour.push_back(Trajectos[i][t].contour);
            if (Trajectos[i].size() > 5)
                drawContours(tmp_output_color, tmp_vector1_contour, -1, color, 1);
            output_color[Trajectos[i][t].frame] = tmp_output_color;
            int countpoint =0;
            if(t>0)
            {
                while(countpoint < 25 && t-countpoint-1 > 0)
                {
                    line(output_color[Trajectos[i][t].frame], Trajectos[i][t-countpoint].center, Trajectos[i][t-countpoint-1].center, color);
                    countpoint += 1;
                }
            }
        }
    }
}
// Draw Trajectories
void DrawTrajectoGray (vector<Mat>& null_gray, vector<vector<Roi> >& Trajectos)
{
    Mat tmp_output;
    for (int i=0; i<Trajectos.size(); i++)
    {
        for (int t=0; t<Trajectos[i].size(); t++)
        {
            tmp_output = null_gray[Trajectos[i][t].frame];
            vector<vector<Point> > tmp_vector1_contour;
            tmp_vector1_contour.push_back(Trajectos[i][t].contour);
            if (Trajectos[i][t].mothertag ==0)
                drawContours(tmp_output, tmp_vector1_contour, -1, Trajectos[i][t].trajtag, -1);
            else
            {
                int giotag = Trajectos[i][t].mothertag;
                drawContours(tmp_output, tmp_vector1_contour, -1, giotag, -1);
            }
            null_gray[Trajectos[i][t].frame] = tmp_output;
        }
    }
    for (int j=0; j<null_gray.size(); j++)
    {
        tmp_output = null_gray[j];
        imshow("Seg&Track", tmp_output);
        waitKey(500);
    }
}
// Save Trajectories
void SaveAllTrajecto2 (FILE* myTrajectofileFrame, FILE* myTrajectofileT, FILE* myTrajectofileX, FILE* myTrajectofileY, vector<vector<Roi> >& Trajectos, string nametrajecto, string namecelltype, double scalefactor, double TimeInterval)
{
    for (int i=0; i<Trajectos.size(); i++)
    {
        for (int j=0; j<Trajectos[i].size(); j++)
        {
            if (j==0)
            {
                string nametrajectoi=nametrajecto;
                std::string snum;
                std::stringstream snumout;
                snumout << Trajectos[i][0].trajtag-1;
                snum = snumout.str();
                snum = "#" + snum;
                nametrajectoi.append(snum);
                fprintf( myTrajectofileFrame, "%s\t%s\n", nametrajectoi.c_str(), namecelltype.c_str());
                fprintf( myTrajectofileT, "%s\t%s\n", nametrajectoi.c_str(), namecelltype.c_str());
                fprintf( myTrajectofileX, "%s\t%s\n", nametrajectoi.c_str(), namecelltype.c_str());
                fprintf( myTrajectofileY, "%s\t%s\n", nametrajectoi.c_str(), namecelltype.c_str());
            }
            fprintf( myTrajectofileFrame, "%i\t", Trajectos[i][j].frame);
            fprintf( myTrajectofileT, "%.3f\t", Trajectos[i][j].frame*TimeInterval);
            fprintf( myTrajectofileX, "%.3f\t", Trajectos[i][j].center.x*scalefactor);
            fprintf( myTrajectofileY, "%.3f\t", Trajectos[i][j].center.y*scalefactor);
        }
        fprintf( myTrajectofileFrame, "\n");
        fprintf( myTrajectofileT, "\n");
        fprintf( myTrajectofileX, "\n");
        fprintf( myTrajectofileY, "\n");
        fflush(myTrajectofileFrame);
        fflush(myTrajectofileT);
        fflush(myTrajectofileX);
        fflush(myTrajectofileY);
    }
    fflush(myTrajectofileFrame);
    fflush(myTrajectofileT);
    fflush(myTrajectofileX);
    fflush(myTrajectofileY);
}
// Save Trajectories
void SaveAllTrajecto3 (FILE* myTrajectofileFrame, FILE* myTrajectofileT, FILE* myTrajectofileX, FILE* myTrajectofileY, vector<vector<Roi> >& Trajectos, string nametrajecto, string namecelltype, double scalefactor, double TimeInterval, int height)
{
    for (int i=0; i<Trajectos.size(); i++)
    {
        for (int j=0; j<Trajectos[i].size(); j++)
        {
            if (j==0)
            {
                string nametrajectoi=nametrajecto;
                std::string snum;
                std::stringstream snumout;
                snumout << Trajectos[i][0].trajtag-1;
                snum = snumout.str();
                snum = "#" + snum;
                nametrajectoi.append(snum);
                fprintf( myTrajectofileFrame, "%s\t%s\n", nametrajectoi.c_str(), namecelltype.c_str());
                fprintf( myTrajectofileT, "%s\t%s\n", nametrajectoi.c_str(), namecelltype.c_str());
                fprintf( myTrajectofileX, "%s\t%s\n", nametrajectoi.c_str(), namecelltype.c_str());
                fprintf( myTrajectofileY, "%s\t%s\n", nametrajectoi.c_str(), namecelltype.c_str());
            }
            fprintf( myTrajectofileFrame, "%i\t", Trajectos[i][j].frame);
            fprintf( myTrajectofileT, "%.3f\t", Trajectos[i][j].frame*TimeInterval);
            fprintf( myTrajectofileX, "%.3f\t", Trajectos[i][j].center.x*scalefactor);
            fprintf( myTrajectofileY, "%.3f\t", (height-Trajectos[i][j].center.y)*scalefactor);
        }
        fprintf( myTrajectofileFrame, "\n");
        fprintf( myTrajectofileT, "\n");
        fprintf( myTrajectofileX, "\n");
        fprintf( myTrajectofileY, "\n");
        fflush(myTrajectofileFrame);
        fflush(myTrajectofileT);
        fflush(myTrajectofileX);
        fflush(myTrajectofileY);
    }
    fflush(myTrajectofileFrame);
    fflush(myTrajectofileT);
    fflush(myTrajectofileX);
    fflush(myTrajectofileY);
}
// Save Trajectories
void SaveAllTrajecto4 (FILE* myTrajectofileFrame, FILE* myTrajectofileT, FILE* myTrajectofileX, FILE* myTrajectofileY, FILE* myTrajectofileA, FILE* myTrajectofileP, vector<vector<Roi> >& Trajectos, string nametrajecto, string namecelltype, double scalefactor, double TimeInterval, int height)
{
    for (int i=0; i<Trajectos.size(); i++)
    {
        for (int j=0; j<Trajectos[i].size(); j++)
        {
            if (j==0)
            {
                string nametrajectoi=nametrajecto;
                std::string snum;
                std::stringstream snumout;
                snumout << Trajectos[i][0].trajtag-1;
                snum = snumout.str();
                snum = "#" + snum;
                nametrajectoi.append(snum);
                fprintf( myTrajectofileFrame, "%s\t%s\n", nametrajectoi.c_str(), namecelltype.c_str());
                fprintf( myTrajectofileT, "%s\t%s\n", nametrajectoi.c_str(), namecelltype.c_str());
                fprintf( myTrajectofileX, "%s\t%s\n", nametrajectoi.c_str(), namecelltype.c_str());
                fprintf( myTrajectofileY, "%s\t%s\n", nametrajectoi.c_str(), namecelltype.c_str());
                fprintf( myTrajectofileA, "%s\t%s\n", nametrajectoi.c_str(), namecelltype.c_str());
                fprintf( myTrajectofileP, "%s\t%s\n", nametrajectoi.c_str(), namecelltype.c_str());
            }
            fprintf( myTrajectofileFrame, "%i\t", Trajectos[i][j].frame);
            fprintf( myTrajectofileT, "%.3f\t", Trajectos[i][j].frame*TimeInterval);
            fprintf( myTrajectofileX, "%.3f\t", Trajectos[i][j].center.x*scalefactor);
            fprintf( myTrajectofileY, "%.3f\t", (height-Trajectos[i][j].center.y)*scalefactor);
            fprintf( myTrajectofileA, "%.3f\t", Trajectos[i][j].area*gsl_pow_2(scalefactor));
            fprintf( myTrajectofileP, "%.3f\t", Trajectos[i][j].peri*scalefactor);
        }
        fprintf( myTrajectofileFrame, "\n");
        fprintf( myTrajectofileT, "\n");
        fprintf( myTrajectofileX, "\n");
        fprintf( myTrajectofileY, "\n");
        fprintf( myTrajectofileA, "\n");
        fprintf( myTrajectofileP, "\n");
        fflush(myTrajectofileFrame);
        fflush(myTrajectofileT);
        fflush(myTrajectofileX);
        fflush(myTrajectofileY);
        fflush(myTrajectofileA);
        fflush(myTrajectofileP);
    }
    fflush(myTrajectofileFrame);
    fflush(myTrajectofileT);
    fflush(myTrajectofileX);
    fflush(myTrajectofileY);
    fflush(myTrajectofileA);
    fflush(myTrajectofileP);
}
void SaveAllTrajecto5 (FILE* myTrajectofileFrame, FILE* myTrajectofileT, FILE* myTrajectofileX, FILE* myTrajectofileY, FILE* myTrajectofileA, FILE* myTrajectofileP, FILE* myTrajectofileC, vector<vector<Roi> >& Trajectos, string nametrajecto, string namecelltype, double scalefactor, double TimeInterval, int height)
{
    for (int i=0; i<Trajectos.size(); i++)
    {
        for (int j=0; j<Trajectos[i].size(); j++)
        {
            if (j==0)
            {
                string nametrajectoi=nametrajecto;
                std::string snum;
                std::stringstream snumout;
                snumout << Trajectos[i][0].trajtag-1;
                snum = snumout.str();
                snum = "#" + snum;
                nametrajectoi.append(snum);
                fprintf( myTrajectofileFrame, "%s\t%s\n", nametrajectoi.c_str(), namecelltype.c_str());
                fprintf( myTrajectofileT, "%s\t%s\n", nametrajectoi.c_str(), namecelltype.c_str());
                fprintf( myTrajectofileX, "%s\t%s\n", nametrajectoi.c_str(), namecelltype.c_str());
                fprintf( myTrajectofileY, "%s\t%s\n", nametrajectoi.c_str(), namecelltype.c_str());
                fprintf( myTrajectofileA, "%s\t%s\n", nametrajectoi.c_str(), namecelltype.c_str());
                fprintf( myTrajectofileP, "%s\t%s\n", nametrajectoi.c_str(), namecelltype.c_str());
                fprintf( myTrajectofileC, "%s\t%s\n", nametrajectoi.c_str(), namecelltype.c_str());
            }
            fprintf( myTrajectofileFrame, "%i\t", Trajectos[i][j].frame);
            fprintf( myTrajectofileT, "%.3f\t", Trajectos[i][j].frame*TimeInterval);
            fprintf( myTrajectofileX, "%.3f\t", Trajectos[i][j].center.x*scalefactor);
            fprintf( myTrajectofileY, "%.3f\t", (height-Trajectos[i][j].center.y)*scalefactor);
            fprintf( myTrajectofileA, "%.3f\t", Trajectos[i][j].area*gsl_pow_2(scalefactor));
            fprintf( myTrajectofileP, "%.3f\t", Trajectos[i][j].peri*scalefactor);
            fprintf( myTrajectofileC, "%.3f\t", Trajectos[i][j].convexity);
        }
        fprintf( myTrajectofileFrame, "\n");
        fprintf( myTrajectofileT, "\n");
        fprintf( myTrajectofileX, "\n");
        fprintf( myTrajectofileY, "\n");
        fprintf( myTrajectofileA, "\n");
        fprintf( myTrajectofileP, "\n");
        fprintf( myTrajectofileC, "\n");
        fflush(myTrajectofileFrame);
        fflush(myTrajectofileT);
        fflush(myTrajectofileX);
        fflush(myTrajectofileY);
        fflush(myTrajectofileA);
        fflush(myTrajectofileP);
        fflush(myTrajectofileC);
    }
    fflush(myTrajectofileFrame);
    fflush(myTrajectofileT);
    fflush(myTrajectofileX);
    fflush(myTrajectofileY);
    fflush(myTrajectofileA);
    fflush(myTrajectofileP);
    fflush(myTrajectofileC);
}
// Draw Objects RoiFrame
void DrawRoisFrame(vector<Mat>& output_color, int ntot, vector<vector<Roi> >& RoisFrames)
{
    Mat tmp_output_color;
    Scalar color(255, 255, 255);
    for (int i=0; i<RoisFrames.size(); i++)
    {
        for (int j=0; j<RoisFrames[i].size(); j++)
        {
            tmp_output_color = output_color[RoisFrames[i][j].frame];
            vector<vector<Point> > tmp_vector1_contour;
            tmp_vector1_contour.push_back(RoisFrames[i][j].contour);
            drawContours(tmp_output_color, tmp_vector1_contour, -1, color, 1);
            output_color[RoisFrames[i][j].frame] = tmp_output_color;
        }
    }
    for (int j=0; j<output_color.size(); j++)
    {
        tmp_output_color = output_color[j];
        /*imshow("Seg&Track", tmp_output_color);
        waitKey(50);*/
    }
}
/////////////////////////////
//////////////////////////////////// Velocity
void Fastest100Trai (double& lengthlimit, vector<Roi>& Trajecton, int& minidx, int& mintime, double& effdist)
{
    int nbroi=Trajecton.size();
    double cumdist[nbroi];
    cumdist[0]=0;
    double cumdistvalue=0;
    int fastestidx;
    for (int i=1; i<nbroi; i++)
    {
        if (Trajecton[i].center.x > Trajecton[i-1].center.x)
            cumdistvalue += RoiDistanceX(Trajecton[i-1], Trajecton[i]);
        else
            cumdistvalue -= RoiDistanceX(Trajecton[i-1], Trajecton[i]);
        cumdist[i] = cumdistvalue;
    }
    int mindiftimelength = INT_MAX;
    for (int i=0; i<nbroi-1; i++)
    {
        for (int j=i+1; j<nbroi; j++)
        {
            if (fabs(cumdist[j]-cumdist[i]) >= lengthlimit)
            {
                if (mindiftimelength > j-i)
                {
                    mindiftimelength = Trajecton[j].frame - Trajecton[i].frame;
                    fastestidx = i;
                    effdist = fabs(cumdist[j]-cumdist[i]);
                }
                break;
            }
        }
    }
    minidx = fastestidx; // index to start
    mintime = mindiftimelength; //  time
}
void extractfastest(vector<Mat>& input_frames, vector<Mat>& outputfastest, vector<Roi>& Trajecto, int& minidx, int& mintime, double lengthlimit, int visout)
{
    //cout << "Start extract!" << endl;
    vector <Point> FullContour;
    int i=minidx;
    while (Trajecto[i].frame-Trajecto[minidx].frame <= mintime)
    {
        for (int j=0; j<Trajecto[i].contour.size(); j++)
        {
            FullContour.push_back(Trajecto[i].contour[j]);
        }
        i++;
        if (i > Trajecto.size()-1)
            break;
    }
    //printf("FullContour!\n");
    Point2f center, vtx[4];
    float rad;
    RotatedRect box = minAreaRect(Mat(FullContour));
    box.points(vtx);

    //int minx=INT_MAX, miny=INT_MAX, maxx=0, maxy=0; /// Corrected on 18 Feb 2013
    int minx=input_frames[0].cols-1, miny=input_frames[0].rows-1, maxx=0, maxy=0;
    //extend the box!
    for (int i=0; i<4; i++)
    {
        if (vtx[i].x < minx)
            minx = vtx[i].x;
        if (vtx[i].y < miny)
            miny = vtx[i].y;
        if (vtx[i].x > maxx)
            maxx = vtx[i].x;
        if (vtx[i].y > maxy)
            maxy = vtx[i].y;
    }
    minx -= 45;
    miny -= 25;
    maxx += 45;
    maxy += 25;
    if( minx < 0)
        minx=0;
    if( miny < 0)
        miny=0;
    if( maxx > input_frames[0].cols-1)
        maxx=input_frames[0].cols-1;
    if( maxy > input_frames[0].rows-1)
        maxy=input_frames[0].rows-1;
    //printf("minx:%i, miny:%i\t maxx:%i, maxy:%i\n", minx, miny, maxx, maxy);
    Size extboxsize = Size(maxx-minx+1, maxy-miny+1);
    i=minidx;
    while (Trajecto[i].frame-Trajecto[minidx].frame <= mintime)
    {
        Mat in;
        int frame = Trajecto[i].frame;
        in = input_frames[frame];
        //normalize(in, in, 0, 255, NORM_MINMAX, -1);
        //in.convertTo(in, CV_8UC1);
        Mat out = Mat::zeros(extboxsize, CV_8UC3);
        for(int c=minx; c<=maxx; c++)
            for(int r=miny; r<=maxy; r++)
                out.at<Vec3b>(r-miny, c-minx) = in.at<Vec3b>(r, c);
        //cvtColor(out, out, CV_GRAY2RGB);
        //cout << "Draw Circle" << i << " of " << minidx+mintime << endl;
        circle(out, Point(Trajecto[i].center.x-minx,Trajecto[i].center.y-miny), 1, Scalar(0,255,0), -1, 8, 0);
        resize(out, out, Size(0,0), 3, 3, INTER_CUBIC);
        outputfastest.push_back(out);
        if(visout==1)
        {
            imshow("Fastest", out);
            waitKey(100);
        }
        i++;
        if (i > Trajecto.size()-1)
            break;
    }
    printf("End extractfastest\n");
}

void extract_all_i(vector<Mat>& input_frames, vector<Mat>& outputfastest, vector<Roi>& Trajecto)
{
    //cout << "Start extract!" << endl;
    vector <Point> FullContour;
    for(int i=0; i<Trajecto.size(); i++)
    {
        for (int j=0; j<Trajecto[i].contour.size(); j++)
        {
            FullContour.push_back(Trajecto[i].contour[j]);
        }
    }
    //cout << "FullContour!" << endl;
    Point2f center, vtx[4];
    float rad;
    RotatedRect box = minAreaRect(Mat(FullContour));
    box.points(vtx);

    int minx=input_frames[0].cols-1, miny=input_frames[0].rows-1, maxx=0, maxy=0;
    //extend the box!
    for (int i=0; i<4; i++)
    {
        if (vtx[i].x < minx)
            minx = vtx[i].x;
        if (vtx[i].y < miny)
            miny = vtx[i].y;
        if (vtx[i].x > maxx)
            maxx = vtx[i].x;
        if (vtx[i].y > maxy)
            maxy = vtx[i].y;
    }
    minx -= 45;
    miny -= 45;
    maxx += 45;
    maxy += 45;
    if( minx < 0)
        minx=0;
    if( miny < 0)
        miny=0;
    if( maxx > input_frames[0].cols-1)
        maxx=input_frames[0].cols-1;
    if( maxy > input_frames[0].rows-1)
        maxy=input_frames[0].rows-1;
    //printf("minx:%i, miny:%i\t maxx:%i, maxy:%i\n", minx, miny, maxx, maxy);
    Size extboxsize = Size(maxx-minx+1, maxy-miny+1);
    //cout << "Box defined" << endl;
    for(int i=0; i<Trajecto.size(); i++)
    {
        //cout << i << endl;
        Mat in;
        int frame = Trajecto[i].frame;
        in = input_frames[frame];
        //cout << "export frames" << endl;
        //normalize(in, in, 0, 255, NORM_MINMAX, -1);
        //in.convertTo(in, CV_8UC1);
        Mat out = Mat::zeros(extboxsize, CV_8UC3);
        for(int c=minx; c<=maxx; c++)
            for(int r=miny; r<=maxy; r++)
                out.at<Vec3b>(r-miny, c-minx) = in.at<Vec3b>(r, c);
        //cvtColor(out, out, CV_GRAY2RGB);
        //cout << "Draw Circle" << i << " of " << minidx+mintime << endl;
        circle(out, Point(Trajecto[i].center.x-minx,Trajecto[i].center.y-miny), 3, Scalar(0,255,0), -1, 8, 0);
        resize(out, out, Size(0,0), 1.5, 1.5, INTER_CUBIC);
        outputfastest.push_back(out);
    }
    //cout << "End extract" << endl;
}

int ReadConfig(char* argv, double& scalefactor, double& TimeInterval, double& SubBackFactor, int& nclose, int& nopen, double& threshold_Mask, double& threshold_Visual, int& GeomFilterLineSizeX, int& GeomFilterLineSizeY, int& GeomFilterChanSizeX, int& GeomFilterChanSizeY, int& GeomFilter2DSizeX, int& GeomFilter2DSizeY, double& maxdisplacement, double& dYdX, double& dYdX_Line, int& disapearancetime, double& costBirthDeath, int& temperaturedecrease, double& normalizedIntensityWeight, int& VisualOutput, int& ObjSize, int& SaveSingleVisualOutputAs)
{
    ifstream infile;
    string templine;
    string line;
    stringstream strstrX, strstrname;
    int countline=0;

    infile.open ( argv, ios::in );
    if(infile.is_open())
    {
        printf("Setting-File: %s -- OK\n", argv);
        while (infile.good())
        {
            getline(infile, templine);
            if (countline == 1)
            {
                float tempfloat;
                sscanf(templine.c_str(), "%*s %*s = %f;", &tempfloat);
                scalefactor=double(tempfloat);
                //printf("%s\t ScaleFactor = %f\n", templine.c_str(), scalefactor);
            }
            else if (countline == 2)
            {
                float tempfloat;
                sscanf(templine.c_str(), "%*s %*s = %f;", &tempfloat);
                TimeInterval=double(tempfloat);
            }
            else if (countline == 4)
            {
                float tempfloat;
                sscanf(templine.c_str(), "%*s %*s = %f;", &tempfloat);
                SubBackFactor=double(tempfloat);
            }
            else if (countline == 5)
            {
                sscanf(templine.c_str(), "%*s %*s = %i;", &nclose);
            }
            else if (countline == 6)
            {
                sscanf(templine.c_str(), "%*s %*s = %i;", &nopen);
            }
            else if (countline == 7)
            {
                float tempfloat;
                sscanf(templine.c_str(), "%*s %*s = %f;", &tempfloat);
                threshold_Mask=double(tempfloat);
            }
            else if (countline == 8)
            {
                float tempfloat;
                sscanf(templine.c_str(), "%*s %*s = %f;", &tempfloat);
                threshold_Visual=double(tempfloat);
            }
            else if (countline == 9)
            {
                sscanf(templine.c_str(), "%*s %*s = %i;", &GeomFilterLineSizeX);
            }
            else if (countline == 10)
            {
                sscanf(templine.c_str(), "%*s %*s = %i;", &GeomFilterLineSizeY);
            }
            else if (countline == 11)
            {
                sscanf(templine.c_str(), "%*s %*s = %i;", &GeomFilterChanSizeX);
            }
            else if (countline == 12)
            {
                sscanf(templine.c_str(), "%*s %*s = %i;", &GeomFilterChanSizeY);
            }
            else if (countline == 13)
            {
                sscanf(templine.c_str(), "%*s %*s = %i;", &GeomFilter2DSizeX);
            }
            else if (countline == 14)
            {
                sscanf(templine.c_str(), "%*s %*s = %i;", &GeomFilter2DSizeY);
            }
            else if (countline == 16)
            {
                float tempfloat;
                sscanf(templine.c_str(), "%*s %*s = %f;", &tempfloat);
                maxdisplacement=double(tempfloat);
            }
            else if (countline == 17)
            {
                float tempfloat;
                sscanf(templine.c_str(), "%*s %*s = %f;", &tempfloat);
                dYdX=double(tempfloat);
            }
            else if (countline == 18)
            {
                float tempfloat;
                sscanf(templine.c_str(), "%*s %*s = %f;", &tempfloat);
                dYdX_Line=double(tempfloat);
            }
            else if (countline == 19)
            {
                sscanf(templine.c_str(), "%*s %*s = %i;", &disapearancetime);
            }
            else if (countline == 20)
            {
                float tempfloat;
                sscanf(templine.c_str(), "%*s %*s = %f;", &tempfloat);
                costBirthDeath=double(tempfloat);
            }
            else if (countline == 21)
            {
                sscanf(templine.c_str(), "%*s %*s = %i;", &temperaturedecrease);
            }
            else if (countline == 22)
            {
                float tempfloat;
                sscanf(templine.c_str(), "%*s %*s = %f;", &tempfloat);
                normalizedIntensityWeight=double(tempfloat);
            }
            else if (countline == 24)
            {
                sscanf(templine.c_str(), "%*s %*s = %i;", &VisualOutput);
            }
            else if (countline == 25)
            {
                sscanf(templine.c_str(), "%*s %*s = %i;", &ObjSize);
            }
            else if (countline == 26)
            {
                sscanf(templine.c_str(), "%*s %*s = %i;", &SaveSingleVisualOutputAs);
            }
            countline += 1;
        }

        printf("Scale = %.3f [um/px]\n", scalefactor);
        printf("Time Interval = %.3f [min]\n", TimeInterval);
        printf("Enhance Background Substraction = %.3f\n", SubBackFactor);
        printf("N close = %i\n", nclose);
        printf("N open = %i\n", nopen);
        printf("Threshold Mask = %.3f\n", threshold_Mask);
        printf("Threshold Visual = %.3f\n", threshold_Visual);
        printf("Show Segmented [0=NO, 1=YES] = %i\n", VisualOutput);
        printf("Object Size = %i\n", ObjSize);
        printf("Save Single Track Movies As [0=avi, 1=tiff] = %i\n", SaveSingleVisualOutputAs);
        return 0;
    }

    else
    {
        printf("Default Setting\n");
//Deafault
//Data
        scalefactor = 1.3;
        TimeInterval = 3;
//Segmentation
        SubBackFactor = 3;
        nclose = 1;
        nopen = 1;
        threshold_Mask = 3;
        threshold_Visual = 0.05;
        GeomFilterLineSizeX = 5;
        GeomFilterLineSizeY = 3;
        GeomFilterChanSizeX = 5;
        GeomFilterChanSizeY = 3;
        GeomFilter2DSizeX = 5;
        GeomFilter2DSizeY = 5;
//Tracking
        maxdisplacement = 77;
        dYdX = 1;
        dYdX_Line = 12;
        disapearancetime = 3;
        costBirthDeath = 1;
        temperaturedecrease = 5;
        normalizedIntensityWeight = 0.1;
//Output
        VisualOutput = 0;
        ObjSize = 0;

        printf("Scale = %.3f [um/px]\n", scalefactor);
        printf("Time Interval = %.3f [min]\n", TimeInterval);
        printf("Enhance Background Substraction = %.3f\n", SubBackFactor);
        printf("N close = %i\n", nclose);
        printf("N open = %i\n", nopen);
        printf("Threshold Mask = %.3f\n", threshold_Mask);
        printf("Threshold Visual = %.3f\n", threshold_Visual);
        printf("Show Segmented = %i\n", VisualOutput);
        printf("Object Size = %i\n", ObjSize);

        return 1;
    }


}

int ReadConfig2Channels(char* argv, double& scalefactor, double& TimeInterval, double& SubBackFactor, int& nclose, int& nopen, double& threshold_Mask, double& threshold_Visual, int& GeomFilterLineSizeX, int& GeomFilterLineSizeY, int& GeomFilterChanSizeX, int& GeomFilterChanSizeY, int& GeomFilter2DSizeX, int& GeomFilter2DSizeY, double& maxdisplacement, double& dYdX, double& dYdX_Line, int& disapearancetime, double& costBirthDeath, int& temperaturedecrease, double& normalizedIntensityWeight, int& VisualOutput, int& ObjSize, double& SubBackFactorB, int& ncloseB, int& nopenB, double& threshold_MaskB, int& GeomFilterLineSizeXB, int& GeomFilterLineSizeYB, int& GeomFilterChanSizeXB, int& GeomFilterChanSizeYB, int& GeomFilter2DSizeXB, int& GeomFilter2DSizeYB, int& ObjSizeB)
{
    ifstream infile;
    string templine;
    string line;
    stringstream strstrX, strstrname;
    int countline=0;

    infile.open ( argv, ios::in );
    if(infile.is_open())
    {
        printf("Setting-File: %s -- OK\n", argv);
        while (infile.good())
        {
            getline(infile, templine);
            if (countline == 1)
            {
                float tempfloat;
                sscanf(templine.c_str(), "%*s %*s = %f;", &tempfloat);
                scalefactor=double(tempfloat);
                //printf("%s\t ScaleFactor = %f\n", templine.c_str(), scalefactor);
            }
            else if (countline == 2)
            {
                float tempfloat;
                sscanf(templine.c_str(), "%*s %*s = %f;", &tempfloat);
                TimeInterval=double(tempfloat);
            }
            else if (countline == 4)
            {
                float tempfloat;
                sscanf(templine.c_str(), "%*s %*s = %f;", &tempfloat);
                SubBackFactor=double(tempfloat);
            }
            else if (countline == 5)
            {
                sscanf(templine.c_str(), "%*s %*s = %i;", &nclose);
            }
            else if (countline == 6)
            {
                sscanf(templine.c_str(), "%*s %*s = %i;", &nopen);
            }
            else if (countline == 7)
            {
                float tempfloat;
                sscanf(templine.c_str(), "%*s %*s = %f;", &tempfloat);
                threshold_Mask=double(tempfloat);
            }
            else if (countline == 8)
            {
                float tempfloat;
                sscanf(templine.c_str(), "%*s %*s = %f;", &tempfloat);
                threshold_Visual=double(tempfloat);
            }
            else if (countline == 9)
            {
                sscanf(templine.c_str(), "%*s %*s = %i;", &GeomFilterLineSizeX);
            }
            else if (countline == 10)
            {
                sscanf(templine.c_str(), "%*s %*s = %i;", &GeomFilterLineSizeY);
            }
            else if (countline == 11)
            {
                sscanf(templine.c_str(), "%*s %*s = %i;", &GeomFilterChanSizeX);
            }
            else if (countline == 12)
            {
                sscanf(templine.c_str(), "%*s %*s = %i;", &GeomFilterChanSizeY);
            }
            else if (countline == 13)
            {
                sscanf(templine.c_str(), "%*s %*s = %i;", &GeomFilter2DSizeX);
            }
            else if (countline == 14)
            {
                sscanf(templine.c_str(), "%*s %*s = %i;", &GeomFilter2DSizeY);
            }
            else if (countline == 16)
            {
                float tempfloat;
                sscanf(templine.c_str(), "%*s %*s = %f;", &tempfloat);
                maxdisplacement=double(tempfloat);
            }
            else if (countline == 17)
            {
                float tempfloat;
                sscanf(templine.c_str(), "%*s %*s = %f;", &tempfloat);
                dYdX=double(tempfloat);
            }
            else if (countline == 18)
            {
                float tempfloat;
                sscanf(templine.c_str(), "%*s %*s = %f;", &tempfloat);
                dYdX_Line=double(tempfloat);
            }
            else if (countline == 19)
            {
                sscanf(templine.c_str(), "%*s %*s = %i;", &disapearancetime);
            }
            else if (countline == 20)
            {
                float tempfloat;
                sscanf(templine.c_str(), "%*s %*s = %f;", &tempfloat);
                costBirthDeath=double(tempfloat);
            }
            else if (countline == 21)
            {
                sscanf(templine.c_str(), "%*s %*s = %i;", &temperaturedecrease);
            }
            else if (countline == 22)
            {
                float tempfloat;
                sscanf(templine.c_str(), "%*s %*s = %f;", &tempfloat);
                normalizedIntensityWeight=double(tempfloat);
            }
            else if (countline == 24)
            {
                sscanf(templine.c_str(), "%*s %*s = %i;", &VisualOutput);
            }
            else if (countline == 25)
            {
                sscanf(templine.c_str(), "%*s %*s = %i;", &ObjSize);
            }
            else if (countline == 27)
            {
                float tempfloat;
                sscanf(templine.c_str(), "%*s %*s = %f;", &tempfloat);
                SubBackFactorB=double(tempfloat);
            }
            else if (countline == 28)
            {
                sscanf(templine.c_str(), "%*s %*s = %i;", &ncloseB);
            }
            else if (countline == 29)
            {
                sscanf(templine.c_str(), "%*s %*s = %i;", &nopenB);
            }
            else if (countline == 30)
            {
                float tempfloat;
                sscanf(templine.c_str(), "%*s %*s = %f;", &tempfloat);
                threshold_MaskB=double(tempfloat);
            }
            else if (countline == 31)
            {
                sscanf(templine.c_str(), "%*s %*s = %i;", &GeomFilterLineSizeXB);
            }
            else if (countline == 32)
            {
                sscanf(templine.c_str(), "%*s %*s = %i;", &GeomFilterLineSizeYB);
            }
            else if (countline == 33)
            {
                sscanf(templine.c_str(), "%*s %*s = %i;", &GeomFilterChanSizeXB);
            }
            else if (countline == 34)
            {
                sscanf(templine.c_str(), "%*s %*s = %i;", &GeomFilterChanSizeYB);
            }
            else if (countline == 35)
            {
                sscanf(templine.c_str(), "%*s %*s = %i;", &GeomFilter2DSizeXB);
            }
            else if (countline == 36)
            {
                sscanf(templine.c_str(), "%*s %*s = %i;", &GeomFilter2DSizeYB);
            }
            else if (countline == 37)
            {
                sscanf(templine.c_str(), "%*s %*s = %i;", &ObjSizeB);
            }
            countline += 1;
        }

        printf("Scale = %.3f [um/px]\n", scalefactor);
        printf("Time Interval = %.3f [min]\n", TimeInterval);
        printf("Enhance Background Substraction = %.3f\n", SubBackFactor);
        printf("N close = %i\n", nclose);
        printf("N open = %i\n", nopen);
        printf("Threshold Mask = %.3f\n", threshold_Mask);
        printf("Threshold Visual = %.3f\n", threshold_Visual);
        printf("Show Segmented = %i\n", VisualOutput);
        printf("Object Size = %i\n", ObjSize);
        return 0;
    }

    else
    {
        printf("Default Setting\n");
//Deafault
//Data
        scalefactor = 1.3;
        TimeInterval = 3;
//Segmentation
        SubBackFactor = 3;
        nclose = 1;
        nopen = 1;
        threshold_Mask = 3;
        threshold_Visual = 0.05;
        GeomFilterLineSizeX = 2;
        GeomFilterLineSizeY = 1;
        GeomFilterChanSizeX = 2;
        GeomFilterChanSizeY = 1;
        GeomFilter2DSizeX = 3;
        GeomFilter2DSizeY = 3;
//Tracking
        maxdisplacement = 77;
        dYdX = 1;
        dYdX_Line = 12;
        disapearancetime = 3;
        costBirthDeath = 1;
        temperaturedecrease = 5;
        normalizedIntensityWeight = 0.1;
//Output
        VisualOutput = 0;
        ObjSize = 0;

        printf("Scale = %.3f [um/px]\n", scalefactor);
        printf("Time Interval = %.3f [min]\n", TimeInterval);
        printf("Enhance Background Substraction = %.3f\n", SubBackFactor);
        printf("N close = %i\n", nclose);
        printf("N open = %i\n", nopen);
        printf("Threshold Mask = %.3f\n", threshold_Mask);
        printf("Threshold Visual = %.3f\n", threshold_Visual);
        printf("Show Segmented = %i\n", VisualOutput);
        printf("Object Size = %i\n", ObjSize);

        return 1;
    }


}
