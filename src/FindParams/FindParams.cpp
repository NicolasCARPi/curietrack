#include <opencv2/opencv.hpp>
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
#include "../CellTracking/TrackingFunctions.h"
#include "../Read_Write_Filter_Paolo/PaFiltering.h"
#include "../Read_Write_Filter_Paolo/PaSaveTiff_Avi.h"
#include "../Read_Write_Filter_Paolo/ReadTiffPa.h"

using namespace cv;
using namespace std;

struct Params_Filter_Segment
{
    // Params!
    //Data
    double scalefactor;
    double TimeInterval; //min;
    //Segmentation
    double SubBackFactor;
    double SubBackFactor_Int;

    int nclose;
    int nopen;
    double threshold_Mask;
    double threshold_Mask_Int;
    double threshold_Visual;
    double threshold_Visual_Int;
    int GeomFilterLineSizeX;
    int GeomFilterLineSizeY;
    int GeomFilterChanSizeX;
    int GeomFilterChanSizeY;
    int GeomFilter2DSizeX;
    int GeomFilter2DSizeY;
    //Tracking
    double maxdisplacement; // [px!!!!]
    double dYdX;
    double dYdX_Line;
    int disapearancetime;
    double costBirthDeath;
    int temperaturedecrease;
    double normalizedIntensityWeight;
    //output
    int VisualOutput;

    //Section 2
    int ntot;
    vector<Mat> input_frames;
    //vector<Mat> input_frames_scaled8bit;

    Mat CumHist;
    Mat CumHistInt;
    int MaxAll;
    int StdevAll3;
    int tamagno;


    int width;
    int height;
    char *argv[6];
};


void Filter_Segment(Params_Filter_Segment& params)
{
    int ntot = params.ntot;
    vector<Mat> input_filtered_frames;
    input_filtered_frames=params.input_frames;
    //cout << mean(input_filtered_frames[0])[0] << endl;
    /*imshow("input", input_filtered_frames[0]);
    waitKey();*/
    vector<Mat> input_scaled8bit_enhanced_frames(ntot);
    //input_scaled8bit_enhanced_frames=params.input_frames_scaled8bit;

    Mat CumHist, CumHistInt;
    CumHist=params.CumHist;
    CumHistInt=params.CumHistInt;
    int MaxAll=params.MaxAll;
    int StdevAll3=params.StdevAll3;
    int tamagno = params.tamagno;
    int width= params.width;
    int height=params.height;
    char** argv;
    argv=params.argv;
    double scalefactor=params.scalefactor;
    double TimeInterval=params.TimeInterval; //min;
    //Segmentation
    double SubBackFactor=params.SubBackFactor;
    int nclose=params.nclose;
    int nopen=params.nopen;
    double threshold_Mask=params.threshold_Mask;
    double threshold_Visual=params.threshold_Visual;
    int GeomFilterLineSizeX=params.GeomFilterLineSizeX;
    int GeomFilterLineSizeY=params.GeomFilterLineSizeY;
    int GeomFilterChanSizeX=params.GeomFilterChanSizeX;
    int GeomFilterChanSizeY=params.GeomFilterChanSizeY;
    int GeomFilter2DSizeX=params.GeomFilter2DSizeX;
    int GeomFilter2DSizeY=params.GeomFilter2DSizeY;
    //Tracking
    double maxdisplacement=params.maxdisplacement; // [px!!!!]
    double dYdX=params.dYdX;
    double dYdX_Line=params.dYdX_Line;
    int disapearancetime=params.disapearancetime;
    double costBirthDeath=params.costBirthDeath;
    int temperaturedecrease=params.temperaturedecrease;
    double normalizedIntensityWeight=params.normalizedIntensityWeight;
    //output
    int VisualOutput=params.VisualOutput;

    printf("Scale = %f [um/px]\n", scalefactor);
    printf("Time Interval = %f [min]\n", TimeInterval);
    printf("\tEnhance Background Substraction = %f\n", SubBackFactor);
    printf("\tN close = %i\n", nclose);
    printf("\tN open = %i\n", nopen);
    printf("\tThreshold Mask = %f\n", threshold_Mask);
    printf("\tThreshold Visual = %f\n", threshold_Visual);
    printf("Show Segmented = %i\n", VisualOutput);



    printf("8bit Enhance\n");
    //double SubBackFactor=3.0;
    int central=round(ntot/2);

    Enhance8bit(input_filtered_frames, ntot, input_scaled8bit_enhanced_frames, CumHistInt, tamagno, SubBackFactor);
    if (width>800 || height>800)
    {
        Mat tempshow=Mat(input_filtered_frames[central],Rect(0,0,min(width,800),min(height,800)));
        imshow("Filtered", tempshow);
        //waitKey(1500);
    }
    else
    {
        imshow("Filtered", input_filtered_frames[central]);
        //waitKey(1500);
    }
    //imshow("Filtered", input_filtered_frames[central] );
    CumHist.release();
    CumHistInt.release();
    printf("Cumulative Histogram 8bit 2 of 2\n");
    CumHistogramInt_ThrGausDecay8bit(input_scaled8bit_enhanced_frames, ntot, CumHist, CumHistInt, MaxAll, StdevAll3);
    if (width>800 || height>800)
    {
        Mat tempshow=Mat(input_scaled8bit_enhanced_frames[central],Rect(0,0,min(width,800),min(height,800)));
        imshow("Enhanced", tempshow);
        //waitKey(1500);
    }
    else
    {
        imshow("Enhanced", input_scaled8bit_enhanced_frames[central]);
        //waitKey(1500);
    }
    //imshow("Enhanced", input_scaled8bit_enhanced_frames[central] );
    printf("Segmentation Start\n", ntot);


    vector<vector<Point> > CellContours;

    vector< vector<Roi> > RoisFrames(ntot); // RoisFrames[i] vector of all rois in frame i
    vector< vector<Roi> > Trajectos; // Trajecto[i] trajectpry i as a vector of roi

    vector<Mat> output_segmented(ntot);
    for (int xx=0; xx<ntot; xx++)
        output_segmented[xx] = Mat::zeros(height, width, CV_16UC1);
    //cout << argv[1] << endl;
    //cout << argv[2] << endl;
    ParallelSegmentCells(input_filtered_frames, input_scaled8bit_enhanced_frames, output_segmented, ntot, tamagno, MaxAll, StdevAll3, RoisFrames, argv, 1, nclose, nopen, threshold_Mask, threshold_Visual, GeomFilterLineSizeX, GeomFilterLineSizeY, GeomFilterChanSizeX, GeomFilterChanSizeY, GeomFilter2DSizeX, GeomFilter2DSizeY, 0);

    vector<Mat> output_frames_color(ntot);
    GenerateColorOut(input_filtered_frames, output_frames_color);


    DrawRoisFrame(output_frames_color, ntot, RoisFrames);
    //cout << "Update Images!" << endl;

    if (width>800 || height>800)
    {
        Mat tempshow=Mat(output_segmented[central],Rect(0,0,min(width,800),min(height,800)));
        imshow("Segmented", tempshow);
        //waitKey(1500);
    }
    else
    {
        imshow("Segmented", output_segmented[central]);
        //waitKey(1500);
    }
    //imshow("Segmented", output_segmented[central] );
    //cout << "Images Updated" << endl;
    //waitKey();

    /*FreeVectorMat(input_scaled8bit_enhanced_frames);
    FreeVectorMat(input_filtered_frames);*/
}

void Filter_Segment_SBF(int SBF, void* p)
{
    Params_Filter_Segment *params_temp = (Params_Filter_Segment *)p;
    Params_Filter_Segment params = *params_temp;
    params.SubBackFactor_Int=SBF;
    params.SubBackFactor=double(SBF)/10;
    int ThM = getTrackbarPos("threshold_Mask", "Filtered");
    params.threshold_Mask=double(ThM)/10;
    //int ThV = getTrackbarPos("threshold_Visual", "Filtered");
    //params.threshold_Visual=double(ThV)/100;
    int nc = getTrackbarPos("nclose", "Filtered");
    params.nclose=nc;
    int no = getTrackbarPos("nopen", "Filtered");
    params.nopen=no;
    printf("New SBF = %f\n", params.SubBackFactor);
    Filter_Segment(params);
    p=(void *)&params;
}
void Filter_Segment_ThM(int ThM, void* p)
{
    Params_Filter_Segment *params_temp = (Params_Filter_Segment *)p;
    Params_Filter_Segment params = *params_temp;
    params.threshold_Mask=double(ThM)/10;
    int SBF = getTrackbarPos("SubBackFactor", "Filtered");
    params.SubBackFactor=double(SBF)/10;
    //int ThV = getTrackbarPos("threshold_Visual", "Filtered");
    //params.threshold_Visual=double(ThV)/100;
    int nc = getTrackbarPos("nclose", "Filtered");
    params.nclose=nc;
    int no = getTrackbarPos("nopen", "Filtered");
    params.nopen=no;
    printf("New ThM = %f\n", params.threshold_Mask);
    Filter_Segment(params);
    p=(void *)&params;
}
/*void Filter_Segment_ThV(int ThV, void* p)
{
    Params_Filter_Segment *params_temp = (Params_Filter_Segment *)p;
    Params_Filter_Segment params = *params_temp;
    params.threshold_Visual=double(ThV)/100;
    int SBF = getTrackbarPos("SubBackFactor", "Filtered");
    params.SubBackFactor=double(SBF)/10;
    int ThM = getTrackbarPos("threshold_Mask", "Filtered");
    params.threshold_Mask=double(ThM)/10;
    int nc = getTrackbarPos("nclose", "Filtered");
    params.nclose=nc;
    int no = getTrackbarPos("nopen", "Filtered");
    params.nopen=no;
    //printf("New ThV = %f\n", params.threshold_Visual);
    Filter_Segment(params);
    p=(void *)&params;
}*/
void Filter_Segment_nc(int nc, void* p)
{
    Params_Filter_Segment *params_temp = (Params_Filter_Segment *)p;
    Params_Filter_Segment params = *params_temp;
    params.nclose=nc;
    int SBF = getTrackbarPos("SubBackFactor", "Filtered");
    params.SubBackFactor=double(SBF)/10;
    int ThM = getTrackbarPos("threshold_Mask", "Filtered");
    params.threshold_Mask=double(ThM)/10;
    //int ThV = getTrackbarPos("threshold_Visual", "Filtered");
    //params.threshold_Visual=double(ThV)/100;
    int no = getTrackbarPos("nopen", "Filtered");
    params.nopen=no;
    Filter_Segment(params);
    p=(void *)&params;
}
void Filter_Segment_no(int no, void* p)
{
    Params_Filter_Segment *params_temp = (Params_Filter_Segment *)p;
    Params_Filter_Segment params = *params_temp;
    params.nopen=no;
    int SBF = getTrackbarPos("SubBackFactor", "Filtered");
    params.SubBackFactor=double(SBF)/10;
    int ThM = getTrackbarPos("threshold_Mask", "Filtered");
    params.threshold_Mask=double(ThM)/10;
    //int ThV = getTrackbarPos("threshold_Visual", "Filtered");
    //params.threshold_Visual=double(ThV)/100;
    int nc = getTrackbarPos("nclose", "Filtered");
    params.nclose=nc;
    Filter_Segment(params);
    p=(void *)&params;
}


int main(int argc, char *argv[])
{

    // Params!
    //Data
    double scalefactor = 1.3;
    double TimeInterval = 3; //min;
    //Segmentation
    double SubBackFactor=0;
    int nclose=0;
    int nopen=0;
    double threshold_Mask=0;
    double threshold_Visual=0;
    int GeomFilterLineSizeX=2;
    int GeomFilterLineSizeY=1;
    int GeomFilterChanSizeX=2;
    int GeomFilterChanSizeY=1;
    int GeomFilter2DSizeX=3;
    int GeomFilter2DSizeY=3;
    //Tracking
    double maxdisplacement=77; // [px!!!!]
    double dYdX=1;
    double dYdX_Line=12;
    int disapearancetime=3;
    double costBirthDeath=1;
    int temperaturedecrease = 5;
    double normalizedIntensityWeight = 0.1;
    //output
    int VisualOutput = 1;
    int SaveSingleVisualOutputAs = 0; /// 0 -> avi; 1 -> tif
    int ObjSize =0;
    int readconfigout=ReadConfig(argv[6], scalefactor, TimeInterval, SubBackFactor, nclose, nopen, threshold_Mask, threshold_Visual, GeomFilterLineSizeX, GeomFilterLineSizeY, GeomFilterChanSizeX, GeomFilterChanSizeY, GeomFilter2DSizeX, GeomFilter2DSizeY, maxdisplacement, dYdX, dYdX_Line, disapearancetime, costBirthDeath, temperaturedecrease, normalizedIntensityWeight, VisualOutput, ObjSize, SaveSingleVisualOutputAs);
    int SubBackFactor_Int = round(SubBackFactor*10);
    int threshold_Mask_Int = round(threshold_Mask*10);
    int threshold_Visual_Int = round(threshold_Visual*100);
    Params_Filter_Segment params;
    params.scalefactor=scalefactor;
    params.TimeInterval=TimeInterval;
    params.SubBackFactor=SubBackFactor;
    params.SubBackFactor_Int=SubBackFactor_Int;
    params.nclose=nclose;
    params.nopen=nopen;
    params.threshold_Mask=threshold_Mask;
    params.threshold_Mask_Int=threshold_Mask_Int;
    params.threshold_Visual=threshold_Visual;
    params.threshold_Visual_Int=threshold_Visual_Int;
    params.GeomFilterLineSizeX=GeomFilterLineSizeX;
    params.GeomFilterLineSizeY=GeomFilterLineSizeY;
    params.GeomFilterChanSizeX=GeomFilterChanSizeX;
    params.GeomFilterChanSizeY=GeomFilterChanSizeY;
    params.GeomFilter2DSizeX=GeomFilter2DSizeX;
    params.GeomFilter2DSizeY=GeomFilter2DSizeY;
    params.maxdisplacement=maxdisplacement;
    params.dYdX=dYdX;
    params.dYdX_Line=dYdX_Line;
    params.disapearancetime=disapearancetime;
    params.costBirthDeath=costBirthDeath;
    params.temperaturedecrease=temperaturedecrease;
    params.normalizedIntensityWeight=normalizedIntensityWeight;
    params.VisualOutput=VisualOutput;

    //printf("type: %s\n", string(argv[2]).c_str());
    /*if (string(argv[3]) != "2D" && string(argv[3]) != "LINE" && string(argv[3]) != "CHANNEL")
    {
        printf("Not supported type!\n");
        exit(0);
    }*/
    int ntot, width, height;
    //string CellTypeName=argv[4];
    vector<Mat> input_frames;
    string str_tif (".tif");
    size_t foundtif;
    string dirname=argv[1];
    foundtif=dirname.rfind(str_tif);
    if (foundtif != string::npos)
    {
        printf("Reading the file: %s\n", argv[1]);
        ReadTIFF(argv[1], ntot, width, height, input_frames);
    }
    else
    {
        printf("Import files from directory: %s\n", argv[1]);
        ImportTIFFdir(argv[1], ntot, width, height, input_frames);
        printf("Files: %i\t Imported: %i\n", ntot, input_frames.size());
    }

    vector<Mat> input_frames_sub;
    int subsel=0;
    int newntot=0;
    for(int i=0; i<ntot; i++)
        if (i==subsel)
        {
            Mat framex;
            input_frames[i].copyTo(framex);
            //input_frames_sub.push_back(input_frames[i]);
            input_frames_sub.push_back(framex);
            /*imshow("imput", framex);
            waitKey(500);*/
            subsel += 1;
            newntot += 1;
        }

    FreeVectorMat(input_frames);
    namedWindow("Filtered", 1);
    printf("n: %i\nInitializing\n", newntot);
    vector<Mat> input_scaled8bit_frames(newntot);

    /*vector<Mat> null_gray(ntot);
    for (int xx=0; xx<ntot; xx++)
    	null_gray[xx] = Mat::zeros(height, width, CV_16UC1);*/
    int MaxAll, StdevAll3;
    Mat CumHist, CumHistInt;
    printf("Cumulative Histogram 16bit\n");
    CumHistogramInt_ThrGausDecay16bit(input_frames_sub, newntot, CumHist, CumHistInt, MaxAll, StdevAll3);
    if (width>800 || height>800)
    {
        Mat tempshow=Mat(input_frames_sub[round(newntot/2)],Rect(0,0,min(width,800),min(height,800)));
        imshow("Input", tempshow);
        waitKey(2000);
    }
    else
    {
        imshow("Input", input_frames_sub[round(newntot/2)]);
        waitKey(2000);
    }
    printf("16bit Enhance\n");
    Enhance16bit(input_frames_sub, newntot, input_frames_sub, MaxAll, StdevAll3); //substract at MaxAll value
    if (width>800 || height>800)
    {
        Mat tempshow=Mat(input_frames_sub[round(newntot/2)],Rect(0,0,min(width,800),min(height,800)));
        imshow("Input", tempshow);
        waitKey(2000);
    }
    else
    {
        imshow("Input", input_frames_sub[round(newntot/2)]);
        waitKey(2000);
    }
    int MinStack, MaxStack;
    printf("Intensity Scale and 8bit conversion\n");
    Scale_Convert8bit(input_frames_sub, newntot, input_scaled8bit_frames, MinStack, MaxStack);
    if (width>800 || height>800)
    {
        Mat tempshow=Mat(input_scaled8bit_frames[round(newntot/2)],Rect(0,0,min(width,800),min(height,800)));
        imshow("Input", tempshow);
        waitKey(2000);
    }
    else
    {
        imshow("Input", input_scaled8bit_frames[round(newntot/2)]);
        waitKey(2000);
    }
    FreeVectorMat(input_frames_sub);
    //cout << "Input Frame 16bit clear!" << endl;
    CumHist.release();
    CumHistInt.release();
    printf("Cumulative Histogram and Threshold\n");
    CumHistogramInt_ThrGausDecay8bit(input_scaled8bit_frames, newntot, CumHist, CumHistInt, MaxAll, StdevAll3);

    vector<Mat> input_scaled8bit_enhanced_frames(newntot);
    vector<Mat> input_filtered_frames(newntot);
    int tamagno = ObjSize; //cell dimension in pixel have to be !!! odd !!!
    if (tamagno == 0)
        tamagno = DefiSizeObjsErodeAll3(input_scaled8bit_frames, newntot, MaxAll, StdevAll3);
    else if (GSL_IS_EVEN(tamagno))
        {
            tamagno += 1;
            cout << "Correct Object size: " << tamagno << endl;
        }
    printf("Filter 8bit\n");
    Filter_8bit1(input_scaled8bit_frames, newntot, input_filtered_frames, tamagno);
    if (width>800 || height>800)
    {
        Mat tempshow=Mat(input_filtered_frames[round(newntot/2)],Rect(0,0,min(width,800),min(height,800)));
        imshow("Input", tempshow);
        waitKey(2000);
    }
    else
    {
        imshow("Input", input_filtered_frames[round(newntot/2)]);
        waitKey(2000);
    }
    FreeVectorMat(input_scaled8bit_frames);
    CumHist.release();
    CumHistInt.release();
    printf("Cumulative Histogram 8bit 1 of 2\n");
    CumHistogramInt_ThrGausDecay8bit(input_filtered_frames, newntot, CumHist, CumHistInt, MaxAll, StdevAll3);


    if (width>800 || height>800)
    {
        Mat tempshow=Mat(input_filtered_frames[round(newntot/2)],Rect(0,0,min(width,800),min(height,800)));
        imshow("Filtered", tempshow);
        waitKey(2000);
    }
    else
    {
        imshow("Filtered", input_filtered_frames[round(newntot/2)]);
        waitKey(2000);
    }


    //cout << mean(input_filtered_frames[0])[0] << endl;

    params.input_frames=input_filtered_frames;
    //params.input_frames_scaled8bit=input_scaled8bit_frames;
    params.CumHist=CumHist;
    params.CumHistInt=CumHistInt;
    params.MaxAll=MaxAll;
    params.StdevAll3=StdevAll3;
    params.tamagno=tamagno;

    params.argv[0]=argv[0];
    params.argv[1]=argv[1];
    params.argv[2]=argv[2];
    params.argv[3]=argv[3];
    params.argv[4]=argv[4];
    params.argv[5]=argv[5];
    params.argv[6]=argv[6];
    /*cout << "Test argv" << endl;
    cout << argv[1] << endl;
    cout << params.argv[1]<< endl;*/

    params.ntot=newntot;
    //params.input_frames=input_frames_sub;
    params.width=width;
    params.height=height;
    void *ParamsVoid;
    ParamsVoid=(void *)&params;
    int SBF=SubBackFactor_Int; // -> 30/10
    int ThM=threshold_Mask_Int; // -> 10/10
    int ThV=threshold_Visual_Int; // -> 5/100
    int nc=nclose;
    int no=nopen;
    createTrackbar("SubBackFactor", "Filtered", &SBF, 300, Filter_Segment_SBF, ParamsVoid);
    createTrackbar("threshold_Mask", "Filtered", &ThM, 300, Filter_Segment_ThM, ParamsVoid);
    //createTrackbar("threshold_Visual", "Filtered", &ThV, 10, Filter_Segment_ThV, ParamsVoid);
    createTrackbar("nclose", "Filtered", &nc, 5, Filter_Segment_nc, ParamsVoid);
    createTrackbar("nopen", "Filtered", &no, 5, Filter_Segment_no, ParamsVoid);
    //waitKey();
    //Filter_Segment(params, input_frames_sub, newntot, width, height, argv);
    Filter_Segment_SBF(SBF ,ParamsVoid);
    waitKey();


    /*namedWindow("Filtered", 1);
    namedWindow("Enhanced", 1);
    namedWindow("Segmented", 1);*/
    return 0;


}
