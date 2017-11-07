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
#include "../Read_Write_Filter_Paolo/PaFiltering.h"
#include "../Read_Write_Filter_Paolo/PaSaveTiff_Avi.h"
#include "../Read_Write_Filter_Paolo/ReadTiffPa.h"

using namespace cv;
using namespace std;

void ReadTrajectoh(string filenameX, string filenameY, string filenameFr, vector<double>& X, vector<double>& Y, vector<int>& Fr, int h)
{
    int countscanX=0;
    int countscanY=0;
    int countscanFr=0;
    ifstream infileX, infileY, infileFr;
    string templineX, templineY, templineFr;
    string lineX, lineY, lineFr;
    stringstream strstrX, strstrY, strstrFr;

    infileX.open ( filenameX.c_str(), ios::in );
    if(infileX.is_open())
        printf("File: %s\t -- OK\n", filenameX.c_str());
    while (infileX.good())
    {
        float float_x;
        getline(infileX, templineX);
        if (countscanX == 2*h+1)
        {
            strstrX << templineX;
            while( getline(strstrX, lineX, '\t'))
            {
                sscanf(lineX.c_str(), "%f", &float_x);
                //cout << float_x << "\t";
                X.push_back(double(float_x));
                cout << X.back() << "\t";
            }
            break;
        }
        countscanX += 1;
    }
    infileX.close();
    printf("\n");

    infileY.open ( filenameY.c_str(), ios::in );
    if(infileY.is_open())
        printf("File: %s\t -- OK\n", filenameY.c_str());
    while (infileY.good())
    {
        float float_y;
        getline(infileY, templineY);
        if (countscanY == 2*h+1)
        {
            strstrY << templineY;
            while(getline(strstrY, lineY, '\t'))
            {
                sscanf(lineY.c_str(), "%f", &float_y);
                //cout << double(float_y) << "\t";
                Y.push_back(double(float_y));
                cout << Y.back() << "\t";
            }
            break;
        }
        countscanY += 1;
    }
    infileY.close();
    printf("\n");

    infileFr.open ( filenameFr.c_str(), ios::in );
    if(infileFr.is_open())
        printf("File: %s\t -- OK\n", filenameFr.c_str());
    while (infileFr.good())
    {
        int f;
        getline(infileFr, templineFr);
        if (countscanFr == 2*h+1)
        {
            strstrFr << templineFr;
            while(getline(strstrFr, lineFr, '\t'))
            {
                sscanf(lineFr.c_str(), "%d", &f);
                //cout << f << "\t";
                Fr.push_back(f);
                cout << Fr.back() << "\t";
            }
            break;
        }
        countscanFr += 1;
    }
    infileFr.close();
    printf("\n");
}


// MAIN type: 2D, LINE, CHANNEL
int main(int argc, char *argv[])
{
//// Read Data
    string namefile; //video file *.tif
    if (argv[1] == NULL)
    {
        cout << "File name: ";
        cin >> 	namefile;
        char *temp_argv = new char[namefile.size() + 1];
        std::strcpy ( temp_argv, namefile.c_str() );
        argv[1] = temp_argv;
    }
    else
        namefile = string(argv[1]);
    int ntot, width, height;
    vector<Mat> input_frames;
    vector<Mat> output_frames;
    printf("Reading the file: %s\n", namefile.c_str());

    string str_tif (".tif");
    size_t foundtif;
    string dirname=argv[1];
    foundtif=dirname.rfind(str_tif);
    if (foundtif != string::npos)
    {
        printf("Reading the file: %s\n", argv[1]);
        ReadTIFF(argv[1], ntot, width, height, input_frames);
        printf("n: %i\nInitializing\n", ntot);
    }
    else
    {
        printf("Import files from directory: %s\n", argv[1]);
        ImportTIFFdir(argv[1], ntot, width, height, input_frames);
        printf("Files: %i\t Imported: %i\n", ntot, input_frames.size());
    }
    //ReadTIFF(argv[1], ntot, width, height, input_frames);
    cout << "Number of Frames: " << ntot << endl;

    int trajectoID;
    printf("Trajecto ID: ");
    cin >> 	trajectoID;
    double PXscale;
    cout << "Scale [um/px]: ";
    cin >> 	PXscale;

    string nameoutFrame;
    string nameoutX;
    string nameoutY;
    if ( namefile.compare(namefile.size()-1, 1, "/") == 0 )
        namefile.resize(namefile.size()-1);
    size_t lastnamepos=namefile.rfind("/");
    string path = namefile.substr(0, lastnamepos+1);
    string lastname = namefile.substr(lastnamepos+1, namefile.size());
    string strtif=".tif";
    /*printf("path: %s\n", path.c_str());
    printf("lastname: %s\n", lastname.c_str());*/
    string lastnametxt=lastname;
    nameoutFrame=path;
    nameoutX=path;
    nameoutY=path;
    string stroutFrame = "Frame_Trajectos_";
    string stroutX = "X_Trajectos_";
    string stroutY = "Y_Trajectos_";
    if (foundtif != string::npos)
        lastnametxt.replace(lastname.find(strtif), strtif.length(), ".txt");
    else
        lastnametxt.append(".txt");
    stroutFrame.append(lastnametxt);
    stroutX.append(lastnametxt);
    stroutY.append(lastnametxt);
    nameoutFrame.append(stroutFrame);
    nameoutX.append(stroutX);
    nameoutY.append(stroutY);


    vector <double> X;
    vector <double> Y;
    vector <int> Fr; //frame
    ReadTrajectoh(nameoutX, nameoutY, nameoutFrame, X, Y, Fr, trajectoID);
    printf("Reading --END--\n");
    if (Fr.size() == 0)
    {
        printf("The trajectory is not in the file!\n");
        exit(0);
    }

    for (int i=0; i<Fr.size(); i++)
    {
        printf("Frame:%i\t X[px]: %.1f\t Y[px]: %.1f\n", Fr[i], X[i]/PXscale, height-Y[i]/PXscale);
        Mat img;
        img = input_frames[Fr[i]];
        normalize(img, img, 0, 255, NORM_MINMAX, -1);
        img.convertTo(img, CV_8UC1);
        //cvtColor(img, img, CV_GRAY2RGB);
        Mat track( img.rows, img.cols, CV_8UC1);
        Mat out( img.rows, img.cols, CV_8UC3);
        //circle(img, Point(X[i]/PXscale,height-Y[i]/PXscale), 1, Scalar(0,0,255), -1, 8, 0);
        int countpoint =0;
        if(i>0)
        {
            while(countpoint < 25 && i-countpoint-1 > 0)
            {
                line(track, Point(X[i-countpoint]/PXscale, height-Y[i-countpoint]/PXscale), Point(X[i-countpoint-1]/PXscale, height-Y[i-countpoint-1]/PXscale), 255, 2);
                /*printf("countpoint: %i\n", countpoint);
                fflush(stdout);*/
                countpoint += 1;
            }
        }
        int from_to1[] = { 0,1};
        mixChannels( &img, 1, &out, 1, from_to1, 1 );
        int from_to2[] = { 0,2};
        mixChannels( &track, 1, &out, 1, from_to2, 1 );
        output_frames.push_back(out);
        imshow("Trajecto", out);
        waitKey(100);
    }

    printf("Save output video\n");
    string nameout;
    nameout = namefile;
    std::string snum;
    std::stringstream snumout;
    snumout << trajectoID;
    snum = snumout.str();
    if (foundtif != string::npos)
        nameout.replace(nameout.find(strtif), strtif.length(), "_Trajecto_"+ snum + ".avi");
    else
        nameout.append("_Trajecto_"+ snum + ".avi");
    Save8bitcolor(nameout, output_frames);
    printf("output saved at %s\n", nameout.c_str());

    exit(0);
}

