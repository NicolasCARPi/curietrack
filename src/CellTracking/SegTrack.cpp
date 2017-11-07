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
#include "./TrackingFunctions.h"
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MAIN type: 2D, 1D, LINE, CHANNEL
int main(int argc, char *argv[])
{

    // Params!
    //Data
    double scalefactor = 1;
    double TimeInterval = 3; //min;
    //Segmentation
    double SubBackFactor=3.0;
    int nclose=1;
    int nopen=1;
    double threshold_Mask=1;
    double threshold_Visual=0.05;
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
    double costBirthDeath=0;
    int temperaturedecrease = 5;
    double normalizedIntensityWeight = 0.1;
    //output
    int VisualOutput = 0;
    int SaveSingleVisualOutputAs = 0; /// 0 -> avi; 1 -> tif
    int ObjSize =0;
    // Second Channel
    double SubBackFactorB = 2;
    int ncloseB = 1;
    int nopenB = 1;
    double threshold_MaskB = 3;
    int GeomFilterLineSizeXB = 5;
    int GeomFilterLineSizeYB = 3;
    int GeomFilterChanSizeXB = 3;
    int GeomFilterChanSizeYB = 3;
    int GeomFilter2DSizeXB = 5;
    int GeomFilter2DSizeYB = 5;
    int ObjSizeB = 31;
    int readconfigout;
    //if (string(argv[3]) == "0")
        readconfigout=ReadConfig(argv[6], scalefactor, TimeInterval, SubBackFactor, nclose, nopen, threshold_Mask, threshold_Visual, GeomFilterLineSizeX, GeomFilterLineSizeY, GeomFilterChanSizeX, GeomFilterChanSizeY, GeomFilter2DSizeX, GeomFilter2DSizeY, maxdisplacement, dYdX, dYdX_Line, disapearancetime, costBirthDeath, temperaturedecrease, normalizedIntensityWeight, VisualOutput, ObjSize, SaveSingleVisualOutputAs);
    //else
    //    readconfigout=ReadConfig2Channels(argv[6], scalefactor, TimeInterval, SubBackFactor, nclose, nopen, threshold_Mask, threshold_Visual, GeomFilterLineSizeX, GeomFilterLineSizeY, GeomFilterChanSizeX, GeomFilterChanSizeY, GeomFilter2DSizeX, GeomFilter2DSizeY, maxdisplacement, dYdX, dYdX_Line, disapearancetime, costBirthDeath, temperaturedecrease, normalizedIntensityWeight, VisualOutput, ObjSize, SubBackFactorB, ncloseB, nopenB, threshold_MaskB, GeomFilterLineSizeXB, GeomFilterLineSizeYB, GeomFilterChanSizeXB, GeomFilterChanSizeYB, GeomFilter2DSizeXB, GeomFilter2DSizeYB, ObjSizeB);
    if (string(argv[4]) != "2D" && string(argv[4]) != "LINE" && string(argv[4]) != "1D" && string(argv[4]) != "CHANNEL")
    {
        printf("Not supported type!\n");
        exit(0);
    }
    else
    {
        printf("Tracking Type: %s\n", argv[4]);
    }
    int ntot, width, height;
    string CellTypeName=argv[5];
    vector<Mat> input_frames;

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
        printf("Files: %i\t Imported: %lu\n", ntot, input_frames.size());
    }

    vector<Mat> input_scaled8bit_frames(ntot);
    int MaxAll, StdevAll3;
    Mat CumHist, CumHistInt;
    printf("Cumulative Histogram 16bit\n");
    CumHistogramInt_ThrGausDecay16bit(input_frames, ntot, CumHist, CumHistInt, MaxAll, StdevAll3);
    printf("16bit Enhance\n");
    Enhance16bit(input_frames, ntot, input_frames, MaxAll, StdevAll3); //substract at MaxAll value
    int MinStack, MaxStack;
    printf("Intansity Scale and 8bit conversion\n");
    Scale_Convert8bit(input_frames, ntot, input_scaled8bit_frames, MinStack, MaxStack);
    FreeVectorMat(input_frames);
    //cout << "Input Frame 16bit clear!" << endl;
    CumHist.release();
    CumHistInt.release();
    printf("Cumulative Histogram and Threshold\n");
    CumHistogramInt_ThrGausDecay8bit(input_scaled8bit_frames, ntot, CumHist, CumHistInt, MaxAll, StdevAll3);

    vector<Mat> input_scaled8bit_enhanced_frames(ntot);
    vector<Mat> input_filtered_frames(ntot);
    int tamagno = ObjSize; //cell dimension in pixel have to be !!! odd !!!
    if (tamagno == 0)
        tamagno = DefiSizeObjsErodeAll3(input_scaled8bit_frames, ntot, MaxAll, StdevAll3);
    else if (GSL_IS_EVEN(tamagno))
    {
        tamagno += 1;
        cout << "Correct Object size: " << tamagno << endl;
    }
    printf("Filter 8bit\n");
    Filter_8bit1(input_scaled8bit_frames, ntot, input_filtered_frames, tamagno);
    FreeVectorMat(input_scaled8bit_frames);
    CumHist.release();
    CumHistInt.release();
    printf("Cumulative Histogram 8bit 1 of 2\n");
    CumHistogramInt_ThrGausDecay8bit(input_filtered_frames, ntot, CumHist, CumHistInt, MaxAll, StdevAll3);
    printf("8bit Enhance\n");
    Enhance8bit(input_filtered_frames, ntot, input_scaled8bit_enhanced_frames, CumHistInt, tamagno, SubBackFactor);
    CumHist.release();
    CumHistInt.release();
    printf("Cumulative Histogram 8bit 2 of 2\n");
    CumHistogramInt_ThrGausDecay8bit(input_scaled8bit_enhanced_frames, ntot, CumHist, CumHistInt, MaxAll, StdevAll3);

    printf("Segmentation Start\n");

    vector<vector<Point> > CellContours;

    vector< vector<Roi> > RoisFrames(ntot); // RoisFrames[i] vector of all rois in frame i
    vector< vector<Roi> > Trajectos; // Trajecto[i] trajectpry i as a vector of roi

    vector<Mat> output_segmented(ntot);
    for (int xx=0; xx<ntot; xx++)
        output_segmented[xx] = Mat::zeros(height, width, CV_16UC1);

    ParallelSegmentCells(input_filtered_frames, input_scaled8bit_enhanced_frames, output_segmented, ntot, tamagno, MaxAll, StdevAll3, RoisFrames, argv, 1, nclose, nopen, threshold_Mask, threshold_Visual, GeomFilterLineSizeX, GeomFilterLineSizeY, GeomFilterChanSizeX, GeomFilterChanSizeY, GeomFilter2DSizeX, GeomFilter2DSizeY, VisualOutput);
    FreeVectorMat(input_scaled8bit_enhanced_frames);
    vector<Mat> output_frames_color(ntot);
    if (string(argv[2]) == "0" && string(argv[3]) == "0")
        GenerateColorOut(input_filtered_frames, output_frames_color); /// To have grey nuclei
    else
    {
        GenerateColorOutMixChannels(input_filtered_frames, output_frames_color); /// To have blue nuclei
        /// read transmitted ligth if possible
        vector<Mat> input_tras;
        string nametrans;
        vector<Mat> input_red;
        string namered;
        if (string(argv[2]) != "0" && string(argv[3]) == "0") /// To have blue nuclei and gray trans
        {
            printf("\nExtract also trasmitted: %s\n", argv[2]);
            nametrans = string(argv[2]);
            foundtif=nametrans.rfind(str_tif);
            if (foundtif != string::npos)
            {
                ReadTIFF(argv[2], ntot, width, height, input_tras);
            }
            else
            {
                ImportTIFFdir(argv[2], ntot, width, height, input_tras);
            }
            for (int z=0; z<ntot; z++)
            {
                Mat temp_trans;
                temp_trans = input_tras[z];
                normalize(temp_trans, temp_trans, 0, 255, NORM_MINMAX, -1);
                temp_trans.convertTo(temp_trans, CV_8UC1);
                cvtColor(temp_trans, temp_trans, CV_GRAY2RGB); /// To have grey trans
                //Mat temp_trans_col(temp_trans.rows, temp_trans.cols, CV_8UC3); /// To have red trans
                //int from_to_trans[] = {0,1};
                //mixChannels(&temp_trans, 1, &temp_trans_col, 1, from_to_trans, 1);
                input_tras[z] = temp_trans ; /// To have grey trans
                //input_tras[z] = temp_trans_col ; /// To have red trans
                output_frames_color[z] = output_frames_color[z] + temp_trans;
            }
        }
        else if (string(argv[3]) != "0") /// To have blue nuclei, green trans and red red
        {
            printf("\nExtract also a second channel: %s\n", argv[3]);
            namered = string(argv[3]);
            foundtif=namered.rfind(str_tif);
            if (foundtif != string::npos)
            {
                ReadTIFF(argv[3], ntot, width, height, input_red);
            }
            else
            {
                ImportTIFFdir(argv[3], ntot, width, height, input_red);
            }
            if (string(argv[2]) != "0") /// To have blue nuclei, green trans and red red
            {
                printf("\nExtract also trasmitted\n");
                nametrans = string(argv[2]);
                foundtif=nametrans.rfind(str_tif);
                if (foundtif != string::npos)
                {
                    ReadTIFF(argv[2], ntot, width, height, input_tras);
                }
                else
                {
                    ImportTIFFdir(argv[2], ntot, width, height, input_tras);
                }
            }
            for (int z=0; z<ntot; z++)
            {
                Mat temp_red, temp_trans;
                temp_red = input_red[z];
                normalize(temp_red, temp_red, 0, 255, NORM_MINMAX, -1);
                temp_red.convertTo(temp_red, CV_8UC1);
                Mat temp_red_col, temp_trans_col;
                temp_red_col = Mat::zeros(temp_red.rows, temp_red.cols, CV_8UC3);/// To have red
                temp_trans_col = Mat::zeros(temp_red.rows, temp_red.cols, CV_8UC3);/// To have red
                int from_to_red[] = {0,2};
                mixChannels(&temp_red, 1, &temp_red_col, 1, from_to_red, 1);
                if (string(argv[2]) != "0") /// To have blue nuclei, green trans and red red
                {
                    temp_trans = input_tras[z];
                    normalize(temp_trans, temp_trans, 0, 255, NORM_MINMAX, -1);
                    temp_trans.convertTo(temp_trans, CV_8UC1);
                    int from_to_trans[] = {0,1};
                    mixChannels(&temp_trans, 1, &temp_trans_col, 1, from_to_trans, 1);
                    output_frames_color[z] = output_frames_color[z] + temp_trans_col + temp_red_col;
                }
                else
                    output_frames_color[z] = output_frames_color[z] + temp_red_col;
            }

        }
        FreeVectorMat(input_tras);
        FreeVectorMat(input_red);
    }
    FreeVectorMat(input_filtered_frames);
    DrawRoisFrame(output_frames_color, ntot, RoisFrames);

    if (string(argv[4]) == "LINE" || string(argv[4]) == "CHANNEL" || string(argv[4]) == "1D")
    {
        dYdX=dYdX_Line;
    }
    SerialBuildTrajecto(argv, RoisFrames, Trajectos, ntot, maxdisplacement, dYdX, disapearancetime, costBirthDeath, height, width, temperaturedecrease, normalizedIntensityWeight);
    string name=string(argv[1]);
    foundtif=name.rfind(str_tif);
    if (foundtif == string::npos)
    {
        string newnameout = name;
        if ( newnameout.compare(newnameout.size()-1, 1, "/") == 0 ) newnameout.resize(newnameout.size()-1);
        string delim="/";
        size_t founddelim;
        founddelim=newnameout.rfind(delim);
        if (founddelim != string::npos)
            newnameout = newnameout.substr(founddelim);
        newnameout.append(".tif");
        if ( name.compare(name.size()-1, 1, "/") != 0 ) name.append("/");
        name.append(newnameout);
    }
    size_t lastnamepos=name.rfind("/");
    string path = name.substr(0, lastnamepos+1);
    string lastname = name.substr(lastnamepos+1, name.size());




    if (Trajectos.size()>0)
    {
        string lastnametxt=lastname;
        string nametrajecto=lastname;
        string nameoutFrame=path;
        string nameoutT=path;
        string nameoutX=path;
        string nameoutY=path;
        string nameoutA=path;
        string nameoutP=path;
        string nameoutC=path;
        string nameoutTrajectos=path;
        string stroutFrame = "Frame_Trajectos_";
        string stroutT = "T_Trajectos_";
        string stroutX = "X_Trajectos_";
        string stroutY = "Y_Trajectos_";
        string stroutA = "A_Trajectos_";
        string stroutP = "P_Trajectos_";
        string stroutC = "C_Trajectos_";
        string stroutTrajectos = lastname;

        lastnametxt.replace(lastname.find(str_tif), str_tif.length(), ".txt");
        nametrajecto.replace(lastname.find(str_tif), str_tif.length(), "");
        stroutFrame.append(lastnametxt);
        stroutT.append(lastnametxt);
        stroutX.append(lastnametxt);
        stroutY.append(lastnametxt);
        stroutA.append(lastnametxt);
        stroutP.append(lastnametxt);
        stroutC.append(lastnametxt);
        nameoutFrame.append(stroutFrame);
        nameoutT.append(stroutT);
        nameoutX.append(stroutX);
        nameoutY.append(stroutY);
        nameoutA.append(stroutA);
        nameoutP.append(stroutP);
        nameoutC.append(stroutC);

        stroutTrajectos.replace(lastname.find(str_tif), str_tif.length(), "_Trajectos.txt");
        nameoutTrajectos.append(stroutTrajectos);
        FILE* MyObjectsTrajectosFile;
        MyObjectsTrajectosFile = fopen (nameoutTrajectos.c_str(), "w");
        fprintf(MyObjectsTrajectosFile, "TrajTAG\tFrame\tOrigID\tX\tY\tArea\tPeri\tConv\tIntInt\tMotherTAG\tNeigbors\n"); /// remember remember idx+1
        for (int z=0; z<Trajectos.size(); z++)
            for(int k=0; k<Trajectos[z].size(); k++)
            {
                fprintf(MyObjectsTrajectosFile, "%i\t%i\t%i\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%i\t%i\n", Trajectos[z][k].trajtag, Trajectos[z][k].frame, Trajectos[z][k].origID, Trajectos[z][k].center.x, Trajectos[z][k].center.y, Trajectos[z][k].area, Trajectos[z][k].peri, Trajectos[z][k].convexity, Trajectos[z][k].iI, Trajectos[z][k].mothertag, Trajectos[z][k].neigbors); /// remember remember idx+1
                fflush(MyObjectsTrajectosFile);
            }
        fflush(MyObjectsTrajectosFile);
        fclose(MyObjectsTrajectosFile);
        FILE* TrajectoFileFrame;
        FILE* TrajectoFileT;
        FILE* TrajectoFileX;
        FILE* TrajectoFileY;
        FILE* TrajectoFileA;
        FILE* TrajectoFileP;
        FILE* TrajectoFileC;
        TrajectoFileFrame = fopen (nameoutFrame.c_str(), "w");
        TrajectoFileT = fopen (nameoutT.c_str(), "w");
        TrajectoFileX = fopen (nameoutX.c_str(), "w");
        TrajectoFileY = fopen (nameoutY.c_str(), "w");
        TrajectoFileA = fopen (nameoutA.c_str(), "w");
        TrajectoFileP = fopen (nameoutP.c_str(), "w");
        TrajectoFileC = fopen (nameoutC.c_str(), "w");
        //SaveAllTrajecto3(TrajectoFileFrame, TrajectoFileT, TrajectoFileX, TrajectoFileY, Trajectos, nametrajecto, CellTypeName, scalefactor, TimeInterval, height);
        //SaveAllTrajecto4(TrajectoFileFrame, TrajectoFileT, TrajectoFileX, TrajectoFileY, TrajectoFileA, TrajectoFileP, Trajectos, nametrajecto, CellTypeName, scalefactor, TimeInterval, height);
        SaveAllTrajecto5(TrajectoFileFrame, TrajectoFileT, TrajectoFileX, TrajectoFileY, TrajectoFileA, TrajectoFileP, TrajectoFileC, Trajectos, nametrajecto, CellTypeName, scalefactor, TimeInterval, height);
        fclose(TrajectoFileFrame);
        fclose(TrajectoFileT);
        fclose(TrajectoFileX);
        fclose(TrajectoFileY);
        fclose(TrajectoFileA);
        fclose(TrajectoFileP);
        fclose(TrajectoFileC);
        if(SaveSingleVisualOutputAs==1)
            DrawTrajectoGreen (output_frames_color, Trajectos);
        else
            DrawTrajecto (output_frames_color, Trajectos);

        /// save a movies of each track to test the output
        cout << "Start to save all" << endl;

        for (int i=0; i<Trajectos.size(); i++)
        {

            if (Trajectos[i].size() > 3)
            {
                vector<Mat> output_single_i;
                //if (string(argv[2]) == "0")
                extract_all_i(output_frames_color, output_single_i, Trajectos[i]);
                /*else
                {
                    vector<Mat> input_tras_nuc_merge;
                    for( int tt=0; tt<input_tras.size(); tt++)
                        input_tras_nuc_merge.push_back(input_tras[tt] + output_frames_color[tt]);
                    extract_all_i(input_tras_nuc_merge, output_single_i, Trajectos[i]);
                }*/

                //cout << "extracted: " << i << endl;
                string AllSingleFolder = path + "AllSingle/";
                mkdir(AllSingleFolder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
                string nametrajectoi=nametrajecto;
                std::string snum;
                std::stringstream snumout;
                snumout << Trajectos[i][0].trajtag-1;
                snum = snumout.str();
                snum = "#" + snum;
                nametrajectoi.append(snum);
                if (SaveSingleVisualOutputAs==1)
                {
                    string name_single = AllSingleFolder+nametrajectoi+"_##_"+CellTypeName+".tiff" ;
                    if (output_single_i.size() > 2)
                    {
                        //cout << "name single movie: " << name_single.c_str() << endl;
                        Save8bitcolorTiff(name_single, output_single_i);
                    }
                }
                else
                {
                    string name_single = AllSingleFolder+nametrajectoi+"_##_"+CellTypeName+".avi" ;
                    if (output_single_i.size() > 2)
                    {
                        //cout << "name single movie: " << name_single.c_str() << endl;
                        Save8bitcolor(name_single, output_single_i);
                    }
                }
            }
        }
        //cout << "End to save all" << endl;
        /// end
    }

    printf("Save output video\n");
    string nameout;
    nameout = name;
    nameout.replace(nameout.find(str_tif), str_tif.length(), ".avi");
    Save8bitcolor(nameout, output_frames_color);
    printf("output saved at %s\n", nameout.c_str());
    //string stravi = ".avi";
    //nameout.replace(nameout.find(stravi), stravi.length(), "_Segmented.tif");
    //SaveTiffOutputNameStr8bit(nameout, output_segmented, ntot, width, height);

    int fastesttraj=Trajectos.size()+1;
    int minidx=0;
    int mintime=ntot+1;
    double lengthlimit= 100/scalefactor; //310.6;	//100um in race UCSF // px/um 0.322
    double meanvel=0;

    if (string(argv[4]) == "LINE" || string(argv[4]) == "1D")
    {

        for (int i=0; i<Trajectos.size(); i++)
        {
            int minidxi;
            int minti;
            double effdisti;
            if (Trajectos[i].size() > 1)
            {
                Fastest100Trai(lengthlimit, Trajectos[i], minidxi, minti, effdisti);

                fastesttraj=i;
                minidx= minidxi;
                mintime= minti;
                meanvel = effdisti/(minti*TimeInterval);
                if (mintime < ntot+1)
                {
                    vector<Mat> outputfastest;
                    vector<Mat> outputfastest_tras;
                    printf("\nFastest Traj:%i\t at frame:%i\t time:%i\t average velocity:%.3f\n", fastesttraj, minidx, mintime, meanvel);
                    fflush(stdout);
                    //plotTrajecto(Trajectos[fastesttraj]);
                    extractfastest(output_frames_color, outputfastest, Trajectos[fastesttraj], minidx, mintime, lengthlimit, VisualOutput);
                    std::string snum;
                    std::stringstream snumout;
                    snumout << mintime;
                    snum = snumout.str();
                    std::string snum1;
                    std::stringstream snumout1;
                    snumout1 << fastesttraj;
                    snum1 = snumout1.str();
                    std::string snum2;
                    std::stringstream snumout2;
                    snumout2 << meanvel;
                    snum2 = snumout2.str();
                    string L1Folder = path + "L100/";
                    mkdir(L1Folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
                    string namefastest = L1Folder + "L100um_" + snum + "_" + snum2 + "_" + snum1 + "_" + lastname;
                    namefastest.replace(namefastest.find(str_tif), str_tif.length(), "_fastest.avi");
                    Save8bitcolor(namefastest, outputfastest);
                    /*if (string(argv[2]) != "0")
                    {
                        extractfastest(input_tras, outputfastest_tras, Trajectos[fastesttraj], minidx, mintime, lengthlimit, VisualOutput);
                        string fastavi = "_fastest.avi";
                        namefastest.replace(namefastest.find(fastavi), fastavi.length(), "_Trans_fastest.avi");
                        Save8bitcolor(namefastest, outputfastest_tras);
                    }*/
                }
            }
        }

        lengthlimit = 350/scalefactor;
        fastesttraj=Trajectos.size()+1;
        minidx=0;
        mintime=ntot+1;
        meanvel=0;


        for (int i=0; i<Trajectos.size(); i++)
        {
            int minidxi;
            int minti;
            double effdisti;
            if (Trajectos[i].size() > 1)
            {
                Fastest100Trai(lengthlimit, Trajectos[i], minidxi, minti, effdisti);

                fastesttraj=i;
                minidx= minidxi;
                mintime= minti;
                meanvel = effdisti/(minti*TimeInterval);
                if (mintime < ntot+1)
                {
                    vector<Mat> outputfastest;
                    vector<Mat> outputfastest_tras;
                    printf("\nFastest Traj:%i\t at frame:%i\t time:%i\t average velocity:%.3f\n", fastesttraj, minidx, mintime, meanvel);
                    fflush(stdout);
                    extractfastest(output_frames_color, outputfastest, Trajectos[fastesttraj], minidx, mintime, lengthlimit, VisualOutput);
                    std::string snum;
                    std::stringstream snumout;
                    snumout << mintime;
                    snum = snumout.str();
                    std::string snum1;
                    std::stringstream snumout1;
                    snumout1 << fastesttraj;
                    snum1 = snumout1.str();
                    std::string snum2;
                    std::stringstream snumout2;
                    snumout2 << meanvel;
                    snum2 = snumout2.str();
                    string L1Folder = path + "L350/";
                    mkdir(L1Folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
                    string namefastest = L1Folder + "L350um_" + snum + "_" + snum2 + "_" + snum1 + "_" + lastname;
                    namefastest.replace(namefastest.find(str_tif), str_tif.length(), "_fastest.avi");
                    Save8bitcolor(namefastest, outputfastest);
                    /*if (string(argv[2]) != "0")
                    {
                        printf("\n\t\t\tExtract also trasmitted\n");
                        extractfastest(input_tras, outputfastest_tras, Trajectos[fastesttraj], minidx, mintime, lengthlimit, VisualOutput);
                        string fastavi = "_fastest.avi";
                        namefastest.replace(namefastest.find(fastavi), fastavi.length(), "_Trans_fastest.avi");
                        Save8bitcolor(namefastest, outputfastest_tras);
                    }*/
                }
            }
        }

        lengthlimit = 700/scalefactor;
        fastesttraj=Trajectos.size()+1;
        minidx=0;
        mintime=ntot+1;
        meanvel=0;

        for (int i=0; i<Trajectos.size(); i++)
        {
            int minidxi;
            int minti;
            double effdisti;
            if (Trajectos[i].size() > 1)
            {
                Fastest100Trai(lengthlimit, Trajectos[i], minidxi, minti, effdisti);

                fastesttraj=i;
                minidx= minidxi;
                mintime= minti;
                meanvel = effdisti/(minti*TimeInterval);
                if (mintime < ntot+1)
                {
                    vector<Mat> outputfastest;
                    vector<Mat> outputfastest_tras;
                    printf("\nFastest Traj:%i\t at frame:%i\t time:%i\t average velocity:%.3f\n", fastesttraj, minidx, mintime, meanvel);
                    fflush(stdout);
                    //plotTrajecto(Trajectos[fastesttraj]);
                    extractfastest(output_frames_color, outputfastest, Trajectos[fastesttraj], minidx, mintime, lengthlimit, VisualOutput);
                    std::string snum;
                    std::stringstream snumout;
                    snumout << mintime;
                    snum = snumout.str();
                    std::string snum1;
                    std::stringstream snumout1;
                    snumout1 << fastesttraj;
                    snum1 = snumout1.str();
                    std::string snum2;
                    std::stringstream snumout2;
                    snumout2 << meanvel;
                    snum2 = snumout2.str();
                    string L1Folder = path + "L700/";
                    mkdir(L1Folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
                    string namefastest = L1Folder + "L700um_" + snum + "_" + snum2 + "_" + snum1 + "_" + lastname;
                    namefastest.replace(namefastest.find(str_tif), str_tif.length(), "_fastest.avi");
                    Save8bitcolor(namefastest, outputfastest);
                    /*if (string(argv[2]) != "0")
                    {
                        printf("\n\t\t\tExtract also trasmitted\n");
                        extractfastest(input_tras, outputfastest_tras, Trajectos[fastesttraj], minidx, mintime, lengthlimit, VisualOutput);
                        string fastavi = "_fastest.avi";
                        namefastest.replace(namefastest.find(fastavi), fastavi.length(), "Trans_fastest.avi"); //UCSF
                        Save8bitcolor(namefastest, outputfastest_tras);
                    }*/
                }
            }
        }
        //FreeVectorMat(input_tras);
    }
    fflush(stdout);
    exit(0);
}
