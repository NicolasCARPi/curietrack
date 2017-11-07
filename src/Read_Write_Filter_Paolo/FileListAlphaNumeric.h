#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/video/background_segm.hpp>
#include <opencv2/imgproc/imgproc_c.h>
#include <opencv2/opencv.hpp>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_math.h>
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
#include <errno.h>
#include <tiff.h>
#include <dirent.h>


using namespace cv;
using namespace std;

////////////////////////////////
bool lessthanAlphabetic(string const& a, string const& b)
{
	int offset = 0;
    	while (true)
    		{
        	if (offset == b.size()) // then b < a or b == a.
            		return false;
        	if (offset == a.size()) // then a < b.
            		return true;
        	if (tolower(a[offset]) < tolower(b[offset]))
            		return true; // true that a < b.
        	if (tolower(b[offset]) < tolower(a[offset]))
            		return false; // because b < a.
        	offset++; // the current characters are the same, so check the next.
		}
}

bool AphaNumeric(string const& a, string const& b)
{
	int la=a.size();
	int lb=b.size();
    	int offset = 0;
    	int diff = 0; // first defference
	while ( diff == 0 && offset<min(la,lb) )
		{
		if ( a[offset] == b[offset] )
			offset++;
		else
			diff = 1;
		}
	if (diff==0)
		{
		if ( la<lb )
			return true;
		else
			return false;
		}
	else
		{
		if ( isdigit(a[offset]) && isdigit(b[offset]))
			{
			int ldigita=0;
			int ldigitb=0;
			int offsetdigita=1;
			int offsetdigitb=1;
			while (isdigit(a[offset+offsetdigita]))
				offsetdigita++;
			while (isdigit(b[offset+offsetdigitb]))
				offsetdigitb++;
			string suba, subb;
			suba = a.substr (offset, offsetdigita);
			subb = b.substr (offset, offsetdigitb);
			int inta, intb;
			inta = atoi(suba.c_str());
			intb = atoi(subb.c_str());
			if (inta < intb)
				return true;
			else
				return false;
			}
		else if (isdigit(a[offset]))
			{
			if (offset>0)
				if (isdigit(b[offset-1]))
					return false;
			}
		else if (isdigit(b[offset]))
			{
			if (offset>0)
				if (isdigit(a[offset-1]))
					return true;
			}

		else
			{
			if (a[offset] < b[offset])
				return true;
			else
				return false;
			}
		}
}

void listFileAphaNumeric(string dirname, vector<string>& filelist)
{
	DIR *dir;
  	struct dirent *entry;
	unsigned char isFile =0x8;
  	if ((dir = opendir(dirname.c_str())) == NULL)
    		perror("opendir() error");
  	else
		{
    		//puts("contents of root:");
    		while ((entry = readdir(dir)) != NULL)
			if ( entry->d_type == isFile)
   				{
      				//printf("  %s\n", entry->d_name);
				filelist.push_back(entry->d_name);
				}
    		closedir(dir);
  		}
	//sort(filelist.begin(), filelist.end(), lessthanAlphabetic);
	sort(filelist.begin(), filelist.end(), AphaNumeric);
	/*for (int i=0; i<filelist.size(); i++)
		{
		printf("  %s\n", filelist[i].c_str());
		}*/
}

















