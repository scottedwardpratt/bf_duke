#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>
#include <string>
#include <cstring>
#include <boost/math/special_functions.hpp>

const double PI=4.0*atan(1.0);
const double HBARC=197.3269602;

using namespace std;

int main(int argc,char *argv[]){
	char cline[300];
	double pt1,pt2,uncplus1,uncplus2,uncminus1,uncminus2,spectraplus,spectraminus,crapplus,crapminus;
  FILE *fptrold=fopen("original/spectra_K_original.txt","r");
	FILE *fptrnew=fopen("spectra_K.txt","w");
	do{
		fscanf(fptrold,"%lf %lf %lf %lf %lf %lf %lf %lf ",&pt1,&pt2,&spectraplus,&uncplus1,&uncplus2,&spectraminus,&uncminus1,&uncminus2);
		if(!feof(fptrold))
			fprintf(fptrnew,"%10.5f %10.5f %10.2f\n",pt1,pt2,0.5*(spectraplus+spectraminus));
	}while(!feof(fptrold));
	
	fclose(fptrold);
	fclose(fptrnew);
	
  return 0;
}


