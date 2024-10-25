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
	double pt,spectra,mult,meanpt,pt1,pt2,d3poverE,ptmin,ptmax,unc;
  FILE *fptrin;
	FILE *fptrout=fopen("experimental_info.txt","a");
	
	// for pions
	fptrin=fopen("spectra_pi.txt","r");
	ptmin=0.2;
	ptmax=2.0;
	mult=meanpt=0.0;
	do{
		fscanf(fptrin,"%lf %lf %lf",&pt1,&pt2,&spectra);
		pt=0.5*(pt1+pt2);
		if(pt>ptmin && pt<ptmax){
			d3poverE=PI*(pt2*pt2-pt1*pt1);
			mult+=d3poverE*spectra;
			meanpt+=d3poverE*spectra*pt;
		}
	}while(!feof(fptrin));
	fclose(fptrin);
	meanpt=meanpt/mult;
	unc=0.05*mult;
	fprintf(fptrout,"mult_pi %8.5f %8.5f  0.0\n",mult,unc);
	unc=0.03*meanpt;
	fprintf(fptrout,"meanpt_pi %8.5f %8.5f  0.0\n",meanpt,unc);
	
	// for kaons
	fptrin=fopen("spectra_K.txt","r");
	ptmin=0.2;
	ptmax=2.0;
	mult=meanpt=0.0;
	do{
		fscanf(fptrin,"%lf %lf %lf",&pt1,&pt2,&spectra);
		pt=0.5*(pt1+pt2);
		if(pt>ptmin && pt<ptmax){
			d3poverE=PI*(pt2*pt2-pt1*pt1);
			mult+=d3poverE*spectra;
			meanpt+=d3poverE*spectra*pt;
		}
	}while(!feof(fptrin));
	fclose(fptrin);
	
	meanpt=meanpt/mult;
	unc=0.05*mult;
	fprintf(fptrout,"mult_K %8.5f %8.5f  0.0\n",mult,unc);
	unc=0.03*meanpt;
	fprintf(fptrout,"meanpt_K %8.5f %8.5f  0.0\n",meanpt,unc);
	
	// for protons
	fptrin=fopen("spectra_p.txt","r");
	ptmin=0.5;
	ptmax=2.50;
	mult=meanpt=0.0;
	do{
		fscanf(fptrin,"%lf %lf %lf",&pt1,&pt2,&spectra);
		pt=0.5*(pt1+pt2);
		if(pt>ptmin && pt<ptmax){
			d3poverE=PI*(pt2*pt2-pt1*pt1);
			mult+=d3poverE*spectra;
			meanpt+=d3poverE*spectra*pt;
		}
	}while(!feof(fptrin));
	fclose(fptrin);
	
	meanpt=meanpt/mult;
	unc=0.05*mult;
	fprintf(fptrout,"mult_p %8.5f %8.5f  0.0\n",mult,unc);
	unc=0.03*meanpt;
	fprintf(fptrout,"meanpt_p %8.5f %8.5f  0.0\n",meanpt,unc);

	fclose(fptrout);
	
  return 0;
}


