#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include <string>
#include <cstring>
#include <iostream>
#include <sstream>
#include <filesystem>
#include "msu_commonutils/misc.h"



using namespace std;
using namespace NMSUPratt;

int main(int argc,char *argv[]){
	const double PI=4.0*atan(1.0);
	string qualifier[4]={"alicePbPb_cent0_5","alicePbPb_cent5_10","alicePbPb_cent10_15","alicePbPb_cent15_20"};
	string filename,varname;
	int iqual,irun,isub,nsubruns,ipt,npt=120;
	double unc,read1,read2,read3,read4,v2_pi,v2_K,v2_p;
	double mult_pi=0.0,mult_K=0.0,mult_p=0.0;
	vector<double> v2pi(npt),v2K(npt),v2p(npt),pt(npt);
	vector<double> spectrapi(npt),spectraK(npt),spectrap(npt);
	
	double ptmin_pi=0.2,ptmin_K=0.2,ptmin_p=0.5;
	double ptmax_pi=2.0,ptmax_K=2.0,ptmax_p=2.5;
	double delpt=0.025,d3poverE,pt1,pt2;
	bool exists;
	FILE *fptrspectra,*fptrv2,*fptrout;
	for(iqual=0;iqual<4;iqual++){
		for(irun=0;irun<6;irun++){
			isub=0;
			for(ipt=0;ipt<npt;ipt++){
				v2pi[ipt]=0.0;
				v2K[ipt]=0.0;
				v2p[ipt]=0.0;
			}
			
			filename="modelruns/run"+to_string(irun)+"/"+qualifier[iqual]+"/results_spectrav2/v2.txt";
			fptrv2=fopen(filename.c_str(),"r");
			filename="modelruns/run"+to_string(irun)+"/"+qualifier[iqual]+"/results_spectrav2/spectra.txt";
			fptrspectra=fopen(filename.c_str(),"r");
			v2_pi=v2_K=v2_p=mult_pi=mult_K=mult_p=0.0;
			for(ipt=0;ipt<npt;ipt++){
				fscanf(fptrv2,"%lf %lf %lf %lf\n",&pt[ipt],&v2pi[ipt],&v2K[ipt],&v2p[ipt]);
				fscanf(fptrspectra,"%lf %lf %lf %lf\n",&pt[ipt],&spectrapi[ipt],&spectraK[ipt],&spectrap[ipt]);
				pt1=pt[ipt]-0.5*delpt;
				pt2=pt[ipt]+0.5*delpt;
				d3poverE=PI*(pt2*pt2-pt1*pt1);
				if(pt[ipt]>ptmin_pi && pt[ipt]<ptmax_pi){
					mult_pi+=spectrapi[ipt]*d3poverE;
					v2_pi+=v2pi[ipt]*spectrapi[ipt]*d3poverE*pt[ipt];
				}
				if(pt[ipt]>ptmin_K && pt[ipt]<ptmax_K){
					mult_K+=spectraK[ipt]*spectraK[ipt]*d3poverE;
					v2_K+=v2K[ipt]*spectraK[ipt]*d3poverE*pt[ipt];
				}
				if(pt[ipt]>ptmin_p && pt[ipt]<ptmax_p){
					mult_p+=spectrap[ipt]*d3poverE;
					v2_p+=v2p[ipt]*spectrap[ipt]*d3poverE*pt[ipt];
				}
		
			};
			fclose(fptrv2);
			fclose(fptrspectra);
			v2_pi=v2_pi/mult_pi;
			v2_K=v2_K/mult_K;
			v2_p=v2_p/mult_p;
	
			// Now write results
			filename="modelruns/run"+to_string(irun)+"/obs.txt";
			fptrout=fopen(filename.c_str(),"a");
			//
			unc=0.04*v2_pi;
			varname="v2_pi_"+qualifier[iqual];
			fprintf(fptrout,"%s   %10.5f  %10.5f  0.0\n",varname.c_str(),v2_pi,unc);
			//
			unc=0.04*v2_K;
			varname="v2_K_"+qualifier[iqual];
			fprintf(fptrout,"%s   %10.5f  %10.5f  0.0\n",varname.c_str(),v2_K,unc);
			//
			unc=0.04*v2_p;
			varname="v2_p_"+qualifier[iqual];
			fprintf(fptrout,"%s   %10.5f  %10.5f  0.0\n",varname.c_str(),v2_p,unc);
			//
	
			fclose(fptrout);	
		}
	}
	return 0;
}


