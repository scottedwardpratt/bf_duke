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
	int iqual,nsubruns,ipt,irun,npt=120;
	double meanpt_pi=0.0,meanpt_K=0.0,meanpt_p=0.0,unc,read1,read2,read3,read4;
	double mult_pi=0.0,mult_K=0.0,mult_p=0.0;
	vector<double> specpi(npt),specK(npt),specp(npt);
	double pt[120];
	double ptmin_pi=0.2,ptmin_K=0.2,ptmin_p=0.5;
	double ptmax_pi=2.0,ptmax_K=2.0,ptmax_p=2.5;
	double delpt=0.025,d3poverE,pt1,pt2;
	bool exists;
	
	FILE *fptrin,*fptrout;
	
	for(iqual=0;iqual<4;iqual++){
		for(irun=0;irun<6;irun++){	
			filename="modelruns/run"+to_string(irun)+"/"+qualifier[iqual]+"/results_spectrav2/spectra.txt";
			fptrin=fopen(filename.c_str(),"r");
			mult_pi=mult_K=mult_p=meanpt_pi=meanpt_K=meanpt_p=0.0;
			for(ipt=0;ipt<npt;ipt++){
				fscanf(fptrin,"%lf %lf %lf %lf\n",&pt[ipt],&specpi[ipt],&specK[ipt],&specp[ipt]);
				pt1=pt[ipt]-0.5*delpt;
				pt2=pt[ipt]+0.5*delpt;
				d3poverE=PI*(pt2*pt2-pt1*pt1);
				if(pt[ipt]>ptmin_pi && pt[ipt]<ptmax_pi){
					mult_pi+=specpi[ipt]*d3poverE;
					meanpt_pi+=specpi[ipt]*d3poverE*pt[ipt];
				}
				if(pt[ipt]>ptmin_K && pt[ipt]<ptmax_K){
					mult_K+=specK[ipt]*d3poverE;
					meanpt_K+=specK[ipt]*d3poverE*pt[ipt];
				}
				if(pt[ipt]>ptmin_p && pt[ipt]<ptmax_p){
					mult_p+=specp[ipt]*d3poverE;
					meanpt_p+=specp[ipt]*d3poverE*pt[ipt];
				}
		
			};
			fclose(fptrin);
			
			meanpt_pi=meanpt_pi/mult_pi;
			meanpt_K=meanpt_K/mult_K;
			meanpt_p=meanpt_p/mult_p;
	
			// Now write results
			filename="modelruns/run"+to_string(irun)+"/obs.txt";
			fptrout=fopen(filename.c_str(),"a");
			//
			unc=0.05*mult_pi;
			varname="mult_pi_"+qualifier[iqual];
			fprintf(fptrout,"%s   %10.5f  %10.5f  0.0\n",varname.c_str(),mult_pi,unc);
			unc=0.025*meanpt_pi;
			varname="meanpt_pi_"+qualifier[iqual];
			fprintf(fptrout,"%s   %10.5f  %10.5f  0.0\n",varname.c_str(),meanpt_pi,unc);
			//
			unc=0.05*mult_K;
			varname="mult_K_"+qualifier[iqual];
			fprintf(fptrout,"%s   %10.5f  %10.5f  0.0\n",varname.c_str(),mult_K,unc);
			unc=0.025*meanpt_K;
			varname="meanpt_K_"+qualifier[iqual];
			fprintf(fptrout,"%s   %10.5f  %10.5f  0.0\n",varname.c_str(),meanpt_K,unc);
			//
			unc=0.05*mult_p;
			varname="mult_p_"+qualifier[iqual];
			fprintf(fptrout,"%s   %10.5f  %10.5f  0.0\n",varname.c_str(),mult_p,unc);
			unc=0.025*meanpt_p;
			varname="meanpt_p_"+qualifier[iqual];
			fprintf(fptrout,"%s   %10.5f  %10.5f  0.0\n",varname.c_str(),meanpt_p,unc);
			//
	
			fclose(fptrout);	
		}
	}
	return 0;
}


