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
	double meanpt_pi=0.0,meanpt_K=0.0,meanpt_p=0.0,unc,read1,read2,read3,read4;
	double mult_pi=0.0,mult_K=0.0,mult_p=0.0;
	vector<double> v2pi(npt),v2K(npt),v2p(npt),pt(npt);
	
	double ptmin_pi=0.2,ptmin_K=0.2,ptmin_p=0.5;
	double ptmax_pi=2.0,ptmax_K=2.0,ptmax_p=2.5;
	double delpt=0.025,d3poverE,pt1,pt2;
	bool exists;
	FILE *fptrin,*fptrout;
	for(iqual=0;iqual<5;iqual++){
		for(irun=0;irun<6;irun++){
			isub=0;
			for(ipt=0;ipt<npt;ipt++){
				v2pi[ipt]=0.0;
				v2K[ipt]=0.0;
				v2p[ipt]=0.0;
			}
			do{
				filename="modelruns/run"+to_string(irun)+"/"+qualifier[iqual]+"/results_spectrav2/subruns/subrun"+to_string(isub)+"/v2.txt";
				exists=filesystem::exists(filename);
				if(exists){
					fptrin=fopen(filename.c_str(),"r");
					for(ipt=0;ipt<npt;ipt++){
						fscanf(fptrin,"%lf %lf %lf %lf",&read1,&read2,&read3,&read4);
						pt[ipt]=read1;
						v2pi[ipt]+=read2;
						v2K[ipt]+=read3;
						v2p[ipt]+=read4;
					}
					fclose(fptrin);
					isub+=1;
				}
			}while(exists);
			nsubruns=isub-1;
			for(ipt=0;ipt<npt;ipt++){
				v2pi[ipt]=v2pi[ipt]/double(nsubruns);
				v2K[ipt]=v2K[ipt]/double(nsubruns);
				v2p[ipt]=v2p[ipt]/double(nsubruns);
			}
		
			filename="modelruns/run"+to_string(irun)+"/"+qualifier[iqual]+"/results_spectrav2/v2.txt";
			fptrout=fopen(filename.c_str(),"w");
			for(ipt=0;ipt<npt;ipt++){
				fprintf(fptrout,"%lf %lf %lf %lf\n",pt[ipt],v2pi[ipt],v2K[ipt],v2p[ipt]);
			};
		}
	}
	fclose(fptrout);
	return 0;
}


