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
	int iqual,irun;
	double bfwidth_y_pipi=0.0,bfnorm_y_pipi=0.0;
	double bfwidth_y_KK=0.0,bfnorm_y_KK=0.0;
	double bfwidth_y_pp=0.0,bfnorm_y_pp=0.0;
	double bfwidth_phi_pipi=0.0,bfnorm_phi_pipi=0.0;
	double bfwidth_phi_KK=0.0,bfnorm_phi_KK=0.0;
	double bfwidth_phi_pp=0.0,bfnorm_phi_pp=0.0;
	
	double ymax_pipi=1.4,ymax_KK=1.2,ymax_pp=1.0;
	double dely=0.1,delphi=PI/14.0,y,unc,b,phi;
	FILE *fptrin,*fptrout;
	
	for(irun=0;irun<6;irun++){
		for(iqual=0;iqual<4;iqual++){
			bfnorm_y_pipi=bfnorm_y_KK=bfnorm_y_pp=0.0;
			bfwidth_y_pipi=bfwidth_y_KK=bfwidth_y_pp=0.0;
		

			filename="modelruns/run"+to_string(irun)+"/"+qualifier[iqual]+"/results/pipi/bf_y.txt";
			fptrin=fopen(filename.c_str(),"r");
			do{
				fscanf(fptrin,"%lf %lf %lf",&y,&b,&unc);
				if(y<ymax_pipi){
					bfnorm_y_pipi+=dely*b;
					bfwidth_y_pipi+=dely*b*y;
				}
			}while(y<ymax_pipi);
			fclose(fptrin);
			filename="modelruns/run"+to_string(irun)+"/"+qualifier[iqual]+"/results/pipi/bf_phi.txt";
			fptrin=fopen(filename.c_str(),"r");
			do{
				fscanf(fptrin,"%lf %lf %lf",&phi,&b,&unc);
				if(!feof(fptrin)){
					bfnorm_phi_pipi+=delphi*b;
					bfwidth_phi_pipi+=delphi*b*fabs(phi);
				}
			}while(!feof(fptrin));
			fclose(fptrin);
			
			filename="modelruns/run"+to_string(irun)+"/"+qualifier[iqual]+"/results/KK/bf_y.txt";
			fptrin=fopen(filename.c_str(),"r");
			do{
				fscanf(fptrin,"%lf %lf %lf",&y,&b,&unc);
				if(y<ymax_KK){
					bfnorm_y_KK+=dely*b;
					bfwidth_y_KK+=dely*b*y;
				}
			}while(y<ymax_KK);
			fclose(fptrin);
			
			filename="modelruns/run"+to_string(irun)+"/"+qualifier[iqual]+"/results/KK/bf_phi.txt";
			fptrin=fopen(filename.c_str(),"r");
			do{
				fscanf(fptrin,"%lf %lf %lf",&phi,&b,&unc);
				if(!feof(fptrin)){
					bfnorm_phi_KK+=delphi*b;
					bfwidth_phi_KK+=delphi*b*fabs(phi);
				}
			}while(!feof(fptrin));
			fclose(fptrin);

			filename="modelruns/run"+to_string(irun)+"/"+qualifier[iqual]+"/results/pp/bf_y.txt";
			fptrin=fopen(filename.c_str(),"r");
			do{
				fscanf(fptrin,"%lf %lf %lf",&y,&b,&unc);
				if(y<ymax_pp){
					bfnorm_y_pp+=dely*b;
					bfwidth_y_pp+=dely*b*y;
				}
			}while(y<ymax_pp);
			fclose(fptrin);
			
			filename="modelruns/run"+to_string(irun)+"/"+qualifier[iqual]+"/results/pp/bf_phi.txt";
			fptrin=fopen(filename.c_str(),"r");
			do{
				fscanf(fptrin,"%lf %lf %lf",&phi,&b,&unc);
				if(!feof(fptrin)){
					bfnorm_phi_pp+=delphi*b;
					bfwidth_phi_pp+=delphi*b*fabs(phi);
				}
			}while(!feof(fptrin));
			
			fclose(fptrin);
	
 
			bfwidth_y_pipi=bfwidth_y_pipi/bfnorm_y_pipi;
			bfwidth_y_KK=bfwidth_y_KK/bfnorm_y_KK;
			bfwidth_y_pp=bfwidth_y_pp/bfnorm_y_pp;
	
			// Now write results
			filename="modelruns/run"+to_string(irun)+"/obs.txt";
			printf("writing results to %s\n",filename.c_str());
			fptrout=fopen(filename.c_str(),"a");
			//
			unc=0.05*bfnorm_y_pipi;
			varname="bfnorm_y_pipi_"+qualifier[iqual];
			fprintf(fptrout,"%s   %10.5f  %10.5f  0.0\n",varname.c_str(),bfnorm_y_pipi,unc);
			unc=0.03*bfwidth_y_pipi;
			varname="bfwidth_y_pipi_"+qualifier[iqual];
			fprintf(fptrout,"%s   %10.5f  %10.5f  0.0\n",varname.c_str(),bfwidth_y_pipi,unc);
			unc=0.03*bfwidth_phi_pipi;
			varname="bfwidth_phi_pipi_"+qualifier[iqual];
			fprintf(fptrout,"%s   %10.5f  %10.5f  0.0\n",varname.c_str(),bfwidth_phi_pipi,unc);
			//
			unc=0.05*bfnorm_y_KK;
			varname="bfnorm_y_KK_"+qualifier[iqual];
			fprintf(fptrout,"%s   %10.5f  %10.5f  0.0\n",varname.c_str(),bfnorm_y_KK,unc);
			unc=0.03*bfwidth_y_KK;
			varname="bfwidth_y_KK_"+qualifier[iqual];
			fprintf(fptrout,"%s   %10.5f  %10.5f  0.0\n",varname.c_str(),bfwidth_y_KK,unc);
			unc=0.03*bfwidth_phi_KK;
			varname="bfwidth_phi_KK_"+qualifier[iqual];
			fprintf(fptrout,"%s   %10.5f  %10.5f  0.0\n",varname.c_str(),bfwidth_phi_KK,unc);
			//
			unc=0.05*bfnorm_y_pp;
			varname="bfnorm_y_pp_"+qualifier[iqual];
			fprintf(fptrout,"%s   %10.5f  %10.5f  0.0\n",varname.c_str(),bfnorm_y_pp,unc);
			unc=0.03*bfwidth_y_pp;
			varname="bfwidth_y_pp_"+qualifier[iqual];
			fprintf(fptrout,"%s   %10.5f  %10.5f  0.0\n",varname.c_str(),bfwidth_y_pp,unc);
			unc=0.03*bfwidth_phi_pp;
			varname="bfwidth_phi_pp_"+qualifier[iqual];
			fprintf(fptrout,"%s   %10.5f  %10.5f  0.0\n",varname.c_str(),bfwidth_phi_pp,unc);
			//
			fclose(fptrout);
		}	
	}
	return 0;
}


