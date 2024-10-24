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
	int Nsubruns_max=99999999;
	const int nspecies=3;
	double unc;
	bool exists;
	char dummy[200];
	string filename,varname;
	int isubrun,Nsubruns,ispecies,ipt,irun;
	double meanpt_pi,meanpt_K,meanpt_p,meanv2_pi,meanv2_K,meanv2_p; 
	double pt_pi,pt_K,pt_p,v2_pi,v2_K,v2_p; 
	FILE *fptrin,*fptrout;
	

	vector<double> bfqinvsum,bfoutsum,bfsidesum,bflongsum,cfqinvsum,cfoutsum,cfsidesum,cflongsum;
	if(argc!=4){
		printf("Usage: summer_spectrav2 qualifier I J are going from modelruns/runI to modelruns/runJ\n");
		exit(1);
	}
	exists=filesystem::exists("modelruns");
	if(!exists){
		printf("this needs to be run from parent directory to modelruns/\n");
		exit(1);
	}
	string qualifier=argv[1];
	int I=atoi(argv[2]);
	int J=atoi(argv[3]);
	string speciesnames[nspecies]={"pi","K","p"};

	for(irun=I;irun<=J;irun++){
		meanpt_pi=meanpt_K=meanpt_p=meanv2_pi=meanv2_K=meanv2_p=0.0;

		// Get Nsubruns
		Nsubruns=0;
		do{
			filename="modelruns/run"+to_string(irun)+"/"+qualifier+"/results_spectrav2/subruns/subrun"+to_string(Nsubruns)+"/meanpt_meanv2.txt";
			printf("filename=%s\n",filename.c_str());
			exists=filesystem::exists(filename);
			if(exists){
				printf("%s exists\n",filename.c_str());
				Nsubruns+=1;
				fptrin=fopen(filename.c_str(),"r");
				while(!feof(fptrin)){
					fscanf(fptrin,"%s %lf",dummy,&pt_pi);
					fscanf(fptrin,"%s %lf",dummy,&pt_K);
					fscanf(fptrin,"%s %lf",dummy,&pt_p);
					fscanf(fptrin,"%s %lf",dummy,&v2_pi);
					fscanf(fptrin,"%s %lf",dummy,&v2_K);
					fscanf(fptrin,"%s %lf",dummy,&v2_p);
					if(!feof(fptrin)){
						meanpt_pi+=pt_pi;
						meanpt_K+=pt_K;
						meanpt_p+=pt_p;
						meanv2_pi+=v2_pi;
						meanv2_K+=v2_K;
						meanv2_p+=v2_p;
					}
				}	
			}	
		}while(exists);
		meanpt_pi=meanpt_pi/double(Nsubruns);
		meanpt_K=meanpt_K/double(Nsubruns);
		meanpt_p=meanpt_p/double(Nsubruns);
		meanv2_pi=meanv2_pi/double(Nsubruns);
		meanv2_K=meanv2_K/double(Nsubruns);
		meanv2_p=meanv2_p/double(Nsubruns);
		fclose(fptrin);
		
		filename="modelruns/run"+to_string(irun)+"/obs.txt";
		fptrout=fopen(filename.c_str(),"a");
		//
		varname="meanpt_pi_"+qualifier;
		unc=0.03*meanpt_pi;
		fprintf(fptrout,"%s  %10.5f  %10.5f  0.0\n",varname.c_str(),meanpt_pi,unc);
		//
		varname="meanpt_K_"+qualifier;
		unc=0.03*meanpt_K;
		fprintf(fptrout,"%s  %10.5f  %10.5f  0.0\n",varname.c_str(),meanpt_K,unc);
		//
		varname="meanpt_p_"+qualifier;
		unc=0.03*meanpt_p;
		fprintf(fptrout,"%s  %10.5f  %10.5f  0.0\n",varname.c_str(),meanpt_p,unc);
		//
		varname="meanv2_pi_"+qualifier;
		unc=0.03*meanv2_pi;
		fprintf(fptrout,"%s  %10.5f  %10.5f  0.0\n",varname.c_str(),meanv2_pi,unc);
		//
		varname="meanv2_K_"+qualifier;
		unc=0.03*meanv2_K;
		fprintf(fptrout,"%s  %10.5f  %10.5f  0.0\n",varname.c_str(),meanv2_K,unc);
		//
		varname="meanv2_p_"+qualifier;
		unc=0.03*meanv2_p;
		fprintf(fptrout,"%s  %10.5f  %10.5f  0.0\n",varname.c_str(),meanv2_p,unc);
		
		fclose(fptrout);
	}
	return 0;
}


