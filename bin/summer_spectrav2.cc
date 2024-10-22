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
	const double DELPT_SPECTRA=0.025;
	const double DELPT_V2=0.025;
	unsigned int ipt,Nsubruns=0,Nsubruns_max=99999999;
	bool exists;
	string dirname,filename,command;
	char dummy[100];
	string meanlabel[6];
	unsigned int isubrun,irun;
	vector<string> subrunfilenames1,subrunfilenames2;
	FILE *fptr;
	double pt,specpi,specK,specp,v2pi,v2K,v2p;
	vector<double> meanpt_sum,meanv2_sum;
	vector<vector<double>> spectra_sum,v2_sum;
	if(argc!=4){
		printf("Usage: summer qualifier I J are going from modelruns/runI to modelruns/runJ\n");
		exit(1);
	}
	exists=filesystem::exists("modelruns");
	if(!exists){
		printf("this needs to be run from parent directory to modelruns/\n");
		exit(1);
	}
	string qualifier=argv[1];
	unsigned int I=atoi(argv[2]);
	unsigned int J=atoi(argv[3]);


	for(irun=I;irun<=J;irun++){

		// Get Nsubruns
		Nsubruns=0;
		do{
			filename="modelruns/run"+to_string(irun)+"/"+qualifier+"/results_spectrav2/subruns/subrun"+to_string(Nsubruns)+"/meanpt_meanv2.txt";
			exists=filesystem::exists(filename);
			if(exists){
				subrunfilenames1.push_back(filename);
				Nsubruns+=1;
			}
		}while(exists);
		if(Nsubruns>Nsubruns_max)
			Nsubruns=Nsubruns_max;
	}
	printf("Nsubruns=%u\n",Nsubruns);
	
	//
	
	// Sum over subruns
	
	spectra_sum.resize(3);
	v2_sum.resize(3);
	
		
	for(irun=I;irun<=J;irun++){
		for(int i=0;i<3;i++){
			spectra_sum[i].clear();
			v2_sum[i].clear();
		}
		
		for(isubrun=0;isubrun<Nsubruns;isubrun++){
			dirname="modelruns/run"+to_string(irun)+"/"+qualifier+"/results_spectrav2/subruns/subrun"+to_string(isubrun);
			filename=dirname+"/spectra.txt";
			fptr=fopen(filename.c_str(),"r");
			ipt=0;
			while(!feof(fptr)){
				fscanf(fptr,"%lf %lf %lf %lf",&pt,&specpi,&specK,&specp);
				if(!feof(fptr)){
					if(isubrun==0){
						spectra_sum[0].push_back(specpi);
						spectra_sum[1].push_back(specK);
						spectra_sum[2].push_back(specp);			
					}
					else{
						spectra_sum[0][ipt]+=specpi;
						spectra_sum[1][ipt]+=specK;
						spectra_sum[2][ipt]+=specp;
					}
				}
				ipt+=1;
			}
			fclose(fptr);
		}
		for(isubrun=0;isubrun<Nsubruns;isubrun++){
			dirname="modelruns/run"+to_string(irun)+"/"+qualifier+"/results_spectrav2/subruns/subrun"+to_string(isubrun);
			filename=dirname+"/v2.txt";
			fptr=fopen(filename.c_str(),"r");
			ipt=0;
			while(!feof(fptr)){
				fscanf(fptr,"%lf %lf %lf %lf",&pt,&v2pi,&v2K,&v2p);
				if(!feof(fptr)){
					if(isubrun==0){
						v2_sum[0].push_back(v2pi);
						v2_sum[1].push_back(v2K);
						v2_sum[2].push_back(v2p);			
					}
					else{
						v2_sum[0][ipt]+=v2pi;
						v2_sum[1][ipt]+=v2K;
						v2_sum[2][ipt]+=v2p;
					}
				}
				ipt+=1;
			}
			fclose(fptr);
		}
		
		meanpt_sum.clear();
		meanv2_sum.clear();
		meanpt_sum.resize(3);
		meanv2_sum.resize(3);
		for(isubrun=0;isubrun<Nsubruns;isubrun++){
			dirname="modelruns/run"+to_string(irun)+"/"+qualifier+"/results_spectrav2/subruns/subrun"+to_string(isubrun);
			filename=dirname+"/meanpt_meanv2.txt";
			fptr=fopen(filename.c_str(),"r");
			
			fscanf(fptr,"%s %lf",dummy,&pt);
			meanlabel[0]=dummy;		
			meanpt_sum[0]=meanpt_sum[0]+pt;
			fscanf(fptr,"%s %lf",dummy,&pt);
			meanlabel[1]=dummy;
			meanpt_sum[1]=meanpt_sum[1]+pt;
			fscanf(fptr,"%s %lf",dummy,&pt);
			meanlabel[2]=dummy;
			meanpt_sum[2]=meanpt_sum[2]+pt;
			
			fscanf(fptr,"%s %lf",dummy,&v2pi);
			meanlabel[3]=dummy;
			meanv2_sum[0]=meanv2_sum[0]+v2pi;
			fscanf(fptr,"%s %lf",dummy,&v2K);
			meanlabel[4]=dummy;
			meanv2_sum[1]=meanv2_sum[1]+v2K;
			fscanf(fptr,"%s %lf",dummy,&v2p);
			meanlabel[5]=dummy;
			meanv2_sum[2]=meanv2_sum[2]+v2p;
			
			fclose(fptr);
		}
			
			
		
		// Write summed results

		
		dirname="modelruns/run"+to_string(irun)+"/"+qualifier+"/results_spectrav2";
		filename=dirname+"/spectra.txt";
		fptr=fopen(filename.c_str(),"w");
		for(ipt=0;ipt<spectra_sum[0].size();ipt++){
			pt=(ipt+0.5)*DELPT_SPECTRA;
			fprintf(fptr,"%8.4f %12.5e %12.5e %12.5e\n",pt,
			spectra_sum[0][ipt]/double(Nsubruns),spectra_sum[1][ipt]/double(Nsubruns),spectra_sum[2][ipt]/double(Nsubruns));
		}
		fclose(fptr);
		
		filename=dirname+"/v2.txt";
		fptr=fopen(filename.c_str(),"w");
		for(ipt=0;ipt<v2_sum[0].size();ipt++){
			pt=(ipt+0.5)*DELPT_V2;
			fprintf(fptr,"%8.4f %12.5e %12.5e %12.5e\n",pt,
			v2_sum[0][ipt]/double(Nsubruns),v2_sum[1][ipt]/double(Nsubruns),v2_sum[2][ipt]/double(Nsubruns));
		}
		fclose(fptr);	
		
		dirname="modelruns/run"+to_string(irun)+"/"+qualifier+"/results_spectrav2";
		filename=dirname+"/meanpt_meanv2.txt";
		fptr=fopen(filename.c_str(),"w");
		fprintf(fptr,"%s  %12.5e\n",meanlabel[0].c_str(),meanpt_sum[0]/double(Nsubruns));
		fprintf(fptr,"%s  %12.5e\n",meanlabel[1].c_str(),meanpt_sum[1]/double(Nsubruns));
		fprintf(fptr,"%s  %12.5e\n",meanlabel[2].c_str(),meanpt_sum[2]/double(Nsubruns));
		fprintf(fptr,"%s  %12.5e\n",meanlabel[3].c_str(),meanv2_sum[0]/double(Nsubruns));
		fprintf(fptr,"%s  %12.5e\n",meanlabel[4].c_str(),meanv2_sum[1]/double(Nsubruns));
		fprintf(fptr,"%s  %12.5e\n",meanlabel[5].c_str(),meanv2_sum[2]/double(Nsubruns));
		fclose(fptr);
		
		
	}
	
	
	
		
	return 0;
}


