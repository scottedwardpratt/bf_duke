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
	const int nspecies=10,nbfs=7;
	bool exists,exists1,exists2;
	string fn,fn1,fn2,fn1sum,fnsum,dirname,command;
	int isubrun,Nsubruns,Nsubruns1=0,Nsubruns2=0,ibf,ispecies,ix,nx,itype,irun;
	vector<string> subrunfilenames1,subrunfilenames2;
	FILE *fptr1,*fptr2,*fptrsum;
	double dx,bf,cf;
	vector<double> delx,bfsum,cfsum;
	double bfqinv,bfout,bfside,bflong,cfqinv,cfout,cfside,cflong;
	vector<double> bfqinvsum,bfoutsum,bfsidesum,bflongsum,cfqinvsum,cfoutsum,cfsidesum,cflongsum;
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
	int I=atoi(argv[2]);
	int J=atoi(argv[3]);
	string bfspecies[nspecies]={"KK","Kp","allcharges","allcharges_phi0","allcharges_phi45","allcharges_phi90","piK","pip","pipi","pp"};
	string bffilenames[nbfs]={"bf_qinv.txt","bf_eta.txt","bf_etas.txt","bf_y1.txt",
	"bf_eta1.txt","bf_phi.txt","bf_y.txt"};


	int dotypes=12;
	printf("Enter 1 to sum type I only, 2 to sum type II only, anything else to do both\n");
	scanf("%d",&dotypes);
	if(dotypes!=1 && dotypes!=2)
		dotypes=12;

	for(irun=I;irun<=J;irun++){

		// Get Nsubruns
		if(dotypes!=2){
			Nsubruns1=0;
			do{
				fn="modelruns/run"+to_string(irun)+"/"+qualifier+"/results_type1/subruns/subrun"+to_string(Nsubruns1)+"/denom.txt";
				exists=filesystem::exists(fn);
				if(exists){
					subrunfilenames1.push_back(fn);
					Nsubruns1+=1;
				}
			}while(exists);
			if(Nsubruns1>Nsubruns_max)
				Nsubruns1=Nsubruns_max;
		}
	
		if(dotypes!=1){
			Nsubruns2=0;
			do{
				fn="modelruns/run"+to_string(irun)+"/"+qualifier+"/results_type2/subruns/subrun"+to_string(Nsubruns2)+"/denom.txt";
				exists=filesystem::exists(fn);
				if(exists){
					subrunfilenames2.push_back(fn);
					Nsubruns2+=1;
				}
			}while(exists);
			if(Nsubruns2>Nsubruns_max)
				Nsubruns2=Nsubruns_max;
		}
		//
	
		// Add type 1 BFs, then add type 2 BFs
		int itype_first=1,itype_last=2;
		if(dotypes==1)
			itype_last=1;
		if(dotypes==2)
			itype_first=2;
		for(itype=itype_first;itype<=itype_last;itype++){
			for(ispecies=0;ispecies<nspecies;ispecies++){
				for(ibf=0;ibf<nbfs;ibf++){
					bfsum.clear();
					cfsum.clear();
					bfqinvsum.clear();
					bfoutsum.clear();
					bfsidesum.clear();
					bflongsum.clear();
					cfqinvsum.clear();
					cfoutsum.clear();
					cfsidesum.clear();
					cflongsum.clear();
					delx.clear();
					fn1="modelruns/run"+to_string(irun)+"/"+qualifier+"/results_type"+to_string(itype)+"/subruns/subrun"+to_string(irun)+"/"+bfspecies[ispecies]+"/"+bffilenames[ibf];
					exists1=filesystem::exists(fn1);
					if(exists1){
						dirname="modelruns/run"+to_string(irun)+"/"+qualifier+"/results_type"+to_string(itype)+"/"+bfspecies[ispecies]+"/";
						if(!filesystem::exists(dirname)){
							command="mkdir -p "+dirname;
							system(command.c_str());
						}
						fnsum="modelruns/run"+to_string(irun)+"/"+qualifier+"/results_type"+to_string(itype)+"/"+bfspecies[ispecies]+"/"+bffilenames[ibf];
						fptrsum=fopen(fnsum.c_str(),"w");
						if(itype==1)
							Nsubruns=Nsubruns1;
						else
							Nsubruns=Nsubruns2;
						for(isubrun=0;isubrun<Nsubruns;isubrun++){
							fn1="modelruns/run"+to_string(irun)+"/"+qualifier+"/results_type"+to_string(itype)+"/subruns/subrun"+to_string(isubrun)+"/"+bfspecies[ispecies]+"/"+bffilenames[ibf];
							fptr1=fopen(fn1.c_str(),"r");
							ix=0;
							while(!feof(fptr1)){
								if(bffilenames[ibf]!="bf_qinv.txt"){
									fscanf(fptr1,"%lf %lf %lf",&dx,&bf,&cf);
									if(!feof(fptr1)){
										if(isubrun==0){
											delx.push_back(dx);
											bfsum.push_back(bf);
											cfsum.push_back(cf);
										}
										else{
											bfsum[ix]+=bf;
											cfsum[ix]+=cf;
											ix+=1;
										}
									}
								}
								else{

									fscanf(fptr1,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",
									&dx,&bfqinv,&bfout,&bfside,&bflong,&cfqinv,&cfout,&cfside,&cflong);
									if(!feof(fptr1)){
										if(isubrun==0){
											delx.push_back(dx);
											bfqinvsum.push_back(bfqinv);
											bfoutsum.push_back(bfout);
											bfsidesum.push_back(bfside);
											bflongsum.push_back(bflong);
											cfqinvsum.push_back(cfqinv);
											cfoutsum.push_back(cfout);
											cfsidesum.push_back(cfside);
											cflongsum.push_back(cflong);
										}
										else{
											bfqinvsum[ix]+=bfqinv;
											bfoutsum[ix]+=bfout;
											bfsidesum[ix]+=bfside;
											bflongsum[ix]+=bflong;
											cfqinvsum[ix]+=cfqinv;
											cfoutsum[ix]+=cfout;
											cfsidesum[ix]+=cfside;
											cflongsum[ix]+=cflong;
											ix+=1;
										}
									}
								}
							}
							fclose(fptr1);
							nx=bfqinvsum.size();
						}
						nx=delx.size();
						for(ix=0;ix<nx;ix++){
							if(bffilenames[ibf]!="bf_qinv.txt"){
								fprintf(fptrsum,"%7.2f %11.4e %11.4e\n",
								delx[ix],bfsum[ix]/double(Nsubruns),cfsum[ix]/double(Nsubruns));
							}
							else{
								fprintf(fptrsum,"%7.2f %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e\n",delx[ix],
								bfqinvsum[ix]/double(Nsubruns),bfoutsum[ix]/double(Nsubruns),
								bfsidesum[ix]/double(Nsubruns),bflongsum[ix]/double(Nsubruns),
								cfqinvsum[ix]/double(Nsubruns),cfoutsum[ix]/double(Nsubruns),
								cfsidesum[ix]/double(Nsubruns),cflongsum[ix]/double(Nsubruns));
							}
						}
						fclose(fptrsum);
					}
				}
			}
		}
	
		// Add type 1 to type 2 BFs
		if(dotypes!=1 && dotypes!=2){
			for(ispecies=0;ispecies<nspecies;ispecies++){
				for(ibf=0;ibf<nbfs;ibf++){
					fn1="modelruns/run"+to_string(irun)+"/"+qualifier+"/results_type1/"+bfspecies[ispecies]+"/"+bffilenames[ibf];
					fn2="modelruns/run"+to_string(irun)+"/"+qualifier+"/results_type2/"+bfspecies[ispecies]+"/"+bffilenames[ibf];
					fnsum="modelruns/run"+to_string(irun)+"/"+qualifier+"/results/"+bfspecies[ispecies]+"/"+bffilenames[ibf];
					exists1=filesystem::exists(fn1);
					exists2=filesystem::exists(fn2);
					if(exists1 && exists2){
						dirname="modelruns/run"+to_string(irun)+"/"+qualifier+"/results/"+bfspecies[ispecies]+"/";
						if(!filesystem::exists(dirname)){
							command="mkdir -p "+dirname;
							system(command.c_str());
						}
						bfsum.clear();
						cfsum.clear();
						bfqinvsum.clear();
						bfoutsum.clear();
						bfsidesum.clear();
						bflongsum.clear();
						cfqinvsum.clear();
						cfoutsum.clear();
						cfsidesum.clear();
						cflongsum.clear();
						delx.clear();
			
						if(bffilenames[ibf]!="bf_qinv.txt"){
							fptr1=fopen(fn1.c_str(),"r");
							while(!feof(fptr1)){
								fscanf(fptr1,"%lf %lf %lf",&dx,&bf,&cf);
								if(!feof(fptr1)){
									delx.push_back(dx);
									bfsum.push_back(bf);
									cfsum.push_back(cf);
								}
							}
							fclose(fptr1);
			
							fptr2=fopen(fn2.c_str(),"r");
							nx=delx.size();
							for(ix=0;ix<nx;ix++){
								fscanf(fptr1,"%lf %lf %lf",&dx,&bf,&cf);
								bfsum[ix]+=bf;
								cfsum[ix]+=cf;
							}
							fclose(fptr2);
						}
						else{
							fptr1=fopen(fn1.c_str(),"r");
							while(!feof(fptr1)){
								fscanf(fptr1,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",
								&dx,&bfqinv,&bfout,&bfside,&bflong,&cfqinv,&cfout,&cfside,&cflong);
								if(!feof(fptr1)){
									delx.push_back(dx);
									bfqinvsum.push_back(bfqinv);
									bfoutsum.push_back(bfout);
									bfsidesum.push_back(bfside);
									bflongsum.push_back(bflong);
									cfqinvsum.push_back(cfqinv);
									cfoutsum.push_back(cfout);
									cfsidesum.push_back(cfside);
									cflongsum.push_back(cflong);
								}
							}
							fclose(fptr1);
		
							fptr2=fopen(fn2.c_str(),"r");
							nx=delx.size();
							for(ix=0;ix<nx;ix++){
								fscanf(fptr2,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",
								&dx,&bfqinv,&bfout,&bfside,&bflong,&cfqinv,&cfout,&cfside,&cflong);
								bfqinvsum[ix]+=bfqinv;
								bfoutsum[ix]+=bfout;
								bfsidesum[ix]+=bfside;
								bflongsum[ix]+=bflong;
								cfqinvsum[ix]+=cfqinv;
								cfoutsum[ix]+=cfout;
								cfsidesum[ix]+=cfside;
								cflongsum[ix]+=cflong;
							
							}
							fclose(fptr2);
						}

				
						nx=delx.size();
					
						fptrsum=fopen(fnsum.c_str(),"w");
						for(ix=0;ix<nx;ix++){
							if(bffilenames[ibf]!="bf_qinv.txt"){
								fprintf(fptrsum,"%7.2f %11.4e %11.4e\n",delx[ix],bfsum[ix],cfsum[ix]);
							}
							else{
								fprintf(fptrsum,"%7.2f %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e\n",delx[ix],
								bfqinvsum[ix],bfoutsum[ix],bfsidesum[ix],bflongsum[ix],
								cfqinvsum[ix],cfoutsum[ix],cfsidesum[ix],cflongsum[ix]);
							}
						}
						fclose(fptrsum);
					}
				}
			}
		}
	}
	return 0;
}


