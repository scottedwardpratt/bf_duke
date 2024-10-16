#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include <string>
#include <cstring>
#include <iostream>
#include <sstream>
#include <filesystem>

const double PI=4.0*atan(1.0);
const double HBARC=197.3269602;

using namespace std;

int main(int argc,char *argv[]){
	int Nsubruns_max=99999999;
	const int nspecies=10,nbfs=9;
	bool exists,exists1,exists2;
	string fn,fn1,fn2,fn1sum,fn2sum,fnsum,dirname,command;
	int isubrun,Nsubruns,Nsubruns1,Nsubruns2,ibf,ispecies,ix,nx,itype;
	vector<string> subrunfilenames1,subrunfilenames2;
	FILE *fptr1,*fptr2,*fptrsum;
	vector<double> delx,bfsum,sigmasum;
	double dx,bf,sigma;
	if(argc!=4){
		printf("Usage: summer qualifier I J are going from modelruns/runI to modelruns/runJ");
		exit(1);
	}
	exists=filesystem::exists("modelruns");
	if(!exists){
		printf("this needs to be run from parent directory to modelruns/\n");
		exit(1);
	}
	string qualifier="argv[1]";
	int I=atoi(argv[2]);
	int J=atoi(argv[3]);
	string bfspecies[nspecies]={"KK","Kp","allcharges","allcharges_phi0","allcharges_phi45","allcharges_phi90","piK","pip","pipi","pp"};
	string bffilenames[nbfs]={"bf_Byphi.txt","bf_eta.txt","bf_etas.txt","bf_qinv.txt","bf_y1.txt",
	"bf_Cyphi.txt","bf_eta1.txt","bf_phi.txt","bf_y.txt"};
	
	// Get Nsubruns
	Nsubruns1=0;
	do{
		fn="modelruns/run0/alicePbPb_cent0_5/results_itype1/subruns/subrun"+to_string(Nsubruns1)+"/denom.txt";
		exists=filesystem::exists(fn);
		if(exists){
			subrunfilenames1.push_back(fn);
			Nsubruns1+=1;
		}
	}while(exists);
	if(Nsubruns1>Nsubruns_max)
		Nsubruns1=Nsubruns_max;
	Nsubruns2=0;
	do{
		fn="modelruns/run0/alicePbPb_cent0_5/results_itype2/subruns/subrun"+to_string(Nsubruns2)+"/denom.txt";
		exists=filesystem::exists(fn);
		if(exists){
			subrunfilenames2.push_back(fn);
			Nsubruns2+=1;
		}
	}while(exists);
	if(Nsubruns2>Nsubruns_max)
		Nsubruns2=Nsubruns_max;
	//
	
	// Add up type 1 BFs then add type 2 BFs
	for(itype=1;itype<=2;itype++){
		for(ispecies=0;ispecies<nspecies;ispecies++){
			for(ibf=0;ibf<nbfs;ibf++){
				bfsum.clear();
				sigmasum.clear();
				delx.clear();
				fn1="modelruns/run0/alicePbPb_cent0_5/results_itype"+to_string(itype)+"/subruns/subrun0/"+bfspecies[ispecies]+"/"+bffilenames[ibf];
				exists1=filesystem::exists(fn1);
				if(exists1){
					dirname="modelruns/run0/alicePbPb_cent0_5/results_itype"+to_string(itype)+"/"+bfspecies[ispecies]+"/";
					if(!filesystem::exists(dirname)){
						command="mkdir -p "+dirname;
						system(command.c_str());
					}
					fn1sum="modelruns/run0/alicePbPb_cent0_5/results_itype"+to_string(itype)+"/"+bfspecies[ispecies]+"/"+bffilenames[ibf];
					fptrsum=fopen(fn1sum.c_str(),"w");
					if(itype==1)
						Nsubruns=Nsubruns1;
					else
						Nsubruns=Nsubruns2;
					for(isubrun=0;isubrun<Nsubruns;isubrun++){
						fn1="modelruns/run0/alicePbPb_cent0_5/results_itype"+to_string(itype)+"/subruns/subrun"+to_string(isubrun)+"/"+bfspecies[ispecies]+"/"+bffilenames[ibf];
						fptr1=fopen(fn1.c_str(),"r");
						ix=0;
						while(!feof(fptr1)){
							fscanf(fptr1,"%lf %lf %lf",&dx,&bf,&sigma);
							if(!feof(fptr1)){
								if(isubrun==0){
									delx.push_back(dx);
									bfsum.push_back(bf);
									sigmasum.push_back(sigma);
								}
								else{
									bfsum[ix]+=bf;
									sigmasum[ix]+=sigma;
									ix+=1;
								}
							}
						}
					}
					fclose(fptr1);
					nx=delx.size();
					for(ix=0;ix<nx;ix++){
						fprintf(fptrsum,"%7.2f %11.4e %11.4e\n",delx[ix],bfsum[ix]/double(nx),sigmasum[ix]/double(nx));
					}
				}
			}
		}
	}
	
	// Add type 1 to type 2 BFs
	for(ispecies=0;ispecies<nspecies;ispecies++){
		for(ibf=0;ibf<nbfs;ibf++){
			bfsum.clear();
			sigmasum.clear();
			delx.clear();
			fn1="modelruns/run0/alicePbPb_cent0_5/results_itype1/"+bfspecies[ispecies]+"/"+bffilenames[ibf];
			fn2="modelruns/run0/alicePbPb_cent0_5/results_itype2/"+bfspecies[ispecies]+"/"+bffilenames[ibf];
			fnsum="modelruns/run0/alicePbPb_cent0_5/results/"+bfspecies[ispecies]+"/"+bffilenames[ibf];
			exists1=filesystem::exists(fn1);
			exists2=filesystem::exists(fn2);
			if(exists1 && exists2){
				dirname="modelruns/run0/alicePbPb_cent0_5/results/"+bfspecies[ispecies]+"/";
				if(!filesystem::exists(dirname)){
					command="mkdir -p "+dirname;
					system(command.c_str());
				}
				bfsum.clear();
				sigmasum.clear();
				delx.clear();
				
				fptr1=fopen(fn1.c_str(),"r");
				while(!feof(fptr1)){
					fscanf(fptr1,"%lf %lf %lf",&dx,&bf,&sigma);
					if(!feof(fptr1)){
						delx.push_back(dx);
						bfsum.push_back(bf);
						sigmasum.push_back(sigma);
					}
				}
				fclose(fptr1);
				
				fptr2=fopen(fn2.c_str(),"r");
				nx=delx.size();
				for(ix=0;ix<nx;ix++){
					fscanf(fptr1,"%lf %lf %lf",&dx,&bf,&sigma);
					bfsum[ix]+=bf;
					sigmasum[ix]+=sigma;
				}
			}
			fclose(fptr2);
				
			fptrsum=fopen(fnsum.c_str(),"w");
			
			for(ix=0;ix<nx;ix++){
				fprintf(fptrsum,"%7.2f %11.4e %11.4e\n",delx[ix],bfsum[ix]/double(nx),sigmasum[ix]/double(nx));
			}
		}
		fclose(fptrsum);
	}
	return 0;
}


