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
	string qualifier[4]={"alicePbPb_cent0_5","alicePbPb_cent5_10","alicePbPb_cent10_15","alicePbPb_cent15_20"};
	string denom_string;
	char denom_read[120];
	double dcountread,subrunweight,subruncfweight;
	double denomcount_pi_1,denomcount_K_1,denomcount_p_1,denomcount_allcharges_1;
	double denomcount_allcharges_phi0_1,denomcount_allcharges_phi45_1,denomcount_allcharges_phi90_1;
	double denomcount_pi_2,denomcount_K_2,denomcount_p_2,denomcount_allcharges_2;
	double denomcount_allcharges_phi0_2,denomcount_allcharges_phi45_2,denomcount_allcharges_phi90_2;
	vector<double> w_pi_1,w_K_1,w_p_1,w_allcharges_1;
	vector<double> w_allcharges_phi0_1,w_allcharges_phi45_1,w_allcharges_phi90_1;
	vector<double> w_pi_2,w_K_2,w_p_2,w_allcharges_2;
	vector<double> w_allcharges_phi0_2,w_allcharges_phi45_2,w_allcharges_phi90_2;
	vector<double> wcf_pi_1,wcf_K_1,wcf_p_1,wcf_allcharges_1;
	vector<double> wcf_allcharges_phi0_1,wcf_allcharges_phi45_1,wcf_allcharges_phi90_1;
	vector<double> wcf_pi_2,wcf_K_2,wcf_p_2,wcf_allcharges_2;
	vector<double> wcf_allcharges_phi0_2,wcf_allcharges_phi45_2,wcf_allcharges_phi90_2;
	bool exists,exists1,exists2;
	string fn,fn1,fn2,fn1sum,fnsum,dirname,command;
	int isubrun,Nsubruns,Nsubruns1=0,Nsubruns2=0,ibf,ispecies,ix,nx,itype,irun,iqual,dotypes;
	vector<string> subrunfilenames1,subrunfilenames2;
	FILE *fptr1,*fptr2,*fptrsum,*fptr_denom;
	double dx,bf,cf;
	vector<double> delx,bfsum,cfsum;
	double bfqinv,bfout,bfside,bflong,cfqinv,cfout,cfside,cflong;
	vector<double> bfqinvsum,bfoutsum,bfsidesum,bflongsum,cfqinvsum,cfoutsum,cfsidesum,cflongsum;
	exists=filesystem::exists("modelruns");
	if(!exists){
		printf("this needs to be run from parent directory to modelruns/\n");
		exit(1);
	}
	string bfspecies[nspecies]={"KK","Kp","allcharges","allcharges_phi0","allcharges_phi45","allcharges_phi90","piK","pip","pipi","pp"};
	string bffilenames[nbfs]={"bf_qinv.txt","bf_eta.txt","bf_etas.txt","bf_y1.txt",
	"bf_eta1.txt","bf_phi.txt","bf_y.txt"};


	dotypes=12;
	//printf("Enter 1 to sum type I only, 2 to sum type II only, anything else to do both\n");
	//scanf("%d",&dotypes);
	if(dotypes!=1 && dotypes!=2)
		dotypes=12;

	for(irun=0;irun<6;irun++){
		for(iqual=0;iqual<4;iqual++){
			
			// Get Nsubruns
			if(dotypes!=2){
				Nsubruns1=0;
				denomcount_pi_1=denomcount_K_1=denomcount_p_1=denomcount_allcharges_1=0;
				denomcount_allcharges_phi0_1=denomcount_allcharges_phi45_1=denomcount_allcharges_phi90_1=0.0;
				do{
					fn="modelruns/run"+to_string(irun)+"/"+qualifier[iqual]+"/results_type1/subruns/subrun"+to_string(Nsubruns1)+"/denom.txt";
					exists=filesystem::exists(fn);
					if(exists){
						subrunfilenames1.push_back(fn);
						if(Nsubruns1<Nsubruns_max){
							fptr_denom=fopen(fn.c_str(),"r");
							while(!feof(fptr_denom)){
								fscanf(fptr_denom,"%s %lf",denom_read,&dcountread);
								if(!feof(fptr_denom)){
									denom_string=denom_read;
									if(denom_string=="denom_pi:"){
										w_pi_1.push_back(dcountread);
										wcf_pi_1.push_back(dcountread*dcountread);
										denomcount_pi_1+=dcountread;
									}
									if(denom_string=="denom_K:"){
										w_K_1.push_back(dcountread);
										wcf_K_1.push_back(dcountread*dcountread);
										denomcount_K_1+=dcountread;
									}
									if(denom_string=="denom_p:"){
										w_p_1.push_back(dcountread);
										wcf_p_1.push_back(dcountread*dcountread);
										denomcount_p_1+=dcountread;
									}
									if(denom_string=="denom_allcharges:"){
										w_allcharges_1.push_back(dcountread);
										wcf_allcharges_1.push_back(dcountread*dcountread);
										denomcount_allcharges_1+=dcountread;
									}
									if(denom_string=="denom_allcharges_phi0:"){
										w_allcharges_phi0_1.push_back(dcountread);
										wcf_allcharges_phi0_1.push_back(dcountread*dcountread);
										denomcount_allcharges_phi0_1+=dcountread;
									}
									if(denom_string=="denom_allcharges_phi45:"){
										w_allcharges_phi45_1.push_back(dcountread);
										wcf_allcharges_phi45_1.push_back(dcountread*dcountread);
										denomcount_allcharges_phi45_1+=dcountread;
									}
									if(denom_string=="denom_allcharges_phi90:"){
										w_allcharges_phi90_1.push_back(dcountread);
										wc_allcharges_phi90_1.push_back((dcountread*dcountread));
										denomcount_allcharges_phi90_1+=dcountread;
									}
									else{
										printf("denom_string not recognized = %s\n",denom_read);
									}
									
								}
							}
							fclose(fptr_denom);
						}
						Nsubruns1+=1;
					}
				}while(exists);
				if(Nsubruns1>Nsubruns_max)
					Nsubruns1=Nsubruns_max;
				for(isubrun=0;isubrun<Nsubruns2;isubrun++){
					w_pi_1[isubrun]=w_pi_1[isubrun]/denomcount_pi_1;
					w_K_1[isubrun]=w_K_1[isubrun]/denomcount_K_1;
					w_p_1[isubrun]=w_p_1[isubrun]/denomcount_p_1;
					w_allcharges_1[isubrun]=w_allcharges_1[isubrun]/denomcount_allcharges_1;
					w_allcharges_phi0_1[isubrun]=w_allcharges_phi0_1[isubrun]/denomcount_allcharges_phi0_1;
					w_allcharges_phi45_1[isubrun]=w_allcharges_phi45_1[isubrun]/denomcount_allcharges_phi45_1;
					w_allcharges_phi90_1[isubrun]=w_allcharges_phi90_1[isubrun]/denomcount_allcharges_phi90_1;
					
				}
			}
	
			if(dotypes!=1){
				Nsubruns2=0;
				denomcount_pi_2=denomcount_K_2=denomcount_p_2=denomcount_allcharges_2=0;
				denomcount_allcharges_phi0_2=denomcount_allcharges_phi45_2=denomcount_allcharges_phi90_2=0.0;
				do{
					fn="modelruns/run"+to_string(irun)+"/"+qualifier[iqual]+"/results_type2/subruns/subrun"+to_string(Nsubruns2)+"/denom.txt";
					exists=filesystem::exists(fn);
					if(exists){
						subrunfilenames2.push_back(fn);
						if(Nsubruns2<Nsubruns_max){
							fptr_denom=fopen(fn.c_str(),"r");
							while(!feof(fptr_denom)){
								fscanf(fptr_denom,"%s %lf",denom_read,&dcountread);
								if(!feof(fptr_denom)){
									denom_string=denom_read;
									if(denom_string=="denom_pi:"){
										w_pi_2.push_back(dcountread);
										wcf_pi_2.push_back((dcountread*dcountread));
										denomcount_pi_2+=dcountread;
									}
									if(denom_string=="denom_K:"){
										w_K_2.push_back(dcountread);
										wcf_K_2.push_back((dcountread*dcountread));
										denomcount_K_2+=dcountread;
									}
									if(denom_string=="denom_p:"){
										w_pi_2.push_back(dcountread);
										wcf_pi_2.push_back((dcountread*dcountread));
										denomcount_p_2+=dcountread;
									}
									if(denom_string=="denom_allcharges:"){
										w_allcharges_2.push_back(dcountread);
										wcf_allcharges_2.push_back((dcountread*dcountread));
										denomcount_allcharges_2+=dcountread;
									}
									if(denom_string=="denom_allcharges_phi0:"){
										w_allcharges_phi0_2.push_back(dcountread);
										wcf_allcharges_phi0_2.push_back(dcountread*dcountread);
										denomcount_allcharges_phi0_2+=dcountread;
									}
									if(denom_string=="denom_allcharges_phi45:"){
										w_allcharges_phi45_2.push_back(dcountread);
										wcf_allcharges_phi45_2.push_back(dcountread*dcountread);
										denomcount_allcharges_phi45_2+=dcountread;
									}
									if(denom_string=="denom_allcharges_phi90:"){
										w_allcharges_phi90_2.push_back(dcountread);
										wcf_allcharges_phi90_2.push_back(dcountread*dcountread);
										denomcount_allcharges_phi90_2+=dcountread;
									}
									else{
										printf("denom_string not recognized = %s\n",denom_read);
									}
									
								}
							}
							fclose(fptr_denom);
						}
						Nsubruns2+=1;
					}
				}while(exists);
				if(Nsubruns2>Nsubruns_max)
					Nsubruns2=Nsubruns_max;
				for(isubrun=0;isubrun<Nsubruns2;isubrun++){
					w_pi_2[isubrun]=w_pi_2[isubrun]/denomcount_pi_2;
					w_K_2[isubrun]=w_K_2[isubrun]/denomcount_K_2;
					w_p_2[isubrun]=w_p_2[isubrun]/denomcount_p_2;
					w_allcharges_2[isubrun]=w_allcharges_2[isubrun]/denomcount_allcharges_2;
					w_allcharges_phi0_2[isubrun]=w_allcharges_phi0_2[isubrun]/denomcount_allcharges_phi0_2;
					w_allcharges_phi45_2[isubrun]=w_allcharges_phi45_2[isubrun]/denomcount_allcharges_phi45_2;
					w_allcharges_phi90_2[isubrun]=w_allcharges_phi90_2[isubrun]/denomcount_allcharges_phi90_2;
					wcf_pi_2[isubrun]=wcf_pi_2[isubrun]/(denomcount_pi_2*denomcount_pi_2);
					wcf_K_2[isubrun]=wcf_K_2[isubrun]/(denomcount_K_2*denomcount_K_2);
					wcf_p_2[isubrun]=wcf_p_2[isubrun]/(denomcount_p_2*denomcount_p_2);
					wcf_allcharges_2[isubrun]=wcf_allcharges_2[isubrun]/(denomcount_allcharges_2*denomcount_allcharges_2);
					wcf_allcharges_phi0_2[isubrun]=wcf_allcharges_phi0_2[isubrun]/(denomcount_allcharges_phi0_2*denomcount_allcharges_phi0_2);
					wcf_allcharges_phi45_2[isubrun]=wcf_allcharges_phi45_2[isubrun]/(denomcount_allcharges_phi45_2*denomcount_allcharges_phi45_2);
					wcf_allcharges_phi90_2[isubrun]=wcf_allcharges_phi90_2[isubrun]/(denomcount_allcharges_phi90_2*denomcount_allcharges_phi90_2);
					
				}
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
					
						double norm,width,normbar,widthbar,sigma2norm,sigma2width;
						normbar=widthbar=sigma2norm=sigma2width=0.0;
				
					
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
						fn1="modelruns/run"+to_string(irun)+"/"+qualifier[iqual]+"/results_type"+to_string(itype)+"/subruns/subrun"+to_string(irun)+"/"+bfspecies[ispecies]+"/"+bffilenames[ibf];
						exists1=filesystem::exists(fn1);
						if(exists1){
							dirname="modelruns/run"+to_string(irun)+"/"+qualifier[iqual]+"/results_type"+to_string(itype)+"/"+bfspecies[ispecies]+"/";
							if(!filesystem::exists(dirname)){
								command="mkdir -p "+dirname;
								system(command.c_str());
							}
							fnsum="modelruns/run"+to_string(irun)+"/"+qualifier[iqual]+"/results_type"+to_string(itype)+"/"+bfspecies[ispecies]+"/"+bffilenames[ibf];
							fptrsum=fopen(fnsum.c_str(),"w");
							if(itype==1)
								Nsubruns=Nsubruns1;
							else
								Nsubruns=Nsubruns2;
							for(isubrun=0;isubrun<Nsubruns;isubrun++){
								
								//"KK","Kp","allcharges","allcharges_phi0","allcharges_phi45","allcharges_phi90","piK","pip","pipi","pp"
								if(itype==1){
									if(ispecies==0){
										subrunweight=w_K_1[isubrun];
										subruncfweight=wcf_K_1[isubrun];
									}
									else if(ispecies==1){
										subrunweight=w_p_1[isubrun];
										subruncfweight=wcf_p_1[isubrun];
									}
									else if(ispecies==2){
										subrunweight=w_allcharges_1[isubrun];
										subruncfweight=wcf_allcharges_1[isubrun];
									}
									else if(ispecies==3){
										subrunweight=w_allcharges_phi0_1[isubrun];
										subruncfweight=wcf_allcharges_phi0_1[isubrun];
									}
									else if(ispecies==4){
										subrunweight=w_allcharges_phi45_1[isubrun];
										subruncfweight=wcf_allcharges_phi45_1[isubrun];
									}
									else if(ispecies==5){
										subrunweight=w_allcharges_phi90_1[isubrun];
										subruncfweight=wcf_allcharges_phi90_1[isubrun];
									}
									else if(ispecies==6){
										subrunweight=w_K_1[isubrun];
										subruncfweight=wcf_K_1[isubrun];
									}
									else if(ispecies==7){
										subrunweight=w_p_1[isubrun];
										subruncfweight=wcf_p_1[isubrun];
									}
									else if(ispecies==8){
										subrunweight=w_pi_1[isubrun];
										subruncfweight=wcf_pi_1[isubrun];
									}
									else if(ispecies==9){
										subrunweight=w_p_1[isubrun];
										subruncfweight=wcf_p_1[isubrun];
									}
								}
								if(itype==2){
									if(ispecies==0){
										subrunweight=w_K_2[isubrun];
										subruncfweight=wcf_K_2[isubrun];
									}
									else if(ispecies==1){
										subrunweight=w_p_2[isubrun];
										subcfrunweight=wcf_p_2[isubrun];
									}
									else if(ispecies==2){
										subrunweight=w_allcharges_2[isubrun];
										subruncfweight=wcf_allcharges_2[isubrun];
									}
									else if(ispecies==3){
										subrunweight=w_allcharges_phi0_2[isubrun];
										subruncfweight=wcf_allcharges_phi0_2[isubrun];
									}
									else if(ispecies==4){
										subrunweight=w_allcharges_phi45_2[isubrun];
										subruncfweight=wcf_allcharges_phi45_2[isubrun];
									}
									else if(ispecies==5){
										subrunweight=w_allcharges_phi90_2[isubrun];
										subruncfweight=w_cfallcharges_phi90_2[isubrun];
									}
									else if(ispecies==6){
										subrunweight=w_K_2[isubrun];
										subruncfweight=wcf_K_2[isubrun];
									}
									else if(ispecies==7){
										subrunweight=w_p_2[isubrun];
										subruncfweight=wcf_p_2[isubrun];
									}
									else if(ispecies==8){
										subrunweight=w_pi_2[isubrun];
										subruncfweight=wcf_pi_2[isubrun];
									}
									else if(ispecies==9){
										subrunweight=w_p_2[isubrun];
										subruncfweight=wcf_p_2[isubrun];
									}
								}
								width=norm=0.0;
								fn1="modelruns/run"+to_string(irun)+"/"+qualifier[iqual]+"/results_type"+to_string(itype)+"/subruns/subrun"+to_string(isubrun)+"/"+bfspecies[ispecies]+"/"+bffilenames[ibf];
								fptr1=fopen(fn1.c_str(),"r");
								ix=0;
								while(!feof(fptr1)){
									if(bffilenames[ibf]!="bf_qinv.txt"){
										fscanf(fptr1,"%lf %lf %lf",&dx,&bf,&cf);
										bf*=subrunweight;
										cf*=subrunweight;
										if(!feof(fptr1)){
											width+=fabs(dx)*bf;
											norm+=bf;
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
								if(bffilenames[ibf]!="bf_qinv.txt"){
									width=width/norm;
								
									normbar+=fabs(norm);
									sigma2norm+=norm*norm;
									widthbar+=fabs(width);
									sigma2width+=width*width;
								}
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
							if(bffilenames[ibf]!="bf_qinv.txt"){
								normbar=normbar/double(Nsubruns);
								widthbar=widthbar/double(Nsubruns);
								sigma2norm=sigma2norm/double(Nsubruns)-normbar*normbar;
								sigma2width=sigma2width/double(Nsubruns)-widthbar*widthbar;
								printf("%s,%s: normbar=%g, widthbar=%g, norm unc percentage=%g, width unc percentage=%g\n",
								bfspecies[ispecies].c_str(),bffilenames[ibf].c_str(),
								normbar,widthbar,
								100.0*sqrt(sigma2norm/double(Nsubruns))/normbar,100.0*sqrt(sigma2width/double(Nsubruns))/widthbar);
							}
						}
					}
				}
			}
	
			// Add type 1 to type 2 BFs
			if(dotypes!=1 && dotypes!=2){
				for(ispecies=0;ispecies<nspecies;ispecies++){
					for(ibf=0;ibf<nbfs;ibf++){
						fn1="modelruns/run"+to_string(irun)+"/"+qualifier[iqual]+"/results_type1/"+bfspecies[ispecies]+"/"+bffilenames[ibf];
						fn2="modelruns/run"+to_string(irun)+"/"+qualifier[iqual]+"/results_type2/"+bfspecies[ispecies]+"/"+bffilenames[ibf];
						fnsum="modelruns/run"+to_string(irun)+"/"+qualifier[iqual]+"/results/"+bfspecies[ispecies]+"/"+bffilenames[ibf];
						exists1=filesystem::exists(fn1);
						exists2=filesystem::exists(fn2);
						if(exists1 && exists2){
							dirname="modelruns/run"+to_string(irun)+"/"+qualifier[iqual]+"/results/"+bfspecies[ispecies]+"/";
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
	}
	return 0;
}


