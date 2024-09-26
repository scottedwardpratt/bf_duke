#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>
#include <vector>
#include <Eigen/Dense>
#include "msu_eos/resonances.h"
#include "msu_eos/eos.h"
#include "msu_commonutils/randy.h"
#include "msu_commonutils/constants.h"

using namespace std;
using namespace NMSUPratt;

void GetLatticeEoS(vector<double> &Parray,vector<double> &sarray,vector<double> &epsilonarray);

int main(){
	CparameterMap parmap;
	parmap.set("MSU_SAMPLER_SFDIRNAME",string("../progdata/resinfo/spectralfunctions"));
	parmap.set("RESONANCES_INFO_FILE",string("../progdata/resinfo/pdg-SMASH.dat"));
	int iT,NT=50;
	double delT=1.0,fq;
	bool use_pole_mass=true;
	double Th,Ph,epsilonh,sh,nhadronsh;
	Eigen::Matrix3d sigmah,chih;
	vector<double> densityh;
	CresList *reslist=new CresList(&parmap);
	densityh.resize(reslist->resmap.size());
	delT=0.005;
	for(iT=1;iT<NT;iT++){
		Th=iT*delT;
		fq=1.0;
		MSU_EOS::CalcEoSandTransportCoefficients(Th,reslist,epsilonh,Ph,
		nhadronsh,densityh,chih,sigmah,use_pole_mass,fq,fq,fq);
		sh=(Ph+epsilonh)/Th;
		printf("%6.4f: %9.5f %9.5f %9.5f %9.5f\n",Th,epsilonh,Ph,sh,nhadronsh);
		//cout << chiharray[iT] << endl;
	}
	
	double epsilon_target,dele=0.005,chiuu,chiud,chius,chiss;
	int iepsilon,Nepsilon=1000;
	fq=1.0;
	FILE *fptr=fopen("eosdata/eos_hadron_vs_epsilon.txt","w");
	fprintf(fptr,"# fugacities set to unity\n");
	for(iepsilon=0;iepsilon<=Nepsilon;iepsilon++){
		epsilon_target=iepsilon*dele;
		if(iepsilon>0.0001){
			MSU_EOS::CalcTFromEpsilonFugacity(epsilon_target,fq,fq,fq,reslist,use_pole_mass,Th);
			MSU_EOS::CalcEoSandTransportCoefficients(Th,reslist,epsilonh,Ph,
			nhadronsh,densityh,chih,sigmah,use_pole_mass,fq,fq,fq);
			sh=(Ph+epsilonh)/Th;
			chiuu=chih(0,0);
			chiud=chih(0,1);
			chius=chih(0,2);
			chiss=chih(2,2);
		}
		else{ // For epsilon<0.05, just scale T=50 EoS quantities (hopefully these are not used
			Th=chiuu=chiud=chius=chiss=nhadronsh=sh=0.000000001;	
		}
		fprintf(fptr,"%15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e\n",
		epsilon_target,Th,nhadronsh,sh,chiuu,chiud,chius,chiss);
	}
	fclose(fptr);
	
	
	
	
	
	vector<double> T_vs_F,chiuu_vs_F,chidd_vs_F,chiud_vs_F,chiss_vs_F,chils_vs_F;
	int iF,NF=100;
	T_vs_F.resize(NF+1);
	chiuu_vs_F.resize(NF+1);
	chidd_vs_F.resize(NF+1);
	chiud_vs_F.resize(NF+1);
	chiss_vs_F.resize(NF+1);
	chils_vs_F.resize(NF+1);


	double delF=1.0/double(NF);
	// First find T at which epsilon(T)=Tf;
	// For each fugacity
	double chifactor_uu,chifactor_ud,chifactor_us,chifactor_ss,chiTc_uu,chiTc_ud,chiTc_us,chiTc_ss;
	epsilon_target=0.36;
	MSU_EOS::CalcTFromEpsilonFugacity(epsilon_target,fq,fq,fq,reslist,use_pole_mass,Th);
	MSU_EOS::CalcEoSandTransportCoefficients(Th,reslist,epsilonh,Ph,
	nhadronsh,densityh,chih,sigmah,use_pole_mass,fq,fq,fq);
	chiTc_uu=chih(0,0);
	chiTc_ud=chih(0,1);
	chiTc_us=chih(0,2);
	chiTc_ss=chih(2,2);
	
	/*
	fptr=fopen("eosdata/eos_hadron_vs_fq.txt\n","w");
	fprintf(fptr,"epsilon=%g\n",epsilon_target);
	for(iF=0;iF<=NF;iF++){
		fq=iF*delF;
		if(iF==0){
			epsilonh=epsilon_target;
			Th=1.0;
			sh=0.5;
			nhadronsh=0.0;
			chiuu=chiud=chius=chiss=0.0;
		}
		else{
			MSU_EOS::CalcTFromEpsilonFugacity(epsilon_target,fq,fq,fq,reslist,use_pole_mass,Th);
			MSU_EOS::CalcEoSandTransportCoefficients(Th,reslist,epsilonh,Ph,
			nhadronsh,densityh,chih,sigmah,use_pole_mass,fq,fq,fq);
			sh=(Ph+epsilonh)/Th;
			chiuu=chih(0,0);
			chiud=chih(0,1);
			chius=chih(0,2);
			chiss=chih(2,2);
		}
		chifactor_uu=chiuu/chiTc_uu;
		chifactor_ud=chiud/chiTc_ud;
		chifactor_us=chius/chiTc_us;
		chifactor_ss=chiss/chiTc_ss;
		//printf("fq=%6.3f sh=%g T=%g, chifactors=(%8.4f %8.4f %8.4f %8.4f\n",fq,sh,Th,chifactor_uu,chifactor_ud,chifactor_us,chifactor_ss);
		fprintf(fptr,"%15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e\n",
		fq,Th,nhadronsh,sh,chiuu,chiud,chius,chiss);		
	}
	fclose(fptr);
	*/
	
	printf("-------------------------------\n");
	NF=100;
	delF=1.0/double(NF);
	epsilon_target=0.36;
	fptr=fopen("eosdata/eos_hadron_chifactors.txt","w");
	do{
		fprintf(fptr,"#fugacity  Th    sh  chifactor_uu  chifactor_ud  chifactor_us chifactor_ss\n");
		MSU_EOS::CalcTFromEpsilonFugacity(epsilon_target,fq,fq,fq,reslist,use_pole_mass,Th);
		MSU_EOS::CalcEoSandTransportCoefficients(Th,reslist,epsilonh,Ph,
		nhadronsh,densityh,chih,sigmah,use_pole_mass,fq,fq,fq);
		chiTc_uu=chih(0,0);
		chiTc_ud=chih(0,1);
		chiTc_us=chih(0,2);
		chiTc_ss=chih(2,2);
		
		for(iF=0;iF<=NF;iF++){
			fq=iF*delF;
			if(iF==0){
				epsilonh=epsilon_target;
				Th=1.0;
				sh=0.5;
				nhadronsh=0.0;
				chiuu=chiud=chius=chiss=0.0;
			}
			else{
				MSU_EOS::CalcTFromEpsilonFugacity(epsilon_target,fq,fq,fq,reslist,use_pole_mass,Th);
				MSU_EOS::CalcEoSandTransportCoefficients(Th,reslist,epsilonh,Ph,
				nhadronsh,densityh,chih,sigmah,use_pole_mass,fq,fq,fq);
				sh=(Ph+epsilonh)/Th;
				chiuu=chih(0,0);
				chiud=chih(0,1);
				chius=chih(0,2);
				chiss=chih(2,2);
			}
			chifactor_uu=chiuu/chiTc_uu;
			chifactor_ud=chiud/chiTc_ud;
			chifactor_us=chius/chiTc_us;
			chifactor_ss=chiss/chiTc_ss;
			fprintf(fptr,"%8.4f %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e\n",
			fq,Th,sh,chiuu,chiud,chius,chiss,chifactor_uu,chifactor_ud,chifactor_us,chifactor_ss);
		}
		printf("Enter epsilon_target: (0 to quit)");
		scanf("%lf",&epsilon_target);
	}while(epsilon_target>0.0001);
	fclose(fptr);

	
}
