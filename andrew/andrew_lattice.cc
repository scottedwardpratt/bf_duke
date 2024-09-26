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
	int iT0,iT1,iT;
	char dummy[200];
	int iepsilon,NT;
	double dE,dele=0.005,T0,T1,chi0,chi1;
	const double HBARC3=pow(HBARC,3);
	double Tread,chiread,error,epsilon_target;
	vector<double> Tlattice,chilattice_uu,chilattice_ud,chilattice_us,chilattice_ss;
	vector<double> Tvec,PvsTvec,EpsilonvsTvec,SvsTvec,ChivsTvec;
	FILE *fptr;
	const int NE=1001;
	Tlattice.resize(NE);
	chilattice_uu.resize(NE);
	chilattice_ud.resize(NE);
	chilattice_us.resize(NE);
	chilattice_ss.resize(NE);

	ChivsTvec.clear(); EpsilonvsTvec.clear();
	GetLatticeEoS(PvsTvec,SvsTvec,EpsilonvsTvec);
	NT=PvsTvec.size();
	
	for(iT=0;iT<=NT;iT++){
		SvsTvec[iT]=SvsTvec[iT]/HBARC3;
	}
	
	// Now read in lattice data for chi_uu
	printf("----- chi_uu --------\n");
	fptr=fopen("latticedata_claudia/chi_uu.dat","r");
	Tvec.clear();
	do{
		fscanf(fptr,"%lf %lf %lf",&Tread,&chiread,&error);
		fgets(dummy,200,fptr);
		if(!feof(fptr) && Tread<=400.00001){
			Tvec.push_back(Tread/1000.0);
			chiread=chiread*pow(Tread/HBARC,3);
			ChivsTvec.push_back(chiread);
		}
	}while(!feof(fptr) && Tread<400.001);
	fclose(fptr);
	
	for(iepsilon=0;iepsilon<NE;iepsilon++){
		epsilon_target=dele*iepsilon;
		iT0=-1;
		do{
			iT0+=1;
			iT1=iT0+1;
		}while(EpsilonvsTvec[iT1]<epsilon_target);
		T0=Tvec[iT0];
		T1=Tvec[iT0];
		chi0=ChivsTvec[iT0];
		chi1=ChivsTvec[iT1];
		dE=EpsilonvsTvec[iT1]-EpsilonvsTvec[iT0];
		Tlattice[iepsilon]=(T0*(EpsilonvsTvec[iT1]-epsilon_target)+T1*(epsilon_target-EpsilonvsTvec[iT0]))/dE;
		chilattice_uu[iepsilon]=(chi0*(EpsilonvsTvec[iT1]-epsilon_target)+chi1*(epsilon_target-EpsilonvsTvec[iT0]))/dE;
	}
	
///////////////////////////////////////////////
		// Now read in lattice data for chi_ud
	
	printf("----- chi_ud --------\n");
	fptr=fopen("latticedata_claudia/c2ud.cont","r");
	Tvec.clear();
	do{
		fscanf(fptr,"%lf %lf %lf",&Tread,&chiread,&error);
		fgets(dummy,200,fptr);
		if(!feof(fptr) && Tread<=400.00001){
			Tvec.push_back(Tread/1000.0);
			chiread=chiread*pow(Tread/HBARC,3);
			ChivsTvec.push_back(chiread);
		}
	}while(!feof(fptr) && Tread<400.001);
	fclose(fptr);
	
	
	for(iepsilon=0;iepsilon<NE;iepsilon++){
		epsilon_target=dele*iepsilon;
		iT0=-1;
		do{
			iT0+=1;
			iT1=iT0+1;
		}while(EpsilonvsTvec[iT1]<epsilon_target);
		T0=Tvec[iT0];
		T1=Tvec[iT1];
		chi0=ChivsTvec[iT0];
		chi1=ChivsTvec[iT1];
		dE=EpsilonvsTvec[iT1]-EpsilonvsTvec[iT0];
		Tlattice[iepsilon]=(T0*(EpsilonvsTvec[iT1]-epsilon_target)+T1*(epsilon_target-EpsilonvsTvec[iT0]))/dE;
		chilattice_ud[iepsilon]=(chi0*(EpsilonvsTvec[iT1]-epsilon_target)+chi1*(epsilon_target-EpsilonvsTvec[iT0]))/dE;
	}
	
	
///////////////////////////////////////////////
		// Now read in lattice data for chi_us
	
	printf("----- chi_us --------\n");
	fptr=fopen("latticedata_claudia/c2us.cont","r");
	Tvec.clear();
	do{
		fscanf(fptr,"%lf %lf %lf",&Tread,&chiread,&error);
		fgets(dummy,200,fptr);
		if(!feof(fptr) && Tread<=400.00001){
			Tvec.push_back(Tread/1000.0);
			chiread=chiread*pow(Tread/HBARC,3);
			ChivsTvec.push_back(chiread);
		}
	}while(!feof(fptr) && Tread<400.001);
	fclose(fptr);
	
	for(iepsilon=0;iepsilon<NE;iepsilon++){
		epsilon_target=dele*iepsilon;
		iT0=-1;
		do{
			iT0+=1;
			iT1=iT0+1;
		}while(EpsilonvsTvec[iT1]<epsilon_target);
		T0=Tvec[iT0];
		T1=Tvec[iT1];
		chi0=ChivsTvec[iT0];
		chi1=ChivsTvec[iT1];
		dE=EpsilonvsTvec[iT1]-EpsilonvsTvec[iT0];
		Tlattice[iepsilon]=(T0*(EpsilonvsTvec[iT1]-epsilon_target)+T1*(epsilon_target-EpsilonvsTvec[iT0]))/dE;
		chilattice_us[iepsilon]=(chi0*(EpsilonvsTvec[iT1]-epsilon_target)+chi1*(epsilon_target-EpsilonvsTvec[iT0]))/dE;
	}
	
///////////////////////////////////////////////
		// Now read in lattice data for chi_ss
	
	printf("----- chi_ud --------\n");
	fptr=fopen("latticedata_claudia/c2S.cont","r");
	Tvec.clear();
	do{
		fscanf(fptr,"%lf %lf %lf",&Tread,&chiread,&error);
		fgets(dummy,200,fptr);
		if(!feof(fptr) && Tread<=400.00001){
			Tvec.push_back(Tread/1000.0);
			chiread=chiread*pow(Tread/HBARC,3);
			ChivsTvec.push_back(chiread);
		}
	}while(!feof(fptr) && Tread<400.001);
	fclose(fptr);
	
	for(iepsilon=0;iepsilon<NE;iepsilon++){
		epsilon_target=dele*iepsilon;
		iT0=-1;
		do{
			iT0+=1;
			iT1=iT0+1;
		}while(EpsilonvsTvec[iT1]<epsilon_target);
		T0=Tvec[iT0];
		T1=Tvec[iT1];
		chi0=ChivsTvec[iT0];
		chi1=ChivsTvec[iT1];
		dE=EpsilonvsTvec[iT1]-EpsilonvsTvec[iT0];
		Tlattice[iepsilon]=(T0*(EpsilonvsTvec[iT1]-epsilon_target)+T1*(epsilon_target-EpsilonvsTvec[iT0]))/dE;
		chilattice_ss[iepsilon]=(chi0*(EpsilonvsTvec[iT1]-epsilon_target)+chi1*(epsilon_target-EpsilonvsTvec[iT0]))/dE;
	}
	
	return 0;
}

void GetLatticeEoS(vector<double> &Parray,vector<double> &sarray,vector<double> &epsilonarray){
	const double h0=0.1396,h1=-0.1800,h2=0.0350,f0=2.76,f1=6.79,f2=-5.29,g1=-0.47,g2=1.04;
	double dT=1.0,T,poverT4,t,I,cs2,hbarc3=pow(HBARC,-3);
	const int NT=401;
	sarray.resize(NT);
	epsilonarray.resize(NT);
	Parray.resize(NT);
	FILE *fptr;
	int iT;
	poverT4=0.0;
	Parray[0]=epsilonarray[0]=sarray[0]=0.0;
	fptr=fopen("eos_lattice.dat","w");
	fprintf(fptr,"# T P(GeV/fm^3) epsilon  s(fm^-3) c_s^2\n");
	for(iT=0;iT<NT;iT++){
		T=(iT+0.5)*dT;
		t=T/200.0;
		I=pow(T,4.0)*exp(-(h1/t)-h2/(t*t))*( h0+f0*(tanh(f1*t+f2)+1.0)/(1.0+g1*t+g2*t*t) );
		poverT4+=(dT/T)*I/pow(T,4.0);
		
		T=(iT+1.0)*dT;
		Parray[iT+1]=pow(T,4)*poverT4;
		I=pow(T,4.0)*exp(-(h1/t)-h2/(t*t))*( h0+f0*(tanh(f1*t+f2)+1.0)/(1.0+g1*t+g2*t*t) );
		epsilonarray[iT+1]=I+3.0*Parray[iT+1];
		sarray[iT+1]=(Parray[iT+1]+epsilonarray[iT+1])/T;
		//printf("%5.1f s=%g\n",T,sarray[iT]*hbarc3);
		if(iT>0){
			cs2=(Parray[iT+1]-Parray[iT-1])/(epsilonarray[iT+1]-epsilonarray[iT-1]);
			fprintf(fptr,"%5.1f %8.4f %8.4f %8.4f %6.4f\n",iT*dT,Parray[iT]*hbarc3/1000.0,epsilonarray[iT]*hbarc3/1000.0,sarray[iT]*hbarc3,cs2);
		}
	}
	for(iT=0;iT<NT;iT++){
		T=(iT+0.5)*dT;
		//printf("%6.1f %10.3e %10.3e %10.3e\n",T,
		//epsilonarray[iT]/(pow(T,4)),Parray[iT]/(pow(T,4)),sarray[iT]/(pow(T,3)));
		Parray[iT]=Parray[iT]*hbarc3/1000.0;
		epsilonarray[iT]=epsilonarray[iT]*hbarc3/1000.0;
	}
	fclose(fptr);

}
