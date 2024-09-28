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
	int iT0,iT1,iepsilon;
	char dummy[300];
	double dE,dele,T0,T1,chi0,chi1,dT=1.0,epsilon_target;
	double Tread,chiread,error;
	vector<double> Tlattice,chilattice_uu,chilattice_ud,chilattice_us,chilattice_ss,slattice,elattice;
	vector<double> TvsTvec,Tvec,PvsTvec,EpsilonvsTvec,SvsTvec,ChivsTvec;
	FILE *fptr;
	const int NE=8169;
	dele=0.005;
	Tlattice.resize(NE);
	chilattice_uu.resize(NE);
	chilattice_ud.resize(NE);
	chilattice_us.resize(NE);
	chilattice_ss.resize(NE);
	slattice.resize(NE);
	elattice.resize(NE);
	

	ChivsTvec.clear(); EpsilonvsTvec.clear();
	GetLatticeEoS(PvsTvec,SvsTvec,EpsilonvsTvec);
	
	for(iepsilon=0;iepsilon<NE;iepsilon++){
		elattice[iepsilon]=iepsilon*dele;
		T0=-dT;
		iT0=-1;
		do{
			iT0+=1;
			iT1=iT0+1;
			T0+=dT;
			T1=T0+dT;
		}while(EpsilonvsTvec[iT1]<elattice[iepsilon]);
		dE=EpsilonvsTvec[iT1]-EpsilonvsTvec[iT0];
		Tlattice[iepsilon]=(EpsilonvsTvec[iT1]-elattice[iepsilon])*T0/dE+(elattice[iepsilon]-EpsilonvsTvec[iT0])*T1/dE;
		Tlattice[iepsilon]=Tlattice[iepsilon]/1000.0;
		slattice[iepsilon]=(EpsilonvsTvec[iT1]-elattice[iepsilon])*SvsTvec[iT0]/dE+(elattice[iepsilon]-EpsilonvsTvec[iT0])*SvsTvec[iT1]/dE;
	}
	Tlattice[0]=0.0;

///////////////////////////////////////////////	
	// Now read in lattice data for chi_uu
	printf("----- chi_uu --------\n");
	fptr=fopen("latticedata_claudia/chi_uu.dat","r");
	Tvec.clear();
	do{
		fscanf(fptr,"%lf %lf %lf",&Tread,&chiread,&error);
		fgets(dummy,300,fptr);
		if(!feof(fptr) && Tread<=400.00001){
			Tvec.push_back(Tread/1000.0);
			chiread=chiread*pow(Tread/HBARC,3);
			ChivsTvec.push_back(chiread);
		}
	}while(!feof(fptr) && Tread<400.001);
	fclose(fptr);
	
	for(iepsilon=0;iepsilon<NE;iepsilon++){
		//printf("Ttarget=%g, e=%g\n",Tlattice[iepsilon],iepsilon*dele);
		iT0=-1;
		do{
			iT0+=1;
			iT1=iT0+1;
		}while(Tvec[iT1]<Tlattice[iepsilon]);
		T0=Tvec[iT0];
		T1=Tvec[iT1];
		chi0=ChivsTvec[iT0];
		chi1=ChivsTvec[iT1];
		dT=Tvec[iT1]-Tvec[iT0];
		chilattice_uu[iepsilon]=(chi0*(Tvec[iT1]-Tlattice[iepsilon])+chi1*(Tlattice[iepsilon]-Tvec[iT0]))/dT;
		//printf("epsilon=%g: iT0=%d, iT1=%d, T0=%g, T1=%g, Tlattice=%g, chi0=%g, chi1=%g, chi=%g\n",
		//elattice[iepsilon],iT0,iT1,T0,T1,Tlattice[iepsilon],chi0,chi1,chilattice_uu[iepsilon]);
	}

	
///////////////////////////////////////////////
		// Now read in lattice data for chi_ud
	
	printf("----- chi_ud --------\n");
	fptr=fopen("latticedata_claudia/c2ud.cont","r");
	Tvec.clear(); ChivsTvec.clear();
	do{
		fscanf(fptr,"%lf %lf %lf",&Tread,&chiread,&error);
		fgets(dummy,300,fptr);
		if(!feof(fptr) && Tread<=400.00001){
			Tvec.push_back(Tread/1000.0);
			chiread=chiread*pow(Tread/HBARC,3);
			ChivsTvec.push_back(chiread);
		}
	}while(!feof(fptr) && Tread<400.001);
	fclose(fptr);
	
	for(iepsilon=0;iepsilon<NE;iepsilon++){
		//printf("Ttarget=%g, e=%g\n",Tlattice[iepsilon],iepsilon*dele);
		iT0=-1;
		do{
			iT0+=1;
			iT1=iT0+1;
		}while(Tvec[iT1]<Tlattice[iepsilon]);
		T0=Tvec[iT0];
		T1=Tvec[iT1];
		chi0=ChivsTvec[iT0];
		chi1=ChivsTvec[iT1];
		dT=Tvec[iT1]-Tvec[iT0];
		chilattice_ud[iepsilon]=(chi0*(Tvec[iT1]-Tlattice[iepsilon])+chi1*(Tlattice[iepsilon]-Tvec[iT0]))/dT;
		//printf("epsilon=%g: iT0=%d, iT1=%d, T0=%g, T1=%g, Tlattice=%g, chi0=%g, chi1=%g, chi=%g\n",
		//elattice[iepsilon],iT0,iT1,T0,T1,Tlattice[iepsilon],chi0,chi1,chilattice_ud[iepsilon]);
	}
	
///////////////////////////////////////////////
		// Now read in lattice data for chi_us
	
	printf("----- chi_us --------\n");
	fptr=fopen("latticedata_claudia/c2us.cont","r");
	Tvec.clear(); ChivsTvec.clear();
	do{
		fscanf(fptr,"%lf %lf %lf",&Tread,&chiread,&error);
		fgets(dummy,300,fptr);
		if(!feof(fptr) && Tread<=400.00001){
			Tvec.push_back(Tread/1000.0);
			chiread=chiread*pow(Tread/HBARC,3);
			ChivsTvec.push_back(chiread);
		}
	}while(!feof(fptr) && Tread<400.001);
	fclose(fptr);
	
	for(iepsilon=0;iepsilon<NE;iepsilon++){
		//printf("Ttarget=%g, e=%g\n",Tlattice[iepsilon],iepsilon*dele);
		iT0=-1;
		do{
			iT0+=1;
			iT1=iT0+1;
		}while(Tvec[iT1]<Tlattice[iepsilon]);
		T0=Tvec[iT0];
		T1=Tvec[iT1];
		chi0=ChivsTvec[iT0];
		chi1=ChivsTvec[iT1];
		dT=Tvec[iT1]-Tvec[iT0];
		chilattice_us[iepsilon]=(chi0*(Tvec[iT1]-Tlattice[iepsilon])+chi1*(Tlattice[iepsilon]-Tvec[iT0]))/dT;
		//printf("epsilon=%g: iT0=%d, iT1=%d, T0=%g, T1=%g, Tlattice=%g, chi0=%g, chi1=%g, chi=%g\n",
		//elattice[iepsilon],iT0,iT1,T0,T1,Tlattice[iepsilon],chi0,chi1,chilattice_us[iepsilon]);
	}

	
///////////////////////////////////////////////
		// Now read in lattice data for chi_ss

	printf("----- chi_ss --------\n");
	fptr=fopen("latticedata_claudia/c2s.cont","r");
	Tvec.clear(); ChivsTvec.clear();
	do{
		fscanf(fptr,"%lf %lf %lf",&Tread,&chiread,&error);
		fgets(dummy,300,fptr);
		if(!feof(fptr) && Tread<=400.00001){
			Tvec.push_back(Tread/1000.0);
			chiread=chiread*pow(Tread/HBARC,3);
			ChivsTvec.push_back(chiread);
		}
	}while(!feof(fptr) && Tread<400.001);
	fclose(fptr);
	
	for(iepsilon=0;iepsilon<NE;iepsilon++){
		//printf("Ttarget=%g, e=%g\n",Tlattice[iepsilon],iepsilon*dele);
		iT0=-1;
		do{
			iT0+=1;
			iT1=iT0+1;
		}while(Tvec[iT1]<Tlattice[iepsilon]);
		T0=Tvec[iT0];
		T1=Tvec[iT1];
		chi0=ChivsTvec[iT0];
		chi1=ChivsTvec[iT1];
		dT=Tvec[iT1]-Tvec[iT0];
		chilattice_ss[iepsilon]=(chi0*(Tvec[iT1]-Tlattice[iepsilon])+chi1*(Tlattice[iepsilon]-Tvec[iT0]))/dT;
		//printf("epsilon=%g: iT0=%d, iT1=%d, T0=%g, T1=%g, Tlattice=%g, chi0=%g, chi1=%g, chi=%g\n",
		//elattice[iepsilon],iT0,iT1,T0,T1,Tlattice[iepsilon],chi0,chi1,chilattice_ss[iepsilon]);
	}	
	

///////////////////////////////////////////////
		// Write lattice vs epsilon
	fptr=fopen("LatticeVsEpsilon.txt","w");
	fprintf(fptr,"#epsilon     T           chi_uu          chi_ud             chi_us             chi_ss\n");
	for(iepsilon=0;iepsilon<NE;iepsilon++){
		epsilon_target=dele*iepsilon;
		fprintf(fptr,"%15.8e %15.8e %15.8e %15.8e %15.8e %15.8e\n",epsilon_target,Tlattice[iepsilon],
		chilattice_uu[iepsilon],chilattice_ud[iepsilon],chilattice_us[iepsilon],chilattice_ss[iepsilon]);
		printf("%8.4f %6.4f %8.5f %8.5f %8.5f %8.5f\n",epsilon_target,Tlattice[iepsilon],
		chilattice_uu[iepsilon]/slattice[iepsilon],chilattice_ud[iepsilon]/slattice[iepsilon],
		chilattice_us[iepsilon]/slattice[iepsilon],chilattice_ss[iepsilon]/slattice[iepsilon]);
	}
	fclose(fptr);

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
		sarray[iT]=(Parray[iT]+epsilonarray[iT])*1000.0/T;
	}
	fclose(fptr);

}
