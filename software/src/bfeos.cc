#include "msu_commonutils/log.h"
#include "msu_eos/eos.h"
#include "msu_eos/resonances.h"
#include "bfduke/bfeos.h"
#include "msu_commonutils/parametermap.h"
#include "msu_commonutils/constants.h"
#include "msu_commonutils/commondefs.h"

using namespace std;
using namespace NMSUPratt;


CparameterMap *CHBEoS::parmap=NULL;
int CHBEoS::NE=0;
double CHBEoS::depsilon=0.0;
double CHBEoS::epsilon_h=0.36;
double CHBEoS::epsilon_qgp=0.76;
vector<double> CHBEoS::TvsE;
vector<double> CHBEoS::svsE;
vector<double> CHBEoS::chillvsE,CHBEoS::chiudvsE,CHBEoS::chilsvsE,CHBEoS::chissvsE;
vector<double> CHBEoS::DvsE;

int CHBEoS::Nchifactors=0;
vector<double> CHBEoS::chifactorll,CHBEoS::chifactorud,CHBEoS::chifactorls,CHBEoS::chifactorss,CHBEoS::TnonequilVec;

CHBEoS::CHBEoS(){	
}

void CHBEoS::ReadDiffusionData_Andrew(){
	string filename="eosdata/DvsEpsilon.txt";
	char dummy[200];
	double epsilon,t,d;
	FILE *fptr=fopen(filename.c_str(),"r");	
	fgets(dummy,100,fptr);
	while(!feof(fptr)){
		fscanf(fptr,"%lf %lf %lf",&epsilon,&t,&d);
		if(!feof(fptr))
			DvsE.push_back(d);
	}	
	fclose(fptr);
	
}

void CHBEoS::ReadEoSData_Andrew(){
	double e,t,s,cuu,cud,cus,css;
	char dummy[100];
	FILE *fptr;
	depsilon=0.005;
	fptr=fopen("eosdata/eos_vs_epsilon.txt","r");

	TvsE.clear();
	chillvsE.clear();
	chiudvsE.clear();
	chilsvsE.clear();
	chissvsE.clear();
	fgets(dummy,100,fptr);
	do{
		fscanf(fptr,"%lf %lf %lf %lf %lf %lf %lf",&e,&t,&s,&cuu,&cud,&cus,&css);
		if(!feof(fptr)){
			TvsE.push_back(t);
			svsE.push_back(s);
			chillvsE.push_back(cuu);
			chiudvsE.push_back(cud);
			chilsvsE.push_back(cus);
			chissvsE.push_back(css);
		}
	}while(!feof(fptr));
	fclose(fptr);
	NE=TvsE.size();
}

void CHBEoS::ReadChiReductionFactors(){
	double fq,t,cfuu,cfud,cfus,cfss;
	double chillvsE,chiudvsE,chilsvsE,chissvsE,sdens;
	char dummy[100];
	FILE *fptr;
	fptr=fopen("eosdata/eos_chifactors.txt","r");
	chifactorll.clear();
	chifactorud.clear();
	chifactorls.clear();
	chifactorss.clear();
	fgets(dummy,100,fptr);
	do{
		fscanf(fptr,"%lf %lf %lf %lf %lf %lf %lf  %lf %lf %lf %lf",&fq,&t,&sdens,
		&chillvsE,&chiudvsE,&chilsvsE,&chissvsE,&cfuu,&cfud,&cfus,&cfss);
		if(!feof(fptr)){
			chifactorll.push_back(cfuu);
			chifactorud.push_back(cfud);
			chifactorls.push_back(cfus);
			chifactorss.push_back(cfss);
			TnonequilVec.push_back(t);
		}
		fgets(dummy,100,fptr);
	}while(!feof(fptr));
	fclose(fptr);
	Nchifactors=chifactorll.size();
}

void CHBEoS::SetChi(){
	int ie0,ie1,if0,if1;
	double w,f_q,delf,crll,crud,crls,crss;
	ie0=floorl(epsilon/depsilon);
	ie1=ie0+1;
	if(ie1>=NE){
		ie1=NE-1;
		ie0=ie1-1;
	}
	w=(ie1*depsilon-epsilon)/depsilon;
	chill=w*chillvsE[ie0]+(1.0-w)*chillvsE[ie1];
	chiud=w*chiudvsE[ie0]+(1.0-w)*chiudvsE[ie1];
	chils=w*chilsvsE[ie0]+(1.0-w)*chilsvsE[ie1];
	chiss=w*chissvsE[ie0]+(1.0-w)*chissvsE[ie1];
	
	
	f_q=(f_u+f_d+f_s)/3.0;
	delf=1.0/Nchifactors;
	if0=floorl(f_q/delf);
	if1=if0+1;
	w=(if1*delf-f_q)/delf;
	crll=w*chifactorll[if0]+(1.0-w)*chifactorll[if1];
	crud=w*chifactorud[if0]+(1.0-w)*chifactorud[if1];
	crls=w*chifactorls[if0]+(1.0-w)*chifactorls[if1];
	crss=w*chifactorss[if0]+(1.0-w)*chifactorss[if1];
	
	if(epsilon<epsilon_h){
		chill*=crll;
		chiud*=crud;
		chils*=crls;
		chiss*=crss;
	}
	else if(epsilon<epsilon_qgp){// above epsilon_qgp reduce off-diag chi by hadron value
		// reduce diagonal by weighted hadron vs qgp value, qgp reducution is f_q
		w=(epsilon_qgp-epsilon)/(epsilon_qgp-epsilon_h);
		
		crll=w*crll+(1.0-w)*f_q;
		crss=w*crss+(1.0-w)*f_q;
		chill*=crll;
		chiud*=crud;
		chils*=crls;
		chiss*=crss;
	}
	else{ 
		chill*=f_q;
		chiss*=f_q;
		chiud*=crud;
		chils*=crls;
	}
}

void CHBEoS::SetTnonequil(){
	int if0,if1;
	double w,f_q,delf;

	f_q=(f_u+f_d+f_s)/3.0;
	delf=1.0/Nchifactors;
	if0=floorl(f_q/delf);
	if1=if0+1;
	w=(if1*delf-f_q)/delf;
	Tnonequil=w*TnonequilVec[if0]+(1.0-w)*TnonequilVec[if1];
}

void CHBEoS::SetTs(){  // T,P and s are equilibrium quantities
	int ie0,ie1;
	double w;
	ie0=floorl(epsilon/depsilon);
	ie1=ie0+1;
	if(ie1>=NE){
		ie1=NE-1;
		ie0=ie1-1;
	}
	w=(ie1*depsilon-epsilon)/depsilon;
	T=w*TvsE[ie0]+(1.0-w)*TvsE[ie1];
	s=w*svsE[ie0]+(1.0-w)*svsE[ie1];
	P=T*s-epsilon;
}



double CHBEoS::GetD(double epsilon){
	unsigned int ie0,ie1;
	double D0,D1,D;
	ie0=floorl(epsilon/depsilon);
	ie1=ie0+1;
	if(ie1<DvsE.size()){
		D0=DvsE[ie0]; D1=DvsE[ie1];
		D=(D0*(ie1*depsilon-epsilon)+D1*(epsilon-ie0*depsilon))/depsilon;
	}
	else{
		D=DvsE[DvsE.size()-1];
	}
	return D;
}

