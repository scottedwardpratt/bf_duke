#include "msu_commonutils/log.h"
#include "msu_eos/eos.h"
#include "msu_eos/resonances.h"
#include "bfduke/bfeos.h"
#include "msu_commonutils/parametermap.h"
#include "msu_commonutils/constants.h"
#include "msu_commonutils/commondefs.h"

using namespace std;
using namespace NMSUPratt;

vector<double> CHBEoS::epsilon_PST,CHBEoS::P_PST,CHBEoS::s_PST,CHBEoS::T_PST;
vector<double> CHBEoS::epsilon_claudia,CHBEoS::P_claudia,CHBEoS::s_claudia,CHBEoS::T_claudia;
vector<double> CHBEoS::twopiTD,CHBEoS::Tdiff;
vector<double> CHBEoS::chill_overs_claudia,CHBEoS::chiud_overs_claudia,CHBEoS::chils_overs_claudia,CHBEoS::chiss_overs_claudia;
vector<double> CHBEoS::chill_HSC,CHBEoS::chiud_HSC,CHBEoS::chils_HSC,CHBEoS::chiss_HSC;
vector<double> CHBEoS::dDdT;
double CHBEoS::depsilon=0.005;
mapdi CHBEoS::etmap;

CHBEoS::CHBEoS(){	
}

CHBEoS::CHBEoS(CparameterMap *parmapset){
	parmap=parmapset;
	ReadDiffusionData();
	ReadEosData_Andrew();
	//FillOutdDdT();
};

void CHBEoS::ReadDiffusionData_Andrew(){
	string filename=dirname+"eosdata/DvsEpsilon.txt";
	char dummy[200];
	double epsilon,t,d;
	FILE *fptr=fopen(filename.c_str(),"r");	
	fgets(dummy,100,fptr);
	while(!feof(fptr)){
		fscanf(fptr,"%lf %d %lf %lf %lf",&epsilon,&t,&d);
		if(!feof(fptr))
			DvsE.push_back(d);
	}	
	fclose(fptr);
	
	//for(unsigned int ie=0;ie<DvsE.size();ie++)
	//	printf("%8.5f %10.5f\n",0.005*ie,DvsE[ie]);
	
	
}X

double CHBEoS::GetDfromE(double epsilon){
	int ie0,ie1;
	double D0,D1,D;
	ie0=floorl(epsilon/depsilon);
	ie1=ie0+ie;
	if(ie1<DvsE.size(){
		D0=DvsE[ie0]; D1=DvsE[ie1];
		D=(D0*(ie1*depsilon-epsilon)+D1*(epsilon-ie0*depsilon))/depsilon;
	}
	return D;
}

void CHBEoS::ReadEoSData_Andrew(){
	double e,t,s,cuu,cud,cus,css;
	int NE=0;
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
			chill.push_back(cuu);
			chiud.push_back(cud);
			chils.push_back(cus);
			chiss.push_back(css);
		}
	}while !feof(fptr);
	fclose(fptr);
	NE=TvsS.size();
}

void CHBEoS::ReadChiReductionFactors(){
	double e,t,cfuu,cfud,cfus,cfss;
	int Nchifactors;
	char dummy[100];
	FILE *fptr;
	fptr=fopen("eosdata/eos_chifactors.txt","r");
	chifactorll.clear();
	chifactorud.clear();
	chifactorls.clear();
	chifactorss.clear();
	fgets(dummy,100,fptr);
	do{
		fscanf(fptr,"%lf %lf %lf",&e,&t,&cfuu,&cfud,&cfus,&cfss);
		if(!feof(fptr)){
			chifactorll.push_back(cuu);
			chifactorud.push_back(cud);
			chifactorls.push_back(cus);
			chifactorss.push_back(css);
		}
	}while !feof(fptr);
	fclose(fptr);
	Nchifactors=chifactorll.size();
}

void qCHBEoS::GetChiOverS(double f_u,double f_d,double f_s){
	int ie0,ie1,if0,if1;
	double w,f_q,delf,crll,crud,crls,crss,;
	ie0=floorl(epsilon/depsilon);
	ie1=ie0+1;
	if(ie1>=NE){
		ie1=NE-1;
		ie0=ie1-1;
	}
	f=(ie1*depsilon-epsilon)/depsilon;
	chill=w*chillvsE(ie0)+(1.0-w)*chillvsE(ie1);
	chiud=w*chiudvsE(ie0)+(1.0-w)*chiudvsE(ie1);
	chils=w*chilsvsE(ie0)+(1.0-w)*chilsvsE(ie1);
	chiss=w*chissvsE(ie0)+(1.0-w)*chissvsE(ie1);
	
	
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
