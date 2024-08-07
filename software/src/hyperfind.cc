#include "msu_commonutils/log.h"
#include "bfduke/hydro2uds.h"
#include "msu_sampler/hyper.h"
#include "bfduke/bfeos.h"
#include "msu_commonutils/misc.h"

using namespace std;
using namespace NMSUPratt;

void CHydroBalance::HyperFindEpsilon(){
	char message[CLog::CHARLENGTH];
	int ix,iy;//,a,b;
	double dEdx,dEdy,dEdt,fugacity_l,fugacity_s,gamma_q;
	Chyper hyper;
	Chyper *newhyper;
	bool GGEt,GGEx,GGEy;
	if(!tau0check){
		for(ix=1;ix<mesh->NX-1;ix++){
			for(iy=1;iy<mesh->NY-1;iy++){
				if(ix<0 || iy<0 || ix==CHBHydroMesh::NX || iy==CHBHydroMesh::NY){
					snprintf(message,CLog::CHARLENGTH,"CHydroBalance::HyperFindEpsilon() -- ix=%d, iy=%d\n",ix,iy);
					CLog::Info(message);
				}
				GetGradEpsilon(ix,iy,dEdt,dEdx,dEdy,GGEt,GGEx,GGEy);
				if(GGEt || GGEx || GGEy){
					hyper.tau=0.5*(newmesh->tau+mesh->tau);
					GetXYBar(ix,iy,hyper.r[1],hyper.r[2]);
					GetUxyBar(ix,iy,hyper.u[1],hyper.u[2]);
					hyper.u[3]=0.0;
					hyper.u[0]=sqrt(1.0+hyper.u[1]*hyper.u[1]+hyper.u[2]*hyper.u[2]);
					GetPiTildeBar(ix,iy,hyper.pitilde[1][1],hyper.pitilde[1][2],
					hyper.pitilde[2][2]);
					double T0;
					GetTBar(ix,iy,T0);
					hyper.T0=T0;
					if(GetDOmega(dEdt,dEdx,dEdy,
					hyper.dOmega[0],hyper.dOmega[1],hyper.dOmega[2],GGEt,GGEx,GGEy)){
						GetGammaFQ(hyper.tau,gamma_q,fugacity_l,fugacity_s);
						hyper.fugacity_u=hyper.fugacity_d=fugacity_l;
						hyper.fugacity_s=fugacity_s;

						/*
						eos->T=hyper.T0;
						eos->GetChiOverS(gamma_q);
						hyper.chi(0,0)=hyper.chi(1,1)=eos->chill;
						hyper.chi(0,1)=hyper.chi(1,0)=eos->chiud;
						hyper.chi(0,2)=hyper.chi(1,2)=hyper.chi(2,0)=hyper.chi(2,1)=eos->chils;
						hyper.chi(2,2)=eos->chiss;
						hyper.chiinv=hyper.chi.inverse();
						for(a=0;a<3;a++){
							for(b=0;b<3;b++){
								chitothyper(a,b)+=hyper.udotdOmega*hyper.chi(a,b);
							}
						}
						*/
						newhyper=new Chyper;
						newhyper->Copy(&hyper);
						hyperlist.push_back(newhyper);
						netUdotOmega+=hyper.u[0]*hyper.dOmega[0]-hyper.u[1]*hyper.dOmega[1]-hyper.u[2]*hyper.dOmega[2];
						if(netUdotOmega!=netUdotOmega){
							hyper.Print();
							exit(1);
						}
					}
				}
			}
		}
	}
}

void CHydroBalance::HyperFindT(){
	char message[CLog::CHARLENGTH];
	int ix,iy,a,b;
	double dTdx,dTdy,dTdt;
	Chyper hyper;
	Chyper *newhyper;
	bool GGTt,GGTx,GGTy;
	if(!tau0check){
		for(ix=1;ix<mesh->NX-1;ix++){
			for(iy=1;iy<mesh->NY-1;iy++){
				if(ix<0 || iy<0 || ix==CHBHydroMesh::NX || iy==CHBHydroMesh::NY){
					snprintf(message,CLog::CHARLENGTH,"CHydroBalance::HyperFindT() -- ix=%d, iy=%d\n",ix,iy);
					CLog::Info(message);
				}
				GetGradT(ix,iy,dTdt,dTdx,dTdy,GGTt,GGTx,GGTy);
				if(GGTt || GGTx || GGTy){
					hyper.tau=0.5*(newmesh->tau+mesh->tau);
					GetXYBar(ix,iy,hyper.r[1],hyper.r[2]);
					GetUxyBar(ix,iy,hyper.u[1],hyper.u[2]);
					GetPiTildeBar(ix,iy,hyper.pitilde[1][1],hyper.pitilde[1][2],
					hyper.pitilde[2][2]);
					double T0;
					GetTBar(ix,iy,T0);
					hyper.T0=T0;
					if(GetDOmega(dTdt,dTdx,dTdy,
					hyper.dOmega[0],hyper.dOmega[1],hyper.dOmega[2],GGTt,GGTx,GGTy)){
						for(a=0;a<3;a++){
							for(b=0;b<3;b++){
								chitothyper(a,b)+=hyper.udotdOmega*chif(a,b);
							}
						}
						newhyper=new Chyper;
						newhyper->Copy(&hyper);
						hyperlist.push_back(newhyper);
					}
				}
			}
		}
	}
}

bool CHydroBalance::GetGradEpsilon(int ix,int iy,
double &dEdt,double &dEdx,double &dEdy,bool &GGEt,bool &GGEx,bool &GGEy){
	bool hypercheck=false;
	double Explus,Exminus,Eyplus,Eyminus,Etplus,Etminus;
	mesh->GetDimensions(NX,NY,DX,DY,DELTAU,TAU0,XMIN,XMAX,YMIN,YMAX);
	
	Exminus=0.25*(mesh->epsilon[ix][iy]+newmesh->epsilon[ix][iy]
		+mesh->epsilon[ix][iy+1]+newmesh->epsilon[ix][iy+1]);
	Explus=0.25*(mesh->epsilon[ix+1][iy]+newmesh->epsilon[ix+1][iy]
		+mesh->epsilon[ix+1][iy+1]+newmesh->epsilon[ix+1][iy+1]);
	
	Eyminus=0.25*(mesh->epsilon[ix][iy]+newmesh->epsilon[ix][iy]
		+mesh->epsilon[ix+1][iy]+newmesh->epsilon[ix+1][iy]);
	Eyplus=0.25*(mesh->epsilon[ix][iy+1]+newmesh->epsilon[ix][iy+1]
		+mesh->epsilon[ix+1][iy+1]+newmesh->epsilon[ix+1][iy+1]);
	
	Etminus=0.25*(mesh->epsilon[ix][iy]+mesh->epsilon[ix][iy+1]
		+mesh->epsilon[ix+1][iy]+mesh->epsilon[ix+1][iy+1]);
	Etplus=0.25*(newmesh->epsilon[ix][iy]+newmesh->epsilon[ix][iy+1]
		+newmesh->epsilon[ix+1][iy]+newmesh->epsilon[ix+1][iy+1]);
	
	dEdx=-(Explus-Exminus)/DX;
	dEdy=-(Eyplus-Eyminus)/DY;
	dEdt=(Etplus-Etminus)/DELTAU;
	
	GGEt=GGEx=GGEy=false;
	if((Etplus-epsilon_f)*(Etminus-epsilon_f)<0.0)
		GGEt=true;
	if((Explus-epsilon_f)*(Exminus-epsilon_f)<0.0)
		GGEx=true;
	if((Eyplus-epsilon_f)*(Eyminus-epsilon_f)<0.0)
		GGEy=true;
	if(GGEt || GGEx || GGEy)
		hypercheck=true;
	return hypercheck;
}

bool CHydroBalance::GetGradT(int ix,int iy,
double &dTdt,double &dTdx,double &dTdy,bool &GGTt,bool &GGTx,bool &GGTy){
	bool hypercheck=false;
	double Txplus,Txminus,Typlus,Tyminus,Ttplus,Ttminus;
	mesh->GetDimensions(NX,NY,DX,DY,DELTAU,TAU0,XMIN,XMAX,YMIN,YMAX);
	
	Txminus=0.25*(mesh->T[ix][iy]+newmesh->T[ix][iy]
		+mesh->T[ix][iy+1]+newmesh->T[ix][iy+1]);
	Txplus=0.25*(mesh->T[ix+1][iy]+newmesh->T[ix+1][iy]
		+mesh->T[ix+1][iy+1]+newmesh->T[ix+1][iy+1]);
	
	Tyminus=0.25*(mesh->T[ix][iy]+newmesh->T[ix][iy]
		+mesh->T[ix+1][iy]+newmesh->T[ix+1][iy]);
	Typlus=0.25*(mesh->T[ix][iy+1]+newmesh->T[ix][iy+1]
		+mesh->T[ix+1][iy+1]+newmesh->T[ix+1][iy+1]);
	
	Ttminus=0.25*(mesh->T[ix][iy]+mesh->T[ix][iy+1]
		+mesh->T[ix+1][iy]+mesh->T[ix+1][iy+1]);
	Ttplus=0.25*(newmesh->T[ix][iy]+newmesh->T[ix][iy+1]
		+newmesh->T[ix+1][iy]+newmesh->T[ix+1][iy+1]);
	
	dTdx=-(Txplus-Txminus)/DX;
	dTdy=-(Typlus-Tyminus)/DY;
	dTdt=(Ttplus-Ttminus)/DELTAU;
	
	GGTt=GGTx=GGTy=false;
	if((Ttplus-Tf)*(Ttminus-Tf)<0.0)
		GGTt=true;
	if((Txplus-Tf)*(Txminus-Tf)<0.0)
		GGTx=true;
	if((Typlus-Tf)*(Tyminus-Tf)<0.0)
		GGTy=true;
	if(GGTt || GGTx || GGTy)
		hypercheck=true;
	return hypercheck;
}

bool CHydroBalance::GetDOmega(double dTdt,double dTdx,double dTdy,
double &dOmega0,double &dOmegaX,double &dOmegaY,bool GGTt,bool GGTx,bool GGTy){
	double dV,tau=0.5*(newmesh->tau+mesh->tau);
	double dTdr=sqrt(dTdx*dTdx+dTdy*dTdy);
	bool success=false;
	dOmega0=dOmegaX=dOmegaY=0.0;
	if(fabs(dTdt)>(DX/DELTAU)*dTdr){
	//if(fabs(dTdt)>dTdr){
		if(GGTt){
			dV=tau*DX*DY;
			dOmega0=-dV*dTdt/fabs(dTdt);
			dOmegaX=-dV*dTdx/fabs(dTdt);
			dOmegaY=-dV*dTdy/fabs(dTdt);
			success=true;
			Omega0tot+=dOmega0;
		}
	}
	else if(fabs(dTdx)>fabs(dTdy)){
		if(GGTx){
			dV=tau*DELTAU*DY;
			dOmega0=-dV*dTdt/fabs(dTdx);
			dOmegaX=-dV*dTdx/fabs(dTdx);
			dOmegaY=-dV*dTdy/fabs(dTdx);
			success=true;
			OmegaXtot+=fabs(dOmegaX);
		}
	}
	else{
		if(GGTy){
			dV=tau*DELTAU*DX;
			dOmega0=-dV*dTdt/fabs(dTdy);
			dOmegaX=-dV*dTdx/fabs(dTdy);
			dOmegaY=-dV*dTdy/fabs(dTdy);
			success=true;
			OmegaYtot+=fabs(dOmegaY);
		}
	}
	return success;
}
