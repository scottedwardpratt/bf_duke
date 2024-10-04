#include "msu_commonutils/log.h"
#include "bfduke/hydro2uds.h"
#include "msu_sampler/hyper.h"
#include "bfduke/bfeos.h"
#include "msu_commonutils/misc.h"

using namespace std;
using namespace NMSUPratt;

void CHydroBalance::HyperFind(){
	CHBEoS hypereos;
	char message[CLog::CHARLENGTH];
	int ix,iy;//,a,b;
	double dFdx,dFdy,dFdt,fugacity_l,fugacity_s,gamma_q;
	double RMAX=0.0,UMAX=0.0,TMAX=0.0;
	Chyper hyper;
	Chyper *newhyper;
	int nhyper=0;
	
	bool GGFt,GGFx,GGFy;
	if(!tau0check){
		for(ix=1;ix<mesh->NX-1;ix++){
			for(iy=1;iy<mesh->NY-1;iy++){
				if(ix<0 || iy<0 || ix==CHBHydroMesh::NX || iy==CHBHydroMesh::NY){
					snprintf(message,CLog::CHARLENGTH,"CHydroBalance::HyperFindEpsilon() -- ix=%d, iy=%d\n",ix,iy);
					CLog::Info(message);
				}
				if(GetGradEpsilon(ix,iy,dFdt,dFdx,dFdy,GGFt,GGFx,GGFy)){
					
					hyper.tau=0.5*(newmesh->tau+mesh->tau);
					GetXYBar(ix,iy,hyper.r[1],hyper.r[2]);
					GetUxyBar(ix,iy,hyper.u[1],hyper.u[2]);
					hyper.u[0]=sqrt(1.0+hyper.u[1]*hyper.u[1]+hyper.u[2]*hyper.u[2]);
					GetPiTildeBar(ix,iy,hyper.pitilde[1][1],hyper.pitilde[1][2],
					hyper.pitilde[2][2]);
					double epsilon0;
					GetEpsilonBar(ix,iy,epsilon0);
					hypereos.epsilon=epsilon0;
					hypereos.SetTnonequil();
					hyper.T0=hypereos.Tnonequil;
					hyper.epsilon=epsilon0;
					
					if(GGFt){
						GetDOmega(dFdt,dFdx,dFdy,hyper.dOmega[0],hyper.dOmega[1],hyper.dOmega[2],GGFt,GGFx,GGFy);
						
						nhyper+=1;
						
						GetGammaFQ(hyper.tau,gamma_q,fugacity_l,fugacity_s);
						hyper.fugacity_u=hyper.fugacity_d=fugacity_l;
						hyper.fugacity_s=fugacity_s;
						double r=sqrt(hyper.r[1]*hyper.r[1]+hyper.r[2]*hyper.r[2]);
						if(r>RMAX)
							RMAX=r;
						double umag=sqrt(hyper.u[1]*hyper.u[1]+hyper.u[2]*hyper.u[2]);
						if(umag>UMAX)
							UMAX=umag;
						if(hyper.T0>TMAX)
							TMAX=hyper.T0;
					
						newhyper=new Chyper;
						newhyper->Copy(&hyper);
						hyperlist.push_back(newhyper);
						double UdotdOmega=hyper.u[0]*hyper.dOmega[0]-hyper.u[1]*hyper.dOmega[1]-hyper.u[2]*hyper.dOmega[2];
						netUdotOmega+=UdotdOmega;
						
					}
				}
			}
		}
	}
	if(nhyper>0){
		NHYPER+=nhyper;
		printf("NHYPER=%d, nhyper=%d, tau=%g, RMAX=%g, UMAX=%g, TMAX=%g\n",
		NHYPER,nhyper,0.5*(newmesh->tau+mesh->tau),RMAX,UMAX,TMAX);
	}
}

bool CHydroBalance::GetGradEpsilon(int ix,int iy,
double &dEdt,double &dEdx,double &dEdy,bool &GGEt,bool &GGEx,bool &GGEy){
	bool hypercheck=false;
	double Explus,Exminus,Eyplus,Eyminus,Etplus,Etminus;
	mesh->GetDimensions(NX,NY,DX,DY,DELTAU,TAU0,XMIN,XMAX,YMIN,YMAX);
	
	Etminus=0.25*(mesh->epsilon[ix][iy]+mesh->epsilon[ix][iy+1]
		+mesh->epsilon[ix+1][iy]+mesh->epsilon[ix+1][iy+1]);
	Etplus=0.25*(newmesh->epsilon[ix][iy]+newmesh->epsilon[ix][iy+1]
		+newmesh->epsilon[ix+1][iy]+newmesh->epsilon[ix+1][iy+1]);
	
	Exminus=0.25*(mesh->epsilon[ix][iy]+newmesh->epsilon[ix][iy]
		+mesh->epsilon[ix][iy+1]+newmesh->epsilon[ix][iy+1]);
	Explus=0.25*(mesh->epsilon[ix+1][iy]+newmesh->epsilon[ix+1][iy]
		+mesh->epsilon[ix+1][iy+1]+newmesh->epsilon[ix+1][iy+1]);
	
	Eyminus=0.25*(mesh->epsilon[ix][iy]+newmesh->epsilon[ix][iy]
		+mesh->epsilon[ix+1][iy]+newmesh->epsilon[ix+1][iy]);
	Eyplus=0.25*(mesh->epsilon[ix][iy+1]+newmesh->epsilon[ix][iy+1]
		+mesh->epsilon[ix+1][iy+1]+newmesh->epsilon[ix+1][iy+1]);
	
	
	dEdt=(Etplus-Etminus)/DELTAU;
	dEdx=-(Explus-Exminus)/DX;
	dEdy=-(Eyplus-Eyminus)/DY;
	
	GGEt=GGEx=GGEy=false;
	if((Etplus-epsilon_f)*(Etminus-epsilon_f)<0.0){
		GGEt=true;
	}
	if((Explus-epsilon_f)*(Exminus-epsilon_f)<0.0)
		GGEx=true;
	if((Eyplus-epsilon_f)*(Eyminus-epsilon_f)<0.0)
		GGEy=true;
	if(GGEt || GGEx || GGEy)
		hypercheck=true;
	return hypercheck;
}

bool CHydroBalance::GetDOmega(double dFdt,double dFdx,double dFdy,
double &dOmega0,double &dOmegaX,double &dOmegaY,bool GGFt,bool GGFx,bool GGFy){
	double dV,tau=0.5*(newmesh->tau+mesh->tau);
	(void) GGFx;
	(void) GGFy;
	//double dFdr=sqrt(dFdx*dFdx+dFdy*dFdy);
	bool success=false;
	dOmega0=dOmegaX=dOmegaY=0.0;
	

	if(GGFt){
		dV=tau*DX*DY;
		dOmega0=-dV*dFdt/fabs(dFdt);
		dOmegaX=-dV*dFdx/fabs(dFdt);
		dOmegaY=-dV*dFdy/fabs(dFdt);
		success=true;
		Omega0tot+=fabs(dOmega0);
		OmegaXtot+=fabs(dOmegaX);
		OmegaYtot+=fabs(dOmegaY);
	}

	
	/*
	if(fabs(dFdt*DELTAU)>(DX*dFdr)){
	//if(fabs(dTdt)>dTdr){
	if(GGFt){
	dV=tau*DX*DY;
	dOmega0=-dV*dFdt/fabs(dFdt);
	dOmegaX=-dV*dFdx/fabs(dFdt);
	dOmegaY=-dV*dFdy/fabs(dFdt);
	success=true;
	Omega0tot+=fabs(dOmega0);
	OmegaXtot+=fabs(dOmegaX);
	OmegaYtot+=fabs(dOmegaY);
	}
	}
	else if(fabs(dFdx)>fabs(dFdy)){
	if(GGFx){
	dV=tau*DELTAU*DY;
	dOmega0=-dV*dFdt/fabs(dFdx);
	dOmegaX=-dV*dFdx/fabs(dFdx);
	dOmegaY=-dV*dFdy/fabs(dFdx);
	success=true;
	Omega0tot+=fabs(dOmega0);
	OmegaXtot+=fabs(dOmegaX);
	OmegaYtot+=fabs(dOmegaY);
	}
	}
	else{
	if(GGFy){
	dV=tau*DELTAU*DX;
	dOmega0=-dV*dFdt/fabs(dFdy);
	dOmegaX=-dV*dFdx/fabs(dFdy);
	dOmegaY=-dV*dFdy/fabs(dFdy);
	success=true;
	Omega0tot+=fabs(dOmega0);
	OmegaXtot+=fabs(dOmegaX);
	OmegaYtot+=fabs(dOmegaY);
	}
	}
	*/

	return success;
}
