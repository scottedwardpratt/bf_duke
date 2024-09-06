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
				if(GetGradEpsilon(ix,iy,dEdt,dEdx,dEdy,GGEt,GGEx,GGEy)){
					hyper.tau=0.5*(newmesh->tau+mesh->tau);
					GetXYBar(ix,iy,hyper.r[1],hyper.r[2]);
					GetUxyBar(ix,iy,hyper.u[1],hyper.u[2]);
					hyper.u[3]=0.0;
					hyper.u[0]=sqrt(1.0+hyper.u[1]*hyper.u[1]+hyper.u[2]*hyper.u[2]);
					GetPiTildeBar(ix,iy,hyper.pitilde[1][1],hyper.pitilde[1][2],
					hyper.pitilde[2][2]);
					double T0,epsilon0;
					GetTBar(ix,iy,T0);
					GetEpsilonBar(ix,iy,epsilon0);
					hyper.T0=T0;
					hyper.epsilon=epsilon0;
					if(GetDOmega(dEdt,dEdx,dEdy,
					hyper.dOmega[0],hyper.dOmega[1],hyper.dOmega[2],GGEt,GGEx,GGEy)){
						GetGammaFQ(hyper.tau,gamma_q,fugacity_l,fugacity_s);
						hyper.fugacity_u=hyper.fugacity_d=fugacity_l;
						hyper.fugacity_s=fugacity_s;
						newhyper=new Chyper;
						newhyper->Copy(&hyper);
						hyperlist.push_back(newhyper);
						double UdotdOmega=hyper.u[0]*hyper.dOmega[0]-hyper.u[1]*hyper.dOmega[1]-hyper.u[2]*hyper.dOmega[2];
						netUdotOmega+=UdotdOmega;
						//printf("%4.2f %g: udotdOmega/tau=%g, u[0]=%g\n",hyper.tau,sqrt(hyper.r[1]*hyper.r[1]+hyper.r[2]*hyper.r[2]),
						//UdotdOmega/hyper.tau,hyper.u[0]);
						/*if(iy==NY/2 && hyper.tau>10 && hyper.epsilon<0.35){
							printf("tau=%g, ix=%d, iy=%d, T0=%g, epsilon=%g, NX=%d, NY=%d\n",
							hyper.tau,ix,iy,T0,hyper.epsilon,NX,NY);
							Misc::Pause();
						}*/
						
					}
				}
			}
		}
	}
}

void CHydroBalance::HyperFindF(){
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
				if(GetGradF(ix,iy,dFdt,dFdx,dFdy,GGFt,GGFx,GGFy)){
					
					hyper.tau=0.5*(newmesh->tau+mesh->tau);
					GetXYBar(ix,iy,hyper.r[1],hyper.r[2]);
					GetUxyBar(ix,iy,hyper.u[1],hyper.u[2]);
					hyper.u[0]=sqrt(1.0+hyper.u[1]*hyper.u[1]+hyper.u[2]*hyper.u[2]);
					GetPiTildeBar(ix,iy,hyper.pitilde[1][1],hyper.pitilde[1][2],
					hyper.pitilde[2][2]);
					double T0,epsilon0;
					GetTBar(ix,iy,T0);
					GetEpsilonBar(ix,iy,epsilon0);
					hyper.T0=T0;
					hyper.epsilon=epsilon0;
					
					if(GGFt){
					//if(epsilon0<2.0*epsilon_f && GetDOmega(dFdt,dFdx,dFdy,
					//hyper.dOmega[0],hyper.dOmega[1],hyper.dOmega[2],GGFt,GGFx,GGFy)){
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
						//printf("%4.2f %g: udotdOmega/tau=%g, u[0]=%g\n",hyper.tau,sqrt(hyper.r[1]*hyper.r[1]+hyper.r[2]*hyper.r[2]),
						//UdotdOmega/hyper.tau,hyper.u[0]);
						/*if(iy==NY/2 && hyper.tau>10 && hyper.epsilon<0.35){
							printf("tau=%g, ix=%d, iy=%d, T0=%g, epsilon=%g, NX=%d, NY=%d\n",
							hyper.tau,ix,iy,T0,hyper.epsilon,NX,NY);
							Misc::Pause();
						}*/
						
					}
				}
			}
		}
	}
	if(nhyper>0)
	  printf("nhyper=%d, tau=%g, RMAX=%g, UMAX=%g, TMAX=%g\n",nhyper,0.5*(newmesh->tau+mesh->tau),RMAX,UMAX,TMAX);
}

bool CHydroBalance::GetGradF(int ix,int iy,
double &dFdt,double &dFdx,double &dFdy,bool &GGFt,bool &GGFx,bool &GGFy){
	bool gg=false;
	if(HYPEREPSILON){
		gg=GetGradEpsilon(ix,iy,dFdt,dFdx,dFdy,GGFt,GGFx,GGFy);
	}
	else if(HYPERT){
		gg=GetGradT(ix,iy,dFdt,dFdx,dFdy,GGFt,GGFx,GGFy);
	}
	else{
		CLog::Fatal("HYPEREPSILON or HYPERT must be true\n");
	}
	return gg;
}


/*
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
					double T0,epsilon0;
					GetTBar(ix,iy,T0);
					GetEpsilonBar(ix,iy,epsilon0);
					hyper.T0=T0;
					hyper.epsilon=epsilon0;
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
						double UdotdOmega=hyper.u[0]*hyper.dOmega[0]-hyper.u[1]*hyper.dOmega[1]-hyper.u[2]*hyper.dOmega[2];
						netUdotOmega+=UdotdOmega;
					}
				}
			}
		}
	}
}
*/

bool CHydroBalance::GetGradEpsilon(int ix,int iy,
double &dEdt,double &dEdx,double &dEdy,bool &GGEt,bool &GGEx,bool &GGEy){
	bool hypercheck=false;
	double Explus,Exminus,Eyplus,Eyminus,Etplus,Etminus;
	mesh->GetDimensions(NX,NY,DX,DY,DELTAU,TAU0,XMIN,XMAX,YMIN,YMAX);
	
	/*if(mesh->T[ix][iy]<0.13 && mesh->T[ix+1][iy]<0.13 && mesh->T[ix][iy+1]<0.13 && mesh->T[ix+1][iy+1]<0.13
	&& newmesh->T[ix][iy]<0.13 && newmesh->T[ix+1][iy]<0.13 && newmesh->T[ix][iy+1]<0.13 && newmesh->T[ix+1][iy+1]<0.13){
		return hypercheck;
	}*/
	
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
		//printf("T=%g,%g,%g,%g,%g,%g,%g,%g\n",mesh->T[ix][iy],mesh->T[ix+1][iy],mesh->T[ix][iy+1],mesh->T[ix+1][iy+1],
		//newmesh->T[ix][iy],newmesh->T[ix+1][iy],newmesh->T[ix][iy+1],newmesh->T[ix+1][iy+1]);
	}
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

bool CHydroBalance::GetDOmega(double dFdt,double dFdx,double dFdy,
double &dOmega0,double &dOmegaX,double &dOmegaY,bool GGFt,bool GGFx,bool GGFy){
	double dV,tau=0.5*(newmesh->tau+mesh->tau);
	double dFdr=sqrt(dFdx*dFdx+dFdy*dFdy);
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
