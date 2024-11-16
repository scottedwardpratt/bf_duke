#include "msu_commonutils/log.h"
#include "bfduke/hydro2uds.h"
#include "msu_commonutils/randy.h"
#include "bfduke/bfeos.h"
#include "bfduke/bfcharge.h"
#include "msu_sampler/hyper.h"

using namespace std;
using namespace NMSUPratt;

CHydroBalance *CHBHydroMesh::hb=NULL;

CHydroBalance::CHydroBalance(){
};

CHydroBalance::CHydroBalance(int run_number_set,int subrun_number_set){
	Omega0tot=OmegaXtot=OmegaYtot=netUdotOmega=NHYPER=0.0;
	run_number=run_number_set;
	subrun_number=subrun_number_set;
	string parfilename="modelruns/fixed_parameters.txt";
	parmap.ReadParsFromFile(parfilename);
	parfilename="modelruns/run"+to_string(run_number)+"/mod_parameters.txt";
	parmap.ReadParsFromFile(parfilename);
	Tf=0.001*parmap.getD("HYPER_FREEZEOUT_TEMP",160.0);
	epsilon_f=parmap.getD("HYPER_FREEZEOUT_EPSILON",0.36);
	SIGMA0=parmap.getD("BF_SIGMA0",0.5);
	DiffusionRatio=parmap.getD("BF_DIFFUSION_RATIO",1.0);
	CHBEoS::parmap=&parmap;
	CHBEoS::ReadEoSData_Andrew();
	CHBEoS::ReadDiffusionData_Andrew();
	CHBEoS::ReadChiReductionFactors();
 	
	CHBHydroMesh::DELTAU=parmap.getD("HYDRO_MESH_DELTAU",0.05);
	CHBHydroMesh::TAU0=parmap.getD("HYDRO_MESH_TAU0",0.6);
	CHBHydroMesh::XMIN=parmap.getD("HYDRO_MESH_XMIN",-25.0);
	CHBHydroMesh::XMAX=parmap.getD("HYDRO_MESH_XMAX",25.0);
	CHBHydroMesh::YMIN=parmap.getD("HYDRO_MESH_YMIN",-25.0);
	CHBHydroMesh::YMAX=parmap.getD("HYDRO_MESH_YMAX",25.0);
	CHBHydroMesh::NX=parmap.getI("HYDRO_MESH_NX",250);
	CHBHydroMesh::NY=parmap.getI("HYDRO_MESH_NY",250);
	CHBHydroMesh::DX=(CHBHydroMesh::XMAX-CHBHydroMesh::XMIN)/double(CHBHydroMesh::NX-1);
	CHBHydroMesh::DY=(CHBHydroMesh::YMAX-CHBHydroMesh::YMIN)/double(CHBHydroMesh::NY-1);
	TAU_EQ=parmap.getD("FUGACITY_TAU_EQ",0.5);
	GAMMA_0=parmap.getD("FUGACITY_GAMMA_0",1.0);
	CHBHydroMesh::TAU_EQ=TAU_EQ;
	CHBHydroMesh::GAMMA_0=GAMMA_0;
	
	DELTAU=CHBHydroMesh::DELTAU;
	TAU0=CHBHydroMesh::TAU0;
	XMIN=CHBHydroMesh::XMIN;
	XMAX=CHBHydroMesh::XMAX;
	YMIN=CHBHydroMesh::YMIN;
	YMAX=CHBHydroMesh::YMAX;
	DX=CHBHydroMesh::DX;
	DY=CHBHydroMesh::DY;
	
	CHBHydroMesh::GetDimensions(NX,NY,DX,DY,DELTAU,TAU0,XMIN,XMAX,YMIN,YMAX);
	WRITE_TRAJ=parmap.getB("HYDRO2UDS_WRITE_TRAJ",false);
	NSAMPLE_HYDRO2UDS=parmap.getD("HYDRO2UDS_NSAMPLE",1);
	randy=new Crandy(run_number);
	mesh=newmesh=oldmesh=NULL;
	Ncollisions=0;
	CHBHydroMesh::hb=this;
	CHBCharge::hb=this;
	tau0check=true;
	tau0readcheck=true;
	itauread=0;
	ransum=0.0;
	ranthresh=randy->ran_exp();
	idmax=0;
	MakeMeshes();
	
	hbcharges.resize(30000);
	nhbcharges=0;
	
	source.resize(30);
	for(int itau=0;itau<30;itau++)
		source[itau].setZero();
	ntraj=0;
}

void CHydroBalance::MakeMeshes(){
	mesh=new CHBHydroMesh();
	newmesh=new CHBHydroMesh();
	oldmesh=new CHBHydroMesh();
}

CHydroBalance::~CHydroBalance(){
}

void CHydroBalance::GetUxyBar(int ix,int iy,double &uxbar,double &uybar){
	uxbar=0.125*(mesh->UX[ix][iy]+mesh->UX[ix+1][iy+1]
		+mesh->UX[ix+1][iy]+mesh->UX[ix][iy+1]
			+newmesh->UX[ix][iy]+newmesh->UX[ix+1][iy+1]
				+newmesh->UX[ix+1][iy]+newmesh->UX[ix][iy+1]);
	uybar=0.125*(mesh->UY[ix][iy]+mesh->UY[ix+1][iy+1]
		+mesh->UY[ix+1][iy]+mesh->UY[ix][iy+1]
			+newmesh->UY[ix][iy]+newmesh->UY[ix+1][iy+1]
				+newmesh->UY[ix+1][iy]+newmesh->UY[ix][iy+1]);
}

void CHydroBalance::GetXYBar(int ix,int iy,double &xbar,double &ybar){
	int delix,deliy;
	double x,y;
	xbar=ybar=0.0;
	for(delix=0;delix<2;delix++){
		for(deliy=0;deliy<2;deliy++){
			mesh->GetXY(ix+delix,iy+deliy,x,y);
			xbar+=0.25*x;
			ybar+=0.25*y;
		}	
	}
}

void CHydroBalance::GetEpsilonBar(int ix,int iy,double &EpsilonBar){
	EpsilonBar=0.125*(mesh->epsilon[ix][iy]+mesh->epsilon[ix+1][iy+1]
		+mesh->epsilon[ix+1][iy]+mesh->epsilon[ix][iy+1]);
	EpsilonBar+=0.125*(newmesh->epsilon[ix][iy]+newmesh->epsilon[ix+1][iy+1]
		+newmesh->epsilon[ix+1][iy]+newmesh->epsilon[ix][iy+1]);
}


void CHydroBalance::GetFugacityBar(double &f_lBar,double &f_sBar){
	double oldtau,newtau,oldf_l,newf_l,oldf_s,newf_s,oldgamma,newgamma;
	oldtau=mesh->tau;
	newtau=mesh->tau;
	GetGammaFQ(oldtau,oldgamma,oldf_l,oldf_s);
	GetGammaFQ(newtau,newgamma,newf_l,newf_s);
	f_lBar=0.5*(oldf_l+newf_l);
	f_sBar=0.5*(oldf_s+newf_s);
}

void CHydroBalance::GetPiTildeBar(int ix,int iy,double &pitildexxbar,double &pitildexybar,
double &pitildeyybar){
	pitildexxbar=0.125*(mesh->pitildexx[ix][iy]+mesh->pitildexx[ix+1][iy+1]
		+mesh->pitildexx[ix+1][iy]+mesh->pitildexx[ix][iy+1]
			+newmesh->pitildexx[ix][iy]+newmesh->pitildexx[ix+1][iy+1]
				+newmesh->pitildexx[ix+1][iy]+newmesh->pitildexx[ix][iy+1]);
	pitildexybar=0.125*(mesh->pitildexy[ix][iy]+mesh->pitildexy[ix+1][iy+1]
		+mesh->pitildexy[ix+1][iy]+mesh->pitildexy[ix][iy+1]
			+newmesh->pitildexy[ix][iy]+newmesh->pitildexy[ix+1][iy+1]
				+newmesh->pitildexy[ix+1][iy]+newmesh->pitildexy[ix][iy+1]);
	pitildeyybar=0.125*(mesh->pitildeyy[ix][iy]+mesh->pitildeyy[ix+1][iy+1]
		+mesh->pitildeyy[ix+1][iy]+mesh->pitildeyy[ix][iy+1]
			+newmesh->pitildeyy[ix][iy]+newmesh->pitildeyy[ix+1][iy+1]
				+newmesh->pitildeyy[ix+1][iy]+newmesh->pitildeyy[ix][iy+1]);
}

void CHydroBalance::MakeCharges(){
	int ix,iy,sign,a,b,itau;
	double Ecutoff=0.01;
	double DQll,DQud,DQls,DQss,g1,g2;
	Eigen::Matrix<double,3,3> DQ;
	CHBCharge *charge1,*charge2;
	for(ix=1;ix<mesh->NX-1;ix++){
		for(iy=1;iy<mesh->NY-1;iy++){
			if(mesh->epsilon[ix][iy]>Ecutoff || newmesh->epsilon[ix][iy]>Ecutoff){
				if(tau0check){
					CalcDQ0(ix,iy,DQll,DQud,DQls,DQss);
				}
				else{
					CalcDQ(ix,iy,DQll,DQud,DQls,DQss);
				}
				DQ(0,0)=DQ(1,1)=DQll;
				DQ(0,1)=DQ(1,0)=DQud;
				DQ(0,2)=DQ(1,2)=DQ(2,1)=DQ(2,0)=DQls;
				DQ(2,2)=DQss;
				for(a=0;a<3;a++){
					for(b=0;b<3;b++){
						ransum+=NSAMPLE_HYDRO2UDS*fabs(DQ(a,b));
						while(ransum>ranthresh){
							ranthresh+=randy->ran_exp();
							if(nhbcharges>=hbcharges.size()){
								//hbcharges.resize(hbcharges.size()+10000);
								CLog::Fatal("hbcharges not big enough: ="+to_string(nhbcharges)+"\n");
							}
							charge1=&hbcharges[nhbcharges];
							nhbcharges+=1;
							if(nhbcharges>=hbcharges.size()){
								//hbcharges.resize(hbcharges.size()+10000);
								CLog::Fatal("hbcharges not big enough: ="+to_string(nhbcharges)+"\n");
							}
							charge2=&hbcharges[nhbcharges];
							nhbcharges+=1;
							mesh->GetXY(ix,iy,charge1->x,charge1->y);
							charge1->x=charge2->x=XMIN+(ix+randy->ran())*DX;
							charge1->y=charge2->y=YMIN+(iy+randy->ran())*DY;
							charge1->tau=charge2->tau=mesh->tau;
							charge1->eta=charge2->eta=0.0;
							if(tau0check){
								randy->ran_gauss2(g1,g2);
								charge1->eta=SIGMA0*g1;
								charge2->eta=SIGMA0*g2;
							}
							charge1->active=charge2->active=true;
							charge1->weight=charge2->weight=1.0;
							charge1->q[0]=charge1->q[1]=charge1->q[2]=0;
							charge2->q[0]=charge2->q[1]=charge2->q[2]=0;
							charge1->SetV(mesh->UX[ix][iy],mesh->UY[ix][iy]);
							charge2->SetV(mesh->UX[ix][iy],mesh->UY[ix][iy]);
							sign=1;
							if(randy->ran()<0.5)
								sign=-1;
							charge1->q[a]=sign;
							if(DQ(a,b)>0.0)
								sign=-sign;
							charge2->q[b]=sign;
							cmap.insert(pairic(idmax,charge1));
							idmax+=1;
							cmap.insert(pairic(idmax,charge2));
							if(WRITE_TRAJ && tau0check && randy->ran()<0.2){
								if(charge1->q[2]!=0 && charge2->q[2]!=0){ // write for ss CF
									charge1->trajinfo=new CTrajInfo(ntraj);
									ntraj+=1;
									charge2->trajinfo=new CTrajInfo(ntraj);
									ntraj+=1;
									charge1->addtraj();
									charge2->addtraj();
								}
							}
							idmax+=1;
							itau=floorl(mesh->tau/0.5);
							if(itau<30)
								source[itau](a,b)-=charge1->q[a]*charge2->q[b];

						}
					}
				}
			}
		}
	}
	tau0check=false;
}

void CHydroBalance::PropagateCharges(){
	char message[CLog::CHARLENGTH];
	mapic::iterator it,oldit,its;
	double dEdt,dEdx,dEdy,u0;
	bool GGEt,GGEx,GGEy;
	Chyper *hyper;
	CHBCharge *charge;
	int ix,iy,id;
	double newtau=newmesh->tau;
	it=cmap.begin();
	its=it;
	while(it!=cmap.end()){
		charge=it->second;
		id=it->first;
		charge->Propagate(newtau);
		mesh->GetIXIY_lower(charge->x,charge->y,ix,iy);
		if(ix<=0 || iy<=0 || ix>=CHBHydroMesh::NX-1 || iy>=CHBHydroMesh::NY-1){
			snprintf(message,CLog::CHARLENGTH,"CHydroBalance::PropagateCharges() disaster, ix=%d, iy=%d, NX=%d, NY=%d\n  epsilon=%g, epsilon_f=%g, T=%g\n",ix,iy,NX,NY,mesh->epsilon[ix][iy],epsilon_f,mesh->T[ix][iy]);
			CLog::Fatal(message);
		}
		if(!(charge->active) || (charge->active && newmesh->epsilon[ix][iy]<epsilon_f)){
			emap.insert(pairic(id,charge));
			hyper=&(charge->hyper);
			if(WRITE_TRAJ)
				charge->addtraj();
			GetGradEpsilon(ix,iy,dEdt,dEdx,dEdy,GGEt,GGEx,GGEy);
			GetUxyBar(ix,iy,hyper->u[1],hyper->u[2]);
			GetXYBar(ix,iy,hyper->r[1],hyper->r[2]);
			hyper->tau=newmesh->tau;
			GetPiTildeBar(ix,iy,hyper->pitilde[1][1],hyper->pitilde[1][2],hyper->pitilde[2][2]);
			hyper->dOmega[0]=-dEdt;
			hyper->dOmega[1]=-dEdx;
			hyper->dOmega[2]=-dEdy;
			u0=sqrt(1.0+hyper->u[1]*hyper->u[1]+hyper->u[2]*hyper->u[2]);
			hyper->udotdOmega=(u0*hyper->dOmega[0]
				-hyper->u[1]*hyper->dOmega[1]-hyper->u[2]*hyper->dOmega[2]);
			charge->active=false;
			oldit=it;
			++it;
			cmap.erase(oldit);
		}
		else{
			++it;
		}

		
	}
	it=cmap.begin();
	while(it!=cmap.end()){
		charge=it->second;
		if(!(charge->active)){
			CLog::Fatal("no dead particles should survive\n");
		}
		++it;
	}
}

void CHydroBalance::ScatterCharges(){
	mapic::iterator it;
	double u0,ux,uy,D;
	CHBCharge *charge;
	int ix,iy;
	it=cmap.begin();
	while(it!=cmap.end()){
		charge=it->second;
		mesh->GetIXIY(charge->x,charge->y,ix,iy);
		D=CHBEoS::GetD(newmesh->epsilon[ix][iy]);
		D*=DiffusionRatio;
		ux=newmesh->UX[ix][iy];
		uy=newmesh->UY[ix][iy];
		u0=sqrt(1.0+ux*ux+uy*uy);
		ransum+=DELTAU/(6.0*D*u0);
		while(ransum>ranthresh){
			ranthresh+=randy->ran_exp();
			charge->SetV(ux,uy);
			if(WRITE_TRAJ)
				charge->addtraj();
			Ncollisions+=1;
		}
		++it;
	}
}

void CHydroBalance::CalcDQ(int ix,int iy,double &DQll,
double &DQud,double &DQls,double &DQss){
	double d4x,s0,sx,sy;
	double f_l,f_s,gamma_q; // fugacities for light and strange quarks
	CHBEoS eos222[2][2][2];
	int j0,jx,jy;
	double ux[2][2][2],uy[2][2][2],u0[2][2][2];
	CHBHydroMesh *mptr;
	
	for(j0=0;j0<2;j0++){
		if(j0==0)
			mptr=mesh;
		else
			mptr=newmesh;
		for(jx=0;jx<2;jx++){
			for(jy=0;jy<2;jy++){
				ux[j0][jx][jy]=mptr->UX[ix+jx][iy+jy];
				uy[j0][jx][jy]=mptr->UY[ix+jx][iy+jy];
				u0[j0][jx][jy]=sqrt(1.0
					+ux[j0][jx][jy]*ux[j0][jx][jy]+uy[j0][jx][jy]*uy[j0][jx][jy]);
				eos222[j0][jx][jy].epsilon=mptr->epsilon[ix+jx][iy+jy];
				GetGammaFQ(mptr->tau,gamma_q,f_l,f_s);
				eos222[j0][jx][jy].f_u=f_l;
				eos222[j0][jx][jy].f_d=f_l;
				eos222[j0][jx][jy].f_s=f_s;
				eos222[j0][jx][jy].SetChi();
			}
		}
	}
	
	DQll=DQud=DQls=DQss=0.0;
	for(j0=0;j0<2;j0++){
		if(j0==0){
			s0=-0.25;
			d4x=DELTAU*DX*DY*mesh->tau;
			GetGammaFQ(mesh->tau,gamma_q,f_l,f_s);
		}
		else{
			s0=0.25;
			d4x=DELTAU*DX*DY*newmesh->tau;
		}
		for(jx=0;jx<2;jx++){
			if(jx==0)
				sx=-0.25;
			else
				sx=0.25;
			for(jy=0;jy<2;jy++){
				if(jy==0)
					sy=-0.25;
				else
					sy=0.25;
				
				DQll+=d4x*s0*u0[j0][jx][jy]*eos222[j0][jx][jy].chill/DELTAU;
				DQud+=d4x*s0*u0[j0][jx][jy]*eos222[j0][jx][jy].chiud/DELTAU;
				DQls+=d4x*s0*u0[j0][jx][jy]*eos222[j0][jx][jy].chils/DELTAU;
				DQss+=d4x*s0*u0[j0][jx][jy]*eos222[j0][jx][jy].chiss/DELTAU;
				
				DQll+=d4x*sx*ux[j0][jx][jy]*eos222[j0][jx][jy].chill/DX;
				DQud+=d4x*sx*ux[j0][jx][jy]*eos222[j0][jx][jy].chiud/DX;
				DQls+=d4x*sx*ux[j0][jx][jy]*eos222[j0][jx][jy].chils/DX;
				DQss+=d4x*sx*ux[j0][jx][jy]*eos222[j0][jx][jy].chiss/DX;
				
				DQll+=d4x*sy*uy[j0][jx][jy]*eos222[j0][jx][jy].chill/DY;
				DQud+=d4x*sy*uy[j0][jx][jy]*eos222[j0][jx][jy].chiud/DY;
				DQls+=d4x*sy*uy[j0][jx][jy]*eos222[j0][jx][jy].chils/DY;
				DQss+=d4x*sy*uy[j0][jx][jy]*eos222[j0][jx][jy].chiss/DY;
								
			}
		}
	}
}

void CHydroBalance::CalcDQ0(int ix,int iy,double &DQll,double &DQud,
double &DQls,double &DQss){
	CHBEoS eos;
	double u0,ux,uy,d3x=DX*DY*newmesh->tau,epsilon;
	double gamma_q,f_l,f_s;
	GetUxyBar(ix,iy,ux,uy);
	u0=sqrt(1.0+ux*ux+uy*uy);
	GetEpsilonBar(ix,iy,epsilon);
	if(epsilon>epsilon_f){
		GetGammaFQ(newmesh->tau,gamma_q,f_l,f_s);
		eos.epsilon=epsilon;
		eos.f_u=eos.f_d=f_l;
		eos.f_s=f_s;
		eos.SetChi();
		DQll=d3x*u0*eos.chill;
		DQud=d3x*u0*eos.chiud;
		DQls=d3x*u0*eos.chils;
		DQss=d3x*u0*eos.chiss;
	}
	else{
		DQll=DQud=DQls=DQss=0.0;
	}
}

void CHydroBalance::SwapMeshes(){
	CHBHydroMesh *swap;
	swap=oldmesh;
	oldmesh=mesh;
	mesh=newmesh;
	newmesh=swap;
}

void CHydroBalance::Reset(){
	Ncollisions=0;
	nhbcharges=0;
	itauread=0;
	idmax=0;
	tau0readcheck=tau0check=true;
	oldmesh->tau=mesh->tau=newmesh->tau=CHBHydroMesh::TAU0;
}

void CHydroBalance::GetGammaFQ(double tau,double &gamma_q,double &fugacity_l,double &fugacity_s){
	double fugacity_meson;
	gamma_q=1.0-(1.0-GAMMA_0)*exp((TAU0-tau)/TAU_EQ);
	fugacity_meson=0.85*gamma_q+0.15;
	fugacity_l=fugacity_s=sqrt(fugacity_meson);	
}


