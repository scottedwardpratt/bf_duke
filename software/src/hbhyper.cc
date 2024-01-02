#include "bfduke/bfhyper.h"
#include "msu_sampler/sampler.h"
#include "msu_boltzmann/msupart.h"
#include "msu_commonutils/randy.h"
#include "msu_commonutils/misc.h"

using namespace std;
using namespace NMSUPratt;

//#define __XYREFLECT__

Csampler* CHBHyperElement::sampler=NULL;

CHBHyperElement::CHBHyperElement(){
	//
}

void CHBHyperElement::Copy(CHBHyperElement *oldhyper){
	tau=oldhyper->tau;
	x=oldhyper->x; y=oldhyper->y;
	dOmegaX=oldhyper->dOmegaX;
	dOmegaY=oldhyper->dOmegaY;
	dOmega0=oldhyper->dOmega0;
	ux=oldhyper->ux;
	uy=oldhyper->uy;
	udotdOmega=oldhyper->udotdOmega;
	pitildexx=oldhyper->pitildexx;
	pitildexy=oldhyper->pitildexy;
	pitildeyy=oldhyper->pitildeyy;
	T=oldhyper->T;
}

/*
void CHyperElement::CalcDOmegaMax(){
	double dOmegaVec;
	double u0=sqrt(1.0+ux*ux+uy*uy);
	double dOmega2=dOmega0*dOmega0-dOmegaX*dOmegaX-dOmegaY*dOmegaY;
	udotdOmega=u0*dOmega0-ux*dOmegaX-uy*dOmegaY;
	dOmegaVec=sqrt(-dOmega2+udotdOmega*udotdOmega);
	dOmegaMax=fabs(udotdOmega)+dOmegaVec;
	vOmega=-udotdOmega/dOmegaVec;
}*/

void CHBHyperElement::Print(){
	printf("HyperElement Info:\n");
	printf("tau=%g, T=%g, x=%g, y=%g, ux=%g, uy=%g, dOmega0=%g, dOmegaX=%g, dOmegaY=%g, udotdOmega=%g\n",
	tau,T,x,y,ux,uy,dOmega0,dOmegaX,dOmegaY,udotdOmega);
	printf("pitildexx=%g, pitildeyy=%g, pitildexy=%g\n",
	pitildexx,pitildeyy,pitildexy);
	printf("---------------------------------------------------------------\n");
}

int CHBHyperElement::MakeParts(CMSU_Boltzmann *boltzmann){
	int nparts=0;
	CMSUPart *part;
	double bweight,mass,eta,rapidity,ETAMAX=boltzmann->ETAMAX;
	int ires,nsample=sampler->NSAMPLE;
	FourVector plab;
	double delN,r[3];
	double delNtot=nhadrons*udotdOmega*nsample;
	CresInfoMap *resmap=&(sampler->reslist->resmap);
	CresInfoMap::iterator rpos;
	CresInfo *resinfo;
	Crandy *randy=sampler->randy;
	if(randy->netprob+delNtot > randy->threshold){
		for(rpos=resmap->begin();rpos!=resmap->end();rpos++){
			resinfo=rpos->second;
			if(resinfo->pid==22){
				rpos++;
				resinfo=rpos->second;
			}
			ires=resinfo->ires;
			delN=(*density)[ires]*udotdOmega*nsample;
			randy->netprob+=delN;
			while(randy->netprob>randy->threshold){
				GetP(resinfo,plab,mass,(*maxweight)[ires]);
				part=boltzmann->GetDeadPart();	
#ifdef __XY_REFLECT__	
				printf("HOWDY, I am XY reflecting\n");	
				if(randy->ran()<0.5){
					r[1]=-r[1];
					plab[1]=-plab[1];
				}
				if(randy->ran()<0.5){
					r[2]=-r[2];
					plab[2]=-plab[2];
				}
#endif
				r[0]=tau; r[1]=x; r[2]=y;
				bweight=1.0;
				rapidity=atanh(plab[3]/plab[0]);
				eta=(1.0-2.0*randy->ran())*ETAMAX;
				rapidity+=eta;
				if(fabs(eta)>ETAMAX){
					printf("eta=%g ??, ETAMAX=%g\n",eta,ETAMAX);
					exit(1);
				}
				part->InitBalance(resinfo->pid,r[1],r[2],r[0],eta,plab[1],plab[2],mass,rapidity,bweight,-1);
				/*
				if(part->balanceID<0){
					for(int a=0;a<3;a++){
						for(int b=0;b<3;b++){
							sampler->chitot(a,b)+=resinfo->q[a]*resinfo->q[b];
						}
					}
				}
				*/
				randy->threshold-=log(randy->ran());
				nparts+=1;
			}
		}	
	}
	else
		randy->netprob+=delNtot;
	return nparts;
}
