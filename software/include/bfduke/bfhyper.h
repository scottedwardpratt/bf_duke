#ifndef __BFHYPER_H__
#define __BFHYPER_H__

#include "msu_commonutils/commondefs.h"
#include "bfduke/bfcommon.h"
using namespace std;
using namespace NMSUPratt;

// Info for HyperSurface
class CHBHyperElement{
public:
	CHBHyperElement();
	~CHBHyperElement(){};
	double tau,x,y,dOmegaX,dOmegaY,dOmega0,ux,uy,T;
	vector<double> *density,*maxweight;
	double epsilon,P,h,nhadrons,lambda;
	double pitildexx,pitildexy,pitildeyy;
	double udotdOmega;
	void GetP(CresInfo *resinfo,FourVector &p,double &mass,double maxweight);
	void Copy(CHBHyperElement *oldhyper);
	int MakeParts();
	void Print();
	// Note dOmegaX and dOmegaY = d\Omega_\mu
	// (subscript mu, so p.Omega = p0*dOmega0+px*dOmegaX+py*dOmegaY)
	static Csampler *sampler;
	static CHydroBalance *hb;
};

#endif