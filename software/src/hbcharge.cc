#include "msu_commonutils/log.h"
#include "bfduke/bfcharge.h"
#include "msu_commonutils/misc.h"
#include "bfduke/hydro2uds.h"
#include "msu_commonutils/randy.h"
#include "msu_commonutils/constants.h"
#include "msu_sampler/hyper.h"

using namespace std;
using namespace NMSUPratt;

CHydroBalance *CHBCharge::hb=NULL;

void CHBCharge::Propagate(double newtau){
	double t0,tf,neweta;
	if(active==true){
		t0=tau*cosh(eta);
		neweta=rapidity-asinh((tau/newtau)*sinh(rapidity-eta));
		tf=newtau*cosh(neweta);
		x+=vx*(tf-t0);
		y+=vy*(tf-t0);
		if(x!=x)
			CLog::Fatal("x!=x inside CHBCharge::Propagate\n");
		tau=newtau;
		eta=neweta;
	}
}

void CHBCharge::SetV(double uxmatter,double uymatter){
	double vz,vperp,phi,vmag=1.0;
	FourVector v,umatter;
	vz=vmag*(1.0-2.0*hb->randy->ran());
	phi=2.0*PI*hb->randy->ran();
	vperp=sqrt(vmag*vmag-vz*vz);
	vx=vperp*cos(phi);
	vy=vperp*sin(phi);
	v[0]=1.0;
	v[1]=vx;
	v[2]=vy;
	v[3]=vz;
	umatter[3]=sinh(eta);
	umatter[1]=uxmatter; umatter[2]=uymatter;
	umatter[0]=sqrt(1.0+uxmatter*uxmatter+uymatter*uymatter+umatter[3]*umatter[3]);
	Misc::Boost(umatter,v);
	rapidity=atanh(v[3]/v[0]);
	vx=v[1]/v[0];
	vy=v[2]/v[0];
}	

void CHBCharge::Print(){
	char message[CLog::CHARLENGTH];
	CLog::Info("Charge Info:\n");
	if(active)
		CLog::Info("active\n");
	else
		CLog::Info("dead\n");
	snprintf(message,CLog::CHARLENGTH,"q = (%d,%d,%d)\n",q[0],q[1],q[2]);
	CLog::Info(message);
	snprintf(message,CLog::CHARLENGTH,
	"x = %g,y=%g, eta=%g, tau=%g\nvx=%g, vy=%g, rapidity=%g, v=%g\n",
	x,y,eta,tau,vx,vy,rapidity,sqrt(vx*vx+vy*vy+tanh(rapidity)*tanh(rapidity)));
	CLog::Info(message);
	CLog::Info("------------------------------------------------------------\n");
}

