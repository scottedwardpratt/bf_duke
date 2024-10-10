#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>
#include "bfduke/hydro2uds.h"
#include "msu_sampler/hyper.h"
#include "msu_commonutils/qualifier.h"

using namespace std;
using namespace NMSUPratt;

int main(int argc,char *argv[]){
	if (argc != 3) {
		printf("Usage: hydro2uds run_number subrun_number\n");
		exit(-1);
	}
	char message[CLog::CHARLENGTH];
	bool oscarfile=true;
	int ievent,nevents,run_number=atoi(argv[1]),subrun_number=atoi(argv[1]);
	CHydroBalance hb(run_number,subrun_number);
	nevents=hb.parmap.getI("HYDRO_NEVENTS",10);
	CQualifiers qualifiers;
	qualifiers.Read("qualifiers.txt");
	for(int iqual=0;iqual<qualifiers.nqualifiers;iqual++){
		hb.qualifier=qualifiers.qualifier[iqual]->qualname;
		printf("--------- BEGIN CALC FOR %s ---------\n",hb.qualifier.c_str());
		for(ievent=0;ievent<nevents;ievent++){
			hb.Omega0tot=hb.OmegaXtot=hb.OmegaYtot=0.0;
			hb.Reset();
			oscarfile=hb.ReadDuke(hb.mesh);
			hb.HyperFind();
			oscarfile=hb.ReadDuke(hb.newmesh);
			hb.HyperFind();
			hb.MakeCharges();
			hb.PropagateCharges();
			do{
				hb.SwapMeshes();
				oscarfile=hb.ReadDuke(hb.newmesh);
				hb.HyperFind();
				hb.MakeCharges();
				hb.PropagateCharges();
				hb.ScatterCharges();
				if(fabs(lrint(hb.mesh->tau)-hb.mesh->tau)<0.001){
					snprintf(message,CLog::CHARLENGTH,"tau=%g, cmap.size=%lu, emap.size=%lu\n",hb.mesh->tau,
					hb.cmap.size(),hb.emap.size());
					CLog::Info(message);
					snprintf(message,CLog::CHARLENGTH,"highestEpsilon=%g, biggestU=%g\n",hb.highestEpsilon,hb.biggestU);
					CLog::Info(message);
				}
			}while(oscarfile);
			snprintf(message,CLog::CHARLENGTH,"final tau=%g, cmap.size=%lu, emap.size=%lu\n",hb.mesh->tau,
			hb.cmap.size(),hb.emap.size());
			CLog::Info(message);
			snprintf(message,CLog::CHARLENGTH,"highestEpsilon=%g, biggestU=%g\n",hb.highestEpsilon,hb.biggestU);
			CLog::Info(message);
			printf("netUdotOmega=%g, Omega0tot=%g, OmegaXtot=%g, OmegaYtot=%g\n",hb.netUdotOmega,hb.Omega0tot,hb.OmegaXtot,hb.OmegaYtot);
	
			hb.WriteCharges(ievent);
			if(run_number==0 && ievent==0)
				hb.WriteHyper_Duke_2D();
			hb.ClearCharges();
		}
		
	}
	
	return 0;
}


