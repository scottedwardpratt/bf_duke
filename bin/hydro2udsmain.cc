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
	int ievent,nevents,run_number=atoi(argv[1]),subrun_number=atoi(argv[2]);
	int subrun_number_max=10000;
	CHydroBalance hb(run_number,subrun_number);
	nevents=hb.parmap.getI("HYDRO_NEVENTS",10);
	CQualifiers qualifiers;
	qualifiers.Read("qualifiers.txt");
	string logfilename="logfiles/run"+to_string(run_number)+"_subrun"+to_string(subrun_number)+".txt";
	CLog::Init(logfilename);
	for(int iqual=0;iqual<qualifiers.nqualifiers;iqual++){
		hb.qualifier=qualifiers.qualifier[iqual]->qualname;
		CLog::Info("--------- BEGIN CALC FOR "+hb.qualifier+" ---------\n");
		for(ievent=0;ievent<nevents;ievent++){
			CLog::Info("xxxx begin for ievent= "+to_string(ievent)+"\n");
			hb.Omega0tot=hb.OmegaXtot=hb.OmegaYtot=0.0;
			hb.Reset();
			hb.randy->reset(nevents*subrun_number_max*run_number+nevents*subrun_number+ievent);
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
			snprintf(message,CLog::CHARLENGTH,"netUdotOmega=%g, Omega0tot=%g, OmegaXtot=%g, OmegaYtot=%g\n",hb.netUdotOmega,hb.Omega0tot,hb.OmegaXtot,hb.OmegaYtot);
			CLog::Info(message);
	
			hb.WriteCharges(ievent);
			if(run_number==0 && ievent==0)
				hb.WriteHyper_Duke_2D();
			hb.ClearCharges();
		}
		
	}
	CLog::Info("----- FINISHING HAPPILY -----\n");
	return 0;
}


