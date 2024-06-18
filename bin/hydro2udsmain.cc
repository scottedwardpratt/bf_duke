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
	if (argc != 2) {
		printf("Usage: hydro2uds run_number\n");
		exit(-1);
  }
	char message[CLog::CHARLENGTH];
	bool oscarfile=true;
	int run_number=atoi(argv[1]);
	string udsfilename="uds"+string(argv[1])+".txt";
	CHydroBalance hb("udsdata/udsparameters.txt",run_number);
	hb.parmap.set("CHARGESINFO_FILENAME",udsfilename);
	CQualifiers qualifiers;
	qualifiers.Read("qualifiers.txt");
	for(int iqual=0;iqual<qualifiers.nqualifiers;iqual++){
		hb.qualifier=qualifiers.qualifier[iqual]->qualname;
		hb.Reset();
		printf("--------- BEGIN CALC FOR %s ---------\n",hb.qualifier.c_str());
		oscarfile=hb.ReadDuke(hb.mesh);
		
		hb.HyperFindEpsilon();
		oscarfile=hb.ReadDuke(hb.newmesh);
		hb.HyperFindEpsilon();
		hb.MakeCharges();
		hb.PropagateCharges();
		do{
			hb.SwapMeshes();
			oscarfile=hb.ReadDuke(hb.newmesh);
			hb.HyperFindEpsilon();
		hb.MakeCharges();
		hb.PropagateCharges();
		hb.ScatterCharges();
		if(fabs(lrint(hb.mesh->tau)-hb.mesh->tau)<0.001){
			snprintf(message,CLog::CHARLENGTH,"tau=%g, cmap.size=%lu, emap.size=%lu\n",hb.mesh->tau,
			hb.cmap.size(),hb.emap.size());
			CLog::Info(message);
			snprintf(message,CLog::CHARLENGTH,"highestT=%g, highestEpsilon=%g, biggestU=%g\n",hb.highestT,hb.highestEpsilon,hb.biggestU);
			CLog::Info(message);
		}
	}while(oscarfile);
	
	snprintf(message,CLog::CHARLENGTH,"tau=%g, cmap.size=%lu, emap.size=%lu\n",hb.mesh->tau,
	hb.cmap.size(),hb.emap.size());
	CLog::Info(message);
	snprintf(message,CLog::CHARLENGTH,"highestT=%g, highestEpsilon=%g, biggestU=%g\n",hb.highestT,hb.highestEpsilon,hb.biggestU);
	CLog::Info(message);
	
		hb.WriteCharges();
		if(run_number==0)
			hb.WriteHyper_Duke_2D();
		hb.ClearCharges();
		
	}
	
	return 0;
}


