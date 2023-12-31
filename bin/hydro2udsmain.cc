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
	bool oscarfile=true;
	int run_number=atoi(argv[1]);
	string udsfilename="uds"+string(argv[1])+".txt";
	CHydroBalance hb("udsdata/udsparameters.txt",run_number);
	hb.parmap.set("CHARGESINFO_FILENAME",udsfilename);
	hb.parmap.set("CHARGESINFO_FILENAME",udsfilename);
	CQualifiers qualifiers;
	qualifiers.Read("qualifiers.txt");
	for(int iqual=0;iqual<qualifiers.nqualifiers;iqual++){
		hb.qualifier=qualifiers.qualifier[iqual]->qualname;
		hb.Reset();
		printf("--------- BEGIN CALC FOR %s ---------\n",hb.qualifier.c_str());
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
			if(fabs(lrint(hb.mesh->tau)-hb.mesh->tau)<0.001)
				CLog::Info("tau="+to_string(hb.mesh->tau)+", cmap.size="+to_string(hb.cmap.size())+", emap.size="+to_string(hb.emap.size())+"\n");
		}while(oscarfile);
		hb.WriteCharges();
		if(run_number==0)
			hb.WriteHyper();
		hb.ClearCharges();
		
	}
	
	return 0;
}


