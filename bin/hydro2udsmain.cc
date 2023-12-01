#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>
#include "hydro2uds.h"
#include "hyper.h"
#include "qualifier.h"

using namespace std;

int main(int argc,char *argv[]){
	if (argc != 2) {
		printf("Usage: hydro2uds run_number\n");
		exit(-1);
  }
	bool oscarfile=true;
	int run_number=atoi(argv[1]);
	string udsfilename="uds"+string(argv[1])+".dat";
	CHydroBalance hb("udsdata/udsparameters.dat",run_number);
	hb.parmap.set("CHARGESINFO_FILENAME",udsfilename);
	CQualifiers qualifiers;
	qualifiers.Read("qualifiers.dat");
	for(int iqual=0;iqual<qualifiers.nqualifiers;iqual++){
		hb.qualifier=qualifiers.qualifier[iqual]->qualname;
		hb.Reset();
		printf("--------- BEGIN CALC FOR %s ---------\n",hb.qualifier.c_str());
		oscarfile=hb.ReadOSCAR(hb.mesh);
		hb.HyperFind();
		oscarfile=hb.ReadOSCAR(hb.newmesh);
		hb.HyperFind();
		hb.MakeCharges();
		hb.PropagateCharges();
		do{
			hb.SwapMeshes();
			oscarfile=hb.ReadOSCAR(hb.newmesh);
			hb.HyperFind();
			hb.MakeCharges();
			hb.PropagateCharges();
			hb.ScatterCharges();
			if(fabs(lrint(hb.mesh->tau)-hb.mesh->tau)<0.001)
				printf("tau=%g, cmap.size=%d, emap.size=%d\n",hb.mesh->tau,int(hb.cmap.size()),int(hb.emap.size()));
		}while(oscarfile);
		hb.WriteCharges();
		if(run_number==0)
			hb.WriteHyper();
		hb.ClearCharges();
	}
	
	return 0;
}


