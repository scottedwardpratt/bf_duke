#include "msu_boltzmann/msu_boltzmann.h"
#include "msu_eos/resonances.h"
#include "msu_sampler/sampler.h"
#include "msu_commonutils/log.h"
#include "msu_boltzmann/balancearrays.h"
#include "msu_commonutils/qualifier.h"
#include "msu_commonutils/randy.h"
#include "msu_commonutils/misc.h"

using namespace std;
using namespace NMSUPratt;

int main(int argc, char *argv[]){
	if (argc != 4) {
		CLog::Fatal("Usage: b3d run_number ievent0 ieventf\n");
  }
	CparameterMap parmap;
	parmap.ReadParsFromFile("model_output/fixed_parameters.txt");
	CresList reslist(&parmap);
	
	CBalanceArrays *barray;
	long long int npartstot,nparts0;
	long long int norm;
	int ievent,iqual,nevents,run_number;
	run_number=atoi(argv[1]);
	int ievent0=atoi(argv[2]),ieventf=atoi(argv[3]);
	nevents=1+ieventf-ievent0;
	CMSU_Boltzmann *b3d=new CMSU_Boltzmann(run_number,&reslist);
	CmasterSampler *ms=new CmasterSampler(&parmap);
	ms->ClearHyperList();
	
	
	b3d->InitCascade();
	barray=b3d->balancearrays;

	CQualifiers qualifiers;
	qualifiers.Read("qualifiers.txt");
	
	for(iqual=0;iqual<qualifiers.nqualifiers;iqual++){
		npartstot=0;
		b3d->SetQualifier(qualifiers.qualifier[iqual]->qualname);
		qualifiers.SetPars(&(b3d->parmap),iqual);
		parmap.set("HYPER_INFO_FILE","model_output/run"+to_string(run_number)+"/"+qualifiers.qualifier[iqual]->qualname+"/udsdata/hyper.txt");
		
		CLog::Info("_________________ iqual="+to_string(iqual)+", nevents="+to_string(nevents)+"\n");
		ms->ReadHyper_OSU_2D();
		for(ievent=ievent0;ievent<=ieventf;ievent++){
			CLog::Info("------ beginning, ievent="+to_string(ievent)+" -------\n");
			
			ms->randy->reset(ievent);
			ms->partlist->Clear();
			printf("check a\n");
			nparts0=ms->MakeEvent();
			printf("check nparts0=%lld\n",nparts0);
			
			b3d->Reset();
			b3d->randy->reset(ievent);
			b3d->InputPartList(ms->partlist);
			
			b3d->PerformAllActions();
			CLog::Info("N initial parts = "+to_string(nparts0)+", N final parts = "+to_string(b3d->PartMap.size())+"\n");
			npartstot+=b3d->PartMap.size();
			barray->ProcessPartMap();
			
		}
		norm=nevents*b3d->NSAMPLE;
		CLog::Info("<Nparts>="+to_string(double(npartstot)/norm)+"\n");
		barray->ConstructBFs();
		barray->WriteBFs();
		barray->WriteDenoms();
		barray->WriteGammaP();
	}
	
	return 0;
}
