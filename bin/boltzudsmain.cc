#include "msu_sampler/master.h"
#include "msu_commonutils/qualifier.h"
#include "msu_boltzmann/msu_boltzmann.h"
#include "msu_commonutils/log.h"
#include "msu_commonutils/sf.h"
#include <cstring>
using namespace std;
using namespace NMSUPratt;

int main(int argc, char *argv[]){
	if (argc != 3) {
		CLog::Info("Usage: msuboltz run_number subrun_number\n");
		exit(-1);
  }
	CparameterMap parmap;
	int run_number=atoi(argv[1]),subrun_number=atoi(argv[2]);
	int subrun_number_max=10000;
	//int ievent0=atoi(argv[2]),ieventf=atoi(argv[3]);
	char message[CLog::CHARLENGTH];
	long long int nmerge,nscatter,nannihilate,ncancel_annihilate,nevents,nparts,npartstot,ievent,ndecay;
	//string logfilename="logfiles/run"+to_string(run_number)+"_subrun"+to_string(subrun_number)+".txt";
	//CLog::Init(logfilename);
	parmap.ReadParsFromFile("modelruns/fixed_parameters.txt");
	
	
	CmasterSampler *ms=new CmasterSampler(&parmap);
	CMSU_Boltzmann::mastersampler=ms;
	//CpartList *pl=new CpartList(&parmap,ms->reslist);
	//CpartList *pl=ms->partlist;
	//pl->Clear();
	//ms->partlist=pl;
	
	//CmasterSampler *ms0=new CmasterSampler(&parmap);
	//CmasterSampler *ms;

	CMSU_Boltzmann *msuboltz=new CMSU_Boltzmann(run_number,subrun_number,ms->reslist);
	
	msuboltz->InitCascade();
	CBalanceArrays *barray=msuboltz->balancearrays;
	barray->FROM_UDS=true;
	nevents=msuboltz->parmap.getI("MSU_BOLTZMANN_NEVENTS_TYPE1",10);
	msuboltz->nevents=0;
	CQualifiers qualifiers;
	int iqual;
	qualifiers.Read("qualifiers.txt");
	for(iqual=0;iqual<qualifiers.nqualifiers;iqual++){
		npartstot=0;
		nmerge=nscatter=nannihilate=ncancel_annihilate=ndecay=0;
		msuboltz->SetQualifier(qualifiers.qualifier[iqual]->qualname);
		qualifiers.SetPars(&(msuboltz->parmap),iqual);
		ms->ReadHyper_Duke_2D(run_number,qualifiers.qualifier[iqual]->qualname);
		for(ievent=0;ievent<nevents;ievent++){
			
			CLog::Info("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
			CLog::Info("--- begin for ievent="+to_string(ievent)+" ---\n");
			CLog::Info("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
			ms->randy->reset(nevents*subrun_number_max*run_number+nevents*subrun_number+ievent);
			if(subrun_number>subrun_number_max){
				CLog::Fatal("OH NO!!! subrun_number>subrun_number_max. Increase this in boltzmain.cc");
			}
			msuboltz->Reset();

			nparts=ms->MakeEvent();
			npartstot+=nparts;
			msuboltz->InputPartList(ms->partlist);
			ms->partlist->Reset();		
		

			if(msuboltz->BFCALC && barray->FROM_UDS){
				msuboltz->ReadCharges(ievent);
				msuboltz->GenHadronsFromCharges(); // Generates inter-correlated parts, with bids = (0,1),(2,3)....
			}

			CLog::Info("begin PerformAllActions\n");
			msuboltz->PerformAllActions();
			CLog::Info("All Actions Performed\n");
			msuboltz->IncrementHadronCount();
		
			nmerge+=msuboltz->nmerge;
			nscatter+=msuboltz->nscatter;
			nannihilate+=msuboltz->nannihilate;
			ncancel_annihilate+=msuboltz->ncancel_annihilate;
			ndecay+=msuboltz->ndecay;
			snprintf(message,CLog::CHARLENGTH,"nevents=%lld <nparts>=%lld, nparts/event=%g\n",ms->NEVENTS,nparts,double(npartstot)/double(ms->NEVENTS));
			CLog::Info(message);
			barray->ProcessPartMap();
			if(msuboltz->BFCALC && barray->FROM_UDS){
				CLog::Info("XXXXX processing BF PartMap\n");
				barray->ProcessBFPartMap();
			}
			snprintf(message,CLog::CHARLENGTH,"Npartstot=?%lu, Nactionstot=?%lu\n",msuboltz->PartMap.size()+msuboltz->DeadPartMap.size(),
			msuboltz->ActionMap.size()+msuboltz->DeadActionMap.size());
			CLog::Info(message);
		}
		snprintf(message,CLog::CHARLENGTH,"ndecay/event=%g, nmerge/event=%g, nscatter/event=%g\n",
		double(ndecay)/double(nevents),double(nmerge)/double(nevents),double(nscatter)/double(nevents));
		CLog::Info(message);
		snprintf(message,CLog::CHARLENGTH,"nannihilate/event=%g, ncancel_annihilate/event=%g\n",
		double(nannihilate)/double(nevents),double(ncancel_annihilate)/double(nevents));
		CLog::Info(message);
		//msuboltz->WriteMuTInfo();
		msuboltz->WriteHadronCount();
		barray->ConstructBFs();
		barray->WriteBFs();
		barray->WriteDenoms();
		barray->WriteGammaP();
		ms->DeleteHyperElements();
	}

	CLog::Info("YIPPEE!!!!! We made it all the way through!\n");
	return 0;
}
