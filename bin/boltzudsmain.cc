#include "msu_sampler/master.h"
#include "msu_commonutils/qualifier.h"
#include "msu_boltzmann/msu_boltzmann.h"
#include "msu_commonutils/log.h"
#include <cstring>
using namespace std;
using namespace NMSUPratt;

int main(int argc, char *argv[]){
	if (argc != 2) {
		CLog::Info("Usage: msuboltz run_number\n");
		exit(-1);
  }
	CparameterMap parmap;
	int run_number=atoi(argv[1]),subrun_number=-1;
	//int ievent0=atoi(argv[2]),ieventf=atoi(argv[3]);
	char message[CLog::CHARLENGTH];
	long long int nmerge,nscatter,nannihilate,ncancel_annihilate,nevents,nparts,npartstot,ievent,ndecay;
	//char logfilename[100];
	//sprintf(logfilename,"msuboltz_log.txt");
	//CLog::Init(logfilename);
	CLog::INTERACTIVE=true;
	printf("check aaa\n");
	parmap.ReadParsFromFile("modelruns/fixed_parameters.txt");
	CmasterSampler ms(&parmap);
	CMSU_Boltzmann::mastersampler=&ms;
	CpartList *pl=new CpartList(&parmap,ms.reslist);
	ms.partlist=pl;

printf("check aa\n");
	CMSU_Boltzmann *msuboltz=new CMSU_Boltzmann(run_number,subrun_number,ms.reslist);
	printf("check a\n");
	msuboltz->InitCascade();
	printf("check b\n");
	CBalanceArrays *barray=msuboltz->balancearrays;
	barray->FROM_UDS=true;
	nevents=msuboltz->parmap.getI("MSU_BOLTZMANN_NEVENTSMAX",10);
	
	//msuboltz->ReadMuTInfo();
	msuboltz->nevents=0;
	CQualifiers qualifiers;
	int iqual;
	qualifiers.Read("qualifiers.txt");
	printf("check c\n");
	for(iqual=0;iqual<qualifiers.nqualifiers;iqual++){
		npartstot=0;
		nmerge=nscatter=nannihilate=ncancel_annihilate=ndecay=0;
		msuboltz->SetQualifier(qualifiers.qualifier[iqual]->qualname);
		qualifiers.SetPars(&(msuboltz->parmap),iqual);
		ms.ReadHyper_Duke_2D(run_number,qualifiers.qualifier[iqual]->qualname);
		for(ievent=0;ievent<nevents;ievent++){
			printf("--- begin for ievent=%lld\n",ievent);
			ms.randy->reset(ievent);
			msuboltz->Reset();
			nparts=ms.MakeEvent();
			npartstot+=nparts;
			msuboltz->InputPartList(pl);
			pl->Clear();
		
			if(msuboltz->BFCALC && barray->FROM_UDS){
				msuboltz->ReadCharges(ievent);
				msuboltz->GenHadronsFromCharges(); // Generates inter-correlated parts, with bids = (0,1),(2,3)....
				msuboltz->DeleteCharges();
			}
		
			msuboltz->PerformAllActions();
			msuboltz->IncrementHadronCount();
		
			nmerge+=msuboltz->nmerge;
			nscatter+=msuboltz->nscatter;
			nannihilate+=msuboltz->nannihilate;
			ncancel_annihilate+=msuboltz->ncancel_annihilate;
			ndecay+=msuboltz->ndecay;
			snprintf(message,CLog::CHARLENGTH,"ievent=%lld nparts=%lld, nparts/event=%g\n",ms.NEVENTS,nparts,double(npartstot)/double(ms.NEVENTS));
			CLog::Info(message);
			barray->ProcessPartMap();
			if(msuboltz->BFCALC && barray->FROM_UDS){
				printf("xXXXX processing BF PartMap\n");
				barray->ProcessBFPartMap();
			}
			msuboltz->KillAllParts();
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
	}

	CLog::Info("YIPPEE!!!!! We made it all the way through!\n");
	return 0;
}
