#include "msu_sampler/master.h"
#include "msu_commonutils/qualifier.h"
#include "msu_boltzmann/msu_boltzmann.h"
#include "msu_commonutils/log.h"
#include <cstring>
using namespace std;
using namespace NMSUPratt;

int main(int argc, char *argv[]){
	if (argc != 3) {
		CLog::Fatal("Usage: boltz run_number subrun_number\n");
  }
	CparameterMap parmap;
	int run_number=atoi(argv[1]),subrun_number=atoi(argv[2]);
	int subrun_number_max=10000;
	int nevents;
	char message[CLog::CHARLENGTH];
	long long int nmerge,nscatter,nannihilate,ncancel_annihilate,nparts,npartstot,ievent,ndecay;
	string logfilename="logfiles/run"+to_string(run_number)+"_subrun"+to_string(subrun_number)+".txt";
	CLog::Init(logfilename);
	parmap.ReadParsFromFile("modelruns/fixed_parameters.txt");
	nevents=parmap.getI("MSU_BOLTZMANN_NEVENTS_TYPE2",1);
	CmasterSampler ms(&parmap);
	CMSU_Boltzmann::mastersampler=&ms;
	CpartList *pl=new CpartList(&parmap,ms.reslist);
	ms.partlist=pl;
	CMSU_Boltzmann *msuboltz=new CMSU_Boltzmann(run_number,subrun_number,ms.reslist);
	msuboltz->InitCascade();
	CBalanceArrays *barray;
	if(msuboltz->BFCALC){
		barray=msuboltz->balancearrays;
		barray->FROM_UDS=false;
	}
	//msuboltz->ReadMuTInfo();
	CQualifiers qualifiers;
	int iqual;
	qualifiers.Read("qualifiers.txt");
	for(iqual=0;iqual<qualifiers.nqualifiers;iqual++){
		npartstot=0;
		nmerge=nscatter=nannihilate=ncancel_annihilate=ndecay=0;
		msuboltz->SetQualifier(qualifiers.qualifier[iqual]->qualname);
		qualifiers.SetPars(&(msuboltz->parmap),iqual);
		ms.ReadHyper_Duke_2D(run_number,qualifiers.qualifier[iqual]->qualname);

		for(ievent=0;ievent<nevents;ievent++){
			CLog::Info("--- begin for ievent="+to_string(ievent)+"\n");
			int idummy=1;
			ms.randy->reset(idummy*nevents*1000+nevents*subrun_number_max*run_number+nevents*subrun_number+ievent);
			if(subrun_number>subrun_number_max){
				CLog::Fatal("OH NO!!! subrun_number>subrun_number_max. Increase this in boltzmain.cc");
			}
			msuboltz->Reset();
			nparts=ms.MakeEvent();
			CLog::Info("created "+to_string(nparts)+" initial parts\n");
			npartstot+=nparts;
			msuboltz->InputPartList(pl);
			pl->Reset();
			CLog::Info("---- begin PerformAllActions\n");
			msuboltz->PerformAllActions();
			CLog::Info("---- actions performed\n");
			msuboltz->IncrementSpectraV2();
			CLog::Info("Nparts final="+to_string(msuboltz->PartMap.size())+"\n");
		
			nmerge+=msuboltz->nmerge;
			nscatter+=msuboltz->nscatter;
			nannihilate+=msuboltz->nannihilate;
			ncancel_annihilate+=msuboltz->ncancel_annihilate;
			ndecay+=msuboltz->ndecay;
			snprintf(message,CLog::CHARLENGTH,"---- ievent=%lld nparts=%lld, nparts/event=%g\n",ievent,nparts,double(npartstot)/double(ievent+1));
			CLog::Info(message);
			if(msuboltz->BFCALC){
				barray->ProcessPartMap();
				CLog::Info("----- partmap processed\n");
			}
			msuboltz->KillAllParts();
			msuboltz->CheckPartMap();
			msuboltz->CheckDeadPartMap();
		}
		snprintf(message,CLog::CHARLENGTH,"ndecay/event=%g, nmerge/event=%g, nscatter/event=%g\n",
		double(ndecay)/double(ievent+1),double(nmerge)/double(ievent+1),
		double(nscatter)/double(ievent+1));
		CLog::Info(message);
		snprintf(message,CLog::CHARLENGTH,"nannihilate/event=%g, ncancel_annihilate/event=%g\n",
		double(nannihilate)/double(ievent+1),double(ncancel_annihilate)/double(ievent+1));
		CLog::Info(message);
		//msuboltz->WriteMuTInfo();
		msuboltz->WriteSpectraV2();
		if(msuboltz->BFCALC){
			barray->ConstructBFs();
			barray->WriteBFs();
			barray->WriteDenoms();
			//barray->WriteGammaP();
		}
		ms.DeleteHyperElements();
	}

	CLog::Info("YIPPEE!!!!! We made it all the way through!\n");
	return 0;
}
