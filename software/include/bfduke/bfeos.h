#ifndef __BFEoS_H__
#define __BFEoS_H__

#include "msu_commonutils/commondefs.h"
#include "bfduke/bfcommon.h"

using namespace std;
using namespace NMSUPratt;

typedef pair<double,int> pairdi;

namespace NMSUPratt{

	class CHBEoS{
	public:
		double T,P,epsilon,s; // T,P,s are all quantities at equilibrium given epsilon (i.e. with fugacities=1)
		double Tnonequil;
		double chill,chiud,chils,chiss;
		double f_l,f_u,f_d,f_s,gamma_q;
		CHBEoS();
		
		void Print();
		void PrintChi();
		void SetChi();
		void SetTs();
		void SetTnonequil();
	
		static void ReadChiReductionFactors();
		static void ReadEoSData_Andrew();
		static void ReadDiffusionData_Andrew();
		static CparameterMap *parmap;
		static int NE;
		static double depsilon,epsilon_h,epsilon_qgp;
		static vector<double> TvsE,svsE;
		static vector<double> chillvsE,chiudvsE,chilsvsE,chissvsE;
		static vector<double> DvsE;
		static double GetD(double epsilon);
		
		static int Nchifactors;
		static vector<double> chifactorll,chifactorud,chifactorls,chifactorss,TnonequilVec;


	};

}


#endif
