#ifndef __BFEoS_H__
#define __BFEoS_H__

#include "msu_commonutils/commondefs.h"
#include "bfduke/bfcommon.h"

using namespace std;
using namespace NMSUPratt;

typedef multimap<double,int> mapdi;
typedef pair<double,int> pairdi;

namespace NMSUPratt{

	class CHBEoS{
	public:
		CparameterMap *parmap;
		double T,P,epsilon,s;
		double chill,chiud,chils,chiss;
		double f_l,f_u,f_d,f_s,gamma_q;
		CHBEoS(CparameterMap *parmapset);
		CHBEoS();
		void ReadDiffusionData();
		static double GetD(double epsilon);

		void ReadEoSData_Andrew();
		void ReadChiReductionFactors();
		void Print();
		void PrintChi();
	
		//private:
		static int NE;
		static double depsilon,epsilon_h,epsilon_qgp;
		static vector<double> TvsE,svsE;
		static vector<double> chillvsE,chiudvsE,chilsvsE,chissvsE;
		static vector<double> DvsE;
		
		static int Nchifactors;
		static vector<double> chifactorll,chifactorud,chifactorus,chifactorss;


	};

}


#endif
