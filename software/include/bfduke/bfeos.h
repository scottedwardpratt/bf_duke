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
		void Print();
		void PrintChi();
	
		//private:
		static int NE;
		static double depsilon;
		static vector<double> TvsE,svsE;
		static vector<double> chillvsE,chiudvsEchilsvsE,chissvsE;
		static vector<double> twopiTD,Tdiff;
		static vector<double> twopiTD,Tdiff
		
		
		static vector<double> epsilon_PST,P_PST,s_PST,T_PST;
		static vector<double> epsilon_claudia,P_claudia,s_claudia,T_claudia;
		static vector<double> twopiTD,Tdiff;
		static vector<double> chill_overs_claudia,chiud_overs_claudia,chils_overs_claudia,chiss_overs_claudia;
		static vector<double> chill_HSC,chiud_HSC,chils_HSC,chiss_HSC;
		static vector<double> dDdT;
		static mapdi etmap;
		static CresList *reslist;

	};

}


#endif
