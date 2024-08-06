#ifndef __HYDROBALANCE_H__
#define __HYDROBALANCE_H__

//#include <boost/math/special_functions/bessel.hpp>
#include "bfduke/bfcommon.h"
#include "bfduke/bfeos.h"
#include "msu_sampler/hyper.h"
#include "msu_commonutils/commondefs.h"
#include "msu_commonutils/parametermap.h"

using namespace std;
using namespace NMSUPratt;

namespace NMSUPratt{

	class CHydroBalance{
	public:
		double Omega0tot,OmegaXtot,OmegaYtot,netUdotOmega;
		CHBEoS *eos;
		bool HYPERT;  // Use T to define hypersurface
		bool HYPEREPSILON; // Use epsilon to define hypersurface
		double DELTAU,TAU0,XMIN,XMAX,YMIN,YMAX,DX,DY;
		int NX,NY,ntraj;
		double NSAMPLE_HYDRO2UDS;
		double biggestomega;
		double highestEpsilon,highestT,biggestU;
		vector<CHBCharge> hbcharges;
		unsigned int nhbcharges;
		int run_number;
		string qualifier;
		double SIGMA0;  //initial movement of one particle at birth
		Crandy *randy;
		CparameterMap parmap;
		double ransum,ranthresh;
		int Ncollisions,idmax,tau0check;
		double GAMMA_0,TAU_EQ;
		double Tf,epsilon_f;
		double DiffusionRatio;
		Eigen::Matrix3d chif,chifinv;
		CHydroBalance();
		CHydroBalance(int run_number);
		~CHydroBalance();
		void SetPars();
		void MakeMeshes();
		void SwapMeshes();
		void MakeCharges();
		void HyperFindT();
		void HyperFindEpsilon();
		void GetUxyBar(int ix,int iy,double &uxbar,double &uybar); // ux,uy for cell ix,iy between mesh and newmesh
		void GetXYBar(int ix,int iy,double &xbar,double &ybar);
		void GetTBar(int ix,int iy,double &Tbar); // ux,uy for cell ix,iy between mesh and newmesh
		void GetEpsilonBar(int ix,int iy,double &EpsilonBar); // ux,uy for cell ix,iy between mesh and newmesh
		void GetPiTildeBar(int ix,int iy,double &pitildexxbar,double &pitildexybar,
		double &pitildeyybar);
		bool GetGradT(int ix,int iy,double &dTdt,double &dTdx,double &dTdy,
		bool &GGTt,bool &GGTx,bool &GGTy); // GradT for same cell
		bool GetGradEpsilon(int ix,int iy,double &dEdt,double &dEdx,double &dEdy,
		bool &GGEt,bool &GGEx,bool &GGEy);
		// returns true if cell has hypersurface, if(forcecalc=false) only calculates if on hypersurface
		bool GetDOmega(double dTdt,double dTdx,double dTdy,double &dOmega0,double &dOmegaX,double &dOmegaY,bool GGTt,bool GGTx,bool GGTy); //Hyper surface element for times between mesh and newmesh
		bool WRITE_TRAJ;
		void PropagateCharges();
		void ScatterCharges();
		void CalcDQ(int ix,int iy,double &DQll,double &DQud,double &DQls,double &DQss);
		void CalcDQ0(int ix,int iy,double &DQll,double &DQud,double &DQls,double &DQss);
		void GetGammaFQ(double tau,double &gamma_q,double &fugacity_l,double &fugacity_s);
		CHBHydroMesh *oldmesh;
		CHBHydroMesh *mesh;
		CHBHydroMesh *newmesh;
		mapic cmap; // active particles
		mapic emap; // emitted particles
		//list<CHBHyperElement *> hyperlist;   // list of hyper-elements
		list<Chyper *> hyperlist;   // list of hyper-elements
		Eigen::Matrix3d chitothyper;
		void WriteCharges(int ichargefile);
		void ClearCharges();
		void WriteSource();
		void WriteHyper();
		void WriteHyper_Duke_2D();
		void WriteFinalCF();
		void Reset();
		bool ReadOSCAR(CHBHydroMesh *hydromesh);
		bool ReadDuke(CHBHydroMesh *hydromesh);
		bool FakeReadOSCAR(CHBHydroMesh *hydromesh);
		double SpectraFromHyper(double mass,double px,double py);
		FILE *fptr_oscar,*fptr_duke;
		string oscar_filename,duke_filename;
		double tau0readcheck;
		int itauread;
		vector<Eigen::Matrix3d> source;
	};

	// Info for one tau
	class CHBHydroMesh{
	public:
		CHBHydroMesh();
		double **T,**DQ,**UX,**UY,**epsilon;
		double **pitildexx,**pitildexy,**pitildeyy;
		double tau;
		int itau;
		void GetXY(int ix,int iy,double &x,double &y);
		void GetIXIY(double x,double y,int &ix,int &iy);
		void GetIXIY_lower(double x,double y,int &ix,int &iy);
		static double GAMMA_0,TAU_EQ;
		static CHydroBalance *hb;
		static int NX,NY;
		static double DX,DY,DELTAU,TAU0;
		static double XMIN,XMAX,YMIN,YMAX;
		static void GetDimensions(int &NXs,int &NYs,double &DXs,double &DYs,double &DELTAUs,double &TAU0s,double &XMINs,double &XMAXs,double &YMINs,double &YMAXs){
			NXs=NX; NYs=NY;
			DXs=DX; DYs=DY; DELTAUs=DELTAU; TAU0s=TAU0;
			XMINs=XMIN; XMAXs=XMAX; YMINs=YMIN; YMAXs=YMAX;
		}
	};

}

#endif
