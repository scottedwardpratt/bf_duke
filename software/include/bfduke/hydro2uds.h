#ifndef __HYDROBALANCE_H__
#define __HYDROBALANCE_H__

//#include <boost/math/special_functions/bessel.hpp>
#include "bfduke/bfcommon.h"
#include "bfduke/bfeos.h"
#include "msu_commonutils/commondefs.h"
#include "msu_commonutils/parametermap.h"

using namespace std;
using namespace NMSUPratt;

class CHydroBalance{
public:
	CHBEoS *eos;
	double DELTAU,TAU0,XMIN,XMAX,YMIN,YMAX,DX,DY;
	int NX,NY,ntraj;
	double NSAMPLE_HYDRO2UDS;
	double biggestomega;
	string qualifier;
	double SIGMA0;  //initial movement of one particle at birth
	Crandy *randy;
	CparameterMap parmap;
	double ransum,ranthresh;
	int Ncollisions,idmax,tau0check;
	double Tf;
	double DiffusionRatio;
	Eigen::Matrix3d chif,chifinv;
	CHydroBalance();
	CHydroBalance(string parfilename,int ranseed);
	~CHydroBalance();
	void SetPars();
	void MakeMeshes();
	void SwapMeshes();
	void MakeCharges();
	void HyperFind();
	void GetUxyBar(int ix,int iy,double &uxbar,double &uybar); // ux,uy for cell ix,iy between mesh and newmesh
	void GetXYBar(int ix,int iy,double &xbar,double &ybar);
	void GetTBar(int ix,int iy,double &Tbar); // ux,uy for cell ix,iy between mesh and newmesh
	void GetPiTildeBar(int ix,int iy,double &pitildexxbar,double &pitildexybar,
	double &pitildeyybar);
	bool GetGradT(int ix,int iy,double &dTdt,double &dTdx,double &dTdy,
	bool &GGTt,bool &GGTx,bool &GGTy); // GradT for same cell
	// returns true if cell has hypersurface, if(forcecalc=false) only calculates if on hypersurface
	bool GetDOmega(double dTdt,double dTdx,double dTdy,double &dOmega0,double &dOmegaX,double &dOmegaY,bool GGTt,bool GGTx,bool GGTy); //Hyper surface element for times between mesh and newmesh
	bool WRITE_TRAJ;
	void PropagateCharges();
	void ScatterCharges();
	void CalcDQ(int ix,int iy,double &DQll,double &DQud,double &DQls,double &DQss);
	void CalcDQ0(int ix,int iy,double &DQll,double &DQud,double &DQls,double &DQss);
	CHBHydroMesh *oldmesh;
	CHBHydroMesh *mesh;
	CHBHydroMesh *newmesh;
	mapic cmap; // active particles
	mapic emap; // emitted particles
	list<CHBHyperElement *> hyperlist;   // list of hyper-elements
	Eigen::Matrix3d chitot,chitothyper;
	void WriteCharges();
	void ClearCharges();
	void WriteSource();
	void WriteHyper();
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
	double **T,**DQ,**UX,**UY;
	double **pitildexx,**pitildexy,**pitildeyy;
	double tau;
	int itau;
	void GetXY(int ix,int iy,double &x,double &y);
	void GetIXIY(double x,double y,int &ix,int &iy);
	void GetIXIY_lower(double x,double y,int &ix,int &iy);
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

#endif
