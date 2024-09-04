#include "bfduke/bfcommon.h"
#include "bfduke/hydro2uds.h"
#include "msu_sampler/hyper.h"


using namespace std;
using namespace NMSUPratt;

int CHBHydroMesh::NX=0;
int CHBHydroMesh::NY=0;
double CHBHydroMesh::DX=0.0;
double CHBHydroMesh::DY=0.0;
double CHBHydroMesh::XMIN=0.0;
double CHBHydroMesh::XMAX=0.0;
double CHBHydroMesh::YMIN=0.0;
double CHBHydroMesh::YMAX=0.0;
double CHBHydroMesh::DELTAU=0.0;
double CHBHydroMesh::TAU0=0.0;
double CHBHydroMesh::TAU_EQ=0.0;
double CHBHydroMesh::GAMMA_0=0.0;

CHBHydroMesh::CHBHydroMesh(){
	int ix,iy;
	T=new double*[NX];
	epsilon=new double *[NX];
	DQ=new double*[NX];
	UX=new double *[NX];
	UY=new double *[NX];
	pitildexx=new double *[NX];
	pitildexy=new double *[NX];
	pitildeyy=new double *[NX];
	for(ix=0;ix<NX;ix++){
		T[ix]=new double[NY];
		epsilon[ix]=new double[NY];
		DQ[ix]=new double[NY];
		UX[ix]=new double[NY];
		UY[ix]=new double[NY];
		pitildexx[ix]=new double[NY];
		pitildexy[ix]=new double[NY];
		pitildeyy[ix]=new double[NY];
		for(iy=0;iy<NY;iy++){
			T[ix][iy]=DQ[ix][iy]=UX[ix][iy]=UY[ix][iy]
				=pitildexx[ix][iy]=pitildexy[ix][iy]=pitildeyy[ix][iy]=0.0;
		}
	}
	tau=TAU0;
}

void CHBHydroMesh::GetXY(int ix,int iy,double &x,double &y){
	x=XMIN+ix*DX;
	y=YMIN+iy*DY;
}

void CHBHydroMesh::GetIXIY(double x,double y,int &ix,int &iy){
	ix=lrint((x-XMIN)/DX);
	iy=lrint((y-YMIN)/DY);
}

void CHBHydroMesh::GetIXIY_lower(double x,double y,int &ix,int &iy){
	ix=floorl((x-XMIN)/DX);
	iy=floorl((y-YMIN)/DY);
}