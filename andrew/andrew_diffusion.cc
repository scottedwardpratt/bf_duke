#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>
#include <vector>
#include <Eigen/Dense>

using namespace std;

int main(){
	
	string filename="latticedata_diffusion/diffusion.txt";
	char dummy[100];
	char voldummy[100];
	int ntaudummy;
	unsigned int ie,it=0,NT;
	double errsysdummy,errsysstatdummy,t,d,Ttarget,T0,T1,dvse,f,epsilon;
	double HBARC_GeV=0.1973269718,PI=4.0*atan(1.0);
	
	vector<double> DvsT,TvsT,DvsE;
	DvsT.push_back(0.0);
	TvsT.push_back(0.0);
	
	FILE *fptr=fopen(filename.c_str(),"r");
	fgets(dummy,100,fptr);
	fscanf(fptr,"%s",voldummy);
	while(!feof(fptr)){
		fscanf(fptr,"%lf %d %lf %lf %lf",&t,&ntaudummy,&d,&errsysdummy,&errsysstatdummy);
		t=t*0.001;
		d=d*HBARC_GeV/(2.0*PI*t);
		TvsT.push_back(t);
		DvsT.push_back(d);
		fscanf(fptr,"%s",voldummy);
	}	
	fclose(fptr);
	NT=TvsT.size();
	DvsT[0]=DvsT[1];
	// add a couple of points
	TvsT.resize(TvsT.size()+2);
	DvsT.resize(DvsT.size()+2);
	TvsT[NT]=2.0*TvsT[NT-1]-TvsT[NT-2];
	DvsT[NT]=2.0*DvsT[NT-1]-DvsT[NT-2];
	NT+=1;
	TvsT[NT]=2.0*TvsT[NT-1]-TvsT[NT-2];
	DvsT[NT]=2.0*DvsT[NT-1]-DvsT[NT-2];
	NT+=1;
	
	fptr=fopen("eosdata/eos_vs_epsilon.txt","r");
	FILE *fptr1=fopen("eosdata/DvsEpsilon.txt","w");
	fprintf(fptr1,"# epsilon      T      D\n");
	fgets(dummy,100,fptr);
	ie=0;
	while(!feof(fptr)){
		fscanf(fptr,"%lf %lf",&epsilon,&t);
		if(!feof(fptr)){
			fgets(dummy,100,fptr);
			Ttarget=t;
			while(TvsT[it]>Ttarget || TvsT[it+1]<Ttarget){
				if(TvsT[it]>Ttarget)
					it-=1;
				else
					it+=1;
			}
		}
		
		T0=TvsT[it];
		T1=TvsT[it+1];
		f=(TvsT[it+1]-Ttarget)/(TvsT[it+1]-TvsT[it]);
		dvse=DvsT[it]*f+(1.0-f)*DvsT[it+1];
		DvsE.push_back(dvse);
		fprintf(fptr1,"%8.4f %15.8e %15.8e\n",epsilon,Ttarget,dvse);
		ie+=1;
	}
	fclose(fptr);
	fclose(fptr1);
	
}	
