#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>
#include <vector>

using namespace std;

int main(){
	const double epsilon_h=0.36,epsilon_qgp=0.72,de=0.005;
	FILE *fptr;
	char dummy[500];
	vector<double> T_lattice,T_hadron,s_lattice,s_hadron;
	vector<double> chiuu_lattice,chiud_lattice,chius_lattice,chiss_lattice;
	vector<double> chiuu_hadron,chiud_hadron,chius_hadron,chiss_hadron;
	
	int ie;
	double epsilon,eread,Tread,sread,chiuuread,chiudread,chiusread,chissread;
	
	fptr=fopen("eosdata/eos_hadron_vs_epsilon.txt","r");
	fgets(dummy,300,fptr);
	ie=0;
	do{
		fscanf(fptr,"%lf %lf %lf %lf %lf %lf %lf",
		&eread,&Tread,&sread,&chiuuread,&chiudread,&chiusread,&chissread);
		fgets(dummy,300,fptr);
		if(!feof(fptr)){
			epsilon=ie*de;
			if(fabs(epsilon-eread)>1.0E-9){
				printf("epsilon doesn't match!!!\n");
				exit(1);
			}
			T_hadron.push_back(Tread);
			s_hadron.push_back(sread);
			chiuu_hadron.push_back(chiuuread);
			chiud_hadron.push_back(chiudread);
			chius_hadron.push_back(chiusread);
			chiss_hadron.push_back(chissread);
			ie+=1;
		}
		//printf("--- %g %g %g\n",eread,Tread,chiuuread);
	}while(!feof(fptr));
	fclose(fptr);
!	
	fptr=fopen("eosdata/eos_lattice_vs_epsilon.txt","r");
	fgets(dummy,300,fptr);
	ie=0;
	do{
		fscanf(fptr,"%lf %lf %lf %lf %lf %lf %lf",
		&eread,&Tread,&sread,&chiuuread,&chiudread,&chiusread,&chissread);
		fgets(dummy,300,fptr);
		if(!feof(fptr)){
			epsilon=ie*de;
			if(fabs(epsilon-eread)>1.0E-9){
				printf("epsilon doesn't match!!!\n");
				exit(1);
			}
			T_lattice.push_back(Tread);
			s_lattice.push_back(sread);
			chiuu_lattice.push_back(chiuuread);
			chiud_lattice.push_back(chiudread);
			chius_lattice.push_back(chiusread);
			chiss_lattice.push_back(chissread);
			ie+=1;
		}
	}while(!feof(fptr));
	fclose(fptr);
	
	double T,s,chiuu,chiud,chius,chiss,f;
	unsigned int NE=T_lattice.size();
	fptr=fopen("eosdata/eos_vs_epsilon.txt","w");
	fprintf(fptr,"#epsilon     T       sdens   chi_uu          chi_ud             chi_us             chi_ss\n");
	for(ie=0;ie<NE;ie++){
		epsilon=ie*de;
		if(epsilon<epsilon_h){
			T=T_hadron[ie];
			s=s_hadron[ie];
			chiuu=chiuu_hadron[ie];
			chiud=chiud_hadron[ie];
			chius=chius_hadron[ie];
			chiss=chiss_hadron[ie];
		}
		else if(epsilon>epsilon_qgp){
			T=T_lattice[ie];
			s=s_lattice[ie];
			chiuu=chiuu_lattice[ie];
			chiud=chiud_lattice[ie];
			chius=chius_lattice[ie];
			chiss=chiss_lattice[ie];
		}
		else{
			f=(epsilon_qgp-epsilon)/(epsilon_qgp-epsilon_h);
			T=f*T_hadron[ie]+(1.0-f)*T_lattice[ie];
			s=f*s_hadron[ie]+(1.0-f)*s_lattice[ie];
			chiuu=f*chiuu_hadron[ie]+(1.0-f)*chiuu_lattice[ie];
			chiud=f*chiud_hadron[ie]+(1.0-f)*chiud_lattice[ie];
			chius=f*chius_hadron[ie]+(1.0-f)*chius_lattice[ie];
			chiss=f*chiss_hadron[ie]+(1.0-f)*chiss_lattice[ie];
		}
		fprintf(fptr,"%15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e\n",
		epsilon,T,s,chiuu,chiud,chius,chiss);
	}
	fclose(fptr);
}
