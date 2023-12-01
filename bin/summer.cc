#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include <string>
#include <cstring>

const double PI=4.0*atan(1.0);
const double HBARC=197.3269602;

using namespace std;

int main(int argc,char *argv[]){
	string qualifier="alice_cent0_10";
	string acceptance="results_cheap";
	string dirprefix="default";
	char qs[200];
	const int NTYPES=10;
	int Nqcut=4;
	int NJDIR;
	FILE *fptr_read;
	FILE *fptr_write;
	string oldfilename,newfilename,newdirname;
	vector<double> ysum[4],xsum;
	double x,y[4],bfnorm,dx;
	int jdir,itype,ifn,ix;
	
	printf("Enter qualifier: ");
	scanf("%s",qs);
	printf("Enter NJDIR: (the number of default_xxx dirs) ");
	scanf("%d",&NJDIR);
	
	qualifier=qs;
	string olddatafilename[6]={"bf_phi.dat","bf_eta.dat","bf_y.dat","bf_eta1.dat","bf_y1.dat","bf_qinv.dat"};
	string newdatafilename[6]={"bf_phi.dat","bf_eta.dat","bf_y.dat","bf_eta1.dat","bf_y1.dat","bf_qinv.dat"};
	string dirname[NTYPES]={"KK","Kp","piK","pip","pipi","pp","allcharges","allcharges_phi0","allcharges_phi45","allcharges_phi90"};
	string command;
	for(itype=0;itype<NTYPES;itype++){
		command="mkdir -p "+dirprefix+"_sum/"+qualifier+"/"+acceptance+"/"+dirname[itype];
		system(command.c_str());
	}
	
	// Read denom sum
	double denom_pi=0,denom_K=0,denom_p=0,denom_allcharges=0;
	double denom_allcharges_phi0=0,denom_allcharges_phi45=0,denom_allcharges_phi90=0;
	char dummy[100];
	string gpstring,dumbo;
	double nplus,nminus;
	printf("NJDIR=%d\n",NJDIR);
	for(jdir=0;jdir<NJDIR;jdir++){
		oldfilename=dirprefix+"_"+to_string(jdir)+"/"+qualifier+"/"+acceptance+"/denom.dat";
		printf("oldfilename=%s\n",oldfilename.c_str());
		fptr_read=fopen(oldfilename.c_str(),"r");
		fscanf(fptr_read,"%s %lf %lf",dummy,&nplus,&nminus);
		denom_pi+=nplus+nminus;
		fscanf(fptr_read,"%s %lf %lf",dummy,&nplus,&nminus);
		denom_K+=nplus+nminus;
		fscanf(fptr_read,"%s %lf %lf",dummy,&nplus,&nminus);
		denom_p+=nplus+nminus;
		fscanf(fptr_read,"%s %lf %lf",dummy,&nplus,&nminus);
		denom_allcharges+=nplus+nminus;
		fscanf(fptr_read,"%s %lf %lf",dummy,&nplus,&nminus);
		denom_allcharges_phi0+=nplus+nminus;
		fscanf(fptr_read,"%s %lf %lf",dummy,&nplus,&nminus);
		denom_allcharges_phi45+=nplus+nminus;
		fscanf(fptr_read,"%s %lf %lf",dummy,&nplus,&nminus);
		denom_allcharges_phi90+=nplus+nminus;
		fclose(fptr_read);
		printf("finished jdir=%d\n",jdir);
	}
	newfilename="default_sum/"+qualifier+"/"+acceptance+"/denom.dat";
	fptr_write=fopen(newfilename.c_str(),"w");
	fprintf(fptr_write,"denom_pi: %g\ndenom_K: %g\ndenom_p: %g\ndenom_allcharges: %g\ndenom_allcharges_phi0: %g\ndenom_allcharges_phi45: %g\ndenom_allcharges_phi90 %g\n",
	denom_pi,denom_K,denom_p,denom_allcharges,denom_allcharges_phi0,denom_allcharges_phi45,denom_allcharges_phi90);
	fclose(fptr_write);
	
	
	// Now average balance funtions
	for(itype=0;itype<NTYPES;itype++){
		if(itype<6)
			Nqcut=6;
		else
			Nqcut=5;
		for(ifn=0;ifn<Nqcut;ifn++){
			dx=0.0;
			if(ifn==0)
				dx=360.0/28.0;
			if(ifn==1)
				dx=0.1;
			if(ifn==2)
				dx=0.1;
			if(ifn==3)
				dx=0.1;
			if(ifn==4)
				dx=0.1;
			if(ifn==5)
				dx=10.0;
			newdirname="default_sum/"+qualifier+"/"+acceptance+"/"+dirname[itype];
			command="mkdir -p "+newdirname;
			system(command.c_str());
			newfilename=newdirname+"/"+newdatafilename[ifn];
			printf("writing to %s\n",newfilename.c_str());
			fptr_write=fopen(newfilename.c_str(),"w");
			xsum.clear();
			ysum[0].clear();
			ysum[1].clear();
			ysum[2].clear();
			ysum[3].clear();
			for(jdir=0;jdir<NJDIR;jdir++){
				oldfilename=dirprefix+"_"+to_string(jdir)+"/"+qualifier+"/"+acceptance+"/"+dirname[itype]+"/"+olddatafilename[ifn];
				fptr_read=fopen(oldfilename.c_str(),"r");
				ix=0;
				fscanf(fptr_read,"%lf %lf",&x,&y[0]);
				if(ifn==5)
					fscanf(fptr_read,"%lf %lf %lf",&y[1],&y[2],&y[3]);
				fgets(dummy,100,fptr_read);

				while(!feof(fptr_read) && ix<30){
					if(jdir==0){
						xsum.push_back(x);
						ysum[0].push_back(y[0]);
						if(ifn>3){
							ysum[1].push_back(y[1]);
							ysum[2].push_back(y[2]);
							ysum[3].push_back(y[3]);
						}
					}
					else{
						ysum[0][ix]+=y[0];
						if(ifn>3){
							ysum[1][ix]+=y[1];
							ysum[2][ix]+=y[2];
							ysum[3][ix]+=y[3];
						}
					}
					ix+=1;
					fscanf(fptr_read,"%lf %lf",&x,&y[0]);
					if(ifn==5)
						fscanf(fptr_read,"%lf %lf %lf",&y[1],&y[2],&y[3]);
					fgets(dummy,100,fptr_read);
				}
				fclose(fptr_read);
			}
			int nx=xsum.size();
			bfnorm=0.0;
			for(ix=0;ix<nx;ix++){
				y[0]=ysum[0][ix];
				if(ifn==0)
					y[0]=0.5*(ysum[0][ix]+ysum[0][nx-ix-1]);
				y[0]=y[0]/double(NJDIR);
				if(ifn==5){
					y[1]=ysum[1][ix]/double(NJDIR);
					y[2]=ysum[2][ix]/double(NJDIR);
					y[3]=ysum[3][ix]/double(NJDIR);
				}
				fprintf(fptr_write,"%5.2f %10.3e",xsum[ix],y[0]);
				if(ifn==5)
					fprintf(fptr_write," %10.3e %10.3e %10.3e",y[1],y[2],y[3]);
				fprintf(fptr_write,"\n");
				bfnorm+=y[0]*dx;
				if(itype==4 && ifn==0){
					printf("%2d %g\n",ix,y[0]);
				}
			}
			fclose(fptr_write);
			printf("=============== bfnorm=%g ===============\n",bfnorm);
		}
	}
	
	
	
	return 0;
}


