#include "msu_commonutils/log.h"
#include "bfduke/hydro2uds.h"
#include "bfduke/bfcharge.h"
#include "msu_sampler/hyper.h"
#include "msu_commonutils/misc.h"

using namespace std;
using namespace NMSUPratt;

bool CHydroBalance::ReadDuke(CHBHydroMesh *hydromesh){
	char message[CLog::CHARLENGTH];
	double r,x,y,rmax=0.0,ur,tau;
	//double xbar=0.0,ybar=0.0,norm=0.0;
	bool keepgoing=true;
	int ix,iy,alpha;
	char dummy[300];
	double **pi,**pitilde;
	FourVector u;
	double s,e,p,t,vx,vy,vz,pi00,pi01,pi02,pi11,pi12,pi22,pi33,Pi;
	pi=new double*[4];
	pitilde=new double*[4];
	for(alpha=0;alpha<4;alpha++){
		pi[alpha]=new double[4];
		pitilde[alpha]=new double[4];
	}
	if(tau0readcheck){
		duke_filename="../hydrodata/"+qualifier+"/"+parmap.getS("HYDRODATA_FILENAME","evolution_xyeta.txt");
		snprintf(message,CLog::CHARLENGTH,"filename=%s\n",duke_filename.c_str());
		CLog::Info(message);
		fptr_duke=fopen(duke_filename.c_str(),"r");
		//for(iline=0;iline<14;iline++){
			//fgets(dummy,300,fptr_duke);
		//}
	}
	tau0readcheck=false;
	hydromesh->tau=TAU0+itauread*DELTAU;
	
	
	for(ix=0;ix<NX;ix++){
		for(iy=0;iy<NY;iy++){
			hydromesh->T[ix][iy]=0.0;
			hydromesh->epsilon[ix][iy]=0.0;
		}
	}
	//if(feof(fptr_duke)){
	//	keepgoing=false;
	//	fclose(fptr_duke);
	//}
	//else{
	
	if(feof(fptr_duke)){
		keepgoing=false;
	}
	else{
		keepgoing=true;
		fscanf(fptr_duke,"%lf",&t);
		if(!feof(fptr_duke)){
			biggestU=highestEpsilon=0.0;
			for(iy=0;iy<NY;iy++){
				for(ix=0;ix<NX;ix++){
					if(ix!=0 && iy!=0)	
						fscanf(fptr_duke,"%lf",&t);
					fscanf(fptr_duke,"%lf %lf %lf %lf %lf %lf %lf",&e,&s,&vx,&vy,&vz,&tau,&p);
					
					tau=TAU0+itauread*DELTAU;
					//if(fabs(tau-hydromesh->tau)>0.00001){
					//CLog::Fatal("reading in tau0="+to_string(tau)+", but hydromesh->tau="+to_string(hydromesh->tau)+"\n");
					//}
					//fscanf(fptr_duke,"%lf %lf %lf %lf %lf %lf %lf %lf",
					//&pi00,&pi01,&pi02,&pi11,&pi12,&pi22,&pi33,&Pi);
					pi00=pi01=pi02=pi11=pi12=pi22=pi33=Pi=0.0;
					pi[0][0]=pi00;
					pi[0][1]=pi[1][0]=pi01;
					pi[0][2]=pi[2][0]=pi02;
					pi[0][3]=pi[3][0]=0.0;
					pi[1][1]=pi11;
					pi[1][2]=pi[2][1]=pi12;
					pi[1][3]=pi[3][1]=0.0;
					pi[2][2]=pi22;
					pi[2][3]=pi[3][2]=0.0;
					pi[3][3]=pi33;
			
					u[0]=1.0/sqrt(1.0-vx*vx-vy*vy);
					u[1]=u[0]*vx;
					u[2]=u[0]*vy;
					u[3]=0.0;
					hydromesh->GetXY(ix,iy,x,y);
					//xbar+=x*e; ybar+=y*e; norm+=e;
					r=sqrt(x*x+y*y);
					if(r>rmax && e>=epsilon_f)
						rmax=r;
					Misc::BoostToCM(u,pi,pitilde);
					hydromesh->pitildexx[ix][iy]=pitilde[1][1];
					hydromesh->pitildexy[ix][iy]=pitilde[1][2];
					hydromesh->pitildeyy[ix][iy]=pitilde[2][2];
		
					fgets(dummy,300,fptr_duke);
					hydromesh->epsilon[ix][iy]=e;
					hydromesh->UX[ix][iy]=u[1];
					hydromesh->UY[ix][iy]=u[2];
					ur=sqrt(u[1]*u[1]+u[2]*u[2]);
					if(ur>biggestU && t>Tf)
						biggestU=ur;
					if(e>highestEpsilon)
						highestEpsilon=e;

				}
			}
		}
		else
			keepgoing=false;

		if(!keepgoing && highestEpsilon>epsilon_f){
			CLog::Info("Highest Epsilon="+to_string(highestEpsilon)+"\n");
			CLog::Fatal("Must increase HYPER_FREEZEOUT_EPSILON, hydro doesn't go long enough to reach this\n");
		}
		if(highestEpsilon<epsilon_f && hydromesh->tau>1.0){
			keepgoing=false;
			fclose(fptr_duke);
		}

		
	}
	for(alpha=0;alpha<4;alpha++){
		delete pi[alpha];
		pi[alpha]=NULL;
		delete pitilde[alpha];
		pitilde[alpha]=NULL;
	}
	delete pi;
	pi=NULL;
	delete pitilde;
	pitilde=NULL;
	
	itauread+=1;	
	return keepgoing;
}

void CHydroBalance::WriteCharges(int ichargefile){
	double f_l,f_s,gamma_q;
	char message[CLog::CHARLENGTH];
	string dirname="model_output/run"+to_string(run_number)+"/"+qualifier+"/udsdata";
	string command="mkdir -p "+dirname;
	system(command.c_str());
	string filename=dirname+"/uds"+to_string(ichargefile)+".txt";
	snprintf(message,CLog::CHARLENGTH,"writing charges to %s\n",filename.c_str());
	mapic::iterator it;
	CHBCharge *charge;
	Chyper *hyper;
	int balanceID;
	unsigned int icharge;
	FILE *fptr=fopen(filename.c_str(),"w");
	fprintf(fptr,"#  id   u  d  s      weight        tau           eta           x             y\n");
	it=emap.begin();
	for(it=emap.begin();it!=emap.end();++it){
		balanceID=it->first;
		charge=it->second;
		hyper=&(charge->hyper);
		GetGammaFQ(hyper->tau,gamma_q,f_l,f_s);
		fprintf(fptr,"%6d %2d %2d %2d %15.9f %15.9f %15.9f %15.9f %15.9f %15.9e %15.9e %15.9e %15.9e %15.9e %15.9e %15.9e %15.9e %15.9e %15.9e %15.9e %15.9e\n",
		balanceID,charge->q[0],charge->q[1],charge->q[2],charge->weight,charge->tau,charge->eta,
		charge->x,charge->y,hyper->T0,hyper->u[1],hyper->u[2],hyper->dOmega[0],hyper->dOmega[1],hyper->dOmega[2],hyper->pitilde[1][1],hyper->pitilde[2][2],hyper->pitilde[1][2],f_l,f_l,f_s);
		if(WRITE_TRAJ){
			if(charge->trajinfo!=NULL){
				for(icharge=0;icharge<charge->trajinfo->x.size();icharge++){
					snprintf(message,CLog::CHARLENGTH,"writing trajectory, icharge=%d\n",icharge);
					CLog::Info(message);
					fprintf(charge->trajinfo->fptr,"%8.5f %8.5f %8.5f %8.5f\n",
					charge->trajinfo->x[icharge],charge->trajinfo->y[icharge],charge->trajinfo->eta[icharge],charge->trajinfo->tau[icharge]);
				}
				fclose(charge->trajinfo->fptr);
			}
		}
		
	}
	fclose(fptr);
	snprintf(message,CLog::CHARLENGTH,"Ncolls/charge=%g\n",2.0*double(Ncollisions)/emap.size());
	CLog::Info(message);
}

void CHydroBalance::ClearCharges(){
	cmap.clear();
	emap.clear();
	nhbcharges=0;
}

void CHydroBalance::WriteSource(){
	// Note: DTAU for source mesh was 0.5 fm/c
	int a,b,itau;
	string dirname="model_output/run"+to_string(run_number)+"/udsdata/"+qualifier;
	string command="mkdir -p "+dirname;
	system(command.c_str());
	string filename=dirname+"/udssource.txt";
	FILE *fptr=fopen(filename.c_str(),"w");
	for(itau=0;itau<30;itau++){
		fprintf(fptr,"%5.2f ",(0.5+itau)*0.5);
		for(a=0;a<3;a++){
			for(b=a;b<3;b++)
				fprintf(fptr,"%7.0f ",0.5*(source[itau](a,b)+source[itau](b,a))/(0.5*NSAMPLE_HYDRO2UDS));
		}
		fprintf(fptr,"\n");
	}
	fclose(fptr);
}

void CHydroBalance::WriteFinalCF(){
	int ieta,a,b,c;
	int Netabins=parmap.getD("CF_NETABINS",50);
	double Deta=parmap.getD("CF_DETA",0.1);
	double eta,Z=NSAMPLE_HYDRO2UDS*Deta;
	Eigen::Matrix3d *cf=new Eigen::Matrix3d[Netabins];
	for(ieta=0;ieta<Netabins;ieta++)
		cf[ieta].setZero(3,3);
	mapic::iterator it;
	CHBCharge *charge1,*charge2;
	it=emap.begin();
	while(it!=emap.end()){
		charge1=it->second;
		++it;
		charge2=it->second;
		++it;
		eta=fabs(charge1->eta-charge2->eta);
		ieta=floorl(eta/Deta);
		if(ieta<Netabins){
			a=b=-1;
			for(c=0;c<3;c++){
				if(charge1->q[c]!=0)
					a=c;
				if(charge2->q[c]!=0)
					b=c;
			}
			cf[ieta](a,b)+=0.5*charge1->q[a]*charge2->q[b];
			cf[ieta](b,a)+=0.5*charge1->q[a]*charge2->q[b];
		}
	}
	
	string dirname="model_output/run"+to_string(run_number)+"/udsdata/"+qualifier;
	string command="mkdir -p "+dirname;
	system(command.c_str());
	string filename=dirname+"/cf_uds.txt";
	FILE *fptr=fopen(filename.c_str(),"w");
	fprintf(fptr,"#Deta       uu            ud          us            ss\n");
	for(ieta=0;ieta<Netabins;ieta++){
		fprintf(fptr,"%4.2f %11.5e %11.5e %11.5e %11.5e\n",
		(0.5+ieta)*Deta,cf[ieta](0,0)/Z,cf[ieta](0,1)/Z,cf[ieta](0,2)/Z,cf[ieta](2,2)/Z);
	}
	fclose(fptr);
	delete [] cf;
	cf=NULL;
}

void CHydroBalance::WriteHyper(){
	char message[CLog::CHARLENGTH];
	string dirname="model_output/run"+to_string(run_number)+"/"+qualifier+"/udsdata";
	string command="mkdir -p "+dirname;
	system(command.c_str());
	string filename=dirname+"/"+parmap.getS("HYPERDATA_FILENAME","hyper.txt");
	snprintf(message,CLog::CHARLENGTH,"writing hyper info to %s\n",filename.c_str());
	Chyper *hyper;
	FILE *fptr=fopen(filename.c_str(),"w");
	fprintf(fptr,"#Tf=%g\n",Tf);
	fprintf(fptr,"#   tau        x            y           Ux              Uy         dOmega0       dOmegaX       dOmegaY\n");
	list<Chyper *>::iterator it;
	for(it=hyperlist.begin();it!=hyperlist.end();++it){
		hyper=*it;
		fprintf(fptr,"%13.7e %13.7e %13.7e %13.7e %13.7e %13.7e %13.7e %13.7e %13.7e %13.7e %13.7e\n",
		hyper->tau,hyper->r[1],hyper->r[2],hyper->u[1],hyper->u[2],hyper->dOmega[0],hyper->dOmega[1],hyper->dOmega[2],hyper->pitilde[1][1],hyper->pitilde[2][2],hyper->pitilde[1][2]);
	}
	snprintf(message,CLog::CHARLENGTH,"Wrote %d hyper-elements\n",int(hyperlist.size()));
	fclose(fptr);
}

void CHydroBalance::WriteHyper_Duke_2D(){
	char message[CLog::CHARLENGTH];
	string dirname="model_output/run"+to_string(run_number)+"/"+qualifier+"/udsdata";
	string command="mkdir -p "+dirname;
	system(command.c_str());
	string filename=dirname+"/"+parmap.getS("HYPERDATA_FILENAME","hyper.txt");
	snprintf(message,CLog::CHARLENGTH,"writing hyper info to %s\n",filename.c_str());
	Chyper *hyper;
	FILE *fptr=fopen(filename.c_str(),"w");
	//fprintf(fptr,"#Tf=%g\n",Tf);
	fprintf(fptr,"#   tau        T        f_u      f_d       f_s      x            y           Ux              Uy         dOmega0       dOmegaX       dOmegaY       pi_xx          pi_yy        pi_xy\n");
	list<Chyper *>::iterator it;
	for(it=hyperlist.begin();it!=hyperlist.end();++it){
		hyper=*it;
		fprintf(fptr,"%13.7e %13.7f %13.7f %13.7f %13.7f %13.7f %13.7e %13.7e %13.7e %13.7e %13.7e %13.7e %13.7e %13.7e %13.7e %13.7e\n",
		hyper->tau,hyper->T0,hyper->epsilon,hyper->fugacity_u,hyper->fugacity_d,hyper->fugacity_s,hyper->r[1],hyper->r[2],hyper->u[1],hyper->u[2],hyper->dOmega[0],hyper->dOmega[1],hyper->dOmega[2],hyper->pitilde[1][1],hyper->pitilde[2][2],hyper->pitilde[1][2]);
	}
	snprintf(message,CLog::CHARLENGTH,"Wrote %d hyper-elements\n",int(hyperlist.size()));
	fclose(fptr);
}
