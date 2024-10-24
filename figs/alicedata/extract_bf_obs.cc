#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>
#include <string>
#include <cstring>
#include <boost/math/special_functions.hpp>

const double PI=4.0*atan(1.0);
const double HBARC=197.3269602;

using namespace std;

int main(int argc,char *argv[]){
	char cline[300];
	double ymax[3]={1.0,1.2,1.4};
	double dy,width,norm,sigma;
	double yread,bf1rread,stat1read,sys1read,bf2read,stat2read,sys2read,bf3read,stat3read,sys3read;
	int iy,ny,itype,ivsv;
	char dummy[200];
	vector<double> y,y1,y2,bf,unc;
  FILE *fptrin;
	FILE *fptrout=fopen("experimental_info.txt","a");
	
	string filename,type[3]={"pipi","KK","pp"};
	string varname,vsv[2]={"vsy","vsphi"};
	
	for(ivsv=0;ivsv<2;ivsv++){
		for(itype=0;itype<3;itype++){
			filename="BF"+vsv[ivsv]+"_"+type[itype]+".txt";
			printf("filename=%s\n",filename.c_str());
			fptrin=fopen(filename.c_str(),"r");
			fgets(dummy,200,fptrin);
			y.clear(); y1.clear(); y2.clear(); bf.clear(); unc.clear();
			do{
				fscanf(fptrin,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
				&yread,&bf1rread,&stat1read,&sys1read,&bf2read,&stat2read,&sys2read,&bf3read,&stat3read,&sys3read);
				if(!feof(fptrin) && yread<ymax[itype]){
					y.push_back(yread);
					//printf("y=%g,bf=%g\n",yread,bf3read);
					bf.push_back(bf3read);
					unc.push_back(sqrt(stat3read*stat3read+sys3read*sys3read));
				}
			}while(!feof(fptrin));
			fclose(fptrin);
			ny=y.size();
			printf("---- ny=%d----\n",ny);
			y1.resize(ny);
			y2.resize(ny);
			norm=width=0.0;
			for(iy=0;iy<ny;iy++){
				printf("-- iy=%d\n",iy);
				if(iy==0)
					y1[iy]=0.0;
				else
					y1[iy]=0.5*(y[iy-1]+y[iy]);
				if(iy<ny-1)
					y2[iy]=0.5*(y[iy]+y[iy+1]);
				else
					y2[iy]=ymax[itype];
				dy=y2[iy]-y1[iy];
				printf("dy=%g\n",dy);
				norm+=dy*bf[iy];
				printf("++++++++  y=%g, y1=%g, y2=%g, bf=%g\n",0.5*(y1[iy]+y2[iy]),y1[iy],y2[iy],bf[iy]);
				width+=dy*bf[iy]*0.5*(y1[iy]+y2[iy]);
			}
			width=width/norm;
			if(ivsv==0){
				sigma=norm*0.05;
				varname="bfnorm_"+type[itype];
				fprintf(fptrout,"%s %9.5f %9.5f 0.0\n",varname.c_str(),norm,sigma);
			}
			sigma=width*0.03;
			varname="bf_"+vsv[ivsv]+"_width_"+type[itype];
			fprintf(fptrout,"%s %9.5f %9.5f 0.0\n",varname.c_str(),width,sigma);
		}
	}

	fclose(fptrout);
	
  return 0;
}


