#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>
#include "bfduke/hydro2uds.h"
#include "msu_sampler/hyper.h"
#include "msu_commonutils/qualifier.h"
#include "msu_commonutils/sf.h"

using namespace std;
using namespace NMSUPratt;

int main(){
	double z,k0;
	z=1.23456;
	k0=Bessel::K0(z);
	printf("k0=%g\n",k0);
	
	return 0;
}


