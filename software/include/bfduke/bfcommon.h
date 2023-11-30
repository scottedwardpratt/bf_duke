#ifndef __HYDROBALANCECOMMON_H__
#define __HYDROBALANCECOMMON_H__

#include "msu_commonutils/commondefs.h"

class CHydroBalance;
class CHBHydroMesh;
class CHBHyperElement;
class CHBEoS;
class CHBCharge;
typedef multimap<int,CHBCharge* > CHBChargeMap;
typedef pair<int,CHBCharge* > CHBChargePair;
typedef multimap<int,CHBCharge* > mapic;
typedef pair<int,CHBCharge* > pairic;
typedef pair<int,CHBCharge* > pairip;

#endif