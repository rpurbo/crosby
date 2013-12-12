#ifndef CROSBY_H
#define CROSBY_H
#include "crosby.h"
#endif

#ifndef RESCONF_H
#define RESCONF_H
#include "resconf.h"
#endif

#ifndef MATH_H
#define MATH_H
#include <math.h>
#endif


struct RESULT score_conflict(struct read_pairs *pair, struct PARAMS *param, int *obs, int i, int start){
	struct RESULT res; 
	int q1,q2;
	float pq1,pq2,Ppq1,Ppq2,chi;	

	q1 = pair->qual1[i] - param->Q_OFFSET;
	q2 = pair->qual2[i-start] - param->Q_OFFSET;

	pq1 = 1- (powf(10, ((float)q1/-10)));
	pq2 = 1- (powf(10,((float)q2/-10)));

	float PA = ((float)obs[base2code[pair->read1[i]]]/pair->read1_len);
	float PC = ((float)obs[base2code[pair->read2[i-start]]]/pair->read1_len);

	// POSTERIOR PROBABILITY OF BASE 1
	Ppq1 = (float)(pq1 * PA) / 0.25;
	Ppq1 *= ((float)((1-pq2)/3) * PA) / 0.25;
	
	// POSTERIOR PROBABILITY OF BASE 2
	Ppq2 = (float)(pq2 * PC) / 0.25;
	Ppq2 *= ((float)((1-pq1)/3) * PC) / 0.25;

	// CALCULATE CHI-SQUARE TEST
	if(Ppq1 > Ppq2){
		chi = ((Ppq1 - Ppq2) * (Ppq1 - Ppq2)) / Ppq2;
		if(chi >= param->CHISQUARE){ //p-value 0.05 
			res.base = pair->read1[i];
			res.qual = pair->qual1[i];
		}else{
			res.base = 'N';
			res.qual = '#';
		}
	}else{
		chi = ((Ppq2 - Ppq1) * (Ppq2 - Ppq1)) / Ppq1;
		if(chi >= param->CHISQUARE){ //p-value 0.05
			res.base = pair->read2[i-start];
			res.qual = pair->qual2[i-start];
		}else{
			res.base = 'N';
			res.qual = '#';
		}	

	}


	return res;
}
