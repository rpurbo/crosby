#ifndef MERGE_KMER_H
#define MERGE_KMER_H
#include "merge_brute.h"
#endif

#ifndef RESCONF_H
#define RESCONF_H
#include "resconf.h"
#endif


#include <math.h>

int merge_stagger(struct read_pairs *pair, struct PARAMS *param ){
	
	int limit,match,mis,c,i,q1,q2;
	float pq1,pq2,Ppq1,Ppq2,chi;
	int maxscore =0;
	int start = -1;
	int obs[4];
	obs[0] = 0;
	obs[1] = 0;
	obs[2] = 0;
	obs[3] = 0;

	// COUNT BASES OBSERVATION
	for (i=0;i<pair->read1_len;i++){
		obs[base2code[pair->read1[i]]]++;
	}

	// == NAIVE BRUTE FORCE - STAGGER MODE==
	for(limit = 0;limit < pair->read2_len - param->MIN_OVERLAP;limit++){
		match =0;
		mis =0;
    		c = 0;

    		for (i=0;i<pair->read1_len - limit;i++){
    
			if(pair->read2[limit+i] == pair->read1[i]){
				match++;
			}else{
				mis++;
			}	

        		if(mis >= (param->ERR_RATIO * (pair->read2_len - limit))){
            			//i = pair->read1_len;
            			exit;
        		}
        
        		c++;
    		}
    
    		if(match > maxscore && mis < (param->ERR_RATIO * (pair->read2_len - limit))){
        		maxscore = match;
        		start = limit;
    		}
	}
	// == NAIVE BRUTE FORCE ==


        // MERGE THE READS
        if(maxscore != 0){
                int total_len = pair->read1_len + pair->read2_len;

                char merge_str[total_len+1] ;
                char merge_qual[total_len+1];
                memset(merge_str, '\0', sizeof( merge_str ));
                memset(merge_qual, '\0', sizeof( merge_qual ));


                // COPY R1 TO MERGE STRING
                strncpy(merge_str, pair->read1,pair->read1_len);
                strncpy(merge_qual, pair->qual1,pair->read1_len);


                // MERGE R2 TO MERGE STRING
                for (i=0;i<pair->read2_len;i++){
                        c = i + start;
                        //printf("%d, %c, %c\n",i,merge_str[c], pair->read2[i]);

                        if(merge_str[c] == '\0'){

                                merge_str[c] = pair->read2[i];
                                merge_qual[c] = pair->qual2[i];

                        }else{

                                if(merge_str[c] == pair->read2[i]){

                                        // IF BASES ARE EQUAL, CHOOSE THE LOWEST QUAL SCORE
                                        q1 = merge_qual[c] - param->Q_OFFSET;
                                        q2 = pair->qual2[i] - param->Q_OFFSET;

                                        if(q1 > q2 ){
                                                merge_qual[c] = pair->qual2[i];
                                        }
                                }else{
                                        // IF BASES ARE DIFFERENT, DO STATISTICAL TEST
                                        struct RESULT res;
                                        res = score_conflict(pair, param, obs, c, start);
                                        merge_str[c] = res.base;
                                        merge_qual[c] = res.qual;


                                }

                        }
                }


		// PRINT TO FILE
                fprintf(param->outFile,"%s\n",pair->header1);
                fprintf(param->outFile,"%s\n", merge_str);
                fprintf(param->outFile,"+\n", merge_str);
                fprintf(param->outFile,"%s\n", merge_qual);

		return 1;
	}else{
	
		return 0;
	}
}
