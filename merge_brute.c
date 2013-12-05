#ifndef MERGE_KMER_H
#define MERGE_KMER_H
#include "merge_brute.h"
#endif

#include <math.h>

void merge_brute(struct read_pairs *pair, struct PARAMS *param ){
	
	int limit,match,mis,c,i,q1,q2;
	float pq1,pq2,Ppq1,Ppq2,chi;
	int maxscore =0;
	int start = -1;
	int obs[4];

	//printf("============================================\n");
	//printf("r1:%s\n",pair->read1 );
        //printf("r2:%s\n",pair->read2 );
	
	/*
	char *reads1_substring[pair->read1_len - param->MIN_OVERLAP];
	char *reads2_substring[pair->read2_len - param->MIN_OVERLAP];

	int len;
	for(i=0;i<pair->read1_len - param->MIN_OVERLAP;i++){
		len = pair->read1_len-i;
		reads1_substring[i] = malloc(sizeof(char) * (len+1));
		memcpy(reads1_substring[i],  &pair->read1[i], sizeof(char) * len);	
		reads1_substring[i][len] = '\0';
		//printf("%s\n",reads1_substring[i]);
		//printf("%p\n",(void *) &pair->read1[i]);

                reads2_substring[i] = malloc(sizeof(char) * (len+1));
                memcpy(reads2_substring[i],  &pair->read2[0], sizeof(char) * len);
                reads2_substring[i][len] = '\0';
		//printf("%s\n",reads2_substring[i]);
                //printf("%p\n",(void *) &pair->read1[i]);

	}
	*/

	// == NAIVE BRUTE FORCE ==
	for(limit = 0;limit < pair->read1_len - param->MIN_OVERLAP;limit++){
		match =0;
		mis =0;
    		c = 0;
		obs[base2code[pair->read1[limit]]]++;       	

    		for (i=0;i<pair->read2_len - limit;i++){
    
			if(pair->read1[limit+i] == pair->read2[i]){
				match++;
			}else{
				mis++;
			}	

        		if(mis >= (param->ERR_RATIO * (pair->read1_len - limit))){
            			i = pair->read1_len;
            			exit;
        		}
        
        		c++;
    		}
    
    		if(match > maxscore && mis < (param->ERR_RATIO * (pair->read1_len - limit))){
        		maxscore = match;
        		start = limit;
    		}
	}
	// == NAIVE BRUTE FORCE ==


        // MERGE THE READS
        if(maxscore != 0){
                int merge_len = pair->read1_len + pair->read2_len - (pair->read1_len - start);
                char merge_str[merge_len+1];
                char merge_qual[merge_len+1];
                memset(merge_str,0,sizeof(char) * (merge_len+1));
                memset(merge_qual,0,sizeof(char) * (merge_len+1));


                for (i=0;i<start;i++){
                        merge_str[i] = pair->read1[i];
                        merge_qual[i] = pair->qual1[i];
                }

                // MIDDLE PART
                for (i=start;i<pair->read1_len;i++){
                        if(pair->read1[i] == pair->read2[i-start]){
                                merge_str[i] = pair->read1[i];
                                //merge_qual[i] = pair->qual1[i];

				// IF BASES ARE EQUAL, CHOOSE THE LOWEST QUAL SCORE
				q1 = pair->qual1[i] - param->Q_OFFSET;
				q2 = pair->qual2[i-start] - param->Q_OFFSET;
			
				if(q1 > q2 ){
					merge_qual[i] = pair->qual2[i-start];
				}else{
					merge_qual[i] = pair->qual1[i];
				}

                        }else{
				// PLACEHOLDER METHOD
                                //merge_str[i] = 'N';
                                //merge_qual[i] = '#';
				

				q1 = pair->qual1[i] - param->Q_OFFSET;
				q2 = pair->qual2[i-start] - param->Q_OFFSET;

				pq1 = pow(10, (q1/-10));				
				pq2 = pow(10,(q2/-10));	
				
				// POSTERIOR PROBABILITY OF BASE 1
				Ppq1 = (pq1 * (obs[base2code[pair->read1[i]]]/pair->read1_len)) / 0.25; 						
				Ppq1 *= (((1-pq2)/3) * (obs[base2code[pair->read1[i]]]/pair->read1_len)) / 0.25;
				
				// POSTERIOR PROBABILITY OF BASE 2
                                Ppq2 = (pq2 * (obs[base2code[pair->read2[i-start]]]/pair->read1_len)) / 0.25;
                                Ppq2 *= (((1-pq1)/3) * (obs[base2code[pair->read2[i-start]]]/pair->read1_len)) / 0.25;

				// CALCULATE CHI-SQUARE TEST
				if(Ppq1 > Ppq2){
					chi = ((Ppq1 - Ppq2) * (Ppq1 - Ppq2)) / Ppq2;			
					if(chi >= 3.84){ //p-value 0.05 
						// different enough
						merge_str[i] = pair->read1[i];
						merge_qual[i] = pair->qual1[i];
					}else{
						merge_str[i] = 'N';
						merge_qual[i] = '#';
					}
				}else{
					chi = ((Ppq2 - Ppq1) * (Ppq2 - Ppq1)) / Ppq1;
                                       if(chi >= 3.84){ //p-value 0.05
                                                // different enough
                                                merge_str[i] = pair->read2[i-start];
                                                merge_qual[i] = pair->qual2[i-start];
                                        }else{
                                                merge_str[i] = 'N';
                                                merge_qual[i] = '#';
                                        }
				}
                        }
                }

                // LAST PART
                for(i= pair->read1_len ;i <= merge_len;i++){
                        merge_str[i] = pair->read2[i-start];
                        merge_qual[i] = pair->qual2[i-start];

                }
                merge_str[merge_len+1] = '\0';
                merge_qual[merge_len+1] = '\0';

		//printf("m:%s\n",merge_str );

		// PRINT TO FILE
                fprintf(param->outFile,"%s\n",pair->header1);
                fprintf(param->outFile,"%s\n", merge_str);
                fprintf(param->outFile,"+\n", merge_str);
                fprintf(param->outFile,"%s\n", merge_qual);

	}
}
