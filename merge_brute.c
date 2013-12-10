#ifndef MERGE_KMER_H
#define MERGE_KMER_H
#include "merge_brute.h"
#endif

#ifndef RESCONF_H
#define RESCONF_H
#include "resconf.h"
#endif

#include <math.h>

int merge_brute(struct read_pairs *pair, struct PARAMS *param ){
	
	int limit,match,mis,c,i,q1,q2;
	float pq1,pq2,Ppq1,Ppq2,chi;
	int maxscore =0;
	int start = -1;
	int obs[4];
	obs[0] = 0;
	obs[1] = 0;
	obs[2] = 0;
	obs[3] = 0;

	/* PROTOTYPE FOR VECTORIZATON
		
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

	// COUNT BASES OBSERVATION
	for (i=0;i<pair->read1_len;i++){
		obs[base2code[pair->read1[i]]]++;
	}


	// == NAIVE BRUTE FORCE ==
	for(limit = 0;limit < pair->read1_len - param->MIN_OVERLAP;limit++){
		match =0;
		mis =0;
    		c = 0;

    		for (i=0;i<pair->read2_len - limit;i++){

			if(pair->read1[limit+i] == pair->read2[i]){
				match++;
			}else{
				mis++;
			}	

        		if(mis > (param->ERR_RATIO * (pair->read1_len - limit))){
            			exit;
        		}
        
        		c++;
    		}
    
    		if(match > maxscore && mis < (param->ERR_RATIO * (pair->read1_len - limit))){
        		maxscore = match;
        		start = limit;
    		}
	}


	//printf("R1: %s\n",pair->read1);
	//printf("R2: %s\n",pair->read2);	

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
		

		/*	
                int merge_len = pair->read1_len + pair->read2_len - (pair->read1_len - start);
                char merge_str[merge_len+1] ;
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

				// IF BASES ARE EQUAL, CHOOSE THE LOWEST QUAL SCORE
				q1 = pair->qual1[i] - param->Q_OFFSET;
				q2 = pair->qual2[i-start] - param->Q_OFFSET;
			
				if(q1 < q2 ){
					merge_qual[i] = pair->qual2[i-start];
				}else{
					merge_qual[i] = pair->qual1[i];
				}

                        }else{
				// IF BASES ARE DIFFERENT, DO STATISTICAL TEST
				struct RESULT res;
				res = score_conflict(pair, param, obs, i, start);
				merge_str[i] = res.base;
				merge_qual[i] = res.qual;
                        
			}

                }
		
                // LAST PART
                for(i= pair->read1_len ;i <= merge_len;i++){
                        merge_str[i] = pair->read2[i-start];
                        merge_qual[i] = pair->qual2[i-start];
			
                }
                merge_str[merge_len+1] = '\0';
                merge_qual[merge_len+1] = '\0';
		*/
		

		// PRINT TO FILE
                fprintf(param->outFile,"%s\n",pair->header1);
                fprintf(param->outFile,"%s\n", merge_str);
                fprintf(param->outFile,"+\n");
                fprintf(param->outFile,"%s\n", merge_qual);

		return 1;
	}else{
	
		return 0;
	}
}
