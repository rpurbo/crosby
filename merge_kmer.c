#ifndef MERGE_KMER_H
#define MERGE_KMER_H
#include "merge_kmer.h"
#endif


void merge_kmer(struct read_pairs *pair, struct PARAMS *param ){
        // build kmer db
        int i,j,k,l,m;
        int idx =0 ;
        int SEED = (param->ERR_RATIO * pair->read1_len) + 1;


        struct kmer_loc* kmer_db[pair->read1_len];   //head of the location linked list
        int kmer_index[4][4][4][4][4] = {0};                // use to locate the head for a 3mer


        for(i=0;i<4;i++){
                for(j=0;j<4;j++){
                        for(k=0;k<4;k++){
                                for(l=0;l<4;l++){
                                        for(m=0;m<4;m++){
                                                kmer_index[i][j][k][l][m] = -1;
                                        }
                                }
                        }
                }
        }


        for (i=0;i<pair->read1_len- 5;i=i+1){
                char ch1 = pair->read1[i];
                char ch2 = pair->read1[i+1];
                char ch3 = pair->read1[i+2];
                char ch4 = pair->read1[i+3];
                char ch5 = pair->read1[i+4];

                if(kmer_index[base2code[ch1]][base2code[ch2]][base2code[ch3]][base2code[ch4]][base2code[ch5]] == -1){

                        kmer_index[base2code[ch1]][base2code[ch2]][base2code[ch3]][base2code[ch4]][base2code[ch5]] = idx;
                        struct kmer_loc *head = (struct kmer_loc*)malloc(sizeof(struct kmer_loc));
                        head->loc = i;
                        head->next = NULL;
                        kmer_db[idx] = head;
                        idx++;

                }else{
                        int id = kmer_index[base2code[ch1]][base2code[ch2]][base2code[ch3]][base2code[ch4]][base2code[ch5]];
                        struct kmer_loc *element = (struct kmer_loc*)malloc(sizeof(struct kmer_loc));
                        struct kmer_loc *head = kmer_db[id];
                        element->loc = i;
                        element->next = head;
                        kmer_db[id] = element;

                }


        }
        
        // TSCAN READ2 and FILL MATRIX and OCC table

        bool diag_matrix[pair->read1_len][SEED];
        int diag_occ[pair->read1_len];
        memset(diag_matrix,0,sizeof(bool) * (pair->read1_len * SEED));
        memset(diag_occ,0,sizeof(int) * pair->read1_len);


        for (i=0;i<SEED;i++){
                char ch1 = pair->read2[i];
                char ch2 = pair->read2[i+1];
                char ch3 = pair->read2[i+2];
                char ch4 = pair->read2[i+3];
                char ch5 = pair->read2[i+4];

                int idx = kmer_index[base2code[ch1]][base2code[ch2]][base2code[ch3]][base2code[ch4]][base2code[ch5]];
                //printf("FOUND %c%c%c[%d]: ",ch1,ch2,ch3,idx);

                if(idx != -1){
                        struct kmer_loc *ptr = kmer_db[idx];
                        while(ptr != NULL){

                                int location = ptr->loc;
                                if((location-i) >= 0){
                                        diag_matrix[location-i][i] = true;
                                        diag_occ[location-i]++;

                                        /*
                                        if((location-i) == 0){
                                                printf("nil:%d-%d-%d\n",location,i,location-i);
                                        }
                                        */
                                }

                                ptr = ptr->next;
                        }
                        //printf("\n");
                }
                //printf("\n");
        }

        // GET THE BEST ALIGNMENT
        int best_score = 1;
        int startx;
        float match;
        float mismatch;
        float score;
        for(i=0;i<pair->read1_len - param->MIN_OVERLAP;i++){
                if(diag_occ[i] > 0){
                        int start = i;
                        match = 0;
                        mismatch = 0;

                        for(j=start;j<pair->read1_len;j++){
                                if(pair->read1[j] == pair->read2[j-start]){
                                        match++;
                                }else{
                                        mismatch++;
                                }

                        }
                        score = mismatch / (match + mismatch);

                        if(score < best_score && score < param->ERR_RATIO){
                                best_score = score;
                                startx = i;
                        }

                }
        }


        // MERGE THE READS
        if(best_score != 1){
                int merge_len = pair->read1_len + pair->read2_len - (pair->read1_len - startx);
                char merge_str[merge_len+1];
                char merge_qual[merge_len+1];
                memset(merge_str,0,sizeof(char) * (merge_len+1));
                memset(merge_qual,0,sizeof(char) * (merge_len+1));


                for (i=0;i<startx;i++){
                        merge_str[i] = pair->read1[i];
                        merge_qual[i] = pair->qual1[i];
                }

                // MIDDLE PART
                for (i=startx;i<pair->read1_len;i++){
                        if(pair->read1[i] == pair->read2[i-startx]){
                                merge_str[i] = pair->read1[i];
                                merge_qual[i] = pair->qual1[i];
                        }else{
                                merge_str[i] = 'N';
                                merge_qual[i] = '#';
                        }
                }

                // LAST PART
                for(i= pair->read1_len ;i <= merge_len;i++){
                        merge_str[i] = pair->read2[i-startx];
                        merge_qual[i] = pair->qual2[i-startx];

                }
                merge_str[merge_len+1] = '\0';
                merge_qual[merge_len+1] = '\0';


                // PRINT TO FILE
                fprintf(param->outFile,"%s\n",pair->header1);
                fprintf(param->outFile,"%s\n", merge_str);
                fprintf(param->outFile,"+\n", merge_str);
                fprintf(param->outFile,"%s\n", merge_qual);

        }	
}


