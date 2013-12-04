/*
 * CROSBY - efficient overlapping paired-end reads merging tool
 *
 * CROSBY is a tool to merge paired-end reads which was designed to originate from a fragment that is smaller than 2xread length. 
 * The design produces two single-end reads that shares several basepairs in the 3'-ends. This enables experimentalist to utilize 
 * long contiguous reads instead of two separate short reads.
 *
 * CROSBY is designed to provide the following advanced features:
 * 	-  K-mer based extension identification. This feature is faster than brute-force algorithm in that it will only test smaller number
 * 	   of extension point indicated by k-mer location.
 * 	-  Automatic quality trimming. CROSBY provides option to automatically trim bases in 3'-end that have quality score lower than 2.
 * 	-  Statistics test for base conflict resolution. CROSBY employs statistics test to decide which bases to retain in the case of
 * 	   mismatch between pair.
 * 	-  Repeat detection and resolution. CROSBY provides a method to detect short repeats in the fragment.
 * 	-  Rescue mode. CROSBY provides a method to link two reads which shared less than 10bps in their extension point. The method
 * 	   examines all possible extension within the 10bps and compares it against a dynamic k-mer database.
 *
 * CROSBY source code can be downloaded from https://github.com/rpurbo/crosby
 *
 * Rikky W. Purbojati
 */

#ifndef CROSBY_H
#define CROSBY_H
#include "crosby.h"
#endif

#ifndef MERGE_KMER_H
#define MERGE_KMER_H
#include "merge_kmer.h"
#endif




void reverse(char*, size_t);
void reverse_complement(char*, size_t);
char complement(char);
	
int main(int argc, char **argv){
        struct PARAMS param;
        param.KMER = 5;
        param.ERR_RATIO = 0.25;
        param.MIN_OVERLAP = 10;

	// GETTING THE ARGUMENTS

	extern char *optarg; 
	extern int optind; 
	int c, err = 0; 
	int mflag=0, pflag=0, fflag=0; 
	char *sname = "default_sname", *fname; 
	static char usage[] = "usage: %s [-x max error ratio] [-L minimum overlap] read1 read2 output_prefix\n"; 
	while ((c = getopt(argc, argv, "k:x:L:")) != -1) 
		switch (c) { 
			case 'k':
				fflag = 1; 
				//param.KMER = atoi(optarg);
				break; 
			case 'x': 
				mflag = 1; 
				param.ERR_RATIO = atof(optarg);
				break; 
			case 'L': 
				pflag = 1; 
				param.MIN_OVERLAP = atoi(optarg);
				break; 
			case '?': 
				err = 1; 
				break; 
		}

	if (err) { 
		fprintf(stderr, usage, argv[0]); 
		exit(1); 
	}

	if((argc-optind) < 3){
		fprintf(stderr, usage, argv[0]);
		exit(1);
	}

	if(param.ERR_RATIO < 0 || param.ERR_RATIO > 1){
		fprintf(stderr, "\nError ratio must be between 0 and 1\n\n");
		fprintf(stderr, usage, argv[0]);
                exit(1);
	}

	size_t len = strlen(argv[optind]);
	char input1[len+1];
	strncpy(input1,argv[optind],len);
	input1[len] = '\0';
	len = strlen(argv[optind+1]);
        char input2[len+1];
        strncpy(input2,argv[optind+1],len);	
	input2[len] = '\0';
	len = strlen(argv[optind+2]);
        char output[len+1];
        strncpy(output,argv[optind+2],len);
	output[len] = '\0';

	// OPENING THE FILES

        FILE * oFile;
        oFile = fopen (output,"w");
	param.outFile= oFile;
	
	FILE* fq1;
	FILE* fq2;
	
        fq1 = fopen(input1,"r");
        fq2 = fopen(input2,"r");	
	
	char * line1 = NULL;
	char * line2 = NULL;
	size_t len1,len2 = 0;
	ssize_t num1,num2;
	
	// ITERATE THE READS	

	while (((num1 = getline(&line1, &len1, fq1)) != -1) && ((num2 = getline(&line2, &len2, fq2)) != -1)) {
		
		char header1[num1];
		char header2[num2];
		strncpy(header1, line1, num1);
		strncpy(header2, line2, num2);
		header1[num1-1] = '\0';
		header2[num2-1] = '\0';
		
		num1 = getline(&line1, &len1, fq1);
		num2 = getline(&line2, &len2, fq2);
		
		char seq1[num1];
		char seq2[num2];
		strncpy(seq1, line1, num1);
		strncpy(seq2, line2, num2);
		seq1[num1-1] = '\0';
		seq2[num2-1] = '\0';
		
		reverse_complement(seq2, num2-1);
			
		num1 = getline(&line1, &len1, fq1);
		num2 = getline(&line2, &len2, fq2);
		
		num1 = getline(&line1, &len1, fq1);
		num2 = getline(&line2, &len2, fq2);
		
		char qual1[num1];
		char qual2[num2];
		strncpy(qual1, line1, num1);
		strncpy(qual2, line2, num2);
		qual1[num1-1] = '\0';
		qual2[num2-1] = '\0';
		
		reverse(qual2,num2-1);
		
		struct read_pairs pair;
		pair.header1 = header1;
		pair.header2 = header2;
		pair.read1 = seq1;
		pair.read2 = seq2;
		pair.qual1 = qual1;
		pair.qual2 = qual2;
		pair.header1_len = strlen(header1);
		pair.header2_len = strlen(header2);
		pair.read1_len = strlen(seq1);
		pair.read2_len = strlen(seq2);
		
		merge_kmer(&pair, &param);
		
	}
	
	fclose(fq1);
	fclose(fq2);
	fclose(oFile);
	
}
/*
void merge_kmer(struct read_pairs pair, struct PARAMS param ){
	// build kmer db
	int i,j,k,l,m;
	int idx =0 ;
	int SEED = (param.ERR_RATIO * pair.read1_len) + 1;
	
	
	struct kmer_loc* kmer_db[pair.read1_len];   //head of the location linked list
	int kmer_index[4][4][4][4][4] = {0};		    // use to locate the head for a 3mer	
	
	
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
	
	
	for (i=0;i<pair.read1_len- 5;i=i+1){
		char ch1 = pair.read1[i];
		char ch2 = pair.read1[i+1];
		char ch3 = pair.read1[i+2];
		char ch4 = pair.read1[i+3];
		char ch5 = pair.read1[i+4];		

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
	
	
	 // TRAVERSE THE DB
	//printf("%s\n", pair.read1);
	//for (i=0;i<pair.read1_len-2;i=i+3){
	//	char ch1 = pair.read1[i];
	//	char ch2 = pair.read1[i+1];
	//	char ch3 = pair.read1[i+2];
	//	int id = kmer_index[base2code[ch1]][base2code[ch2]][base2code[ch3]];
	//	struct kmer_loc *ptr = kmer_db[id];
	//		
	//	printf("%c%c%c:",ch1,ch2,ch3);
	//	while(ptr != NULL){
	//		printf("%d-",ptr->loc);
	//		ptr = ptr->next;
	//	}
	//		
	//	printf("\n");		
	//	
	//}
	
	
	//printf("%s\n",pair.read1);
	//printf("%s\n",pair.read2);
	
	// TSCAN READ2 and FILL MATRIX and OCC table
	
	bool diag_matrix[pair.read1_len][SEED];
	int diag_occ[pair.read1_len];
	memset(diag_matrix,0,sizeof(bool) * (pair.read1_len * SEED));
	memset(diag_occ,0,sizeof(int) * pair.read1_len);
	
	
	//for (i=0;i<pair.read2_len-2 - param.MIN_OVERLAP;i++){
	for (i=0;i<SEED;i++){
		char ch1 = pair.read2[i];
		char ch2 = pair.read2[i+1];
		char ch3 = pair.read2[i+2];
		char ch4 = pair.read2[i+3];
		char ch5 = pair.read2[i+4];
		
		int idx = kmer_index[base2code[ch1]][base2code[ch2]][base2code[ch3]][base2code[ch4]][base2code[ch5]];
		//printf("FOUND %c%c%c[%d]: ",ch1,ch2,ch3,idx);
		
		if(idx != -1){
			struct kmer_loc *ptr = kmer_db[idx];
			while(ptr != NULL){
				
				int location = ptr->loc;
				if((location-i) >= 0){
					diag_matrix[location-i][i] = true;
					diag_occ[location-i]++;
					
					
					//if((location-i) == 0){
					//	printf("nil:%d-%d-%d\n",location,i,location-i);	
					//}
					
				}
				
				ptr = ptr->next;
			}
			//printf("\n");
		}
		//printf("\n");
	}
	
	
	//for(i=0;i<pair.read2_len;i++){
	//	printf("i:%d - cnt: %d\n",i,diag_occ[i]);
	//}
	//
	//printf("========\n");
	//for(i=0;i<pair.read1_len;i++){
	//	for(j=0;j<SEED;j++){
	//		printf("%d",diag_matrix[i][j]);
	//	}
	//	printf("\n");
	//}
	//printf("\n");
	
	
	
	// GET THE BEST ALIGNMENT
	int best_score = 1;
	int startx;
	float match;
	float mismatch;
	float score;
	for(i=0;i<pair.read1_len - param.MIN_OVERLAP;i++){
		if(diag_occ[i] > 0){
			int start = i;
			match = 0;
			mismatch = 0;
			
			for(j=start;j<pair.read1_len;j++){
				if(pair.read1[j] == pair.read2[j-start]){
					match++;
				}else{
					mismatch++;
				}
				
			}
			score = mismatch / (match + mismatch);
			
			if(score < best_score && score < param.ERR_RATIO){
				best_score = score;
				startx = i;
			}
			
		}
	}
	
	
	// MERGE THE READS
	if(best_score != 1){
		int merge_len = pair.read1_len + pair.read2_len - (pair.read1_len - startx);
		char merge_str[merge_len+1];
		char merge_qual[merge_len+1];
		memset(merge_str,0,sizeof(char) * (merge_len+1));
		memset(merge_qual,0,sizeof(char) * (merge_len+1));
	
			
		for (i=0;i<startx;i++){
			merge_str[i] = pair.read1[i];
			merge_qual[i] = pair.qual1[i];
		}	

		// MIDDLE PART
		for (i=startx;i<pair.read1_len;i++){
			if(pair.read1[i] == pair.read2[i-startx]){
				merge_str[i] = pair.read1[i]; 
				merge_qual[i] = pair.qual1[i];
			}else{
				merge_str[i] = 'N'; 
				merge_qual[i] = '#';
			}
		}

		// LAST PART
		for(i= pair.read1_len ;i <= merge_len;i++){
			merge_str[i] = pair.read2[i-startx];
			merge_qual[i] = pair.qual2[i-startx];
			
		}
		merge_str[merge_len+1] = '\0';
		merge_qual[merge_len+1] = '\0';
		
		
		// PRINT TO FILE
		fprintf(param.outFile,"%s\n",pair.header1);
		fprintf(param.outFile,"%s\n", merge_str);
		fprintf(param.outFile,"+\n", merge_str);
		fprintf(param.outFile,"%s\n", merge_qual);

				

	}
	
	
	
}
*/

void reverse(char *p, size_t len)
{
	char *pp = p + len - 1;
	while (p < pp) {
		char tmp = *p;
		*p = *pp;
		*pp = tmp;
		p++;
		pp--;
	}
}

void reverse_complement(char *p, size_t len)
{
	char *pp = p + len - 1;
	while (p <= pp) {
		char tmp = *p;
		*p = complement(*pp);
		*pp = complement(tmp);
		p++;
		pp--;
	}
}

char complement(char c)
{
	return complement_tab[(unsigned char)c];
}
