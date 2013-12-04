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

#ifndef MERGE_BRUTE_H
#define MERGE_BRUTE_H
#include "merge_brute.h"
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
		
		//merge_kmer(&pair, &param);
		merge_brute(&pair, &param);		

	}
	
	fclose(fq1);
	fclose(fq2);
	fclose(oFile);
	
}





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
