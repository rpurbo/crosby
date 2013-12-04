#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

static char complement_tab[] = {
        ['A'] = 'T',
        ['C'] = 'G',
        ['G'] = 'C',
        ['T'] = 'A',
        ['N'] = 'N',
};


static int base2code[] = {
        ['A'] = 0,
        ['C'] = 1,
        ['G'] = 2,
        ['T'] = 3,
        ['N'] = 4,
};

struct read_entry {
        char *header;
        char *read;
        char *qual;
        int header_len;
        int read_len;
};

struct read_pairs {
        char *header1;
        char *read1;
        char *qual1;
        int header1_len;
        int read1_len;
        char *header2;
        char *read2;
        char *qual2;
        int header2_len;
        int read2_len;
};

struct PARAMS{
        int FRAGMENT_LEN;
        int MIN_OVERLAP;
        int MAX_OVERLAP;
        float ERR_RATIO;
        int KMER;
        char *PREFIX;
        FILE *outFile;
};

struct kmer_loc{
        int loc;
        struct kmer_loc* next;
};
