gcc -O3 -ftree-vectorize -ftree-vectorizer-verbose=1 -march=native -mtune=native -lm -o crosby crosby.c merge_kmer.c merge_brute.c merge_stagger.c
