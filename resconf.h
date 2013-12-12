#ifndef CROSBY_H
#define CROSBY_H
#include "crosby.h"
#endif

struct RESULT {
        char base;
        char qual;
};

struct RESULT score_conflict(struct read_pairs*, struct PARAMS*, int*, int, int);
