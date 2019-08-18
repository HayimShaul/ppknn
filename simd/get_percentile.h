#ifndef ___ONE_MEAN___
#define ___ONE_MEAN___

#include <vector>
#include <helib_number.h>

#include "point.h"

extern long L;
extern long p;
extern long r;

extern int thread_num;
extern int keySize;

extern std::vector<Point<int> > rawDiscreteData;
extern std::vector<int> rawDataClasses;
extern bool measureAccuracy;

extern Point<int> discreteQuery;

void initialize(int argc, char **argv);

#endif
