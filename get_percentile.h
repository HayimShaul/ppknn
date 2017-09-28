#ifndef ___ONE_MEAN___
#define ___ONE_MEAN___

#include <vector>
#include <helib_number.h>

extern int minValue;
extern int maxValue;

extern long L;
extern long p;
extern long r;
extern std::vector<int> rawData;
extern bool measureAccuracy;

void initialize(int argc, char **argv);

#endif
