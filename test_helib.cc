#include <iostream>
#include <string.h>

#include <carray_iterator.h>
#include <helib_number.h>
#include <eq.h>
#include <cmp.h>
#include <binomial_tournament.h>
#include <unsigned_word.h>
#include <average.h>

#include "point.h"
#include "one_mean_framework.h"
#include "one_mean.h"

typedef UnsignedWord<16, HelibNumber> HelibBits;


int main(int argc, char**argv) {
	measureAccuracy = false;

	initialize(argc, argv);

	HelibKeys keys;

	long R = 1;
	long r = 1;
	long d = 1;
	long c = 2;
	long k = 80;
	long s = 0;
	long chosen_m = 0;
	Vec<long> gens;
	Vec<long> ords;

	keys.initKeys(s, R, p, r, d, c, k, 64, L, chosen_m, gens, ords);
	HelibNumber::set_global_keys(&keys);


	Point<double> avgPoint;
	double avgNorm;

	one_mean<HelibNumber, HelibBits>(avgPoint, avgNorm, rawPoints, minPointManhattan, resolution);
}

