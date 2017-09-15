#include <zp.h>
#include <helib_number.h>

#include <eq.h>
#include "polynomial.h" // Before it gets into Liphe
#include <cmp.h>
#include <binomial_tournament.h>
#include <first_non_zero.h>
#include <unsigned_word.h>
#include <average.h>

#include "get_percentile_framework.h"
#include "get_percentile.h"

#define SIMD_FACTOR 2
typedef ZP<SIMD_FACTOR> MyZP;
typedef UnsignedWord<16, MyZP> MyZPBits;

int main(int argc, char **argv) {
	initialize(argc, argv);
	MyZP::set_global_p(p);

	get_percentile<MyZP, MyZPBits>(rawData, 15.8);
}

