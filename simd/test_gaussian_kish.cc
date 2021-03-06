#include <liphe/zp.h>
#include <liphe/helib_number.h>

#include <liphe/eq.h>
#include <liphe/polynomial.h>
#include <liphe/cmp.h>
#include <liphe/binomial_tournament.h>
#include <liphe/first_non_zero.h>
#include <liphe/unsigned_word.h>
#include <liphe/average.h>

#include "get_percentile_framework.h"
#include "get_percentile.h"
#include "special_polynomials.h"

#define SIMD_FACTOR 600
typedef ZP MyZP;
typedef UnsignedWord<18, MyZP> MyZPBits;

template<class Number>
Polynomial<Number> SpecialPolynomials<Number>::square_msd_polynomial;

template<class Number>
Polynomial<Number> SpecialPolynomials<Number>::sqrt_msd_polynomial;

template<class Number>
Polynomial<Number> SpecialPolynomials<Number>::sqrt_polynomial;

template<class Number>
Polynomial<Number> SpecialPolynomials<Number>::is_positive_polynomial;

//template<class Number>
//std::vector<Polynomial<Number> > SpecialPolynomials<Number>::convert_to_bit;

int main(int argc, char **argv) {
	ZP::set_global_simd_factor(SIMD_FACTOR);

	RETRIES = 5;
	MAX_CANDIDATES = -1;

//	srand(time(NULL));

	initialize(argc, argv);
	MyZP::set_global_p(p, r);

	SpecialPolynomials<MyZP>::init_polynomials(p);

//	secure_geo_search<MyZP, MyZPBits>(rawDiscreteData, discreteQuery);

//	secure_knn_classifier<MyZP, MyZPBits>(rawDiscreteData, rawDataClasses, discreteQuery);

	test_secure_knn_classifier<MyZP, MyZPBits>(rawDiscreteData, rawDataClasses);
	return 0;
}


