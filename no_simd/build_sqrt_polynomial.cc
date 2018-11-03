#include <assert.h>

#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_p.h>

#include <eq.h>
#include <primes.h>

using namespace NTL;

int p;
int phi_p;

const ZZ_pX poly_mul(const ZZ_pX &p1, const ZZ_pX &p2) {
	ZZ_pX ret = p1 * p2;

	for (int i = p; i <= deg(ret); ++i) {
		int c = i % phi_p;
		if (c == 0)
			c = phi_p;
//std::cerr << "folding " << i << " onto " << c << std::endl;
		ret[c] = ret[c] + ret[i];
		ret[i] = 0;
	}
	ret.normalize();
//std::cerr << "deg: " << deg(p1) << " * " << deg(p2) << " = " << deg(ret) << std::endl;
	return ret;
}

const ZZ_pX poly_power(const ZZ_pX &poly, int e) {
	if (e == 1)
		return poly;

	if ((e % 2) == 0)
		return poly_power(poly_mul(poly,poly), e/2);

	return poly_mul(poly, poly_power(poly, e-1));
}

int main(int argc, char **argv) {
	if (argc < 1) {
		std::cerr << "Usage: " << argv[0] << " <mod>" << std::endl;
		exit(1);
	}

	p = atoi(argv[1]);
	phi_p = ::phi(p);
	NTL::ZZ_p::init(NTL::ZZ(p));

	NTL::ZZ_pX sqrt;

	int x = 1;
	int x2 = 1;
	while (x2 < p) {
		NTL::ZZ_pX val(NTL::INIT_MONO, 0, 1);
		while (((0.5+x)*(0.5+x) > x2) && (x2 < p)) {
std::cerr << "x2 == " << x2 << "\r";
			val *= NTL::ZZ_pX(NTL::INIT_MONO, 0, x2) + NTL::ZZ_pX(NTL::INIT_MONO, 1, -1);  // val *= x2-x
			++x2;
		}
		val = poly_power(val, phi_p);
		val = (-val + 1) * x;
		sqrt += val;

std::cerr << "deg = " << deg(sqrt) << std::endl;

		++x;
	}

	std::cout << sqrt << std::endl;
}
