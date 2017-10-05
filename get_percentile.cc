#include <string.h>
#include <string>
#include <fstream>
#include <istream>
#include <iostream>
#include <sstream>
#include <algorithm>

#include <helib_number.h>
#include <unsigned_word.h>
#include <print.h>
#include <primes.h>
#include "get_percentile.h"

bool measureAccuracy = true;


// bounding box for Manhattan
int minValue;
int maxValue;


std::vector<int> read_data(const char *fname) {
	std::vector<int> points;
	int x;

	std::istream *in;

	if (strcmp(fname, "-") == 0) {
		in = &(std::cin);
	} else {
		in = new std::ifstream(fname);
	}

	while ((*in) >> x) {
		points.push_back(x);
	}

	if (in != &(std::cin))
		delete in;

	return points;
}


int mylog2(int x) {
	int ret = 1;
	while ((1 << ret) < x)
		++ret;
	return ret;
}


long L = 0;
long p = 0;
long r = 0;
int resolution = 0;
bool flex = false;
std::vector<int> rawData;

void initialize(int argc, char **argv) {
	char defaultFname[] = "-";
	char *fname = defaultFname;

	for (int argc_i = 0; argc_i < argc; ++argc_i) {
		if (memcmp(argv[argc_i], "--p=", 4) == 0)
			p = atoi(argv[argc_i] + 4);
		if (memcmp(argv[argc_i], "--L=", 4) == 0)
			L = atoi(argv[argc_i] + 4);
		if (memcmp(argv[argc_i], "--in=", 5) == 0)
			fname = argv[argc_i] + 5;

		if (strcmp(argv[argc_i], "--help") == 0) {
			std::cout << "   --L=" << std::endl;
			std::cout << "   --size=" << std::endl;
			std::cout << "   --in=" << std::endl;
		}
	}

	rawData = read_data(fname);
	if (rawData.size() == 0) {
		std::cout << "No data read" << std::endl;
		exit(0);
	}

	minValue = *(std::min_element(rawData.begin(), rawData.end()));
	maxValue = *(std::max_element(rawData.begin(), rawData.end()));

	std::cout << "min value " << minValue << std::endl;
	std::cout << "min value " << maxValue << std::endl;

	if (p == 0)
		p = Primes::find_prime_bigger_than(maxValue);

//	p = 2;
	r = 1;
//	while ((1 << r) < 2 * maxValue * maxValue)
//		++r;
//	++r;

	int expected_L = mylog2(mylog2(maxValue)) + mylog2(mylog2(rawData.size()));
	if (L == 0)
		L = expected_L;

	std::cout << "using p = " << p << std::endl;
	std::cout << "using r = " << r << std::endl;
	std::cout << "using L = " << L << std::endl;
}

