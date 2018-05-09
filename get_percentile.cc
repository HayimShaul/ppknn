#include <string.h>
#include <string>
#include <fstream>
#include <istream>
#include <iostream>
#include <sstream>
#include <algorithm>

#include "json.hpp"
using nlohmann::json;

#include <helib_number.h>
#include <unsigned_word.h>
#include <print.h>
#include <primes.h>

#include "thread_pool.h"
#include "get_percentile.h"


bool measureAccuracy = true;

std::vector<Point2D<float> > rawData;
std::vector<Point2D<int> > rawDiscreteData;

Point2D<float> query;
Point2D<int> discreteQuery;


// bounding box for Boston
Point2D<float> minPoint;
Point2D<float> maxPoint;

int maxX = 100;
int maxY = 100;

std::vector<Point2D<float> > read_data(const char *fname) {
	std::vector<Point2D<float> > points;

	std::istream *in;

	if (strcmp(fname, "-") == 0) {
		in = &(std::cin);
	} else {
		in = new std::ifstream(fname);
	}

	json j;
	(*in) >> j;

	if (j.find("HotelListResponse") != j.end()) {
		points.resize( j["HotelListResponse"]["HotelList"]["HotelSummary"].size() );
		j = j["HotelListResponse"]["HotelList"]["HotelSummary"];
		int i = 0;
		for (json::iterator it = j.begin(); it != j.end(); ++it) {
			Point2D<float> p( (*it)["longitude"], (*it)["latitude"] );

			std::cout << "point[" << i << "] = " << p << std::endl;
			points[i] = p;
			++i;
		}
	} else if (j.find("2Ddata") != j.end()) {
		points.resize( j["2Ddata"].size() );
		j = j["2Ddata"];
		int i = 0;
		for (json::iterator it = j.begin(); it != j.end(); ++it) {
			Point2D<float> p( (*it)["x"], (*it)["y"] );

			std::cout << "point[" << i << "] = " << p << std::endl;
			points[i] = p;
			++i;
		}
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
bool flex = false;
int thread_num = 1;


void initialize(int argc, char **argv) {
	char defaultFname[] = "-";
	char *fname = defaultFname;
	bool no_query = true;
	int n = -1;
	bool has_fname = false;

	for (int argc_i = 0; argc_i < argc; ++argc_i) {
		if (memcmp(argv[argc_i], "--p=", 4) == 0) {
			p = atoi(argv[argc_i] + 4);
			maxX = maxY = p/2;
		}
		if (memcmp(argv[argc_i], "--n=", 4) == 0)
			n = atoi(argv[argc_i] + 4);
		if (memcmp(argv[argc_i], "--L=", 4) == 0)
			L = atoi(argv[argc_i] + 4);
		if (memcmp(argv[argc_i], "--in=", 5) == 0) {
			fname = argv[argc_i] + 5;
			has_fname = true;
		}
		if (memcmp(argv[argc_i], "--q=", 4) == 0) {
			float x, y;
			if ((sscanf(argv[argc_i] + 4, "%f,%f", &x, &y) != 2) &&
					(sscanf(argv[argc_i] + 4, "(%f,%f)", &x, &y) != 2)) {
				std::cerr << "Error parsing query point " << (argv[argc_i] + 5) << std::endl;
				exit(1);
			}
			query = Point2D<float>(x,y);
			no_query = false;
		}
		if (memcmp(argv[argc_i], "--t=", 4) == 0)
			ThreadPool::init(atoi(argv[argc_i] + 4));

		if (strcmp(argv[argc_i], "--help") == 0) {
			std::cout << "   --L=" << std::endl;
			std::cout << "   --p=" << std::endl;
			std::cout << "   --n= how many points to generate (can't go with --in)" << std::endl;
			std::cout << "   --in=" << std::endl;
			std::cout << "   --q= query point" << std::endl;
			std::cout << "   --t= thread number" << std::endl;
		}
	}

	if (has_fname && ((n != -1) || (p > 0))) {
		std::cout << "Can't have --in together with --p or --n" << std::endl;
		exit(1);
	}

	if ((n > -1) && (p == 0)) {
		std::cout << "--p should be given when --n is given" << std::endl;
		exit(1);
	}


	if (has_fname) {
		rawData = read_data(fname);

		maxPoint = rawData[0];
		minPoint = rawData[0];
		for (auto i = rawData.begin(); i != rawData.end(); ++i) {
			maxPoint = max(maxPoint, *i);
			minPoint = min(minPoint, *i);
		}

		std::cout << "min value " << minPoint << std::endl;
		std::cout << "min value " << maxPoint << std::endl;

		discreteBase = minPoint;
		discreteResolution = Point2D<float>(maxX, maxY) / (maxPoint - minPoint);

		std::cout << "Discrete data:" << std::endl;
		rawDiscreteData.resize(0);
		for (auto i = rawData.begin(); i != rawData.end(); ++i) {
			Point2D<int> p = discretify(*i);
			rawDiscreteData.push_back(p);
			std::cout << "  " << p << std::endl;
		}
	} else {
		std::default_random_engine _generator(clock());
		std::uniform_int_distribution<int> xDistribution(0, maxX);
		std::uniform_int_distribution<int> yDistribution(0, maxY);
		for (int i = 0; i < n; ++i) {
			int x = xDistribution(_generator);
			int y = yDistribution(_generator);
			Point2D<float> raw_p( x, y );
			Point2D<int> p( x, y );

			std::cout << "point[" << i << "] = " << p << std::endl;
			rawData.push_back(raw_p);
			rawDiscreteData.push_back(p);
		}

		maxPoint = Point2D<float>(maxX, maxY);
		minPoint = Point2D<float>(0, 0);

		std::cout << "min value " << minPoint << std::endl;
		std::cout << "min value " << maxPoint << std::endl;

		discreteBase = minPoint;
		discreteResolution = Point2D<float>(maxX, maxY) / (maxPoint - minPoint);

		std::cout << "Discrete data:" << std::endl;
		for (auto i = rawDiscreteData.begin(); i != rawDiscreteData.end(); ++i) {
			std::cout << "  " << (*i) << std::endl;
		}
	}




	if (rawData.size() == 0) {
		std::cout << "No data read" << std::endl;
		exit(1);
	}


	if (no_query) {
		std::default_random_engine generator(clock());
		std::uniform_int_distribution<int> randomX(0, maxX);
		std::uniform_int_distribution<int> randomY(0, maxY);

		discreteQuery = Point2D<int>(randomX(generator), randomY(generator));
		query = undiscretify(discreteQuery);
	} else {
		discreteQuery = discretify(query);
	}
	


	std::cout << "Discrete query: "  << discreteQuery << std::endl;

	std::cout << "Discrete distances:" << std::endl;
	for (auto i = rawDiscreteData.begin(); i != rawDiscreteData.end(); ++i) {
		int dist = abs(i->x - discreteQuery.x) + abs(i->y - discreteQuery.y);
		std::cout << dist << std::endl;
	}

	int avg = 0;
	for (auto i = rawDiscreteData.begin(); i != rawDiscreteData.end(); ++i) {
		int dist = abs(i->x - discreteQuery.x) + abs(i->y - discreteQuery.y);
		avg += dist;
	}
	avg /= rawDiscreteData.size();

	int sigma = 0;
	for (auto i = rawDiscreteData.begin(); i != rawDiscreteData.end(); ++i) {
		int dist = abs(i->x - discreteQuery.x) + abs(i->y - discreteQuery.y);
		sigma += (dist - avg) * (dist - avg);
	}
	sigma /= rawDiscreteData.size();
	sigma = sqrt(sigma);

	std::cout << "avg = " << avg << std::endl;
	std::cout << "sigma = " << sigma << std::endl;

	std::vector<int> distribution_test(20);
	for (auto i = distribution_test.begin(); i != distribution_test.end(); ++i)
		*i = 0;

	for (auto i = rawDiscreteData.begin(); i != rawDiscreteData.end(); ++i) {
		int dist = abs(i->x - discreteQuery.x) + abs(i->y - discreteQuery.y);
		int bucket = (dist - avg) / ((sigma+2)/3) + 10;

		if (bucket > 18)
			bucket = 18;
		if (bucket < 0)
			bucket = 0;

		++distribution_test[bucket];
	}

	int bucket = 0;
	for (auto i = distribution_test.begin(); i != distribution_test.end(); ++i) {
		std::cout << ((bucket - 10) * ((sigma+2)/3) + avg)  << ": " << *i << std::endl;
		++bucket;
	}

	p = Primes::find_prime_bigger_than(maxX + maxY);

//	p = 2;
	r = 1;
//	while ((1 << r) < 2 * maxValue * maxValue)
//		++r;
//	++r;

//	int expected_L = mylog2(mylog2()) + mylog2(mylog2(rawData.size()));
	int expected_L = 15;
	if (L == 0)
		L = expected_L;

	std::cout << "using p = " << p << std::endl;
	std::cout << "using r = " << r << std::endl;
	std::cout << "using L = " << L << std::endl;
}

