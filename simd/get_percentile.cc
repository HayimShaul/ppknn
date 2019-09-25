#include <string.h>
#include <string>
#include <fstream>
#include <istream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <random>

#include "json.hpp"
using nlohmann::json;

#include <helib_number.h>
#include <unsigned_word.h>
#include <print.h>
#include <primes.h>

#include "thread_pool.h"
#include "get_percentile.h"
#include "configuration.h"


bool measureAccuracy = true;

std::vector<Point<float> > rawData;
std::vector<Point<int> > rawDiscreteData;
std::vector<int> rawDataClasses;

Point<float> query;
Point<int> discreteQuery;


// bounding box for Boston
Point<float> minPoint;
Point<float> maxPoint;

Point<int> resolutionInt;
Point<float> resolutionFloat;

int keySize = 80;

void stringToType(float &n, const std::string &s) { n = std::stof(s); }
void stringToType(int &n, const std::string &s) { n = std::stoi(s); }

template<class Number>
std::vector<Number> parse_vector(const char *line) {
	std::vector<Number> v;
	std::string str(line);
	std::string delimiter(",");

	v.resize(0);
	size_t pos = 0;
	std::string token;
	while ((pos = str.find(delimiter)) != std::string::npos) {
		token = str.substr(0, pos);
		Number n;
		stringToType(n, token);
		v.insert(v.begin(), n);
		str.erase(0, pos + delimiter.length());
	}
	v.insert(v.begin(), std::stof(str));

	return v;
}

bool read_csv_line(std::istream *in, std::vector<float> &v, int &cls) {
	std::string str;
	std::string delimiter(",");

	try {
		v.resize(0);
		if (!getline(*in, str))
			return false;

		size_t pos = 0;
		std::string token;
		while ((pos = str.find(delimiter)) != std::string::npos) {
			token = str.substr(0, pos);
			v.insert(v.begin(), std::stof(token));
			str.erase(0, pos + delimiter.length());
		}

		cls = std::stoi(str);
	} catch (std::exception &e) {
		return false;
	}

	return true;
}

// read a csv file given as x,y,c   where c=0,1 is the class of a point and (x,y) is the coordinates
void read_classifier_data(const char *fname) {
	std::istream *in;

	if (strcmp(fname, "-") == 0) {
		in = &(std::cin);
	} else {
		in = new std::ifstream(fname);
	}

	std::vector<float> v;
	int c;
	while (read_csv_line(in, v, c)) {
		std::cout << "read ";
		for (unsigned int i = 0; i < v.size(); ++i)
			std::cout << v[i] << ", ";
		std::cout << c << std::endl;
		Point<float> p(v);
		rawData.push_back(p);
		rawDataClasses.push_back(c);
	}

	if (in != &(std::cin))
		delete in;
}

// read data from expedia json
//void read_hotel_data(const char *fname) {
//	std::istream *in;
//
//	if (strcmp(fname, "-") == 0) {
//		in = &(std::cin);
//	} else {
//		in = new std::ifstream(fname);
//	}
//
//	json j;
//	(*in) >> j;
//
//	if (j.find("HotelListResponse") != j.end()) {
//		rawData.resize( j["HotelListResponse"]["HotelList"]["HotelSummary"].size() );
//		j = j["HotelListResponse"]["HotelList"]["HotelSummary"];
//		int i = 0;
//		for (json::iterator it = j.begin(); it != j.end(); ++it) {
//			Point<float> p( (*it)["longitude"], (*it)["latitude"] );
//
//			std::cout << "point[" << i << "] = " << p << std::endl;
//			rawData[i] = p;
//			++i;
//		}
//	} else if (j.find("2Ddata") != j.end()) {
//		rawData.resize( j["2Ddata"].size() );
//		j = j["2Ddata"];
//		int i = 0;
//		for (json::iterator it = j.begin(); it != j.end(); ++it) {
//			Point<float> p( (*it)["x"], (*it)["y"] );
//
//			std::cout << "point[" << i << "] = " << p << std::endl;
//			rawData[i] = p;
//			++i;
//		}
//	}
//
//	if (in != &(std::cin))
//		delete in;
//}


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

enum SetResolutionMethod {
	ByP,
	Default,
	ByUser,
};

void initialize(int argc, char **argv) {
	char defaultFname[] = "-";
	char *fname = defaultFname;
	bool no_query = true;
	int n = -1;
	bool has_fname = false;
	SetResolutionMethod setResolutionMethod = SetResolutionMethod::Default;

	for (int argc_i = 0; argc_i < argc; ++argc_i) {
		if (memcmp(argv[argc_i], "--p=", 4) == 0) {
			p = atoi(argv[argc_i] + 4);
			setResolutionMethod = SetResolutionMethod::ByP;
		}
		if (memcmp(argv[argc_i], "--key=", 6) == 0)
			keySize = atoi(argv[argc_i] + 6);
		if (memcmp(argv[argc_i], "--res=", 6) == 0) {
			resolutionInt = Point<int>(parse_vector<int>(argv[argc_i] + 6));
			resolutionFloat = Point<float>(parse_vector<float>(argv[argc_i] + 6));
			setResolutionMethod = SetResolutionMethod::ByUser;
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
			query = Point<float>(parse_vector<float>(argv[argc_i] + 4));
			no_query = false;
		}
		if (memcmp(argv[argc_i], "--t=", 4) == 0)
			ThreadPool::init(atoi(argv[argc_i] + 4));
		if (memcmp(argv[argc_i], "--cv=", 5) == 0) {
			Configuration::cv_k_fold = atoi(argv[argc_i] + 5);
		}

		if (strcmp(argv[argc_i], "--help") == 0) {
			std::cout << "   --key= key size (default 80)" << std::endl;
			std::cout << "   --res=maxX,maxY,...   set the resolution of the grid" << std::endl;
			std::cout << "   --L=" << std::endl;
			std::cout << "   --p=" << std::endl;
			std::cout << "   --n= how many points to generate (can't go with --in)" << std::endl;
			std::cout << "   --in=  the csv file name of the input" << std::endl;
			std::cout << "   --q= query point" << std::endl;
			std::cout << "   --t= thread number" << std::endl;
			std::cout << "   --cv=1  how many samples to remove when testing. k-fold cv" << std::endl;
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


	if (!has_fname) {
		std::cerr << "Must provide an input filename" << std::endl;
		exit(1);
	}

	read_classifier_data(fname);

	auto maxit = std::max_element(rawDataClasses.begin(), rawDataClasses.end());
	if (maxit == rawDataClasses.end()) {
		std::cerr << "cannot find max element" << std::endl;
		exit(1);
	}
	Configuration::classNumber = *maxit + 1;

	if (setResolutionMethod == SetResolutionMethod::ByP) {
		resolutionInt.dim( rawData[0].dim() );
		resolutionFloat.dim( rawData[0].dim() );
		for (unsigned int i = 0; i < rawData[0].dim(); ++i) {
			resolutionFloat[i] = resolutionInt[i] = p/rawData[0].dim();
		}
	} else if (setResolutionMethod == SetResolutionMethod::Default) {
		resolutionInt.dim( rawData[0].dim() );
		resolutionFloat.dim( rawData[0].dim() );
		for (unsigned int i = 0; i < rawData[0].dim(); ++i) {
			resolutionFloat[i] = resolutionInt[i] = 100;
		}
	} else if (setResolutionMethod == SetResolutionMethod::ByUser) {
		if (resolutionInt.dim() == 1) {
			int res = resolutionInt[0];
			resolutionInt.dim(rawData[0].dim());
			resolutionFloat.dim(rawData[0].dim());
			for (unsigned int i = 0; i < resolutionFloat.dim(); ++i) {
				resolutionFloat[i] = res;
				resolutionInt[i] = res;
			}
		}
	}


	// Find the maximal and minimal x and y coordinates
	maxPoint = rawData[0];
	minPoint = rawData[0];
	for (auto i = rawData.begin(); i != rawData.end(); ++i) {
		maxPoint = max(maxPoint, *i);
		minPoint = min(minPoint, *i);
	}

	std::cout << "min value " << minPoint << std::endl;
	std::cout << "min value " << maxPoint << std::endl;

	discreteBase = minPoint;
	discreteResolution = resolutionFloat / (maxPoint - minPoint);

	std::cout << "Discrete data:" << std::endl;
	rawDiscreteData.resize(0);
	for (auto i = rawData.begin(); i != rawData.end(); ++i) {
		Point<int> p = discretify(*i);
		rawDiscreteData.push_back(p);
		std::cout << "  " << p << std::endl;
	}
//	} else {
//		std::cerr << "must provide input filename" << std::endl;
//		std::default_random_engine _generator(clock());
//		std::uniform_int_distribution<int> xDistribution(0, maxX);
//		std::uniform_int_distribution<int> yDistribution(0, maxY);
//		for (int i = 0; i < n; ++i) {
//			int x = xDistribution(_generator);
//			int y = yDistribution(_generator);
//			Point<float> raw_p( x, y );
//			Point<int> p( x, y );
//
//			std::cout << "point[" << i << "] = " << p << std::endl;
//			rawData.push_back(raw_p);
//			rawDiscreteData.push_back(p);
//		}
//
//		maxPoint = Point<float>(maxX, maxY);
//		minPoint = Point<float>(0, 0);
//
//		std::cout << "min value " << minPoint << std::endl;
//		std::cout << "min value " << maxPoint << std::endl;
//
//		discreteBase = minPoint;
//		discreteResolution = Point<float>(maxX, maxY) / (maxPoint - minPoint);
//
//		std::cout << "Discrete data:" << std::endl;
//		for (auto i = rawDiscreteData.begin(); i != rawDiscreteData.end(); ++i) {
//			std::cout << "  " << (*i) << std::endl;
//		}
//	}




	if (rawData.size() == 0) {
		std::cout << "No data read" << std::endl;
		exit(1);
	}


	if (no_query) {
		std::default_random_engine generator(clock());
		discreteQuery.dim(resolutionInt.dim());
		for (unsigned int i = 0; i < resolutionInt.dim(); ++i) {
			std::uniform_int_distribution<int> randomX(0, resolutionInt[i]);
			discreteQuery[i] = randomX(generator);
		}
		query = undiscretify(discreteQuery);
	} else {
		discreteQuery = discretify(query);
	}
	


	std::cout << "Discrete query: "  << discreteQuery << std::endl;

	std::cout << "Discrete distances:" << std::endl;
	for (auto i = rawDiscreteData.begin(); i != rawDiscreteData.end(); ++i) {
		int dist = (*i - discreteQuery).normL1();
		std::cout << dist << std::endl;
	}

	int avg = 0;
	for (auto i = rawDiscreteData.begin(); i != rawDiscreteData.end(); ++i) {
		int dist = (*i - discreteQuery).normL1();
		avg += dist;
	}
	avg /= rawDiscreteData.size();

	int sigma = 0;
	for (auto i = rawDiscreteData.begin(); i != rawDiscreteData.end(); ++i) {
		int dist = (*i - discreteQuery).normL1();
		sigma += (dist - avg) * (dist - avg);
	}
	sigma /= rawDiscreteData.size();
	sigma = sqrt(sigma);

	std::cout << "avg = " << avg << std::endl;
	std::cout << "sigma = " << sigma << std::endl;

	std::cout << "Histogram:" << std::endl;
	std::vector<int> distribution_test(20);
	for (auto i = distribution_test.begin(); i != distribution_test.end(); ++i)
		*i = 0;

	for (auto i = rawDiscreteData.begin(); i != rawDiscreteData.end(); ++i) {
		int dist = (*i - discreteQuery).normL1();
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

	std::cout << "Resolution = " << resolutionInt << " with norm " << resolutionInt.normL1() << std::endl;
	p = Primes::find_prime_bigger_than(2 * resolutionInt.normL1());

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

