#ifndef ___GET_PERCENTILE_FRAMEWORK___
#define ___GET_PERCENTILE_FRAMEWORK___

#include <sys/resource.h>
#include <sys/types.h>
#include <sys/time.h>
#include <vector>
#include <thread>

#include <time_measurements.h>

#include "polynomial.h"
#include "get_percentile.h"
#include "special_polynomials.h"
#include "input_database.h"


TakeTimes global_timer;


void print_stat_prefix(int i) {
	std::cout << "After ";

	if (i == -1)
		std::cout << "all";
	else
		std::cout << i;

	std::cout << " ";
}

int get_mem() {
	struct rusage rus;
	int res;

	res = getrusage(RUSAGE_SELF, &rus);
	if (res == -1) {
		perror("rusage");
		return 0;
	}

	return rus.ru_maxrss / 1024;
}

int realAvg;
int realAvgSqr;
int n;

void print_stat(int i) {
	print_stat_prefix(i);
	std::cout << global_timer.stats("Everything") << std::endl;

	print_stat_prefix(i);
	std::cout << "Used " << get_mem() << " MegaBytes" << std::endl;
}

template<class NUMBER>
void print_detailed_stats(int step, NUMBER &avgTemp, NUMBER &avgSqrMsdTemp, NUMBER &avgSqrLsdTemp) {
	int avg = avgTemp.to_int();
	int avgSqrLsd = avgSqrLsdTemp.to_int();   // this is x mod p
	int avgSqrMsd = avgSqrMsdTemp.to_int();   // this is x/2^a

	int avgSqr = avgSqrMsd *p + avgSqrLsd;

	print_stat_prefix(step);
	std::cout << std::endl;
	std::cout << "Avg " << avg << " / " << (realAvg / n) << " = " << ((double)avg * n / realAvg) << std::endl;
	std::cout << "AvgSqrLsd = " << avgSqrLsd << std::endl;
	std::cout << "AvgSqrMsd = " << avgSqrMsd << std::endl;
	std::cout << "AvgSqr " << avgSqr << " / " << (realAvgSqr / n) << " = " << ((double)avgSqr * n / realAvgSqr) << std::endl;

}

inline int sqr(int x) { return x*x; }

unsigned int statistic_rate = 10;


template<class Number>
Number fold(Number n) {
	AddBinomialTournament<Number> ret;

	int end = n.simd_factor();
	int simd_factor = n.simd_factor();
	int batch = 1;
	while (simd_factor != 0) {
		if ((simd_factor & 1) != 0) {
			ret.add_to_tournament( n.rotate_left(end - batch) );
			end -= batch;
		}

		n += n.rotate_left(batch);

		simd_factor >>= 1;
		batch <<= 1;
	}

	return ret.unite_all();
}

template<class Number, class NumberBits>
void multithreaded_average(Distances<Number, NumberBits> &_x, std::function<int(int)> f_m, Number &output) {
	std::thread *thread = new std::thread[thread_num];
	AverageLiphe<Number, NumberBits, CompareNative<NumberBits>, TruncConversion> avg;
	std::mutex mutex;

	avg.set_n(_x.size());
	avg.set_max_sample(output.p());
	avg.compute_resample_constant(0.1, 0.05);

	if ((bool)f_m)
		avg.set_f_m( f_m );


	for (int thread_i = 0; thread_i < thread_num; ++thread_i) {
		thread[thread_i] = std::thread( [thread_i, &avg, &_x, &mutex](){
			typename Distances<Number, NumberBits>::Iterator input = _x.begin();
			input += thread_i;

			while (input != _x.end()) {
				NumberBits xi = *input;
				avg.add_simd(xi, &mutex);

				++input;
			}
		} );
	}

	output = 0;
	for (int thread_i = 0; thread_i < thread_num; ++thread_i) {
		thread[thread_i].join();
	}
	output = avg.getAverage();

	delete[] thread;
}

template<class Number, class NumberBits>
void multithreaded_averages(Distances<Number, NumberBits> &_x, Number &avg, Number &sqrMsd, Number &sqrLsd) {
	std::thread *thread = new std::thread[thread_num];
	AverageLiphe<Number, NumberBits, CompareNative<NumberBits>, TruncConversion> avgAvg;
	AverageLiphe<Number, NumberBits, CompareNative<NumberBits>, TruncConversion> avgAvgSqrMsd;
	AverageLiphe<Number, NumberBits, CompareNative<NumberBits>, TruncConversion> avgAvgSqrLsd;

	std::mutex mutexAvg;
	std::mutex mutexAvgSqrMsd;
	std::mutex mutexAvgSqrLsd;

	avgAvg.set_n(_x.size());
	avgAvg.set_max_sample(avg.p());
	avgAvg.compute_resample_constant();

	avgAvgSqrMsd.set_n(_x.size());
	avgAvgSqrMsd.set_max_sample(sqrMsd.p());
	avgAvgSqrMsd.set_f_m( [](int m)->int{ return ::sqrt(m*p); } );
	avgAvgSqrMsd.compute_resample_constant();

	avgAvgSqrLsd.set_n(_x.size());
	avgAvgSqrLsd.set_max_sample(sqrLsd.p());
	avgAvgSqrLsd.set_f_m( [](int m)->int{ return ::sqrt(m); } );
	avgAvgSqrLsd.compute_resample_constant();


	ThreadPool threads;

	{
		AutoTakeTimes tt("populating cache");

		typename Distances<Number, NumberBits>::Iterator input = _x.begin();
		input.set_thread_pool(&threads);

		// populate the cache
		while (input != _x.end()) {
			while (!threads.has_free_cpu())
				threads.process_jobs(1);

			std::shared_ptr<NumberBits> xi(new NumberBits(*input));

			++input;
		}
		threads.process_jobs();
	}


	{
		AutoTakeTimes tt("computing averages");

		typename Distances<Number, NumberBits>::Iterator input = _x.begin();
		while (input != _x.end()) {
			while (!threads.has_free_cpu())
			threads.process_jobs(1);

			std::shared_ptr<NumberBits> xi(new NumberBits(*input));

			avgAvg.add_simd(xi, &mutexAvg, &threads);
			avgAvgSqrMsd.add_simd(xi, &mutexAvgSqrMsd, &threads);
			avgAvgSqrLsd.add_simd(xi, &mutexAvgSqrLsd, &threads);

			++input;
		}
		threads.process_jobs();
	}


	{
		AutoTakeTimes tt("getting averages");

		avg = avgAvg.getAverage();
		sqrMsd = avgAvgSqrMsd.getAverage();
		sqrLsd = avgAvgSqrLsd.getAverage();
	}

	delete[] thread;
}


#define MULTI_THREADED

inline float sqr(float x) { return x*x; }

void print_histogram(const std::vector<Point2D<int> > &sites, const std::vector<int> &classes, const Point2D<int> &query) {
	int avg = 0;
	for (auto i = sites.begin(); i != sites.end(); ++i) {
		int dist = abs(i->x - query.x) + abs(i->y - query.y);
		avg += dist;
	}
	avg /= sites.size();

	int sigma = 0;
	for (auto i = sites.begin(); i != sites.end(); ++i) {
		int dist = abs(i->x - query.x) + abs(i->y - query.y);
		sigma += (dist - avg) * (dist - avg);
	}
	sigma /= sites.size();
	sigma = sqrt(sigma);

	std::cout << "avg = " << avg << std::endl;
	std::cout << "sigma = " << sigma << std::endl;

	std::cout << "Histogram:" << std::endl;
	std::vector<int> distribution_test(20);
	std::vector<int> distribution_test_zero(20);
	std::vector<int> distribution_test_one(20);

	for (auto i = distribution_test.begin(); i != distribution_test.end(); ++i)
		*i = 0;
	distribution_test_zero = distribution_test;
	distribution_test_one = distribution_test;

	int i_site = 0;
	for (auto i = sites.begin(); i != sites.end(); ++i, ++i_site) {
		int dist = abs(i->x - query.x) + abs(i->y - query.y);
		int bucket = (dist - avg) / ((sigma+2)/3) + 10;

		if (bucket > 18)
			bucket = 18;
		if (bucket < 0)
			bucket = 0;

		++distribution_test[bucket];
		if (classes[i_site] == 0)
			++distribution_test_zero[bucket];
		else
			++distribution_test_one[bucket];
	}

	for (unsigned int bucket = 0; bucket < distribution_test.size(); ++bucket) {
		std::cout << ((bucket - 10) * ((sigma+2)/3) + avg)  << ": " << distribution_test[bucket] <<
		" = " << distribution_test_zero[bucket] << " + " << distribution_test_one[bucket] << std::endl;
	}
}

int real_knn_classifier(const std::vector<Point2D<int> > &sites, const std::vector<int> &classes, const Point2D<int> &query) {
	print_histogram(sites, classes, query);

	std::vector<float> distances;
	for (auto i_sites = sites.begin(); i_sites != sites.end(); ++i_sites) {
		float dist = sqr((*i_sites).x - query.x) + sqr((*i_sites).y - query.y);
		distances.push_back(dist);
	}
	std::sort(distances.begin(), distances.end());

	float threshold = distances[distances.size() * 0.02];

	int classOne = 0;
	int classZero = 0;

	for (unsigned int i_sites = 0; i_sites < sites.size(); ++i_sites) {
		float dist = sqr(sites[i_sites].x - query.x) + sqr(sites[i_sites].y - query.y);
		if (dist < threshold) {
			if (classes[i_sites] == 0)
				++classZero;
			else if (classes[i_sites] == 1)
				++classOne;
			else {
				std::cerr << "Error: classes should be 0s of 1s\n";
				exit(1);
			}
		}
	}

	std::cout << "real count of class 0: " << classZero << std::endl;
	std::cout << "real count of class 1: " << classOne << std::endl;

	return (classZero > classOne) ? 0 : 1;
}


template<class Number, class NumberBits>
void secure_knn_classifier_gaussian(const std::vector<Point2D<int> > &sites, const std::vector<int> &classes, const Point2D<int> &query, std::vector<int> &classZeroCountVector, std::vector<int> &classOneCountVector) {

	std::cout << "Starting classifier KNN" << std::endl;

	Distances<Number, NumberBits> distances(sites, classes, query);

	Number avgValue;
	Number avgSqrValue;

	n = sites.size();

	////////////////////////////////////
	// Start computing average and average of squares
	////////////////////////////////////

	Number avgEnc;
	Number avgSqrMsdEnc;
	Number avgSqrLsdEnc;

	{
		global_timer.start();
		multithreaded_averages<Number, NumberBits>(distances, avgEnc, avgSqrMsdEnc, avgSqrLsdEnc);
		std::cout << global_timer.end("computin averages");
	}


	{
		AutoTakeTimes tt("folding vectors");

		global_timer.start();
		avgEnc = fold(avgEnc);
		avgSqrMsdEnc = fold(avgSqrMsdEnc);
		avgSqrLsdEnc = fold(avgSqrLsdEnc);
		std::cout << global_timer.end("folding");
	}


	////////////////////////////////////
	// End computing average and average of squares
	////////////////////////////////////



	////////////////////////////////////
	// start computing threshold
	////////////////////////////////////

	std::vector<NumberBits> thresholdCandidatesBits;
	{
		ThreadPool threads;

		global_timer.start();

		Number avgLsdEnc = avgEnc * avgEnc;
		Number avgMsdEnc;
		{
			AutoTakeTimes tt("computing avg sqaure");
			SpecialPolynomials<Number>::square_msd_polynomial.compute(avgMsdEnc, avgEnc, &threads);
		}

		Number sigma;
		{
			AutoTakeTimes tt("computing sigma");
			sigma = SpecialPolynomials<Number>::sqrt_msd_polynomial.compute( avgSqrMsdEnc - avgMsdEnc, &threads );
			sigma += SpecialPolynomials<Number>::sqrt_polynomial.compute( avgSqrLsdEnc - avgLsdEnc, &threads );
		}

		int inv2 = power_mod(2, phi(Number::get_global_ring_size()) - 1, Number::get_global_ring_size());
		std::cout << "2^{-1} mod " << Number::get_global_ring_size() << " = " << inv2 << std::endl;

		{
			std::vector<Number> thresholdCandidates;
			{
				Number threshold = avgEnc;
				thresholdCandidates.push_back(threshold);
			}
			{
				Number threshold = avgEnc - sigma;
				thresholdCandidates.push_back(threshold);
			}
			{
				Number threshold = avgEnc - sigma*2;
				thresholdCandidates.push_back(threshold);
			}
			{
				Number threshold = avgEnc - sigma*inv2;
				thresholdCandidates.push_back(threshold);
			}
			{
				Number threshold = avgEnc - (sigma+1)*inv2;
				thresholdCandidates.push_back(threshold);
			}

			{
				AutoTakeTimes tt("converting threshold to bits");
				thresholdCandidatesBits.resize(thresholdCandidates.size());

				for (unsigned int i_candidate = 0; i_candidate < thresholdCandidates.size(); ++i_candidate) {
					convert_to_bits<NumberBits, Number>(thresholdCandidatesBits[i_candidate], thresholdCandidates[i_candidate], &threads);
				}
				threads.process_jobs();
			}
		}


		std::cout << global_timer.end("compute threasholds");





		int realAvg = 0;
		int realAvgSqr = 0;
		for (auto i = sites.begin(); i != sites.end(); ++i) {
			int dist = abs((*i).x - query.x) + abs((*i).y - query.y);
			realAvg += dist;
			realAvgSqr += dist*dist;
		}
		realAvg /= sites.size();
		realAvgSqr /= sites.size();

//		std::cout << "real avg = " << realAvg << std::endl;
//
//		std::cout << "real avg^2 = " << (realAvg*realAvg) << std::endl;
//		std::cout << "real avg^2 / p = " << (realAvg*realAvg / p) << std::endl;
//		std::cout << "real avg^2 % p = " << (realAvg*realAvg % p) << std::endl;
//
//		std::cout << "real avgSqr = " << realAvgSqr << std::endl;
//		std::cout << "real avgSqr / p = " << (realAvgSqr / p) << std::endl;
//		std::cout << "real avgSqr % p = " << (realAvgSqr % p) << std::endl;

		std::cout << "avg = " << avgEnc.to_int()  << "         (real = " << realAvg << ")" << std::endl;

		std::cout << "avg^2 = " << (avgMsdEnc.to_int()*p + avgLsdEnc.to_int()) << "        (real = " << (realAvg*realAvg) << ")" << std::endl;
		std::cout << "avg^2 / p = " << avgMsdEnc.to_int() <<  "            (real = " << (realAvg*realAvg / p) << ")" << std::endl;
		std::cout << "avg^2 % p = " << avgLsdEnc.to_int() << "             (real = " << (realAvg*realAvg % p) << ")" << std::endl;

		std::cout << "avgSqr = " << (avgSqrMsdEnc.to_int()*p + avgSqrLsdEnc.to_int()) << "          (real = " << realAvgSqr << ")" << std::endl;
		std::cout << "avgSqr / p = " << avgSqrMsdEnc.to_int() << "        (real = " << (realAvgSqr / p) << ")" << std::endl;
		std::cout << "avgSqr % p = " <<  avgSqrLsdEnc.to_int() << "        (real = " << (realAvgSqr % p) << ")" << std::endl;


		std::cout << "sigma = sqrt_msd(" << (avgSqrMsdEnc - avgMsdEnc).to_int() << ") + " << "sqrt_lsd(" << (avgSqrLsdEnc - avgLsdEnc).to_int() << ")" << std::endl;
		std::cout << "sigma = " << sigma.to_int() << "            (real = " << ::sqrt(-realAvg*realAvg + realAvgSqr) << ")" << std::endl;
		for (unsigned int i_candidate = 0; i_candidate < thresholdCandidatesBits.size(); ++i_candidate)
			std::cout << "Threshold candidate " << i_candidate << " = " << thresholdCandidatesBits[i_candidate].to_int() << std::endl;

	}

	////////////////////////////////////
	// end computing threshold
	////////////////////////////////////


	print_stat(-1);
	if (measureAccuracy) {
		print_detailed_stats(-1, avgEnc, avgSqrMsdEnc, avgSqrLsdEnc);
	}


	////////////////////////////////////
	// start counting classes
	////////////////////////////////////

	{
		auto dist = distances.begin();
		while (dist != distances.end()) {
			NumberBits xi = *dist;
			++dist;
		}
	}
	
	////////////////////////////////////
	// end counting classes
	////////////////////////////////////





	////////////////////////////////////
	// start reporting points
	////////////////////////////////////

	std::vector<Number> classOneCountEnc(thresholdCandidatesBits.size());
	std::vector<Number> classZeroCountEnc(thresholdCandidatesBits.size());

	for (unsigned int i_candidate = 0; i_candidate < thresholdCandidatesBits.size(); ++i_candidate) {
		classOneCountEnc[i_candidate] = 0;
		classZeroCountEnc[i_candidate] = 0;
	}

	auto dist = distances.begin();
	while (dist != distances.end()) {
		NumberBits xi = dist.getDistances();
		std::vector<int> classes = dist.getClasses();
		++dist;

		// classOne and classZero are bitmaps with 1 for samples with 
		std::vector<long int> classOne(classes.size());
		std::vector<long int> classZero(classes.size());
		for (unsigned int i_class = 0; i_class < classes.size(); ++i_class) {
			classZero[i_class] = (classes[i_class] == 0) ? 1 : 0;
			classOne[i_class] = (classes[i_class] == 1) ? 1 : 0;
		}

		global_timer.start();
		{
			for (unsigned int i_candidate = 0; i_candidate < thresholdCandidatesBits.size(); ++i_candidate) {
				Number knnEnc = xi < thresholdCandidatesBits[i_candidate];
				classOneCountEnc[i_candidate] += knnEnc * classOne;
				classZeroCountEnc[i_candidate] += knnEnc * classZero;

				std::cout << "The KNN indicator vector for candidate " << thresholdCandidatesBits[i_candidate].to_int() << std::endl;
				std::vector<long int> knn = knnEnc.to_vector();
				unsigned int i = 0;
				for (auto ri = knn.begin(); ri != knn.end(); ++ri) {
					if (i < sites.size()) {
						if ((*ri != 0) && (*ri != 1))
							std::cout << "Error: ri = " << (*ri) << " which is not binary" << std::endl;
						std::cout << i << ") " << "x = " << *ri << "   dist = " << distances.getPlaintextDistance(i) << "    class = " << classes[i] << std::endl;
						++i;
					}
				}
			}
		}
	}

	std::cout << global_timer.end("reporting points");
//	std::cout << "depth = mul " << classOneCountEnc0.mul_depth() << " add " << classOneCountEnc0.add_depth() << std::endl;

//	if (OK)
//		std::cout << "test is ok" << std::endl;

	
	for (unsigned int i_candidate = 0; i_candidate < thresholdCandidatesBits.size(); ++i_candidate) {
		std::vector<long int> classOneCountV = classOneCountEnc[i_candidate].to_vector();
		std::vector<long int> classZeroCountV = classZeroCountEnc[i_candidate].to_vector();
		int classOneCount = 0;
		for (int n : classOneCountV) classOneCount += n;
		int classZeroCount = 0;
		for (int n : classZeroCountV) classZeroCount += n;

		std::cout << "candidate " << i_candidate << " secure count of class 0: " << classZeroCount << std::endl;
		std::cout << "candidate " << i_candidate << " secure count of class 1: " << classOneCount << std::endl;
		classZeroCountVector.push_back(classZeroCount);
		classOneCountVector.push_back(classOneCount);
	}


//	int secureKnnClassifier = (classZeroCount > classOneCount) ? 0 : 1;

	////////////////////////////////////
	// end reporting points
	////////////////////////////////////

	print_stat(-1);

//	return secureKnnClassifier;
}

template<class Number, class NumberBits>
int secure_knn_classifier_blackhole(const std::vector<Point2D<int> > &sites, const std::vector<int> &classes, const Point2D<int> &query, int &classZeroCount, int &classOneCount) {

	std::cout << "Starting classifier blackhole" << std::endl;

	Distances<Number, NumberBits> distances(sites, classes, query);
	unsigned int i = 0;

	Number avgValue;
	Number avgSqrValue;

	n = sites.size();


	////////////////////////////////////
	// start counting classes
	////////////////////////////////////

	{
		auto dist = distances.begin();
		while (dist != distances.end()) {
			NumberBits xi = *dist;
			++dist;
		}
	}
	
	////////////////////////////////////
	// end counting classes
	////////////////////////////////////


	std::vector<long int> thres_vector;
	for (unsigned int i_threshold = 0; i_threshold < Number::simd_factor(); ++i_threshold) 
		thres_vector.push_back(1);


	////////////////////////////////////
	// start reporting points
	////////////////////////////////////

	bool OK = true;
	i = 0;
	Number classOneCountEnc(0);
	Number classZeroCountEnc(0);

	NumberBits thresholdBits(thres_vector);
	
	auto dist = distances.begin();
	while (dist != distances.end()) {
		NumberBits xi = dist.getDistances();
		std::vector<int> classes = dist.getClasses();
		++dist;

		// classOne and classZero are bitmaps with 1 for samples with 
		std::vector<long int> classOne(classes.size());
		std::vector<long int> classZero(classes.size());
		for (unsigned int i_class = 0; i_class < classes.size(); ++i_class) {
			classZero[i_class] = (classes[i_class] == 0) ? 1 : 0;
			classOne[i_class] = (classes[i_class] == 1) ? 1 : 0;
		}

		Number classOneEnc(classOne);
		Number classZeroEnc(classZero);

		global_timer.start();
		Number knnEnc = xi < thresholdBits;
		classOneCountEnc += knnEnc * classOneEnc;
		classZeroCountEnc += knnEnc * classZeroEnc;
		std::cout << global_timer.end("reporting points");

		std::cout << "depth = mul " << classOneCountEnc.mul_depth() << " add " << classOneCountEnc.add_depth() << std::endl;

		std::vector<long int> knn = knnEnc.to_vector();
		for (auto ri = knn.begin(); ri != knn.end(); ++ri) {
			if (i < sites.size()) {
				if ((*ri != 0) && (*ri != 1))
					OK = false;
				std::cout << i << ") " << "x = " << *ri << std::endl;
				++i;
			}


		}
	}

	if (OK)
		std::cout << "test is ok" << std::endl;

	std::vector<long int> classOneCountVector = classOneCountEnc.to_vector();
	std::vector<long int> classZeroCountVector = classZeroCountEnc.to_vector();
	classOneCount = 0;
	for (int n : classOneCountVector) classOneCount += n;
	classZeroCount = 0;
	for (int n : classZeroCountVector) classZeroCount += n;

	std::cout << "secure count of class 0 (blackhole): " << classZeroCount << std::endl;
	std::cout << "secure count of class 1 (blackhole): " << classOneCount << std::endl;


	int secureKnnClassifier = (classZeroCount > classOneCount) ? 0 : 1;

	////////////////////////////////////
	// end reporting points
	////////////////////////////////////

	print_stat(-1);

	return secureKnnClassifier;
}

int avgIterations = 0;

template<class Number, class NumberBits>
int secure_knn_classifier(const std::vector<Point2D<int> > &sites, const std::vector<int> &classes, const Point2D<int> &query) {
	std::vector<int> classZeroCount;
	std::vector<int> classOneCount;

	for (int i = 0; i < 5; ++i) {
		++avgIterations;
		classZeroCount.resize(0);
		classOneCount.resize(0);
		secure_knn_classifier_gaussian<Number, NumberBits>(sites, classes, query, classZeroCount, classOneCount);

		for (unsigned int sigmaFactor = 0; sigmaFactor < classOneCount.size(); ++sigmaFactor) {
			int one = classOneCount[sigmaFactor];
			int zero = classZeroCount[sigmaFactor];
			std::cout << "Iteration " << i << " sigmaFactor " << sigmaFactor << " one=" << one << "  zero=" << zero << std::endl;

			if (one + zero < 0.02 * sites.size()) {
				std::cout << "Too little neighbors" << std::endl;
				continue;
			}
			// Deal with the case where we have a small count of 1 but huge amount of 0
			if (one < 0.05 * sites.size()) {
				std::cout << "enough neighbors, and 1 count is very small. Classifying as 0" << std::endl;
				return 0;
			}
			if (zero < 0.05 * sites.size()) {
				std::cout << "enough neighbors, and 0 count is very small. Classifying as 1" << std::endl;
				return 1;
			}
			if (one + zero > 0.1 * sites.size()) {
				std::cout << "Too many neighbors" << std::endl;
				continue;
			}
			std::cout << "enough neighbors. Classifying by majority" << std::endl;
			return (one > zero) ? 1 : 0;
		}
	}
	std::cout << "ITERATIONS: did not get a small neighborhood - trying blackhole" << std::endl;

//	ret = secure_knn_classifier_blackhole<Number, NumberBits>(sites, classes, query, classZeroCount, classOneCount);
//	if (((classZeroCount + classOneCount) < sites.size() * 0.04) ||
//			((classZeroCount > sites.size() * 0.1) && (classOneCount > sites.size() * 0.1))) {
//		std::cout << "ITERATIONS: blackhole failed to find a good neighborhood: zero=" << classZeroCount << " one=" << classOneCount << std::endl;
//		return -1;
//	}

	return -1;
}

template<class Number, class NumberBits>
void test_secure_knn_classifier(const std::vector<Point2D<int> > &sites, const std::vector<int> &classes) {
	int match = 0;
	int mismatch = 0;
	int secClassificationFailed = 0;

	int secureCorrect = 0;
	int secureIncorrect = 0;

	int realCorrect = 0;
	int realIncorrect = 0;

	for (unsigned int i_query = 0; i_query < sites.size(); ++i_query) {

		std::vector<Point2D<int> > sub_sites;
		std::vector<int> sub_classes;
		for (unsigned int i_copy = 0; i_copy < sites.size(); ++i_copy) {
			if (i_copy != i_query) {
				sub_classes.push_back(classes[i_copy]);
				sub_sites.push_back(sites[i_copy]);
			}
		}

		int secKnnClass = secure_knn_classifier<Number, NumberBits>(sub_sites, classes, sites[i_query]);
		int realKnnClass = real_knn_classifier(sub_sites, sub_classes, sites[i_query]);

		std::cout << "Secure KNN classifier classified as: " << secKnnClass << std::endl;
		std::cout << "Original KNN classifier classified as: " << realKnnClass << std::endl;

		if (secKnnClass == -1) {
			++secClassificationFailed;
			std::cout << "SecKnn failed to classify" << std::endl;
		} else {
			if (realKnnClass == secKnnClass) {
				++match;
			} else {
				std::cout << "real KNN and secure KNN mismatch" << std::endl;
				++mismatch;
			}
			if (realKnnClass == classes[i_query]) { ++realCorrect; } else { ++realIncorrect; }
			if (secKnnClass == classes[i_query]) { ++secureCorrect; } else { ++secureIncorrect; }
		}

		if (match + mismatch > 0) {
			std::cout << "OUTPUT [1]: " << "matched: " << match << " out of " << (match+mismatch) << " = " << ((int)100*match/(match+mismatch)) << "%" << std::endl;
			std::cout << "OUTPUT [2]: " << "correct secure: " << secureCorrect << " out of " << (secureCorrect+secureIncorrect) << " = " << ((int)100*secureCorrect/(secureCorrect+secureIncorrect)) << "%" << std::endl;
			std::cout << "OUTPUT [3]: " << "correct real: " << realCorrect << " out of " << (realCorrect+realIncorrect) << " = " << ((int)100*realCorrect/(realCorrect+realIncorrect)) << "%" << std::endl;
		}
		std::cout << "OUTPUT [4]: " << "classification failed: " << secClassificationFailed << std::endl;

		std::cout << "ITERATIONS: average iterations needed: " << (avgIterations / (i_query+1)) << std::endl;
	}

}


#endif
