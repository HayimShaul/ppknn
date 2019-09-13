#ifndef ___GET_PERCENTILE_FRAMEWORK___
#define ___GET_PERCENTILE_FRAMEWORK___

#include <sys/resource.h>
#include <sys/types.h>
#include <sys/time.h>
#include <vector>
#include <thread>

#include <time_measurements.h>

#define USE_SIMD

#include "polynomial.h"
#include "get_percentile.h"
#include "special_polynomials.h"
#include "input_database.h"
#include "configuration.h"


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
Number fold(Number n, unsigned int count) {
	if (count == 0)
		count = n.simd_factor();
 
#ifdef RECURSIVE_FOLDING
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
#else
	Number ret = n;

	ThreadPool threads;
	std::mutex access_ret;

	if (n.simd_factor() < count*2) {
		for (unsigned int i = 1; i < n.simd_factor(); ++i) {
			while (!threads.has_free_cpu())
				threads.process_jobs(1);

			threads.submit_job(std::function<void(void)>([&n, &ret, &access_ret, i](){
				Number a = n.rotate_left(i);
				access_ret.lock();
				ret += a;
				access_ret.unlock();
			}));
		}
	} else {
		for (unsigned int i = 1; i < count; ++i) {
			while (!threads.has_free_cpu())
				threads.process_jobs(1);

			threads.submit_job(std::function<void(void)>([&n, &ret, &access_ret, i](){
				Number a = n.shift_left(i);
				access_ret.lock();
				ret += a;
				access_ret.unlock();

				a = n.shift_right(i);
				access_ret.lock();
				ret += a;
				access_ret.unlock();
			}));
		}
	}

	threads.process_jobs();

	return ret;
#endif
}

int popcnt(int i) {
	i = ((i >> 1) & 0x5555) + (i & 0x5555);
	i = ((i >> 2) & 0x3333) + (i & 0x3333);
	i = ((i >> 4) & 0x0f0f) + (i & 0x0f0f);
	i = ((i >> 8) & 0x00ff) + (i & 0x00ff);
	return i;
}


int argmax(const std::vector<int> &v) {
	int a = 0;
	for (unsigned int i = 1; i < v.size(); ++i)
		if (v[i] > v[a])
			a = i;
	return a;
}

//template<class Number>
//Number fold_power_2(Number n, int simd_factor) {
//#	ifdef DONT_USE_SIMD
//	return n;
//#	endif
//
//	assert(popcnt(simd_factor) == 1);
//
//	Number ret;
//
//	Number *a = &n;
//	Number *b = &ret;
//
//	int iter = simd_factor;
//	while (iter > 1) {
//		std::vector<long int> mask_left;
//		std::vector<long int> mask_right;
//
//		for (int i = 0; i < simd_factor; ++i) {
//			if ((i / (iter/2)) & 1) {
//				mask_left[i] = 0;
//				mask_rightt[i] = 1;
//			} else
//				mask_left[i] = 1;
//				mask_rightt[i] = 0;
//			}
//		}
//
//		*b = (a->shift_left(iter/2) * Number::from_vector(mask_right)) + (a->shift_right(iter/2) * Number::from_vector(mask_left));
//
//		iter /= 2;
//
//		Number *temp = a;
//		a = b;
//		b = temp;
//	}
//
//	return *b;
//}

template<class Number, class NumberBits>
void multithreaded_average(Distances<Number, NumberBits> &_x, std::function<int(int)> f_m, Number &output) {
	std::thread *thread = new std::thread[thread_num];
	AverageLiphe<Number, Number, ComparePoly<Number>, NoConversion<Number>> avg;
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
#				ifdef USE_SIMD
				avg.add_simd(xi, &mutex);
#				else
				avg.add(xi, &mutex);
#				endif

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
	AverageLiphe<Number, Number, ComparePoly<Number>, NoConversion<Number> > avgAvg;
	AverageLiphe<Number, Number, ComparePoly<Number>, NoConversion<Number> > avgAvgSqrMsd;
	AverageLiphe<Number, Number, ComparePoly<Number>, NoConversion<Number> > avgAvgSqrLsd;

	std::mutex mutexAvg;
	std::mutex mutexAvgSqrMsd;
	std::mutex mutexAvgSqrLsd;

	avgAvg.set_n(_x.size());
	avgAvg.set_max_sample(avg.p() / 2);
	avgAvg.compute_resample_constant();

	avgAvgSqrMsd.set_n(_x.size());
	avgAvgSqrMsd.set_max_sample(sqrMsd.p() / 2);
	avgAvgSqrMsd.set_f_m( [](int m)->int{ return ::sqrt(m*p); } );
	avgAvgSqrMsd.compute_resample_constant();

	avgAvgSqrLsd.set_n(_x.size());
	avgAvgSqrLsd.set_max_sample(sqrLsd.p() / 2);
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

			std::shared_ptr<Number> xi(new Number(*input));

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
}


#define MULTI_THREADED

inline float sqr(float x) { return x*x; }

void print_histogram(const std::vector<Point<int> > &sites, const std::vector<int> &classes, const Point<int> &query) {
	int avg = 0;
	for (auto i = sites.begin(); i != sites.end(); ++i) {
		int dist = ((*i) - query).normL1();
		avg += dist;
	}
	avg /= sites.size();

	int sigma = 0;
	for (auto i = sites.begin(); i != sites.end(); ++i) {
		int dist = ((*i) - query).normL1();
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
		int dist = ((*i) - query).normL1();
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

int real_knn_classifier(const std::vector<Point<int> > &sites, const std::vector<int> &classes, const Point<int> &query, int k = -1) {
	print_histogram(sites, classes, query);

	std::vector<float> distances;
	for (auto i_sites = sites.begin(); i_sites != sites.end(); ++i_sites) {
		float dist = (*i_sites - query).normL2sqr();
		distances.push_back(dist);
	}
	std::sort(distances.begin(), distances.end());

	if (k == -1)
		k = distances.size() * 0.02;
	float threshold = distances[k];

	std::vector<int> classCount(Configuration::classNumber, 0);

	for (unsigned int i_sites = 0; i_sites < sites.size(); ++i_sites) {
		float dist = (sites[i_sites] - query).normL2sqr();
		if (dist < threshold) {
			++classCount[ classes[i_sites] ];
		}
	}
//	int counted = std::accumulate(classCount.begin(), classCount.end());
	for (unsigned int i_sites = 0; i_sites < sites.size(); ++i_sites) {
		float dist = (sites[i_sites] - query).normL2sqr();
		if ((dist == threshold) && (std::accumulate(classCount.begin(), classCount.end(), 0) < k)) {
			++classCount[ classes[i_sites] ];
		}
	}

	int maxClass = 0;
	for (unsigned int i = 0; i < classCount.size(); ++i) {
		std::cout << "real count of class " << i << ": " << classCount[i] << std::endl;
		if (classCount[maxClass] < classCount[i])
			maxClass = i;
	}

	return maxClass;
}


bool OK = true;
int MAX_CANDIDATES = -1;

int foldMod(const std::vector<long int> &a, int mod) {
	int ret = 0;
	for (unsigned int i = 0; i < a.size(); ++i) {
		ret = (ret + a[i]) % mod;
	}
	return ret;
}

std::vector<long int> random_vector(unsigned int size, long mod) {
	std::vector<long int> ret(size);
	for (unsigned int i = 0; i < size; ++i)
		ret[i] = random() % mod;
	return ret;
}

// This function masks the classCount ciphertext vector, sends the masked version to the client who then decrypts folds, encrypts again and sends back to the server
// The server then removes the mask
// classCountEnc has the structure classCountEnc[i_candidate][i_class] = 1/0 with SIMD slot s, if site s has class i_class for candidate i_candidate
// foldedClassCountEnc has the structure foldedClassCountEnc[i_class] = n with SIMD slot s, then class i_class and candidate s has n neighbors
template<class Number>
void maskFoldRecrypt(std::vector<std::vector<Number> > &classCountEnc, std::vector<Number> &foldedClassCountEnc) {
	unsigned int candidateNum = classCountEnc.size();
	unsigned int classNum = classCountEnc[0].size();

	// mask
	std::vector<std::vector<std::vector<long int> > > mask;
	mask.resize(candidateNum);
	for (unsigned int i_candidate = 0; i_candidate < candidateNum; ++i_candidate) {
		mask[i_candidate].resize(classNum);
		for (unsigned int i_class = 0; i_class < classNum; ++i_class)
			mask[i_candidate][i_class] = random_vector(Number::simd_factor(), Number::get_global_ring_size());
	}

	for (unsigned int i_candidate = 0; i_candidate < candidateNum; ++i_candidate) {
		for (unsigned int i_class = 0; i_class < classNum; ++i_class)
			classCountEnc[i_candidate][i_class] += mask[i_candidate][i_class];
	}

	//////////////////////////////////////////////
	// send the masked classCountEnc to the client
	//////////////////////////////////////////////

	// decrypt fold and encrypt the classCount vector
	std::vector<std::vector<long int> > classCount;
	classCount.resize(candidateNum);
	for (unsigned int i_candidate = 0; i_candidate < candidateNum; ++i_candidate) {
		classCount[i_candidate].resize(classNum);
		for (unsigned int i_class = 0; i_class < classNum; ++i_class) {
			// decrypt
			std::vector<long int> v = classCountEnc[i_candidate][i_class].to_vector();
			// fold
			classCount[i_candidate][i_class] = foldMod(v, Number::get_global_ring_size());
		}
	}
	// encrypt
	foldedClassCountEnc.resize(classNum);
	for (unsigned int i_class = 0; i_class < classNum; ++i_class) {
		std::vector<long int> v(candidateNum);
		for (unsigned int i_candidate; i_candidate < candidateNum; ++i_candidate)
			v[i_candidate] = classCount[i_candidate][i_class];
		foldedClassCountEnc[i_class].from_vector(v);
	}

	//////////////////////////////////////////////
	// send the masked classCountEnc to the client
	//////////////////////////////////////////////

	// fold and transpose the mask matrix
	std::vector<std::vector<long int> > foldedMask;
	foldedMask.resize(classNum);
	for (unsigned int i_class = 0; i_class < classNum; ++i_class) {
		foldedMask[i_class].resize(candidateNum, 0);
		for (unsigned int i_candidate = 0; i_candidate < candidateNum; ++i_candidate) {
			foldedMask[i_class][i_candidate] = foldMod(mask[i_candidate][i_class], Number::get_global_ring_size());
		}
	}

	// unmask
	for (unsigned int i_class = 0; i_class < classNum; ++i_class) {
		foldedClassCountEnc[i_class] -= foldedMask[i_class];
	}
}

// classesCountVector[i_candidate][i_class] = how many neighbors of class i_class we found for candidate i_candidate
template<class Number, class NumberBits>
void old_secure_knn_classifier_gaussian(const std::vector<Point<int> > &sites, const std::vector<int> &classes, const Point<int> &query, std::vector<std::vector<int> > &classesCountVector) {

	std::cout << "Starting classifier KNN" << std::endl;

	Distances<Number, NumberBits> distances(sites, classes, query);

	Number avgValue;
	Number avgSqrValue;

	////////////////////////////////////
	// Start computing average and average of squares
	////////////////////////////////////

	Number avgEnc;
	Number avgSqrMsdEnc;
	Number avgSqrLsdEnc;

	{
		AutoTakeTimes tt("computing averages");

		global_timer.start();
		multithreaded_averages<Number, NumberBits>(distances, avgEnc, avgSqrMsdEnc, avgSqrLsdEnc);
		std::cout << global_timer.end("computing averages");
	}


	{
		AutoTakeTimes tt("folding vectors");

		global_timer.start();
		avgEnc = fold(avgEnc, sites.size());
		avgSqrMsdEnc = fold(avgSqrMsdEnc, sites.size());
		avgSqrLsdEnc = fold(avgSqrLsdEnc, sites.size());
		std::cout << global_timer.end("folding");
	}


	////////////////////////////////////
	// End computing average and average of squares
	////////////////////////////////////



	////////////////////////////////////
	// start computing threshold
	////////////////////////////////////

	std::vector<Number> thresholdCandidates;
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

		long inv2 = power_mod(2, phi(Number::get_global_ring_size()) - 1, Number::get_global_ring_size());
		std::cout << "2^{-1} mod " << Number::get_global_ring_size() << " = " << inv2 << std::endl;

		{
			{
				Number threshold = avgEnc;
				thresholdCandidates.push_back(threshold);
			}
			{
				Number threshold = avgEnc - sigma;
				thresholdCandidates.push_back(threshold);
			}
			{
//				Number threshold = avgEnc - sigma*2;
				Number threshold = avgEnc - sigma - sigma;
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

		}


		std::cout << global_timer.end("compute threashold candidates");
		if (MAX_CANDIDATES > 0)
			thresholdCandidates.resize(MAX_CANDIDATES);





		int realAvg = 0;
		int realAvgSqr = 0;
		for (auto i = sites.begin(); i != sites.end(); ++i) {
			int dist = ((*i) - query).normL1();
			realAvg += dist;
			realAvgSqr += dist*dist;
		}
		realAvg /= sites.size();
		realAvgSqr /= sites.size();

		std::cout << "avg = " << avgEnc.to_int()  << "         (real = " << realAvg << ")" << std::endl;

		std::cout << "avg^2 = " << (avgMsdEnc.to_int()*p + avgLsdEnc.to_int()) << "        (real = " << (realAvg*realAvg) << ")" << std::endl;
		std::cout << "avg^2 / p = " << avgMsdEnc.to_int() <<  "            (real = " << (realAvg*realAvg / p) << ")" << std::endl;
		std::cout << "avg^2 % p = " << avgLsdEnc.to_int() << "             (real = " << (realAvg*realAvg % p) << ")" << std::endl;

		std::cout << "avgSqr = " << (avgSqrMsdEnc.to_int()*p + avgSqrLsdEnc.to_int()) << "          (real = " << realAvgSqr << ")" << std::endl;
		std::cout << "avgSqr / p = " << avgSqrMsdEnc.to_int() << "        (real = " << (realAvgSqr / p) << ")" << std::endl;
		std::cout << "avgSqr % p = " <<  avgSqrLsdEnc.to_int() << "        (real = " << (realAvgSqr % p) << ")" << std::endl;


		std::cout << "sigma = sqrt_msd(" << (avgSqrMsdEnc - avgMsdEnc).to_int() << ") + " << "sqrt_lsd(" << (avgSqrLsdEnc - avgLsdEnc).to_int() << ")" << std::endl;
		std::cout << "sigma = " << sigma.to_int() << "            (real = " << ::sqrt(-realAvg*realAvg + realAvgSqr) << ")" << std::endl;
		for (unsigned int i_candidate = 0; i_candidate < thresholdCandidates.size(); ++i_candidate)
			std::cout << "Threshold candidate " << i_candidate << " = " << thresholdCandidates[i_candidate].to_int() << std::endl;

	}

	////////////////////////////////////
	// end computing threshold
	////////////////////////////////////



	////////////////////////////////////
	// start counting classes
	////////////////////////////////////


	std::vector< std::vector<Number> > classCountEnc(thresholdCandidates.size());
	for (unsigned int i_candidate = 0; i_candidate < thresholdCandidates.size(); ++i_candidate) {
		classCountEnc[i_candidate].resize(Configuration::classNumber, 0);
	}

	auto dist = distances.begin();
	while (dist != distances.end()) {
		Number xi = dist.getDistances();

		{
			for (unsigned int i_candidate = 0; i_candidate < thresholdCandidates.size(); ++i_candidate) {
				global_timer.start();

				Number knnEnc = ComparePoly<Number>(xi) < thresholdCandidates[i_candidate];

				for (unsigned int i_class = 0; i_class < Configuration::classNumber; ++i_class) {
					Number knnClassEnc = knnEnc * dist.getClass(i_class);
					classCountEnc[i_candidate][i_class] += knnClassEnc;
				}

				std::cout << global_timer.end("count classes for one candidate for one batch of points");

				// debugging
//				std::cout << "The KNN indicator vector for candidate " << thresholdCandidates[i_candidate].to_int() << std::endl;
//				std::vector<long int> knn = knnEnc.to_vector();
//				std::vector<std::vector<long int> > knnClasses(Configuration::classNumber);
//				for (unsigned int i_class = 0; i_class < Configuration::classNumber; ++i_class)
//					knnClasses[i_class] = knnClassEnc[i_class].to_vector();
//				std::vector<long int> plaintextDistances = dist.getPlaintextDistances();
//				std::vector<long int> plaintextClasses = dist.getPlaintextClasses();
//				for (unsigned int i = 0; (i < knn.size()) && (i < sites.size()); ++i) {
//					std::cout << dist.loc() << ", " << i << ") " << "x = " << knn[i] << "   dist = " << plaintextDistances[i] << "    class = " << plaintextClasses[i];
//					int numberOfClasses = 0;
//					for (unsigned int i_class = 0; i_class < Configuration::classNumber; ++i_class) {
//						std::cout << " class " << i_class << ": " << knnClass[i_candidate][i];
//						numberOfClasses += knnClass[i_class][i];
//					}
//					std::cout << std::endl;
//					if (numberOfClasses != 1)
//						OK = false;
//				}
			}
		}

		++dist;
	}

//	std::cout << "depth = mul " << classOneCountEnc0.mul_depth() << " add " << classOneCountEnc0.add_depth() << std::endl;

	if (OK)
		std::cout << "test is OK" << std::endl;


	classesCountVector.resize(thresholdCandidates.size());
	for (unsigned int i_candidate = 0; i_candidate < thresholdCandidates.size(); ++i_candidate) {
		std::cout << "Candidate " << i_candidate << std::endl;

		classesCountVector[i_candidate].resize(Configuration::classNumber);
		for (unsigned int i_class = 0; i_class < Configuration::classNumber; ++i_class) {
			classesCountVector[i_candidate][i_class] = 0;
			std::vector<long int> classCountV = classCountEnc[i_candidate][i_class].to_vector();

			std::cout << "   candidate " << i_candidate << " class " << i_class << " neighbors:" << std::endl;
			for (unsigned int n = 0; n < classCountV.size(); ++n) {
				std::cout << "       " << classCountV[n] << std::endl;
				classesCountVector[i_candidate][i_class] += classCountV[n];
			}
		}
	}


	std::cout << " Total of:" << std::endl;
	for (unsigned int i_candidate = 0; i_candidate < thresholdCandidates.size(); ++i_candidate) {
		for (unsigned int i_class = 0; i_class < Configuration::classNumber; ++i_class) {
			for (unsigned int i_class = 0; i_class < Configuration::classNumber; ++i_class) {
				std::cout << "candidate " << i_candidate << " class " << i_class << ": " << classesCountVector[i_candidate][i_class] << " neighbors" << std::endl;
			}
		}
	}

	print_stat(-1);
//	return secureKnnClassifier;
}

// classesCountVector[i_candidate][i_class] = how many neighbors of class i_class we found for candidate i_candidate
template<class Number, class NumberBits>
int secure_knn_classifier_gaussian(const std::vector<Point<int> > &sites, const std::vector<int> &classes, int k, const Point<int> &query) {

	std::cout << "Starting classifier KNN" << std::endl;

	Distances<Number, NumberBits> distances(sites, classes, query);

	Number avgValue;
	Number avgSqrValue;

	////////////////////////////////////
	// Start computing average and average of squares
	////////////////////////////////////

	Number avgEnc;
	Number avgSqrMsdEnc;
	Number avgSqrLsdEnc;

	{
		AutoTakeTimes tt("computing averages");

		global_timer.start();
		multithreaded_averages<Number, NumberBits>(distances, avgEnc, avgSqrMsdEnc, avgSqrLsdEnc);
		std::cout << global_timer.end("computing averages");
	}


	{
		AutoTakeTimes tt("folding vectors");

		global_timer.start();
		avgEnc = fold(avgEnc, sites.size());
		avgSqrMsdEnc = fold(avgSqrMsdEnc, sites.size());
		avgSqrLsdEnc = fold(avgSqrLsdEnc, sites.size());
		std::cout << global_timer.end("folding");
	}


	////////////////////////////////////
	// End computing average and average of squares
	////////////////////////////////////



	////////////////////////////////////
	// start computing threshold
	////////////////////////////////////

	std::vector<Number> thresholdCandidates;
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

		long inv2 = power_mod(2, phi(Number::get_global_ring_size()) - 1, Number::get_global_ring_size());
		std::cout << "2^{-1} mod " << Number::get_global_ring_size() << " = " << inv2 << std::endl;

		{
			{
				Number threshold = avgEnc;
				thresholdCandidates.push_back(threshold);
			}
			{
				Number threshold = avgEnc - sigma;
				thresholdCandidates.push_back(threshold);
			}
			{
//				Number threshold = avgEnc - sigma*2;
				Number threshold = avgEnc - sigma - sigma;
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

		}


		std::cout << global_timer.end("compute threashold candidates");
		if (MAX_CANDIDATES > 0)
			thresholdCandidates.resize(MAX_CANDIDATES);





		int realAvg = 0;
		int realAvgSqr = 0;
		for (auto i = sites.begin(); i != sites.end(); ++i) {
			int dist = ((*i) - query).normL1();
			realAvg += dist;
			realAvgSqr += dist*dist;
		}
		realAvg /= sites.size();
		realAvgSqr /= sites.size();

		std::cout << "avg = " << avgEnc.to_int()  << "         (real = " << realAvg << ")" << std::endl;

		std::cout << "avg^2 = " << (avgMsdEnc.to_int()*p + avgLsdEnc.to_int()) << "        (real = " << (realAvg*realAvg) << ")" << std::endl;
		std::cout << "avg^2 / p = " << avgMsdEnc.to_int() <<  "            (real = " << (realAvg*realAvg / p) << ")" << std::endl;
		std::cout << "avg^2 % p = " << avgLsdEnc.to_int() << "             (real = " << (realAvg*realAvg % p) << ")" << std::endl;

		std::cout << "avgSqr = " << (avgSqrMsdEnc.to_int()*p + avgSqrLsdEnc.to_int()) << "          (real = " << realAvgSqr << ")" << std::endl;
		std::cout << "avgSqr / p = " << avgSqrMsdEnc.to_int() << "        (real = " << (realAvgSqr / p) << ")" << std::endl;
		std::cout << "avgSqr % p = " <<  avgSqrLsdEnc.to_int() << "        (real = " << (realAvgSqr % p) << ")" << std::endl;


		std::cout << "sigma = sqrt_msd(" << (avgSqrMsdEnc - avgMsdEnc).to_int() << ") + " << "sqrt_lsd(" << (avgSqrLsdEnc - avgLsdEnc).to_int() << ")" << std::endl;
		std::cout << "sigma = " << sigma.to_int() << "            (real = " << ::sqrt(-realAvg*realAvg + realAvgSqr) << ")" << std::endl;
		for (unsigned int i_candidate = 0; i_candidate < thresholdCandidates.size(); ++i_candidate)
			std::cout << "Threshold candidate " << i_candidate << " = " << thresholdCandidates[i_candidate].to_int() << std::endl;

	}

	////////////////////////////////////
	// end computing threshold
	////////////////////////////////////



	////////////////////////////////////
	// start counting classes
	////////////////////////////////////


	std::vector< std::vector<Number> > classCountEnc(thresholdCandidates.size());
	for (unsigned int i_candidate = 0; i_candidate < thresholdCandidates.size(); ++i_candidate) {
		classCountEnc[i_candidate].resize(Configuration::classNumber, 0);
	}

	auto dist = distances.begin();
	while (dist != distances.end()) {
		Number xi = dist.getDistances();

		{
			for (unsigned int i_candidate = 0; i_candidate < thresholdCandidates.size(); ++i_candidate) {
				global_timer.start();

				Number knnEnc = ComparePoly<Number>(xi) < thresholdCandidates[i_candidate];

				for (unsigned int i_class = 0; i_class < Configuration::classNumber; ++i_class) {
					Number knnClassEnc = knnEnc * dist.getClass(i_class);
					classCountEnc[i_candidate][i_class] += knnClassEnc;
				}

				std::cout << global_timer.end("count classes for one candidate for one batch of points");

				// debugging
//				std::cout << "The KNN indicator vector for candidate " << thresholdCandidates[i_candidate].to_int() << std::endl;
//				std::vector<long int> knn = knnEnc.to_vector();
//				std::vector<std::vector<long int> > knnClasses(Configuration::classNumber);
//				for (unsigned int i_class = 0; i_class < Configuration::classNumber; ++i_class)
//					knnClasses[i_class] = knnClassEnc[i_class].to_vector();
//				std::vector<long int> plaintextDistances = dist.getPlaintextDistances();
//				std::vector<long int> plaintextClasses = dist.getPlaintextClasses();
//				for (unsigned int i = 0; (i < knn.size()) && (i < sites.size()); ++i) {
//					std::cout << dist.loc() << ", " << i << ") " << "x = " << knn[i] << "   dist = " << plaintextDistances[i] << "    class = " << plaintextClasses[i];
//					int numberOfClasses = 0;
//					for (unsigned int i_class = 0; i_class < Configuration::classNumber; ++i_class) {
//						std::cout << " class " << i_class << ": " << knnClass[i_candidate][i];
//						numberOfClasses += knnClass[i_class][i];
//					}
//					std::cout << std::endl;
//					if (numberOfClasses != 1)
//						OK = false;
//				}
			}
		}

		++dist;
	}

	// classCountEnc has the structure classCountEnc[i_candidate][i_class] = 1/0 with SIMD slot s, if site s has class i_class for candidate i_candidate
	// foldedClassCountEnc has the structure foldedClassCountEnc[i_class] = n with SIMD slot s, then class i_class and candidate s has n neighbors
	std::vector<Number> foldedClassCountEnc;
	{
		AutoTakeTimes("mask fold and recrypt");
		maskFoldRecrypt(classCountEnc, foldedClassCountEnc);
	}

	std::vector<long int> vec1(Number::simd_factor(), 1);

	Number classEnc;
	{
		AutoTakeTimes("compute Max");

		Number totalCount = 0;
		for (unsigned int i_class = 0; i_class < foldedClassCountEnc.size(); ++i_class)
			totalCount += foldedClassCountEnc[i_class];

		Number tooFew = ComparePoly<Number>(totalCount) < (k/2);
		Number tooMany = ComparePoly<Number>(totalCount) > (2*k);


		Number tooManyClass = 0;
		 	// else if (totalCount > 2*k) 
		for (unsigned int i_class = 0; i_class < foldedClassCountEnc.size(); ++i_class) {
			tooManyClass += 
				// If we can reduce the "neighbors sphere" until we are left with 2*k s.t. one class has more than treshold neighbors
				(ComparePoly<Number>(totalCount - foldedClassCountEnc[i_class]) < k) * (1+i_class);
		}


		std::vector<Number> classIsMax(foldedClassCountEnc.size());
		for (unsigned int i_class = 0; i_class < foldedClassCountEnc.size(); ++i_class) {
std::cerr << "check that addConstant considers simd" << std::endl;
			classIsMax[i_class] = (-tooFew + 1)*(i_class + 1);
		}

		for (unsigned int i_class = 0; i_class < foldedClassCountEnc.size(); ++i_class) {
			for (unsigned int j_class = i_class + 1; j_class < foldedClassCountEnc.size(); ++j_class) {
				Number isIBiggerJ = ComparePoly<Number>(foldedClassCountEnc[i_class] - foldedClassCountEnc[j_class]) < (2*k);
				classIsMax[i_class] *= isIBiggerJ;
				classIsMax[j_class] *= (-isIBiggerJ + 1);
			}
		}

		Number inRangeClass = 0;
		for (unsigned int i_class = 0; i_class < foldedClassCountEnc.size(); ++i_class) {
			inRangeClass += classIsMax[i_class];
		}

		classEnc = tooMany*tooManyClass + (-tooMany+1)*inRangeClass;
	}

	print_stat(-1);

	// Send to the user who decrypts and figures out the class

	std::vector<long int> classVector = classEnc.to_vector();

	int ret = -1;
	for (unsigned int i_candidate = 0; i_candidate < classVector.size(); ++i_candidate) {
		if ((classVector[i_candidate] != 0) && (ret != -1))
			ret = classVector[i_candidate] - 1;
	}

	// debugging
	if ((ret < -1) || (ret > Configuration::classNumber))
		OK = false;

	return ret;
}


int avgIterations = 0;
int RETRIES = 5;

template<class Number, class NumberBits>
int secure_knn_classifier(const std::vector<Point<int> > &sites, const std::vector<int> &classes, const Point<int> &query) {
//	std::vector<std::vector<int> > classesCountVector;

	int k = 0.05 * sites.size();
	for (int i = 0; i < RETRIES; ++i) {
		++avgIterations;
//		classesCountVector.resize(0);
		int prediction = secure_knn_classifier_gaussian<Number, NumberBits>(sites, classes, k, query);
		if (prediction != -1)
			return prediction;

//		for (unsigned int i_candidate = 0; i_candidate < classesCountVector.size(); ++i_candidate) {
//			std::cout << "Iteration " << i << " candidate " << i_candidate << ":" << std::endl;
//			int totalCount = std::accumulate(classesCountVector[i_candidate].begin(), classesCountVector[i_candidate].end(), 0);
//
//			float threshold = 0.05 * sites.size();
//
//			if (totalCount < 0.5 * threshold) {
//				std::cout << "Too little neighbors" << std::endl;
//				continue;
//			}
//			// Deal with the case where we have a small count of 1 but huge amount of 0
//			if (totalCount > 2 * threshold) {
//				for (unsigned int i_class = 0; i_class < classesCountVector[i_candidate].size(); ++i_class) {
//					// If we can reduce the "neighbors sphere" until we are left with 2*threshold s.t. one class has more than treshold neighbors
//					if (totalCount - classesCountVector[i_candidate][i_class] < threshold) {
//						std::cout << "too many neighbors, but class " << i_class << " is considerably big. Classifying as " << i_class << std::endl;
//						return i_class;
//					}
//				}
//
//				std::cout << "Too many neighbors" << std::endl;
//				continue;
//			}
//			std::cout << "enough neighbors. Classifying by majority" << std::endl;
//			return argmax(classesCountVector[i_candidate]);
//		}
	}

	return 0;
//	return -1;
}

template<class Number, class NumberBits>
void test_secure_knn_classifier(const std::vector<Point<int> > &sites, const std::vector<int> &classes) {
	int match = 0;
	int mismatch = 0;
	int secClassificationFailed = 0;

	int secureCorrect[2] = {0};
	int secureIncorrect[2] = {0};

	int realCorrect[2] = {0};
	int realIncorrect[2] = {0};

	for (unsigned int i_query = 0; i_query < sites.size(); ++i_query) {

		std::vector<Point<int> > sub_sites;
		std::vector<int> sub_classes;
		for (unsigned int i_copy = 0; i_copy < sites.size(); ++i_copy) {
			if (i_copy != i_query) {
				sub_classes.push_back(classes[i_copy]);
				sub_sites.push_back(sites[i_copy]);
			}
		}

		int secKnnClass = secure_knn_classifier<Number, NumberBits>(sub_sites, classes, sites[i_query]);
		int realKnnClass = real_knn_classifier(sub_sites, sub_classes, sites[i_query]);

		if (OK)
			std::cout << "test is OK";

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
			if (realKnnClass == classes[i_query]) { ++realCorrect[classes[i_query]]; } else { ++realIncorrect[classes[i_query]]; }
			if (secKnnClass == classes[i_query]) { ++secureCorrect[classes[i_query]]; } else { ++secureIncorrect[classes[i_query]]; }
		}

		if (match + mismatch > 0)
			std::cout << "OUTPUT [1]: " << "matched: " << match << " out of " << (match+mismatch) << " = " << ((int)100*match/(match+mismatch)) << "%" << std::endl;
		std::cout << "OUTPUT [2]: " << "correct secure: class0 = " << secureCorrect[0] << " class1 = " << secureCorrect[1] << std::endl;
		std::cout << "OUTPUT [3]: " << "incorrect secure: class0 = " << secureIncorrect[0] << " class1 = " << secureIncorrect[1] << std::endl;
		std::cout << "OUTPUT [4]: " << "total: " << (secureCorrect[0]+secureCorrect[1]+secureIncorrect[0]+secureIncorrect[1]) << std::endl;

		if (secureCorrect[0]+secureCorrect[1]+secureIncorrect[0]+secureIncorrect[1] > 0)
			std::cout << "OUTPUT [5]: " << "secure ratio: " << ((int)100*(secureCorrect[0]+secureCorrect[1])/(secureCorrect[0]+secureCorrect[1]+secureIncorrect[0]+secureIncorrect[1])) << "%" << std::endl;
		if (secureCorrect[1] + secureIncorrect[1])
			std::cout << "OUTPUT [6]: " << "secure Dice score1: " << (2.0*secureCorrect[1] / (secureCorrect[1] + secureIncorrect[1])) << std::endl;
		if (secureCorrect[0] + secureIncorrect[0])
			std::cout << "OUTPUT [7]: " << "secure Dice score0: " << (2.0*secureCorrect[0] / (secureCorrect[0] + secureIncorrect[0])) << std::endl;




		std::cout << "OUTPUT [8]: " << "correct real: class0 = " << realCorrect[0] << " class1 = " << realCorrect[1] << std::endl;
		std::cout << "OUTPUT [9]: " << "incorrect real: class0 = " << realIncorrect[0] << " class1 = " << realIncorrect[1] << std::endl;
		std::cout << "OUTPUT [10]: " << "total: " << (realCorrect[0]+realCorrect[1]+realIncorrect[0]+realIncorrect[1]) << std::endl;
		if (realCorrect[0]+realCorrect[1]+realIncorrect[0]+realIncorrect[1] > 0)
			std::cout << "OUTPUT [11]: " << "insecure ratio: " << ((int)100*(realCorrect[0]+realCorrect[1])/(realCorrect[0]+realCorrect[1]+realIncorrect[0]+realIncorrect[1])) << "%" << std::endl;
		if (realCorrect[1] + realIncorrect[1] > 0)
			std::cout << "OUTPUT [12]: " << "insecure Dice score1: " << (2.0*realCorrect[1] / (realCorrect[1] + realIncorrect[1])) << std::endl;
		if (realCorrect[0] + realIncorrect[0] > 0)
			std::cout << "OUTPUT [13]: " << "insecure Dice score0: " << (2.0*realCorrect[0] / (realCorrect[0] + realIncorrect[0])) << std::endl;

		std::cout << "OUTPUT [14]: " << "classification failed: " << secClassificationFailed << std::endl;

		std::cout << "ITERATIONS: average iterations needed: " << (avgIterations / (i_query+1)) << std::endl;
	}

	if (OK)
		std::cout << "test is OK";

}

void test_kish_classifier(const std::vector<Point<int> > &sites, const std::vector<int> &classes) {
	unsigned int MIN_K = sites.size() * 0.02;
	unsigned int MAX_K = MIN_K + 1;

	for (unsigned int k = MIN_K; k < MAX_K; ++k) {

//		int correct[2] = {0};
//		int incorrect[2] = {0};

		std::vector<float> f1Scores(Configuration::classNumber, 0);
		std::vector<unsigned int> truePositives(Configuration::classNumber, 0);
		std::vector<unsigned int> falsePositives(Configuration::classNumber, 0);
		std::vector<unsigned int> trueNegatives(Configuration::classNumber, 0);
		std::vector<unsigned int> falseNegatives(Configuration::classNumber, 0);
		std::vector<unsigned int> weights(Configuration::classNumber, 0);

		for (unsigned int i_query = 0; i_query < sites.size(); ++i_query) {
			std::vector<Point<int> > sub_sites;
			std::vector<int> sub_classes;
			std::set<int> exclude;
			exclude.insert(i_query);
			while (exclude.size() < Configuration::cv_k_fold) {
				exclude.insert(random() % sites.size());
			}

			for (unsigned int i_copy = 0; i_copy < sites.size(); ++i_copy) {
				if (exclude.find(i_copy) == exclude.end()) {
					sub_classes.push_back(classes[i_copy]);
					sub_sites.push_back(sites[i_copy]);
				}
			}

			int knnClass = real_knn_classifier(sub_sites, sub_classes, sites[i_query], k);

			++weights[classes[i_query]];
			for (int i_class = 0; i_class < (int)Configuration::classNumber; ++i_class) {
				if ((knnClass == i_class) && (classes[i_query] == i_class))
					++truePositives[i_class];
				if ((knnClass != i_class) && (classes[i_query] == i_class))
					++falseNegatives[i_class];
				if ((knnClass == i_class) && (classes[i_query] != i_class))
					++falsePositives[i_class];
				if ((knnClass != i_class) && (classes[i_query] != i_class))
					++trueNegatives[i_class];
			}
		}

		std::vector<float> precision(Configuration::classNumber, 0);
		std::vector<float> recall(Configuration::classNumber, 0);

		std::cout << "KNN k = " << k << std::endl;
		for (unsigned int i_class = 0; i_class < Configuration::classNumber; ++i_class) {
			precision[i_class] = 1.0 - (1.0 * falsePositives[i_class]) / (truePositives[i_class] + falsePositives[i_class]);
			recall[i_class] = 1.0 - (1.0 * falseNegatives[i_class]) / (truePositives[i_class] + falseNegatives[i_class]);
			f1Scores[i_class] = 2.0 * precision[i_class] * recall[i_class] / (precision[i_class] + recall[i_class]);
			std::cout << "KISH OUTPUT: " << k << " Dice score " << i_class << ": " << f1Scores[i_class] << std::endl;
		}


	}

	if (OK)
		std::cout << "test is OK";

}

#endif
