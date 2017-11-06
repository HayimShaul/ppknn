#ifndef ___GET_PERCENTILE_FRAMEWORK___
#define ___GET_PERCENTILE_FRAMEWORK___

#include <sys/resource.h>
#include <sys/types.h>
#include <vector>
#include <thread>

#include "polynomial.h"
#include "get_percentile.h"
#include "special_polynomials.h"
#include "input_database.h"



time_t start = 0;
time_t total = 0;

void start_stat_interval() {
	start = clock();
}

void end_stat_interval() {
	total += clock() - start;
}

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
	std::cout << "Took " << total << " micro" << std::endl;

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

	avg = avgAvg.getAverage();
	sqrMsd = avgAvgSqrMsd.getAverage();
	sqrLsd = avgAvgSqrLsd.getAverage();

	delete[] thread;
}


#define MULTI_THREADED


template<class Number, class NumberBits>
void secure_geo_search(const std::vector<Point2D<int> > &sites, const Point2D<int> &query) {
	Distances<Number, NumberBits> distances(sites, query);
	unsigned int i = 0;

	Number avgValue;
	Number avgSqrValue;

	n = sites.size();

	////////////////////////////////////
	// Start computing average and average of squares
	////////////////////////////////////

#ifdef MULTI_THREADED
	Number avgEnc;
	Number avgSqrMsdEnc;
	Number avgSqrLsdEnc;

	multithreaded_averages<Number, NumberBits>(distances, avgEnc, avgSqrMsdEnc, avgSqrLsdEnc);

	avgEnc = fold(avgEnc);
	avgSqrMsdEnc = fold(avgSqrMsdEnc);
	avgSqrLsdEnc = fold(avgSqrLsdEnc);

#else
#	error no more single thread
//	AverageLiphe<Number, NumberBits, CompareNative<NumberBits>, TruncConversion> avg(_x.size());
//	avg.set_max_sample(*(std::max_element(_x.begin(), _x.end())));
//	avg.compute_resample_constant(0.1, 0.05);
//
//	AverageLiphe<Number, NumberBits, CompareNative<NumberBits>, TruncConversion> avgSqrMsd(_x.size());
//	avgSqrMsd.set_max_sample(*(std::max_element(_x.begin(), _x.end())));
//	avgSqrMsd.compute_resample_constant(0.1, 0.05);
//	avgSqrMsd.set_f_m( [](int m)->int{ return ::sqrt(m*p); } );
//
//	AverageLiphe<Number, NumberBits, CompareNative<NumberBits>, TruncConversion> avgSqrLsd(_x.size());
//	avgSqrLsd.set_max_sample(*(std::max_element(_x.begin(), _x.end())));
//	avgSqrLsd.compute_resample_constant(0.1, 0.05);
//	avgSqrLsd.set_f_m( [](int m)->int{ return ::sqrt(m); } );
//
//
//	i = 0;
//	input.set_to_begin();
//	while (!input.is_end()) {
//		NumberBits xi = *input;
//		++input;
//
//		start_stat_interval();
//
//		avg.add_simd(xi);
//		avgSqrLsd.add_simd(xi);
//		avgSqrMsd.add_simd(xi);
//
//
////		std::cout << "\n\n===================\n";
////		std::cout << "avg (raw) = " << avg.get_bits().to_bit_stream() << std::endl;
////		std::cout << "avg (bits) = " << avg.getAverage().to_bit_stream() << std::endl;
////		std::cout << "avg  = " << avg.getAverage().to_int() << std::endl;
////		std::cout << std::endl;
////		std::cout << "avgSqr (raw) = " << avgSqr.get_bits().to_bit_stream() << std::endl;
////		std::cout << "avgSqr (bits) = " << avgSqr.getAverage().to_bit_stream() << std::endl;
//
//
//		unsigned int xii = i * NumberBits::simd_factor();
//		while ((xii < _x.size()) && (xii < (i + 1) * NumberBits::simd_factor())) {
//			realAvg += _x[xii];
//			realAvgSqr += _x[xii] * _x[xii];
//			++xii;
//		}
//
//		end_stat_interval();
//
//		if ((i % statistic_rate) == statistic_rate - 1) {
//			print_stat(i);
//
//			if (measureAccuracy) {
//				Number avgTempEnc = fold(avg.getAverage());
//				Number avgSqrMsdTempEnc = fold(avgSqrMsd.getAverage());
//				Number avgSqrLsdTempEnc = fold(avgSqrLsd.getAverage());
//
//				print_detailed_stats(i, avgTempEnc, avgSqrMsdTempEnc, avgSqrLsdTempEnc);
//			}
//		}
//		++i;
//	}
////
//	start_stat_interval();
//
//	Number avgEnc = fold(avg.getAverage());
//	Number avgSqrMsdEnc = fold(avgSqrMsd.getAverage());
//	Number avgSqrLsdEnc = fold(avgSqrLsd.getAverage());
//
//	end_stat_interval();
//
#endif

	int realAvg = 0;
	int realAvgSqr = 0;
	for (auto i = sites.begin(); i != sites.end(); ++i) {
		int dist = abs((*i).x - query.x) + abs((*i).y - query.y);
		realAvg += dist;
		realAvgSqr += dist*dist;
	}
	realAvg /= sites.size();
	realAvgSqr /= sites.size();

	std::cout << "real avg = " << realAvg << std::endl;
	std::cout << "real avgSqr = " << realAvgSqr << std::endl;
	std::cout << "real avgSqr / p = " << (realAvgSqr / p) << std::endl;
	std::cout << "real avgSqr % p = " << (realAvgSqr % p) << std::endl;

	std::cout << "avg = " << avgEnc.to_int() << std::endl;
	std::cout << "avgSqr / p = " << avgSqrMsdEnc.to_int() << std::endl;
	std::cout << "avgSqr % p = " <<  avgSqrLsdEnc.to_int() << std::endl;

	////////////////////////////////////
	// End computing average and average of squares
	////////////////////////////////////



	////////////////////////////////////
	// start computing threshold
	////////////////////////////////////

	NumberBits thresholdBits;

	{
		ThreadPool threads;

		start_stat_interval();

		Number avgLsdEnc = avgEnc * avgEnc;
		Number avgMsdEnc;
		SpecialPolynomials<Number>::square_msd_polynomial.compute(avgMsdEnc, avgEnc, &threads);

//std::cerr << "square of " << avgEnc.to_int() << " = " << " lsd= " << avgLsdEnc.to_int() << " msd= " << avgMsdEnc.to_int() << std::endl;
//		Number sigma =  simulate_sqrt_polynomial( avgSqrMsdEnc - avgMsdEnc )
//							+ simulate_sqrt_polynomial( avgSqrLsdEnc - avgLsdEnc );
		Number sigma = SpecialPolynomials<Number>::sqrt_msd_polynomial.compute( avgSqrMsdEnc - avgMsdEnc, &threads );
		sigma += SpecialPolynomials<Number>::sqrt_polynomial.compute( avgSqrLsdEnc - avgLsdEnc, &threads );

		Number threshold = avgEnc;
		threshold -= sigma;

		convert_to_bits<NumberBits, Number>(thresholdBits, threshold, &threads);

		threads.process_jobs();

		end_stat_interval();

		std::cout << "sigma = sqrt_msd(" << (avgSqrMsdEnc - avgMsdEnc).to_int() << ") + " << "sqrt_lsd(" << (avgSqrLsdEnc - avgLsdEnc).to_int() << std::endl;
		std::cout << "sigma = " << sigma.to_int() << std::endl;
		std::cout << "real sigma = " <<  ( ::sqrt(realAvgSqr/n - realAvg*realAvg/(n*n)) ) << std::endl;
		std::cout << "Threshold = " << threshold.to_int() << std::endl;
		std::cout << "converting to bits " << threshold.to_int() << " turned into " << thresholdBits.to_int() << std::endl;

	}

	////////////////////////////////////
	// end computing threshold
	////////////////////////////////////


	print_stat(-1);
	if (measureAccuracy) {
		print_detailed_stats(-1, avgEnc, avgSqrMsdEnc, avgSqrLsdEnc);
	}



	////////////////////////////////////
	// start reporting points
	////////////////////////////////////

	i = 0;
	auto dist = distances.begin();
	while (dist != distances.end()) {
		NumberBits xi = *dist;
		++dist;

		start_stat_interval();
		Number reportEnc = xi < thresholdBits;
		end_stat_interval();

		std::vector<long int> report = reportEnc.to_vector();
		for (auto ri = report.begin(); ri != report.end(); ++ri) {
			if (i < sites.size()) {
				std::cout << i << ") " << "x = " << *ri << std::endl;
				++i;
			}
		}
	}


	////////////////////////////////////
	// end reporting points
	////////////////////////////////////


	print_stat(-1);

//	report(sigmaDist, compactRepresentation)
}

#endif
