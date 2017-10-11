#ifndef ___GET_PERCENTILE_FRAMEWORK___
#define ___GET_PERCENTILE_FRAMEWORK___

#include <sys/resource.h>
#include <sys/types.h>
#include <vector>

#include "polynomial.h"
#include "get_percentile.h"
#include "special_polynomials.h"
#include "input_iterator.h"





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







template<class Number, class Bits>
class PercentileIterator {
private:
	const std::vector<int> *_raw_data;
	Bits _threshold;
	unsigned int _location;

	Number _last_value;
	unsigned int _last_value_location;

	void copy(const PercentileIterator &a) { _raw_data = a._raw_data; _threshold = a._threshold; _location = a._location; _last_value_location = a._last_value_location; _last_value = a._last_value; }
public:
	PercentileIterator(const std::vector<int> *arr, const Bits &thr) : _raw_data(arr), _threshold(thr), _location(0), _last_value(-1) {}
	PercentileIterator(const PercentileIterator &a) { copy(a); }

	PercentileIterator &operator=(const PercentileIterator &a) { copy(a); return *this; }

	void begin() { _location = 0; }
	void end() { _location = _raw_data->size(); }

	bool is_end() { return _location == _raw_data->size(); }
	bool is_begin() { return _location == 0; }

	void operator++() { ++_location; }

	bool operator==(const PercentileIterator &a) const {
		return _location == a._location;
	}
	bool operator!=(const PercentileIterator &a) const { return !operator==(a); }

	const Number &operator*() {
		if (_location != _last_value_location) {
			Bits xi = Bits::static_from_int((*_raw_data)[_location]);
			_last_value = Number(1) - (xi < _threshold);

			if (((_last_value.to_int() == 1) && ((xi.to_int() < _threshold.to_int()))) ||
				((_last_value.to_int() == 0) && (!(xi.to_int() < _threshold.to_int())))) {
					std::cout << "Error comparing " << xi.to_int() << " and " << _threshold.to_int() << std::endl;
					_last_value = Number(1) - (xi < _threshold);
			}
			_last_value_location = _location;
		}
		return _last_value;
	}
};

//template<class Number>
//Number simulate_square_msd_polynomial(const Number &x) {
//#warning this is a simulation
//	int _x = x.to_int();
//	_x *= _x;
//	_x /= p;

//	return Number(_x);
//}

//template<class Number>
//Number  simulate_sqrt_polynomial(const Number &x) {
//#warning this is a simulation
//
//	int _x = x.to_int();
//	_x = ::sqrt(_x);
//
//	return Number(_x);
//}

template<class NumberBits, class Number>
NumberBits convert_to_bits(const Number &x) {
	int bits = 0;
	while ((1 << bits) < p)
		++bits;

	NumberBits ret;
	ret.set_bit_length(bits);

	int batch_size = 0;
	Number *powers = SpecialPolynomials<Number>::convert_to_bit[0].compute_powers(x, batch_size);
	for (int i = 0; i < bits; ++i)
		ret.set_bit(i, SpecialPolynomials<Number>::convert_to_bit[i].compute(x, powers, batch_size) );

	delete[] powers;

	return ret;

//#warning this is a simulation
//
//	int _x = x.to_int();
//
//	int bits = 0;
//	while ((1 << bits) < _x)
//		++bits;
//
//	NumberBits ret;
//	ret.set_bit_length(bits);
//
//	for (int i = 0; i < bits; ++i)
//		ret.set_bit(i, Number((_x >> i) & 1 ));
//
//	return ret;
}


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


template<class OutNumber, class TempNumber>
void get_percentile(const std::vector<int> &_x, float percentile) {

	OutNumber avgValue;
	OutNumber avgSqrValue;

	n = _x.size();

	AverageLiphe<OutNumber, TempNumber, CompareNative<TempNumber>, TruncConversion> avg(_x.size());
	avg.set_max_sample(*(std::max_element(_x.begin(), _x.end())));
	avg.compute_resample_constant(0.1, 0.05);

	AverageLiphe<OutNumber, TempNumber, CompareNative<TempNumber>, TruncConversion> avgSqrMsd(_x.size());
	avgSqrMsd.set_max_sample(*(std::max_element(_x.begin(), _x.end())));
	avgSqrMsd.compute_resample_constant(0.1, 0.05);
	avgSqrMsd.set_f_m( [](int m)->int{ return ::sqrt(m)*p; } );

	AverageLiphe<OutNumber, TempNumber, CompareNative<TempNumber>, TruncConversion> avgSqrLsd(_x.size());
	avgSqrLsd.set_max_sample(*(std::max_element(_x.begin(), _x.end())));
	avgSqrLsd.compute_resample_constant(0.1, 0.05);
	avgSqrLsd.set_f_m( [](int m)->int{ return ::sqrt(m); } );




	////////////////////////////////////
	// Start computing average and average of squares
	////////////////////////////////////

//	unsigned int i = 0;
//	while (i < _x.size()) {

//		unsigned int simd_begin = i;
//		unsigned int simd_end = i + TempNumber::simd_factor();
//		if (simd_end > _x.size())
//			simd_end = _x.size();
//
//		std::vector<int> simd_vector(TempNumber::simd_factor());
//		auto to_zero = std::copy(_x.begin() + simd_begin,  _x.begin() + simd_end, simd_vector.begin());
//		for (; to_zero != simd_vector.end(); ++to_zero)
//			*to_zero = 0;
//
//		while (i < simd_end) {
//			realAvg += _x[i];
//			realAvgSqr += _x[i] * _x[i];
//			++i;
//		}
//
//		TempNumber xi(simd_vector);

	InputIterator<TempNumber, std::vector<int>::const_iterator> input(_x.begin(), _x.end());
	unsigned int i = 0;
	while (!input.is_end()) {
		TempNumber xi = *input;
		++input;

		start_stat_interval();

		avg.add_simd(xi);
		avgSqrLsd.add_simd(xi);
		avgSqrMsd.add_simd(xi);


//		std::cout << "\n\n===================\n";
//		std::cout << "avg (raw) = " << avg.get_bits().to_bit_stream() << std::endl;
//		std::cout << "avg (bits) = " << avg.getAverage().to_bit_stream() << std::endl;
//		std::cout << "avg  = " << avg.getAverage().to_int() << std::endl;
//		std::cout << std::endl;
//		std::cout << "avgSqr (raw) = " << avgSqr.get_bits().to_bit_stream() << std::endl;
//		std::cout << "avgSqr (bits) = " << avgSqr.getAverage().to_bit_stream() << std::endl;


		unsigned int xii = i * TempNumber::simd_factor();
		while ((xii < _x.size()) && (xii < (i + 1) * TempNumber::simd_factor())) {
			realAvg += _x[xii];
			realAvgSqr += _x[xii] * _x[xii];
			++xii;
		}

		end_stat_interval();

		if ((i % statistic_rate) == statistic_rate - 1) {
			print_stat(i);

			if (measureAccuracy) {
				OutNumber avgTempEnc = fold(avg.getAverage());
				OutNumber avgSqrMsdTempEnc = fold(avgSqrMsd.getAverage());
				OutNumber avgSqrLsdTempEnc = fold(avgSqrLsd.getAverage());

				print_detailed_stats(i, avgTempEnc, avgSqrMsdTempEnc, avgSqrLsdTempEnc);
			}
		}
		++i;
	}

	////////////////////////////////////
	// End computing average and average of squares
	////////////////////////////////////




	////////////////////////////////////
	// start computing threshold
	////////////////////////////////////

	start_stat_interval();

	OutNumber avgEnc = fold(avg.getAverage());
	OutNumber avgLsdEnc = avgEnc * avgEnc;
//	OutNumber avgMsdEnc = simulate_square_msd_polynomial(avgEnc);
	OutNumber avgMsdEnc = SpecialPolynomials<OutNumber>::square_msd_polynomial.compute(avgEnc);

//std::cerr << "square of " << avgEnc.to_int() << " = " << " lsd= " << avgLsdEnc.to_int() << " msd= " << avgMsdEnc.to_int() << std::endl;
	OutNumber avgSqrMsdEnc = fold(avgSqrMsd.getAverage());
	OutNumber avgSqrLsdEnc = fold(avgSqrLsd.getAverage());
//	OutNumber sigma =  simulate_sqrt_polynomial( avgSqrMsdEnc - avgMsdEnc )
//						+ simulate_sqrt_polynomial( avgSqrLsdEnc - avgLsdEnc );
	OutNumber sigma = SpecialPolynomials<OutNumber>::sqrt_polynomial.compute( avgSqrMsdEnc - avgMsdEnc )
						+ SpecialPolynomials<OutNumber>::sqrt_polynomial.compute( avgSqrLsdEnc - avgLsdEnc );

	OutNumber threshold = avgEnc;
	threshold += sigma;
	threshold += sigma;

	TempNumber thresholdBits = convert_to_bits<TempNumber, OutNumber>(threshold);

	////////////////////////////////////
	// end computing threshold
	////////////////////////////////////

	end_stat_interval();






	print_stat(-1);
	if (measureAccuracy) {
		print_detailed_stats(-1, avgEnc, avgSqrMsdEnc, avgSqrLsdEnc);
	}

//	OutNumber sigma = sqrt_poly.compute(avgSqrEnc - avgEnc);

	std::cout << "avg = " << avgEnc.to_int() << std::endl;
	std::cout << "avgSqr = " << (avgSqrMsdEnc.to_int() * p + avgSqrLsdEnc.to_int()) << std::endl;
	std::cout << "sigma = " << sigma.to_int() << std::endl;
	std::cout << "real sigma = " <<  ( ::sqrt(realAvgSqr/n - realAvg*realAvg/(n*n)) ) << std::endl;
	std::cout << "Threshold = " << threshold.to_int() << std::endl;


std::cout << "converting to bits " << threshold.to_int() << " turned into " << thresholdBits.to_int() << std::endl;

//	if (percentile == 15.8) {
//		// take 1 sigma away from avg
//		threshold -= sigma;
//	} else if (percentile == 2.2) {
//		// take 2 sigma away from avg
//		threshold -= sigma * 2;
//	} else if (percentile == 0.1) {
//		// take 3 sigma away from avg
//		threshold -= sigma * 3;
//	}

//	TempNumber thresholdBits = get_bits<OutNumber, TempNumber>(threshold);

//	PercentileIterator<OutNumber, TempNumber> sigmaDist(&rawData, thresholdBits);






	////////////////////////////////////
	// start reporting points
	////////////////////////////////////

	start_stat_interval();

//	i = 0;
//	sigmaDist.begin();
//	while (!sigmaDist.is_end()) {
//		OutNumber x = *sigmaDist;
//
//		end_stat_interval();
//		std::cout << i << ") ";
//		std::cout << "x = " << x.to_int() << "  " << std::endl;
//		print_stat(i);
//		++sigmaDist;
//		++i;
//		start_stat_interval();
//	}
//	std::cout << std::endl;

	i = 0;
	input.set_to_begin();
	while (!input.is_end()) {
		TempNumber xi = *input;
		++input;

		start_stat_interval();
		OutNumber reportEnc = xi > thresholdBits;
		end_stat_interval();

		std::vector<long int> report = reportEnc.to_vector();
		for (auto ri = report.begin(); ri != report.end(); ++ri) {
			if (i < _x.size()) {
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
