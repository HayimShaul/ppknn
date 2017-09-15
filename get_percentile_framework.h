#ifndef ___GET_PERCENTILE_FRAMEWORK___
#define ___GET_PERCENTILE_FRAMEWORK___

#include <sys/resource.h>
#include <sys/types.h>
#include <vector>

#include "polynomial.h"
#include "get_percentile.h"







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
void print_detailed_stats(int step, NUMBER &avgTemp, NUMBER &avgSqrTemp) {
	int avg = avgTemp.to_int();
	int avgSqr = avgSqrTemp.to_int();

	print_stat_prefix(step);
	std::cout << "Avg " << avg << " / " << (realAvg / n) << " = " << ((double)avg * n / realAvg) << std::endl;
	std::cout << "Avg Sqr " << avgSqr << " / " << (realAvgSqr / n) << " = " << ((double)avgSqr * n / realAvgSqr) << std::endl;

}

inline int sqr(int x) { return x*x; }

int statistic_rate = 10;







template<class Number>
class PercentileIterator {
private:
	const std::vector<int> *_raw_data;
	Number _threshold;
	int _location;

	Number _last_value;
	int _last_value_location;

	void copy(const PercentileIterator &a) { _raw_data = a._raw_data; _threshold = a._threshold; _location = a._location; _last_value_location = a._last_value_location; _last_value = a._last_value; }
public:
	PercentileIterator(const std::vector<int> *arr, const Number &thr) : _raw_data(arr), _threshold(thr), _location(0), _last_value(-1) {}
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

	Number operator*() {
		if (_location != _last_value_location) {
			Number xi = Number::static_from_int((*_raw_data)[_location]);
			_last_value = (xi < _threshold);
			_last_value_location = _location;
		}
		std::cout << "* operator = " << _last_value.to_int() << std::endl;
		return _last_value;
	}
};


template<class Number>
Polynomial<Number> build_sqrt_polynomial(int p) {
	Polynomial<Number> poly(0);

	for (int i = 0; i < p; ++i) {
		std::cerr << "Building sqrt polynomial " << i << " / " << p << "\r";
		std::cerr << std::endl;
		Polynomial<Number> i_x = Polynomial<Number>(i, "-x").set_mod(p);
		i_x ^= phi(p);
		
		poly += (Polynomial<Number>(1) - i_x) * sqrt(i);
//		poly += (Polynomial<Number>(1) - (Polynomial<Number>(i, "-x").set_mod(p))^phi(p)) * sqrt(i);
	}

	return poly;
}


template<class Number>
Number get_bit(const Number &x, int bit) {
	Polynomial<Number> p(0);

	for (int i = 0; i < x.p(); ++i) {
		if ((i & (1 << bit)) == 1) {
			Polynomial<Number> i_x = Polynomial<Number>(i, "-x");
			i_x ^= phi(x.p());
			p += Polynomial<Number>(1) - i_x;

//			p += Polynomial<Number>(1) - Polynomial<Number>(i, "-x")^phi(x.p());
		}
	}

	return p.compute(x);
}


template<class Number, class BitNumber>
BitNumber get_bits(const Number &x) {
	BitNumber ret;

	int bit_size = 1;
	while ((1 << bit_size) < x.p())
		++bit_size;

	ret.set_bit_length(bit_size);

	for (int i = 0; i < bit_size; ++i) {
		ret.set_bit(i, get_bit(x, i));
	}

	return ret;
}


template<class NumberBits>
NumberBits fhe_sqrt_approximation(const NumberBits &x) {
	NumberBits ret;
	ret.set_bit_length((x.bitLength() + 1) / 2);

	for (int i = 0; i < x.bitLength() / 2; ++i) {
		ret[i] = x[2*i] + x[2*i+1] - ( x[2*i] * x[2*i+1] );
	}

	if ((x.bitLength() % 2) == 1) {
		ret[x.bitLength() / 2] = x[x.bitLength() - 1];
	}

	return ret;
}

template<class OutNumber, class TempNumber>
void get_percentile(const std::vector<int> &_x, float percentile) {
	OutNumber avgValue;
	OutNumber avgSqrValue;

	n = _x.size();

	AverageLipheBits<TempNumber, OutNumber, TempNumber, CompareNative<TempNumber>, TruncConversion> avg(_x.size());
	avg.set_bit_number(5); // TODO: read it from the input
	avg.compute_resample_constant(0.1, 0.05, 11);
//	avg.set_resample_constant( TODO ... )

	AverageLipheBits<TempNumber, OutNumber, TempNumber, CompareNative<TempNumber>, TruncConversion> avgSqr(_x.size());
	avgSqr.set_bit_number(9); // TODO: read it from the input
//	avgSqr.set_resample_constant( TODO ... )

	for (int i = 0; i < _x.size(); ++i) {

		TempNumber xi(_x[i]);

		start_stat_interval();

		avg.add(xi);
		OutNumber cost = xi.bits_to_number();
		avgSqr.add_with_cost(xi, cost);

//		std::cout << "\n\n===================\n";
//		std::cout << "avg (raw) = " << avg.get_bits().to_bit_stream() << std::endl;
//		std::cout << "avg (bits) = " << avg.getAverage().to_bit_stream() << std::endl;
//		std::cout << "avg  = " << avg.getAverage().to_int() << std::endl;
//		std::cout << std::endl;
//		std::cout << "avgSqr (raw) = " << avgSqr.get_bits().to_bit_stream() << std::endl;
//		std::cout << "avgSqr (bits) = " << avgSqr.getAverage().to_bit_stream() << std::endl;

		end_stat_interval();

		if ((i % statistic_rate) == statistic_rate - 1) {
			print_stat(i);

			if (measureAccuracy) {
				TempNumber avgTempEnc = avg.getAverage();
				TempNumber avgSqrTempEnc = avgSqr.getAverage();

				print_detailed_stats(i, avgTempEnc, avgSqrTempEnc);
			}
		}
	}

		std::cout << "\n\n===================\n";
		std::cout << "avg (raw) = " << avg.get_bits().to_bit_stream() << std::endl;
		std::cout << "avg (bits) = " << avg.getAverage().to_bit_stream() << std::endl;
		std::cout << "avg  = " << avg.getAverage().to_int() << std::endl;
		std::cout << std::endl;
		std::cout << "avgSqr (raw) = " << avgSqr.get_bits().to_bit_stream() << std::endl;
		std::cout << "avgSqr (bits) = " << avgSqr.getAverage().to_bit_stream() << std::endl;

	std::cout << "Exiting ... \n";
	return;

	TempNumber avgEnc = avg.getAverage();
	TempNumber avgSqrEnc = avgSqr.getAverage();

	print_stat(-1);
	if (measureAccuracy) {
		print_detailed_stats(-1, avgEnc, avgSqrEnc);
	}

//	OutNumber sigma = sqrt_poly.compute(avgSqrEnc - avgEnc);

	std::cout << "avg = " << avgEnc.to_int() << std::endl;
	std::cout << "avgSqr = " << avgSqrEnc.to_int() << std::endl;

	TempNumber xxx = avgSqrEnc - avgEnc;
	TempNumber sigma = fhe_sqrt_approximation(xxx);

	TempNumber threshold = avgEnc;

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

	threshold -= (sigma << 1);

//	TempNumber thresholdBits = get_bits<OutNumber, TempNumber>(threshold);

	PercentileIterator<TempNumber> sigmaDist(&rawData, threshold);


	int i = 0;
	sigmaDist.begin();
	while (!sigmaDist.is_end()) {
		OutNumber x = (*sigmaDist).bits_to_number();
		std::cout << i << ") ";
		std::cout << "x = " << x.to_int() << "  " << std::endl;
		++sigmaDist;
		++i;
	}
	std::cout << std::endl;
//	report(sigmaDist, compactRepresentation)
}


#endif
