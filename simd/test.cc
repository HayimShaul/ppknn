#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <iostream>
#include <vector>
#include <chrono>
#include <random>

int sample_size = 10000;
int avg = 1000;
int sigma = 30;

void make_sample(std::vector<int> &sample) {
	int seed = 1;

	seed = std::chrono::system_clock::now().time_since_epoch().count();

	std::default_random_engine generator (seed);
	std::normal_distribution<float> distribution(avg, sigma);

	avg = random() % 2000;
	sigma = random() % 100;

	sample.resize( sample_size );
	for (int i = 0; i < sample.size(); ++i) {
		sample[i] = distribution(generator);
	}
}

int bit_length(int i) {
	int ret = 0;
	while ((1 << ret) <= i)
		++ret;
	return ret;
}

int get_bit(int x, int i) {
	return (x >> i) & 1;
}

int set_bit(int b, int i) {
	assert ((b == 0) || (b == 1));
	return b << i;
}

int approx_sqr(int x) {
	int l = bit_length(x);
	int ret = 0;

	for (int i = 0; i < 2*l; ++i) {
		ret += set_bit( get_bit(x, i/2) * get_bit(x, (i+1)/2), i);
	}

	return ret;
}


int approx_sqrt(int x) {
	int l = bit_length(x);
	int ret = 0;

	for (int i = 0; i < (l+1)/2; ++i) {
		int a_2i = get_bit(x, 2*i);
		int a_2i1 = get_bit(x, 2*i + 1);

		ret += set_bit( a_2i + a_2i1 - a_2i*a_2i1, i);
	}

	return ret;
}



int get_avg(const std::vector<int> x) {
	int ret = 0;
	for (int i = 0; i < x.size(); ++i)
		ret += x[i];
	return ret / x.size();
}

int get_avg2(const std::vector<int> x) {
	int ret = 0;
	for (int i = 0; i < x.size(); ++i)
		ret += x[i] * x[i];
	return ret / x.size();
}




void test_dan(const std::vector<int> &sample) {
	int avg = get_avg(sample);
	int avg2 = get_avg2(sample);

//	int sigma2 = avg2 - approx_sqr(avg);
	int sigma2 = avg2 - avg*avg;
//	int T = approx_sqr(avg) - 4*sigma2;
	int T = avg*avg - 4*sigma2;
std::cout << "T = " << T << std::endl;

	int count = 0;
	for (int i = 0; i < sample.size(); ++i) {
//std::cout << "sample: " << sample[i] << "   sample^2= " << sample[i]*sample[i] << std::endl;
		if (sample[i]*sample[i] < T)
			++count;
	}

	std::cout << "dan counted " << count << " millionairs" << std::endl;
}

void test_mine(const std::vector<int> &sample) {
	int avg = get_avg(sample);
	int avg2 = get_avg2(sample);

	int sigma = approx_sqrt(avg2 - approx_sqr(avg));
	int T = avg - 2*sigma;

	int count = 0;
	for (int i = 0; i < sample.size(); ++i) {
		if (sample[i] < T)
			++count;
	}

	std::cout << "i   counted " << count << " millionairs" << std::endl;
}

void real_results(const std::vector<int> &sample) {
	int avg = get_avg(sample);
	int avg2 = get_avg2(sample);

	int sigma = sqrt(avg2 - avg*avg);
	int T = avg - 2*sigma;

	int count = 0;
	for (int i = 0; i < sample.size(); ++i) {
		if (sample[i] < T)
			++count;
	}

	std::cout << "real counted " << count << " millionairs" << std::endl;
}


int main(int, char**) {
	std::vector<int> sample;

	make_sample(sample);
	
	real_results(sample);
	test_dan(sample);
	test_mine(sample);
}
