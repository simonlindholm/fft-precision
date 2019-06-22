#include "common.h"
#include <random>

int main() {
	const int iters = 5*3*3*3000000;
	double eps = 1 - nextafter(1, 0); // 2^-53
	default_random_engine rng;
	auto rc = [&rng]() -> C {
		// return C(
				// uniform_real_distribution<double>(-10, 10)(rng),
				// uniform_real_distribution<double>(-10, 10)(rng));
		return C(
				random_range(-10, 10),
				random_range(-10, 10));
	};
	auto test = [&rc, eps](C x) {
		Cd X = x;
		double res = 0;
		rep(it,0,iters) {
			C y = rc();
			Cd Y = y;
			C prod = x * y;
			Cd PROD = X * Y;
			C diff = C((Cd)prod - PROD);
			res = max(res, abs(diff / prod));
		}
		return res / eps;
	};

	auto test2 = [&rc, eps]() {
		double res = 0;
		rep(it,0,iters) {
			C x = rc();
			Cd X = x;
			C y = rc();
			Cd Y = y;
			C prod = x * y;
			Cd PROD = X * Y;
			C diff = C((Cd)prod - PROD);
			res = max(res, abs(diff / prod));
		}
		return res / eps;
	};

	cout << "sqrt 5: " << sqrt(5) << endl;
	cout << "random: " << test(rc()) << endl;
	cout << "random pairs: " << test2() << endl;

	/*
	double res = 0;
	rep(it,0,100) {
		res = max(res, test(rc()));
	}

	cout << "many random: " << res << endl;
	*/

	rep(it,0,20) {
		C x = (C)polar(1.0L, M_PIl / (1 << it));
		cout << "root " << it << ": " << test(x) << endl;
	}
}

