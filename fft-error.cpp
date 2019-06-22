#include "common.h"

int main() {
	const int n = 1 << 16;
	const int maxcoef = 1000;
	vector<C> a(n);
	rep(i,0,n) a[i] = C(rand() % maxcoef, rand() % maxcoef);

	vector<Cd> r1{all(a)};
	fftLd(r1);

	auto r2 = a;
	fftAccurate(r2);

	auto r3 = a;
	fftOld(r3);

	double eps = nextafter(1.0, 2.0) - 1.0; // 2^-52

	auto pr = [&](auto&& vec) {
		ld avg = 0;
		ld maxerr = 0;
		rep(i,0,n) {
			auto diff = Cd(vec[i]) - r1[i];
			ld err = abs(diff) / eps;
			avg += err;
			maxerr = max(maxerr, err);
		}
		avg /= n;
		cout << setprecision(7) << maxerr << ' ' << avg << endl;
	};

	pr(r1);
	pr(r2);
	pr(r3);
}

