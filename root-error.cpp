#include "common.h"

int main() {
	const int logn = 16;
	int n = 1 << logn;

	vector<Cd> rt1(2, 1);
	vector<Cd> rt2(2, 1);
	vector<C> rt3(2, 1);
	vector<C> rt4(2, 1);
	vector<C> rt5(2, 1);
	vector<C> rt6(2, 1);
	vector<C> rt7(2, 1);
	vector<C> rt8(2, 1);
	vector<C> rt9(2, 1);
	rt1.resize(n);
	rt2.resize(n);
	rt3.resize(n);
	rt4.resize(n);
	rt5.resize(n);
	rt6.resize(n);
	rt7.resize(n);
	rt8.resize(n);
	rt9.resize(n);

	for (int k = 2; k < n; k *= 2) {
		rep(i, k, 2 * k) rt1[i] = polar(1.0L, M_PIl / k * (i - k));
	}

	for (int k = 2; k < n; k *= 2) {
		Cd z[] = {1, polar(1.0L, M_PIl / k)};
		rep(i, k, 2 * k) rt2[i] = Cd(rt2[i / 2]) * z[i & 1];
	}

	for (int k = 2; k < n; k *= 2) {
		rep(i, k, 2 * k) rt3[i] = polar(1.0L, M_PIl / k * (i - k));
	}

	for (int k = 2; k < n; k *= 2) {
		rep(i, k, 2 * k) rt4[i] = polar(1.0L, 1.L * M_PI / k * (i - k));
	}

	for (int k = 2; k < n; k *= 2) {
		rep(i, k, 2 * k) rt5[i] = polar(1.0, M_PI / k * (i - k));
	}

	for (int k = 2; k < n; k *= 2) {
		Cd z[] = {1, polar(1.0L, M_PIl / k)};
		rep(i, k, 2 * k) rt6[i] = Cd(rt6[i / 2]) * z[i & 1];
	}

	for (int k = 2; k < n; k *= 2) {
		Cd z[] = {1, polar(1.0L, 1.L * M_PI / k)};
		rep(i, k, 2 * k) rt7[i] = Cd(rt7[i / 2]) * z[i & 1];
	}

	for (int k = 2; k < n; k *= 2) {
		Cd z[] = {1, polar(1.0, M_PI / k)};
		rep(i, k, 2 * k) rt8[i] = Cd(rt8[i / 2]) * z[i & 1];
	}

	for (int k = 2; k < n; k *= 2) {
		C z[] = {1, polar(1.0, M_PI / k)};
		rep(i, k, 2 * k) rt9[i] = C(rt9[i / 2]) * z[i & 1];
	}

	double eps = nextafter(1.0, 2.0) - 1.0; // 2^-52

	auto pr = [&](auto&& vec) {
		ld avg = 0;
		ld maxerr = 0;
		rep(i,0,n) {
			auto diff = Cd(vec[i]) - rt1[i];
			ld err = abs(diff) / eps;
			avg += err;
			maxerr = max(maxerr, err);
		}
		avg /= n;
		cout << setprecision(6) << fixed << maxerr << ' ' << avg << endl;
	};

	pr(rt1);
	pr(rt2);
	pr(rt3);
	pr(rt4);
	pr(rt5);
	pr(rt6);
	pr(rt7);
	pr(rt8);
	pr(rt9);
}

