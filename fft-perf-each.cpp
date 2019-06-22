#include "common.h"

int main(int argc, char** argv) {
	// It is possible to bit reversal halfway through the FFT, to make the
	// rep(j,0,k) loops either all long or all short. This function measures
	// the performance of long vs short loops -- it turns out not matter very
	// much. There's a ~10% slowdown for the first two k values.
	// This measurement is imperfect because more roots are used for later
	// iterations, making them slower.
	const int logn = atoi(argv[1]);
	const int n = 1 << logn;
	const int iters = atoi(argv[2]);
	vector<C> a(n);
	trav(x, a) x = C(
				random_range(-10, 10),
				random_range(-10, 10));
	int L = 31 - __builtin_clz(n);
	static vector<complex<ld>> R(2, 1);
	static vector<C> rt(2, 1);
	for (static int k = 2; k < n; k *= 2) {
		R.resize(n); rt.resize(n);
		auto x = polar(1.0L, M_PIl / k);
		rep(i, k, 2*k) rt[i] = R[i] = i&1 ? R[i/2] * x : R[i/2];
	}
	vi rev(n);
	rep(i,0,n) rev[i] = (rev[i / 2] | (i & 1) << L) / 2;
	vector<clock_t> times(logn);
	rep(its,0,iters) {
		rep(i,0,n) if (i < rev[i]) swap(a[i], a[rev[i]]);
		int ctr = 0;
		for (int k = 1; k < n; k *= 2) {
			clock_t t = clock();
			for (int i = 0; i < n; i += 2 * k) rep(j,0,k) {
				auto x = (double *)&rt[j+k], y = (double *)&a[i+j+k];
				C z(x[0]*y[0] - x[1]*y[1], x[0]*y[1] + x[1]*y[0]);
				a[i + j + k] = a[i + j] - z;
				a[i + j] += z;
			}
			times[ctr++] += clock() - t;
		}
	}
	vol = a.data();
	clock_t sum = 0, sum2 = 0, sum3 = 0;
	rep(i,0,logn) {
		cout << times[i] << endl;
		sum += times[i];
		if (i*2 < logn-1) sum2 += times[i];
		else sum3 += times[i];
		if (i*2 < logn) sum2 += times[i];
		else sum3 += times[i];
	}
	cout << endl;
	cout << sum << ' ' << (sum2 + sum3) / 2 << endl;
	cout << sum2 << ' ' << sum3 << endl;
}
