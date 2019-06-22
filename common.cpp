#include "common.h"

volatile void *vol;

double r() {
	static_assert(RAND_MAX == INT_MAX, "");
	long long max = 1LL << 31;
	long long r = rand();
	r *= max;
	r ^= rand();
	return (double)r / (double)(max * max);
}

void fftLd(vector<Cd>& a) {
	int n = sz(a), L = 31 - __builtin_clz(n);
	vi rev(n);
	rep(i,0,n) rev[i] = (rev[i / 2] | (i & 1) << L) / 2;
	static vector<Cd> rt(2, 1);
	for (static int k = 2; k < n; k *= 2) {
		rt.resize(n);
		Cd z[] = {1, polar(1.0L, M_PIl / k)};
		rep(i, k, 2 * k) rt[i] = rt[i / 2] * z[i & 1];
		// rep(i, k, 2 * k) rt[i] = polar(1.0L, M_PIl / k * (i - k));
	}
	rep(i,0,n) if (i < rev[i]) swap(a[i], a[rev[i]]);
	for (int k = 1; k < n; k *= 2)
		for (int i = 0; i < n; i += 2 * k) rep(j,0,k) {
			auto x = (ld *)&rt[j+k], y = (ld *)&a[i+j+k];
			Cd z(x[0]*y[0] - x[1]*y[1], x[0]*y[1] + x[1]*y[0]);
			a[i + j + k] = a[i + j] - z;
			a[i + j] += z;
		}
}

void fftAccurate(vector<C>& a) {
	int n = sz(a), L = 31 - __builtin_clz(n);
	static vector<complex<ld>> R(2, 1);
	static vector<C> rt(2, 1);   // ^ 10% faster if double
	for (static int k = 2; k < n; k *= 2) {
		R.resize(n); rt.resize(n);
		auto x = polar(1.0L, M_PIl / k); // (lower-case L)
		// rep(i, k, 2 * k) r[i] = R[i] = R[i / 2] * R[!(i & 1)];
		rep(i, k, 2*k) rt[i] = R[i] = i&1 ? R[i/2] * x : R[i/2];
	}
	vi rev(n);
	rep(i,0,n) rev[i] = (rev[i / 2] | (i & 1) << L) / 2;
	rep(i,0,n) if (i < rev[i]) swap(a[i], a[rev[i]]);
	for (int k = 1; k < n; k *= 2) {
		for (int i = 0; i < n; i += 2 * k) rep(j,0,k) {
			// C z = rt[j+k] * a[i+j+k]; // (25% faster if hand-rolled)  /// include-line
			auto x = (double *)&rt[j+k], y = (double *)&a[i+j+k];        /// exclude-line
			C z(x[0]*y[0] - x[1]*y[1], x[0]*y[1] + x[1]*y[0]);           /// exclude-line
			a[i + j + k] = a[i + j] - z;
			a[i + j] += z;
		}
	}
}

void fftOld(vector<C>& a) {
	int n = sz(a), L = 31 - __builtin_clz(n);
	vi rev(n); static vector<C> rt(2, 1);
	rep(i,0,n) rev[i] = (rev[i / 2] | (i & 1) << L) / 2;
	for (static int k = 2; k < n; k *= 2) {
		rt.resize(n);
		Cd z[] = {1, polar(1.0L, M_PIl / k)};
		rep(i, k, 2 * k) rt[i] = Cd(rt[i / 2]) * z[i & 1];
	}
	rep(i,0,n) if (i < rev[i]) swap(a[i], a[rev[i]]);
	for (int k = 1; k < n; k *= 2)
		for (int i = 0; i < n; i += 2 * k) rep(j,0,k) {
			// C z = rt[j+k] * a[i+j+k]; // (25% faster if hand-rolled)  /// include-line
			auto x = (double *)&rt[j+k], y = (double *)&a[i+j+k];        /// exclude-line
			C z(x[0]*y[0] - x[1]*y[1], x[0]*y[1] + x[1]*y[0]);           /// exclude-line
			a[i + j + k] = a[i + j] - z;
			a[i + j] += z;
		}
}
