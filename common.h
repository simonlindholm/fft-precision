#include <vector>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <climits>
#include <complex>
#include <cmath>
#include <cassert>
#include <cstring>
#include <algorithm>
using namespace std;

#define rep(i, a, b) for(int i = a; i < (b); ++i)
#define trav(a, x) for(auto& a : x)
#define all(x) x.begin(), x.end()
#define sz(x) (int)(x).size()
typedef long long ll;
typedef pair<int, int> pii;
typedef vector<int> vi;

extern volatile void *vol;

double r();

inline double random_range(int range_min, int range_max) {
	return r() * (range_max - range_min) + range_min;
}

template<class T>
T fast_mul(T a, T b) {
	return T(real(a) * real(b) - imag(a) * imag(b),
			real(a) * imag(b) + imag(a) * real(b));
}

typedef vector<double> vd;
typedef complex<double> C;
typedef complex<long double> Cd;

void fftLd(vector<Cd>& a);
void fftAccurate(vector<C>& a);
void fftOld(vector<C>& a);

template<class C, class D, class F>
vector<D> convMod(const vi& alo, const vi& ahi, const vi& blo, const vi& bhi, const F& fft) {
	assert(sz(alo) == sz(ahi));
	assert(sz(blo) == sz(bhi));
	assert(!alo.empty());
	assert(!blo.empty());
	vector<D> res((sz(alo) + sz(blo) - 1) * 4);
	int B = 32 - __builtin_clz(sz(res)), n = 1 << B;
	vector<C> L(n), R(n), outs(n), outl(n);
	rep(i,0,sz(alo)) L[i] = C(alo[i], ahi[i]);
	rep(i,0,sz(blo)) R[i] = C(blo[i], bhi[i]);
	fft(L), fft(R);
	rep(i,0,n) {
		int j = -i & (n - 1);
		outl[j] = (L[i] + conj(L[j])) * R[i] / (D(2) * n);
		outs[j] = (L[i] - conj(L[j])) * R[i] / (D(2) * n) / C(0, 1);
	}
	fft(outl), fft(outs);
	rep(i,0,sz(res)/4) {
		res[i*4 + 0] = real(outl[i]);
		res[i*4 + 1] = imag(outl[i]);
		res[i*4 + 2] = real(outs[i]);
		res[i*4 + 3] = imag(outs[i]);
	}
	return res;
}

template<class C, class D, class F>
vector<D> conv(const vd& a, const vd& b, F fft) {
	if (a.empty() || b.empty()) return {};
	vector<D> res(sz(a) + sz(b) - 1);
	int L = 32 - __builtin_clz(sz(res)), n = 1 << L;
	vector<C> in(n), out(n);
	copy(all(a), begin(in));
	rep(i,0,sz(b)) in[i].imag(b[i]);
	fft(in);
	trav(x, in) x *= x;
	rep(i,0,n) out[i] = in[-i & (n - 1)] - conj(in[i]);
	fft(out);
	rep(i,0,sz(res)) res[i] = D(imag(out[i]) / (4 * n));
	return res;
}

template<class C1, class C2, class D, class F1, class F2>
vector<D> conv2(const vd& a, const vd& b, F1 fft1, F2 fft2) {
	if (a.empty() || b.empty()) return {};
	vector<D> res(sz(a) + sz(b) - 1);
	int L = 32 - __builtin_clz(sz(res)), n = 1 << L;
	vector<C1> in(n);
	vector<C2> out(n);
	copy(all(a), begin(in));
	rep(i,0,sz(b)) in[i].imag(b[i]);
	fft1(in);
	trav(x, in) x *= x;
	rep(i,0,n) out[i] = C2(in[-i & (n - 1)] - conj(in[i]));
	fft2(out);
	rep(i,0,sz(res)) res[i] = D(imag(out[i]) / (4 * n));
	return res;
}

template<class C, class D, class F>
vector<D> convNaive(const vd& a, const vd& b, F fft) {
	if (a.empty() || b.empty()) return {};
	vector<D> res(sz(a) + sz(b) - 1);
	int L = 32 - __builtin_clz(sz(res)), n = 1 << L;
	vector<C> in1(n), in2(n), out(n);
	copy(all(a), begin(in1));
	copy(all(b), begin(in2));
	fft(in1);
	fft(in2);
	rep(i,0,n) out[i] = in1[i] * in2[i];
	reverse(1 + all(out));
	fft(out);
	rep(i,0,sz(res)) res[i] = D(real(out[i])) / n;
	return res;
}
