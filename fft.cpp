#include <bits/stdc++.h>
using namespace std;

#define rep(i, a, b) for(int i = a; i < (b); ++i)
#define trav(a, x) for(auto& a : x)
#define all(x) x.begin(), x.end()
#define sz(x) (int)(x).size()
typedef long long ll;
typedef pair<int, int> pii;
typedef vector<int> vi;

volatile void *vol;

template<class T>
T fast_mul(T a, T b) {
	return T(real(a) * real(b) - imag(a) * imag(b),
			real(a) * imag(b) + imag(a) * real(b));
}

typedef complex<double> C;
typedef complex<long double> Cd;

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
			auto x = (long double *)&rt[j+k], y = (long double *)&a[i+j+k];
			Cd z(x[0]*y[0] - x[1]*y[1], x[0]*y[1] + x[1]*y[0]);
			a[i + j + k] = a[i + j] - z;
			a[i + j] += z;
		}
}

void fftAccurate(vector<C>& a) {
	int n = sz(a), L = 31 - __builtin_clz(n);
	static vector<complex<long double>> R(2, 1);
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

void fftCurrent(vector<C>& a) {
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

double fftCountMultipliedMass(vector<C>& a) {
	int n = sz(a), L = 31 - __builtin_clz(n);
	static vector<complex<long double>> R(2, 1);
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
	double mass = 0;
	for (int k = 1; k < n; k *= 2) {
		for (int i = 0; i < n; i += 2 * k) rep(j,0,k) {
			auto x = (double *)&rt[j+k], y = (double *)&a[i+j+k];
			mass += norm(a[i+j+k]);
			C z(x[0]*y[0] - x[1]*y[1], x[0]*y[1] + x[1]*y[0]);
			a[i + j + k] = a[i + j] - z;
			a[i + j] += z;
		}
	}
	return mass;
}

void fftError() {
	const int n = 1 << 16;
	const int maxcoef = 1000;
	vector<C> a(n);
	rep(i,0,n) a[i] = C(rand() % maxcoef, rand() % maxcoef);

	vector<Cd> r1{all(a)};
	fftLd(r1);

	auto r2 = a;
	fftAccurate(r2);

	auto r3 = a;
	fftCurrent(r3);

	double eps = nextafter(1.0, 2.0) - 1.0; // 2^-52

	auto pr = [&](auto&& vec) {
		long double avg = 0;
		long double maxerr = 0;
		rep(i,0,n) {
			auto diff = Cd(vec[i]) - r1[i];
			long double err = abs(diff) / eps;
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

typedef vector<double> vd;

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

void convError(int argc, char** argv) {
	// ./a.out conv-error 16 100000 100 10 -1,1
	//                    logn maxcoef its maxrand mode
	const int logn = atoi(argv[2]);
	const int n = 1 << logn;
	const int its = atoi(argv[4]); // 16
	const int maxcoef = atoi(argv[3]); // 100000 for logn = 16
	const int maxrand = atoi(argv[5]); // 2000;
	const double randratio = maxrand / (double)maxcoef;
	int mode = 0;
	int bisa = 0;
	if (argc > 6) {
		char *modestr = argv[6];
		int parti = 0;
		while (char *tok = strtok(modestr, ",")) {
			int parsed = atoi(tok);
			if (parti == 0) mode = parsed;
			if (parti == 1) bisa = parsed;
			parti++;
			modestr = 0;
		}
	}
	srand(its);

	long double avgs[100] = {};
	long double maxerrs[100] = {};

	// In theory we have a bound |n * (log2(n) - 1) * maxcoef^2 < 9.3e14| where results
	// are guaranteed to round correctly. In practice we can go slightly higher.
	cout << (2*n * log2(n) * maxcoef * maxcoef) << ' ' << 9e14 << endl;

	auto r = []() -> double { return rand() / (RAND_MAX + 1.0); };

	rep(it,0,its) {
		int good[] = {0,1,6,8};
		vd a(n), b(n);
		int curmode = mode == -1 ? good[rand() % 4] : mode;
		int curbisa = bisa == -1 ? rand() % 4 : bisa;
		rep(i,0,n) a[i] = maxcoef - (r() * maxrand);
		rep(i,0,n) b[i] = maxcoef - (r() * maxrand);
		if (curmode == 1) rep(i,0,n) {
			a[i] *= (double)sin(M_PIl * 2 / n * i) * M_SQRT2;
			b[i] *= (double)cos(M_PIl * 2 / n * i) * M_SQRT2;
		}
		if (curmode == 2) rep(i,0,n) {
			a[i] = sqrt(r() * maxcoef * maxcoef * 2);
			b[i] = sqrt(r() * maxcoef * maxcoef * 2);
		}
		if (curmode == 3) {
			rep(i,0,n) if (rand() & 128) a[i] *= -1;
			rep(i,0,n) if (rand() & 128) b[i] *= -1;
		}
		if (curmode == 4) {
			rep(i,0,n) a[i] = r() * (maxcoef + 1) - maxcoef / 2.0;
			rep(i,0,n) b[i] = r() * (maxcoef + 1) - maxcoef / 2.0;
		}
		if (curmode == 5) {
			rep(i,0,n) a[i] = b[i] = 0;
			a[n-1] = maxcoef * sqrt(n) * (1 - r() * randratio);
			b[n-1] = maxcoef * sqrt(n) * (1 - r() * randratio);
		}
		if (curmode == 6) {
			rep(i,0,n) a[i] = b[i] = 0;
			a[n-2] = -maxcoef * sqrt(n/2) * (1 - r() * randratio);
			b[n-2] =  maxcoef * sqrt(n/2) * (1 - r() * randratio);
			a[n-1] =  maxcoef * sqrt(n/2) * (1 - r() * randratio);
			b[n-1] = -maxcoef * sqrt(n/2) * (1 - r() * randratio);
		}
		if (curmode == 7 || curmode == 8) {
			bool preferHigh = (curmode == 7);
			for (vd* ap : {&a, &b}) {
				vd& ar = *ap;
				double ws = 0;
				rep(i,0,n) {
					int pc = (int)__builtin_popcount(i);
					if (!preferHigh) pc = logn - pc;
					double w = pow(10000, (double)pc / logn);
					ar[i] = (double)maxcoef * maxcoef * (1 - r() * randratio) * w;
					ws += w;
				}
				ws = n/ws;
				rep(i,0,n) ar[i] *= ws;
				rep(i,0,n) ar[i] = sqrt(ar[i]);
			}
		}
		if (curmode == 9 || curmode == 10) {
			// random vectors with more weight on higher popcounts, with same
			// total L2 mass as the flat distribution
			// this seems to behave almost exactly the same as random...
			bool preferHigh = (curmode == 9);
			for (vd* ap : {&a, &b}) {
				vd& ar = *ap;
				rep(i,0,n) ar[i] += r() * maxrand * maxrand;
				rep(_,0,n*2) {
					int i = 0, ib = 0;
					if (!preferHigh) i = n-1, ib = logn;
					rep(it,0,2) {
						int i2 = rand() % n, i2b = __builtin_popcount(i2);
						if ( preferHigh && i2b > ib) ib = i2b, i = i2;
						if (!preferHigh && i2b < ib) ib = i2b, i = i2;
					}
					ar[i] += (1 - r() * randratio) / 2 * maxcoef * maxcoef;
				}
				rep(i,0,n) ar[i] = sqrt(ar[i]);
			}
		}
		if (curmode == 11) rep(i,0,n) {
			a[i] *= (double)sin(M_PIl * 2 / n * i) * M_SQRT2;
			b[i] *= (double)cos(M_PIl * 2 / n * i) * M_SQRT2;
			a[i] += (r() - 0.5) * maxrand;
			b[i] += (r() - 0.5) * maxrand;
		}

		if (it == 0) {
			double mass = 0;
			rep(i,0,n) mass += a[i] * a[i];
			cout << "L2 norm compared to base " << mass / maxcoef / maxcoef / n << endl;

			vector<C> vec(2*n);
			copy(all(a), begin(vec));
			double m = fftCountMultipliedMass(vec);
			cout << "root-multipled L2 mass " << m << endl;
		}

		if (curbisa == 1) b = a;
		if (curbisa == 2) {
			b = a;
			reverse(all(b));
		}
		if (curbisa == 3) rep(i,0,n) if (rand() & 128) {
			a[i] *= -1;
			b[n-1 - i] *= -1;
		}

		vector<long double> r1 = conv<Cd, long double>(a, b, fftLd);
		vd r2 = conv<C, double>(a, b, fftAccurate);
		vd r3 = conv<C, double>(a, b, fftCurrent);
		vd r4 = conv2<Cd, C, double>(a, b, fftLd, fftAccurate);
		vd r5 = conv2<C, Cd, double>(a, b, fftAccurate, fftLd);
		vd r6 = convNaive<C, double>(a, b, fftAccurate);

		auto pr = [&](auto&& vec, int ind) {
			long double avg = 0;
			long double maxerr = 0;
			rep(i,0,n) {
				auto diff = vec[i] - r1[i];
				long double err = abs(diff);
				avg += err;
				maxerr = max(maxerr, err);
			}
			avg /= n;
			avgs[ind] += avg;
			maxerrs[ind] = max(maxerrs[ind], maxerr);
		};

		pr(r2, 1);
		pr(r3, 2);
		pr(r4, 3);
		pr(r5, 4);
		pr(r6, 5);
	}
	rep(i,1,6) {
		cout << setprecision(5) << maxerrs[i] << ' ' << avgs[i] / its << endl;
	}
}

void convPerf(int argc, char** argv) {
	int which = atoi(argv[2]);
	const int iters = 1;
	const int n = 1 << 23;
	const int maxcoef = 1000;
	for (int it = 0; it < iters; it++) {
		vd a(n), b(n);
		rep(i,0,n) a[i] = rand() % maxcoef;
		rep(i,0,n) b[i] = rand() % maxcoef;

		if (which == 1) {
			vd r1 = conv<Cd, double>(a, b, fftLd);
			vol = r1.data();
		} else if (which == 2) {
			vd r2 = conv<C, double>(a, b, fftAccurate);
			vol = r2.data();
		} else if (which == 3) {
			vd r3 = conv<C, double>(a, b, fftCurrent);
			vol = r3.data();
		}
	}
}

void rootPerf(int argc, char** argv) {
	int which = atoi(argv[2]);
	const int logn = 24;
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

	if (which == 1) {
	rt1.resize(n);
	for (int k = 2; k < n; k *= 2) {
		rep(i, k, 2 * k) rt1[i] = polar(1.0L, M_PIl / k * (i - k));
	}

	} else if (which == 2) {
	rt2.resize(n);
	for (int k = 2; k < n; k *= 2) {
		Cd z[] = {1, polar(1.0L, M_PIl / k)};
		rep(i, k, 2 * k) rt2[i] = rt2[i / 2] * z[i & 1];
	}

	} else if (which == 21) {
	rt2.resize(n);
	for (int k = 2; k < n; k *= 2) {
		rt2[0] = polar(1.0L, M_PIl / k);
		rep(i, k, 2 * k) rt2[i] = rt2[i / 2] * rt2[(i & 1) ^ 1];
	}

	} else if (which == 22) {
	rt2.resize(n);
	for (int k = 2; k < n; k *= 2) {
		auto x = polar(1.0L, M_PIl / k);
		rep(i, k, 2*k) rt2[i] = i&1 ? rt2[i >> 1] * x : rt2[i >> 1];
	}

	} else if (which == 3) {
	rt3.resize(n);
	for (int k = 2; k < n; k *= 2) {
		rep(i, k, 2 * k) rt3[i] = polar(1.0L, M_PIl / k * (i - k));
	}

	} else if (which == 4) {
	rt4.resize(n);
	for (int k = 2; k < n; k *= 2) {
		rep(i, k, 2 * k) rt4[i] = polar(1.0L, 1.L * M_PI / k * (i - k));
	}

	} else if (which == 5) {
	rt5.resize(n);
	for (int k = 2; k < n; k *= 2) {
		rep(i, k, 2 * k) rt5[i] = polar(1.0, M_PI / k * (i - k));
	}

	} else if (which == 6) {
	rt6.resize(n);
	for (int k = 2; k < n; k *= 2) {
		Cd z[] = {1, polar(1.0L, M_PIl / k)};
		rep(i, k, 2 * k) rt6[i] = Cd(rt6[i / 2]) * z[i & 1];
	}

	} else if (which == 7) {
	rt7.resize(n);
	for (int k = 2; k < n; k *= 2) {
		Cd z[] = {1, polar(1.0L, 1.L * M_PI / k)};
		rep(i, k, 2 * k) rt7[i] = Cd(rt7[i / 2]) * z[i & 1];
	}

	} else if (which == 8) {
	rt8.resize(n);
	for (int k = 2; k < n; k *= 2) {
		Cd z[] = {1, polar(1.0, M_PI / k)};
		rep(i, k, 2 * k) rt8[i] = Cd(rt8[i / 2]) * z[i & 1];
	}

	} else if (which == 9) {
	rt9.resize(n);
	for (int k = 2; k < n; k *= 2) {
		C z[] = {1, polar(1.0, M_PI / k)};
		rep(i, k, 2 * k) rt9[i] = C(rt9[i / 2]) * z[i & 1];
	}
	}
}

void rootError() {
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
		long double avg = 0;
		long double maxerr = 0;
		rep(i,0,n) {
			auto diff = Cd(vec[i]) - rt1[i];
			long double err = abs(diff) / eps;
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

void rootAccError() {
	const int lim = 25;
	int n = 1 << lim;
	static vector<complex<long double>> R(2, 1);
	static vector<C> rt(2, 1);
	for (static int k = 2; k < n; k *= 2) {
		R.resize(n); rt.resize(n);
		auto x = polar(1.0L, M_PIl / k);
		rep(i, k, 2*k) rt[i] = R[i] = i&1 ? R[i/2] * x : R[i/2];
	}
	vector<long double> errors(n, 0);
	double eps = (nextafter(1.0, 2.0) - 1.0) / 2; // 2^-53
	int iter = 0;
	for (int k = 1; k < n; k *= 2) {
		iter++;
		for (int i = 0; i < n; i += 2 * k) rep(j,0,k) {
			long double err = abs((Cd)rt[j + k] - R[j + k]);
			errors[i + j] = errors[i + j + k] =
				err + max(errors[i + j], errors[i + j + k]);
		}
		long double maxerr = 0;
		rep(i,0,n) {
			maxerr = max(maxerr, errors[i]);
		}
		cout << iter << ": " << maxerr / eps << endl; // ~0.45 * iter - 0.8
	}
}

double random_range(int range_min, int range_max)
{
	long double r =
		(long double)((unsigned long long)(rand() & 0x7fff) |
		 ((unsigned long long)(rand() & 0x7fff) << 15) |
		 ((unsigned long long)(rand() & 0x7fff) << 30) |
		 ((unsigned long long)(rand() & 0x7fff) << 45));
	r /= (long double)(1ULL << 60);
	return (double)(range_min + r * (long double)(range_max - range_min));
}

void multError() {
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
	auto test = [&rc, eps, iters](C x) {
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

	auto test2 = [&rc, eps, iters]() {
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

void fftPerfEach(int argc, char** argv) {
	// It is possible to bit reversal halfway through the FFT, to make the
	// rep(j,0,k) loops either all long or all short. This function measures
	// the performance of long vs short loops -- it turns out not matter very
	// much. There's a ~10% slowdown for the first two k values.
	// This measurement is imperfect because more roots are used for later
	// iterations, making them slower.
	const int logn = atoi(argv[2]);
	const int n = 1 << logn;
	const int iters = atoi(argv[3]);
	vector<C> a(n);
	trav(x, a) x = C(
				random_range(-10, 10),
				random_range(-10, 10));
	int L = 31 - __builtin_clz(n);
	static vector<complex<long double>> R(2, 1);
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

int main(int argc, char** argv) {
	assert(argc >= 2);
	string mode = argv[1];
	if (mode == "root-error") rootError();
	else if (mode == "root-acc-error") rootAccError();
	else if (mode == "root-perf") rootPerf(argc, argv);
	else if (mode == "mult-error") multError();
	else if (mode == "fft-error") fftError();
	else if (mode == "conv-error") convError(argc, argv);
	else if (mode == "conv-perf") convPerf(argc, argv);
	else if (mode == "fft-perf-each") fftPerfEach(argc, argv);
	else cout << "no such mode" << endl;
}
