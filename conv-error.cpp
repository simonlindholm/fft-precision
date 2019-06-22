#include "common.h"

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

int main(int argc, char** argv) {
	// ./a.out 16 100000 100 10 -1,1
	// ./a.out 8 2000000 50000 100 0,1
	//         logn maxcoef its maxrand mode
	const int logn = atoi(argv[1]);
	const int n = 1 << logn;
	const int maxcoef = atoi(argv[2]); // e.g. 100000 for logn = 16
	const int its = atoi(argv[3]); // e.g. 10
	const int maxrand = atoi(argv[4]); // e.g. 100
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

	// In theory we have a bound |2n * log2(n) * maxcoef^2 < 9.3e14| where results
	// are guaranteed to round correctly. In practice we can go slightly higher.
	cout << (2*n * log2(n) * maxcoef * maxcoef) << ' ' << 9e14 << endl;

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
					double w = pow(100, (double)pc / logn);
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
			a[i] *= sin(M_PI * 2 / n * i) * M_SQRT2;
			b[i] *= cos(M_PI * 2 / n * i) * M_SQRT2;
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
			cout << "root-multiplied L2 mass " << m << endl;
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
