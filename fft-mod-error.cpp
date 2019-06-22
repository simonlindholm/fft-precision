#include "common.h"

int main(int argc, char** argv) {
	// ./a.out 16 100000 100 10 -1,1
	// ./a.out 8 2000000 50000 100 0,1
	//         logn maxcoef its maxrand mode
	const int logn = atoi(argv[1]);
	const int n = 1 << logn;
	const int maxcoef = atoi(argv[2]); // e.g. 100000
	const int its = atoi(argv[3]); // e.g. 10
	const int maxrand = atoi(argv[4]); // e.g. 100
	int mode = 0;
	int bisa = 0;
	if (argc > 6) {
		char *modestr = argv[5];
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

	// In theory we have a bound |2n * log2(n) * maxcoef^2 < 8.6e14| where results
	// are guaranteed to round correctly. In practice we can go slightly higher.
	cout << (2*n * log2(n) * maxcoef * maxcoef) << ' ' << 8.6e14 << endl;

	rep(it,0,its) {
		vi a(n), b(n), c(n), d(n);
		int curmode = mode;
		int curbisa = bisa;
		rep(i,0,n) a[i] = maxcoef - rand() % maxrand;
		rep(i,0,n) b[i] = maxcoef - rand() % maxrand;
		rep(i,0,n) c[i] = maxcoef - rand() % maxrand;
		rep(i,0,n) d[i] = maxcoef - rand() % maxrand;
		if (curmode == 1) {
			rep(i,0,n) a[i] = rand() % maxcoef;
			rep(i,0,n) b[i] = rand() % maxcoef;
			rep(i,0,n) c[i] = rand() % maxcoef;
			rep(i,0,n) d[i] = rand() % maxcoef;
		}

		if (curbisa == 1) b = c = d = a;
		if (curbisa == 2) {
			b = c = d = a;
			reverse(all(c));
			reverse(all(d));
		}

		vector<long double> r1 = convMod<Cd, long double>(a, b, c, d, fftLd);
		vd r2 = convMod<C, double>(a, b, c, d, fftAccurate);

		auto pr = [&](auto&& vec, int ind) {
			long double avg = 0;
			long double maxerr = 0;
			rep(i,0,sz(vec)) {
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
	}
	rep(i,1,2) {
		cout << setprecision(5) << maxerrs[i] << ' ' << avgs[i] / its << endl;
	}
}

