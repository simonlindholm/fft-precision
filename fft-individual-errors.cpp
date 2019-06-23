#include "common.h"

struct Stats {
	vector<array<ld, 4>> partial;
	ld total = 0;
	ld actual = 0;
	ld weighted = 0;

	Stats& operator+=(const Stats& other) {
		assert(sz(partial) == sz(other.partial));
		rep(i,0,sz(partial)) rep(j,0,4) partial[i][j] += other.partial[i][j];
		total += other.total;
		actual += other.actual;
		weighted += other.weighted;
		return *this;
	}

	Stats& operator/=(ld c) {
		rep(i,0,sz(partial)) rep(j,0,4) partial[i][j] /= c;
		total /= c;
		actual /= c;
		weighted /= c;
		return *this;
	}
};

ostream& operator<<(ostream& os, Stats& s) {
	rep(i,0,sz(s.partial)) {
		auto ar = s.partial[i];
		os << "level " << i << ":\t" << ar[0] << '\t' << ar[1] << '\t' << ar[2] << '\t' << ar[3] << endl;
	}
	os << "bound on L2 error " << s.total << endl;
	os << "possible average " << s.total / sqrt(sz(s.partial)) << endl;
	os << "actual L2 error " << s.actual << endl;
	os << "weighted L1 error " << s.weighted << endl;
	return os;
}

Stats myfft(vector<C>& a) {
	int n = sz(a), L = 31 - __builtin_clz(n);
	static vector<complex<ld>> R(2, 1);
	static vector<C> rt(2, 1);   // ^ 10% faster if double
	for (static int k = 2; k < n; k *= 2) {
		R.resize(n); rt.resize(n);
		auto x = polar(1.0L, M_PIl / k); // (lower-case L)
		rep(i, k, 2*k) rt[i] = R[i] = i&1 ? R[i/2] * x : R[i/2];
	}
	vi rev(n);
	rep(i,0,n) rev[i] = (rev[i / 2] | (i & 1) << L) / 2;
	rep(i,0,n) if (i < rev[i]) swap(a[i], a[rev[i]]);
	int kit = 0;
	ld totalErr = 0;
	Stats s;
	s.partial.resize(L);
	for (int k = 1; k < n; k *= 2, kit++) {
		double pw2 = n / k / 2;
		ld multErrs = 0, addErrs = 0, subErrs = 0, sumErrs;
		for (int i = 0; i < n; i += 2 * k) rep(j,0,k) {
			auto x = (double *)&rt[j+k], y = (double *)&a[i+j+k];
			auto X = (ld *)&R[j+k];
			C z(x[0]*y[0] - x[1]*y[1], x[0]*y[1] + x[1]*y[0]);
			Cd Z(X[0]*y[0] - X[1]*y[1], X[0]*y[1] + X[1]*y[0]);
			Cd Z2((ld)x[0]*y[0] - (ld)x[1]*y[1], (ld)x[0]*y[1] + (ld)x[1]*y[0]);
			ld multErr = norm(Z - (Cd)z) * 2 * pw2;
			// cout << "complex mult error: " << log2(norm(Z - (Cd)z) * 2 * pw2) << ' ' << kit << endl;
			// if (multErr != 0.0) cout << log2(multErr) << ' ' << log2(norm(Z2 - (Cd)z) * 2 * pw2) << ' ' << kit << endl;
			auto sub = a[i + j] - z;
			auto add = a[i + j] + z;
			auto SUB = (Cd)a[i + j] - (Cd)z;
			auto ADD = (Cd)a[i + j] + (Cd)z;
			ld addErr = norm(ADD - (Cd)add) * pw2;
			ld subErr = norm(SUB - (Cd)sub) * pw2;
			multErrs += multErr;
			addErrs += addErr;
			subErrs += subErr;
			// cout << "addition error: " << addErr << endl;
			// cout << "subtraction error: " << norm(SUB - (Cd)sub) * pw2 << endl;
			a[i + j + k] = sub;
			a[i + j] = add;
		}
		s.partial[kit][0] = multErrs = sqrt(multErrs);
		s.partial[kit][1] = sumErrs = sqrt(addErrs + subErrs);
		s.partial[kit][2] = addErrs = sqrt(addErrs);
		s.partial[kit][3] = subErrs = sqrt(subErrs);
		totalErr += multErrs;
		totalErr += sumErrs;
	}
	s.total = totalErr;
	return s;
}

int main(int argc, char** argv) {
	// ./a.out 16 10000000000 2000 10 0
	//         logn maxcoef maxrand its mode
	const int logn = atoi(argv[1]);
	const int n = 1 << logn;
	static_assert(sizeof(long) == sizeof(ll));
	const double maxcoef = (double)atol(argv[2]);
	const double maxrand = (double)atol(argv[3]);
	const int its = atoi(argv[4]);
	const double randratio = maxrand / maxcoef;
	int mode = atoi(argv[5]);

	Stats sum;
	sum.partial.resize(logn+1);
	rep(it,0,its) {
		vd a(n), b(n);
		rep(i,0,n) a[i] = maxcoef - (r() * maxrand);
		rep(i,0,n) b[i] = maxcoef - (r() * maxrand);

		if (mode == 1 || mode == 2) {
			bool preferHigh = (mode == 1);
			vd& ar = a;
			double ws = 0;
			rep(i,0,n) {
				int pc = (int)__builtin_popcount(i);
				if (!preferHigh) pc = logn - pc;
				double w = pow(5, pc);
				ar[i] = maxcoef * maxcoef * (1 - r() * randratio) * w;
				ws += w;
			}
			ws = n/ws;
			rep(i,0,n) ar[i] *= ws;
			rep(i,0,n) ar[i] = sqrt(ar[i]);
		}
		if (mode == 3) {
			rep(i,0,n) a[i] *= r();
		}
		if (mode == 4) {
			rep(i,0,n) a[i] *= (r() - 0.5) * 2;
		}

		if (it == 0) {
			double mass = 0;
			rep(i,0,n) mass += a[i] * a[i];
			cout << "L2 norm compared to base " << mass / maxcoef / maxcoef / n << endl;
		}

		vector<C> ac(2*n);
		copy(all(a), begin(ac));
		vector<Cd> AC(2*n);
		copy(all(a), begin(AC));
		Stats s = myfft(ac);
		fftLd(AC);
		ld err = 0;
		ld weighted = 0;
		rep(i,0,2*n) {
			Cd E = (Cd)ac[i] - AC[i];
			err += norm(E);
			weighted += abs(E * AC[i]);
		}
		s.actual = sqrt(err);
		s.weighted = weighted;
		sum += s;
	}
	sum /= its;
	cout << sum;
}
