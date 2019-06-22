#include "common.h"

int main() {
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

