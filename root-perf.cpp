#include "common.h"

int main(int argc, char** argv) {
	int which = atoi(argv[1]);
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

