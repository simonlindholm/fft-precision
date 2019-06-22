#include "common.h"

int main(int argc, char** argv) {
	int which = atoi(argv[1]);
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

