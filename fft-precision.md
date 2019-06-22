# Notes on the accuracy bounds of KACTL's (2-in-1) FFT and FFT-mod

This write-up discusses the numerical accuracy of:

* content/numerical/FastFourierTransform.h
* content/numerical/FastFourierTransformMod.h

It mainly covers theoretical bounds, but see the bottom for practical ones.

## Basic bound

Our analysis will derive from the following paper, with some modifications:
http://www.daemonology.net/papers/fft.pdf

The paper considers the following method of convolution:

```
complex[] conv(complex[] a, complex[] b):
  A = fft(a)
  B = fft(b)
  complex C[N]
  for i in 0..N-1:
    C[i] = A[i] * B[i] / N
  c = invfft(C)
  return c
```

and shows the following bound on the precision:

`|c - c'|_∞ ≤ |a|₂ * |b|₂ * ((1 + ε)^3n * (1 + ε√5)^(3n+1) * (1 + εβ)^3n - 1)`

where:
- c' is the computed version of c
- `|a|_∞ = max(|a_i|)`
- `|a|₂ = ∑a_i²`
- ε = 2^-53
- εβ is the maximum absolute error of any computed root of unity
  (1/√2 assuming accurate computation)
- N is the FFT input size
- n = log N

Note that `((1 + ε)^3n * (1 + ε√5)^(3n+1) * (1 + εβ)^3n - 1) ≈ 3nε(1 + √5 + 1/√2) + ε√5`,
so for this to give accurate rounding for integer inputs we want approximately:

`3nε(1 + √5 + 1/√2)|a|₂|b|₂ < 1/2`

or

`nKM² < 1/2 / ε / 3 / (1 + √5 + 1/√2) ~ 3.8e14`

if a and b have lengths ≤ K, and entries which are all ≤ M in absolute value.

## Proof of basic bound

The proof of that bound goes roughly as follows:

Additions/subtractions introduce relative errors of at most ε.
Complex multiplications introduce relative errors of at most ε√5 (the paper's
proof is wrong, but it's shown true in a follow-up paper).

For the forward FFTs, we will keep an error bound on |A|₂ (and |B|₂,
symmetrically). This makes sense since it is the quantity that is preserved
across the FFT (up to a factor of N), by Parseval's theorem.

Each FFT step replaces two values (a, b) by (a + ωb, a - ωb) – this replacement
acts linearly on errors, and preserves the L2 norm of the pair.
The multiplication ωb introduces a relative error of at most (1 + ε√5)(1 + εβ),
and the addition an additional (1 + ε) factor in L2 measure across the pair.
In the end, after the FFT steps, we obtain |A' - A|₂ ≤ |A|₂ * (((1 + ε)(1 + ε√5)(1 + εβ))^n - 1).

C is the component-wise product of A and B, so we can use Cauchy-Schwarz
combined with triangle inequalities/linearity to obtain:

`∑|C_i - C'_i| ≤ |A|₂ * |B|₂ * (((1 + ε)(1 + ε√5)(1 + εβ))^2n * (1 + ε√5) - 1) / N`, i.e.
`|C - C'|₁ ≤ |a|₂ * |b|₂ * (((1 + ε)(1 + ε√5)(1 + εβ))^2n * (1 + ε√5) - 1)`

Note the extra complex multiplication which also adds relative error to each
term, and that the division by N adds no error since N is a power of 2.

When we do the inverse FFT, every output term will be the sum of N terms, each
rescaled by some complex root of unity. Thus we can follow each output term's
addition tree backward and look at the relative error in L1 measure, and get:

`|c - c'|_∞ ≤ |a|₂ * |b|₂ * ((1 + ε)^3n * (1 + ε√5)^(3n+1) * (1 + εβ)^3n - 1)`

as wanted.

## Application to the KACTL code

The above analysis does not apply directly to KACTL, which has two versions of (numerical) FFTs:

* a 2-in-1 FFT for convolving doubles
* an FFT-mod for convolutions of integers (modulo other integers, generally)

The two convolution methods look as follows:

```
double[] conv2in1(double[] a, double[] b):
  complex in[N], C[N]
  for i in 0..N-1:
    in[i].real = a[i] if i < len(a) else 0
    in[i].imag = b[i] if i < len(b) else 0
  AB = fft(in)
  for i in 0..N-1:
    j = -i & (N - 1)
    C[i] = (AB[j] * AB[j] - conj(AB[i] * AB[i])) / (4 * N)
    C[i] = (AB[j] - conj(AB[i])) * (AB[j] + conj(AB[i])) / (4 * N)
  out = fft(C)
  double ret[N]
  for i in 0..N-1:
    ret[i] = out[i].imag
  return ret

int[] fftmod(int[] a, int[] b):
  complex ina[N], inb[N], D[N], C[N];
  cut = ceil(sqrt(mod))
  for i in 0..N-1:
    assert 0 <= a[i] < mod
    assert 0 <= b[i] < mod
    ina[i].real = a[i] / cut
    ina[i].imag = a[i] % cut
    inb[i].real = b[i] / cut
    inb[i].imag = b[i] % cut
  L = fft(ina)
  R = fft(inb)
  for i in 0..N-1:
    j = -i & (N - 1)
    C[j] = (L[i] + conj(L[j])) * R[i] / (2 * N)
    D[j] = (L[i] - conj(L[j])) * R[i] / (2 * N) / (0+1i)
  c = fft(C)
  d = fft(D)
  int res[N]
  for i in 0..N-1:
    hi = round(c[i].real)
    lo = round(d[i].imag)
    mid = round(c[i].imag) + round(d[i].real)
    res[i] = (hi * cut * cut + mid * cut + lo) % mod
  return res
```

We look into these individually below.

In addition, the value β = 1/√2 assumes perfectly rounded roots, which is not a priori clear.
KACTL *does* achieve this bound, by computing roots using long doubles:

```cpp
vector<complex<double>> rt(n);
vector<complex<long double>> R(n);
rt[1] = R[1] = 1;
for (int k = 2; k < n; k *= 2) {
  auto x = polar(1.0L, M_PIl / k);
  rep(i, k, 2*k) rt[i] = R[i] = i&1 ? R[i/2] * x : R[i/2];
}
```

Note that the more obvious root computation method:

```cpp
for (int k = 2; k < n; k *= 2) {
  rep(i, k, 2*k) rt[i] = polar(1.0, M_PI / k * (i - k));
}
```

does not fare as well – its worst-case error is 5 times as bad, and average error 2.6 times.

### 2-in-1 FFT

In the 2-in-1 FFT, we first create a vector `in` with a L2 norm of at most `√(|a|₂² + |b|₂²)`.
Then, we FFT it, introducing `(1 + ε)^n * (1 + ε√5)^n * (1 + εβ)^n` relative error as before.

If we compute C as C[i] = AB[j] * AB[j], we would as before get a bound of

`(|a|₂² + |b|₂²) * (((1 + ε)(1 + ε√5)(1 + εβ))^2n * (1 + ε√5) - 1)`

on the L1 norm of C's error. The same holds for conj(AB[i] * AB[i]), and the errors add up
additively together with a (1 + ε) factor to give

`|C - C'|₁ ≤ 1/2 * (|a|₂² + |b|₂²) * (((1 + ε)(1 + ε√5)(1 + εβ))^2n * (1 + ε)(1 + ε√5) - 1)`

after dividing by 4N. This error propagates to the final output, together with
some error-on-error which we ignore for simplicity.

Additionally, we have some error on the inverse FFT, proportional to the L1
norm of C. We can simplify C algebraically to

`C = FFT(a) * FFT(b) / N`,

so by Cauchy-Schwarz we get an additional error of

`|a|₂|b|₂ (((1 + ε)(1 + ε√5)(1 + εβ))^n - 1)`.

for a total of approximately

`(|a|₂² + |a|₂|b|₂ + |b|₂²) * n * ε(1 + √5 + β)`,

which by AM-GM is at most

`(|a|₂² + |b|₂²) * n * ε(1 + √5 + β) * 3/2`.

Thus, for rounding we get an accuracy bound of:

`n(|a|₂² + |b|₂²) < 1/2 / ε / (1 + √5 + β) * 2 / 3 ~ 7.6e14`

which is the same bound `nKM² < 3.8e14` that we had before in the case when the norms are equal.

### FFT-MOD

This works similarly to the 2-in-1 FFT, except for the bound on C. We have L2
norm bounds on `ina` and `inb` as

`|ina|₂ ≤ √2 √mod √len(a)`
`|inb|₂ ≤ √2 √mod √len(b)`

and thus on the FFT errors as

`|L - L'|₂ ≤ √N √2 √mod √len(a) (((1 + ε)(1 + ε√5)(1 + εβ))^n - 1)`
`|R - R'|₂ ≤ √N √2 √mod √len(b) (((1 + ε)(1 + ε√5)(1 + εβ))^n - 1)`.

`(L[i] + conj(L[-i])) / 2` simplifies algebraically to `FFT(Re(ina))`, and
`(L[i] - conj(L[-i])) / 2` to `FFT(Im(ina))`.

Using Cauchy-Schwarz and the triangle inequality we can bound

|C - C'|₁ ≈ |AB - A'B'|₁
          ≤ |A * (B - B')| + |(A - A') * B| + |(A - A') * (B - B')|
          ≈ |A * (B - B')| + |(A - A') * B|
          ≤ |A|₂|B-B'|₂ + |A-A'|₂|B|₂`

with

A = (L[i] + conj(L[-i])) / 2√N,
B = R / √N,

where we have cheated a bit by ignoring the error in the multiplication and the
error-on-error term. This yields approximately

`|C - C'|₁ ≤ (√2 + √2 * √2) * mod √len(a) √len(b) * εn * (1 + √5 + β)`

and similarly for D. C simplifies algebraically to FFT(Re(ina)) * FFT(inb) / N,
and D similarly but with Im. This gives an additional error of at most roughly

`mod √2 √len(a) √len(b) * εn * (1 + √5 + β)`

for the inverse transform. Summing these yields

`|c - c'|_∞ ≤ ~2 * (1 + √2) n mod √(len(a)len(b)) ε(1 + √5 + β)`

and for rounding accurately:

`n mod √(len(a)len(b)) < 1/2 / ε / (4 + √2) / (1 + √5 + 1/√2) ~ 2.36e14`

or `nKM < 2.36e14`.

## Better bounds

It turns out that these bounds can be strengthened slightly, particularly for FFT-MOD.

A first observation we can make is that the first two steps of each FFT involve
multiplication by the unit roots 1 and i, which involves no numerical error.
Thus, an error factor like (1 + ε√5)⁶(1 + εβ)⁶ will vanish from each bound.
If the inputs are integers, an additional factor (1 + ε)⁴ will vanish, because
reasonably-sized integers can be exactly represented by doubles.
Introducing this into the final bound makes it look uglier (it would replace n by
`n - <small constant>`), so we mostly ignore this, but it does mean we can
justify throwing away single factors of (1 + ε) or (1 + ε√5) like we did above.

A second observation is that we reduce the (1 + εβ) factor somewhat. Every one
of the N outputs of the FFT is computed as the sum of N inputs, in a tree-like
fashion, with multiplications by roots of unity along the way. Thus, there are
only N^2 "paths" that multiplications can take from input to output, and none
of the paths is likely to contain only worst-case errors. In fact we compute
the worst-case error along all paths for some reasonable maximum n:

```cpp
int N = 1 << 25;
vector<complex<long double>> rtErr = /* difference between exact and computed roots */;
vector<long double> errors(N, 0);
double eps = (nextafter(1.0, 2.0) - 1.0) / 2; // 2^-53
for (int k = 1; k < N; k *= 2) {
  for (int i = 0; i < N; i += 2 * k) rep(j,0,k) {
    int a = i + j, b = i + j + k;
    errors[a] = errors[b] = abs(rtErr[a]) + max(errors[a], errors[b]);
  }
  auto maxerr = *max_element(begin(errors), end(errors));
  cout << ++iter << ": " << maxerr / eps << endl;
}
```

and it turns out that (at least for n ≤ 25) this grows as ~0.45 * n - 0.84,
rather than the naive ~0.707 * n. This seems to extrapolate to larger n as
well, see fft-root-error.png, although it is hard to see how to prove this
formally. (0.84 comes from the first two roots being particularly nice.)
Thus, we can replace (1 + √5 + 1/√2) ~ 3.94 by (1 + √5 + 0.45) ~ 3.69.
This improves the 2-in-1 bound to `n(|a|₂² + |b|₂²) < 8.1e14`.

A third observation is that only half of all values get multiplied by roots
in each step, so we should be able to reduce the `√5 + β` term by a factor 2.
We do this in different ways for the forward and backward FFTs.

For the forward FFT, we split the array into groups which have been combined
with each other – initially there are N group, and in each step we have a bound
on the L2 norm of each group. In each FFT step, groups are paired up with each
other, and one group in each pair gets its entries multiplied by a root of unity.
Then, the groups are replaced by ± each other. Thus, the relative error that is
introduced introduced is only (1 + ε)(1 + ½(√5 + β)ε).

For the backward FFT, we no longer have bounds on individual coefficients, so
we cannot do the above. However, we can use that the output vector is real to
see that `C[x] = conj(C[-x])` and thus `|C[x]| = |C[-x]|`. x and -x have almost
inverse binary representations of each other: to get the binary representation
of -x, you can take the representation of x, flip every bit, find the last zero
(at dummy position "-1" if x = 0), then flip that zero and all the ones that
come after it. As a consequence, it holds that w(x) + w(-x) ≤ n + 1 where w(x)
is the number of ones in x.

When performing the FFT, the number of times a value will be multiplied by a
root is the same as the number of bits in its binary representation. Hence,
if we are lazy and count only absolute error, the rounding error coming from
indices x and -x together will be at most
(|C[x]| + |C[-x]|) * nε + max(|C[x]|, |C[-x]|) * (n+1)(√5 + β)ε,
which equals |C[x]| * (2n + (n+1)(√5 + β))ε.

We can disregard the +1 by using the first observation, and get an approximate
bound that replaces (1 + √5 + β) ~ 3.69 by (1 + ½(√5 + β)) ~ 2.34.
This improves the 2-in-1 bound to `nM²(len(a) + len(b)) < ~1.3e15` when all
`|a_i|, |b_i| ≤ M`.

The bound of the form `n(|a|₂² + |b|₂²)` only improves slightly, because the
forward FFT might have all its weight on worst-case coefficients. One can
work around that by reversing the arrays before convolving if in a bad case,
but it's not worth the cost in our case. More precisely, the bound becomes
`n(|a|₂² + |b|₂²) < 1/2 / ε / ((1 + √5 + β) + ½(1 + ½(√5 + β))) ~ 9.3e14`.

A fourth observation is that if all input values are non-negative, as in KACTL's
FFT-MOD, a decent amount of the Fourier mass will be on the first coefficient –
the one that represents the sum of the inputs. However, computing the sum of
inputs involves no floating point error at all at long as N * √mod < 2^53. Thus,
we can improve our error bound by handling it specially. We don't care too much
about this for the regular FFT, since it allows negative inputs, but we will
use this to derive a slightly better error bound for the FFT-MOD.

Write p+qi for 1/√mod times the average of the first len(a) values of ina, and
r+si similarly for inb. Then

`N * len(a) * mod * (p + q) ≥ |L|₂² = |L_0|₂² + |L_(≠0)|₂² = len(a)² * mod * (p² + q²) + |L_(≠0)|₂²`,

`|ina_(≠0)|₂² ≤ len(a) * mod * (p + q - len(a)/N * (p² + q²))`

When bounding the error of the FFT of ina, we split all terms into two
categories: those only consisting of sums of inputs ("sum-terms"), and the
rest. Each FFT step will consist of taking two terms (a, b) and creating two
new terms (a + ωb, a - ωb). Sum terms will only be combined with other sum
terms, resulting in one sum term and one normal term. Thus, the number of
sum terms will be halved in each step, and in the end only one sum term will
remain.

The process of replacing a pair of terms with another maintains the L2 norm
of that pair, up to a normalizing factor of √2. Hence, the (normalized)
combined L2 norm of all sum terms will decrease over time, as terms get removed
from the set. The final sum term will have squared L2 norm len(a)² * mod * (p² + q²).
During each step, we will therefore introduce an error which post-FFT gets
a squared L2 norm of at most

`N * len(a) * mod * (p + q - len(a)/N * (p² + q²)) * ε² (1 + ½(√5 + β))²`

disregarding tiny error-on-error effects for simplicity. The non-squared L2
norm of the final error becomes:

`|L' - L|₂ ≤ n √N √len(a) √mod ε √(p + q - len(a)/N * (p² + q²)) * (1 + ½(√5 + β))`.

by the triangle inequality. Cauchy-Schwarz bounds the error on `C = A * B` as
(ignoring cumulative error):

`|C' - C|₁ ≤ |A' - A|₂|B|₂ + |A|₂|B' - B|₂`

In our case, B's L2 error is given by the above bound, A's is twice as much
(due to adding the conjugate), while

`|B|₂ ≤ √len(b) √N √mod √(r + s)`
`|A|₂ = 2 |FFT(Re(ina))| ≤ 2 √len(a) √N √mod √p`

Hence, dividing by 2N,

`|C' - C|₁ ≤ n √len(a) √len(b) mod ε E`

where

`E = √(p + q - len(a)/N * (p² + q²)) (1 + ½(√5 + β)) √(r + s) + √(r + s - len(b)/N * (r² + s²)) (1 + ½(√5 + β)) √p`.

If N grows, E does as well. In the worst case we will have N = 2 * (len(a) + len(b)).
Let x = len(a) / (len(a) + len(b)). Then we can have

`E ≤ √(p + q - x/2 * (p² + q²)) (1 + ½(√5 + β)) √(r + s) + √(r + s - (1 - x)/2 * (r² + s²)) (1 + ½(√5 + β)) √p`.

This error is added up together with numerical error from the inverse FFT,
which is at most

`|C|₁ n ε (1 + ½(√5 + β)) ≤ |A|₂|B|₂/2N n ε (1 + ½(√5 + β)) ≤ n √len(a) √len(b) mod ε F`

where

`F = √p √(r + s) (1 + ½(√5 + β))`

Now it remains to maximize `E + F`, subject to `0 ≤ p,q,r,s,x ≤ 1`.

Doing this numerically, with β = 0.45, results in

`E + F ≤ 10.33`

with the maximum being achieved by `x = 0`, `p = q = r = s = 1`.
This is a 20% win compared to the value achieved before. Our final bound is

`|res|_∞ ≤ n mod √len(a) √len(b) * 1.15e-15`

or for accurate rounding

`n mod √len(a) √len(b) < 4.3e14`.

(For simplicity in the bound KACTL's description strengthens this into

`n mod (len(a) + len(b)) < 4.3e14 * 2`.

Note that √mod * max(len(a), len(b)) < 2^53 under these conditions.)

Instead of doing all this work analyzing the effects of non-negative
coefficients, a saner strategy from a numerical perspective would be to reduce
all input values into [-√mod/2, √mod/2), which saves a factor 4 over the original
estimate. However, that is more code, and unlikely to be worth doing in KACTL.

## Practical bounds

The above has all concerned worst-case theoretical bounds. In practice, the
precision is a fair bit better than those, even on worst-case-like input.
This comes from two factors: floating point errors becoming random and starting
to cancel, and from numerical errors of roots, multiplications and
additions/subtractions being less than their worst cases. According to the
literature, the first aspect converts the factor n in the bound into √n, but we
will stick with n (≈ 15) to compare bounds more easily.

We will *not* be considering input that tries to trigger specific numerical
error: that's a hard task. It requires controling the error of O(N log N)
operations given only O(N) inputs, and the inputs of those operations
depend on each other in complex ways. (However, the number of bits in the input
that can be changed without substantially changing the L2 norm is large enough
that it *might* be possible in theory.)

TODO

## Blinding

We have seen how practical bounds beat theoretical worst-case ones by a decent
margin. Can we somehow force the practical ones to happen? For FFT-MOD it turns
out that we can, through randomization/blinding.

In an FFT-MOD, we want to compute the outputs of a convolution `a * b` modulo
some large integer M. If we take random numbers x and y that are relatively
prime to M (such numbers are plentiful), and multiply every coefficient of a
by x and of b by y, all mod M, then we can multiply by `(xy)^-1` afterwards
to get back our result. However, a and b will look random to the FFT, and we
should be able to use the practical bounds for ina/inb. In fact, if we use
balanced representations for `% M` and `% cut`, coefficients of ina/inb will
point in random directions, have average squared L2 norms of 1/12 times their
usual value, and there will be much cancellation in the end results and thus
in the middle pointwise multiplications, shaving off perhaps a factor √N from
the error.

The existing structure of a and b that the FFT could possibly see and become
worse by is if there are very many of some value. In particular, the ina/inb
coefficients or absolute values thereof may not be concentrated around their
mean if this happens. We could work around this by not (only?) multiplying
everything by a constant, but instead/also by adding a random bias term for
each coefficient. However, correcting for that in the end may be difficult, and
require a second FFT, at which point we would instead get larger wins by
switching the FFT's underlying float type to something more precise.
As an optimization, it should be possible to use the same bias term for larger
chunks of coefficients.

This sort of blinding process is complicated enough that it's not worth doing
in KACTL, but it's fun to think about.
