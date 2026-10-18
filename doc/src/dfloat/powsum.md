---
title: Sieved power sums with dfloat balls
functions: gr_powsum_sieved dfloat_powsum_sieved _dfloat_powsum_sieved _dfloat_powsum_sieved_ball
---

These functions compute the partial sums
$$S = \sum_{n=1}^{N} n^{-s}$$
that make up the main sum of the Riemann–Siegel formula for $\zeta(s)$
at large height $t = \operatorname{Im}(s)$. `gr_powsum_sieved` is a
generic implementation over three rings given at run time, and
`dfloat_powsum_sieved` instantiates it with `dfloat` balls and plain
expansions. `acb_dirichlet_zeta_rs_r` calls `dfloat_powsum_sieved` for
$N \ge 16$ and falls back to `acb_dirichlet_powsum_sieved` when it
returns 0; this makes the Riemann–Siegel evaluations used to isolate
zeros 17 to 18 times faster at the heights of the $10^{15}$-th to
$10^{17}$-th zeros.

## The algorithm {#algorithm}

The algorithm is that of `acb_dirichlet_powsum_sieved`. The function
$f(n) = n^{-s}$ is completely multiplicative, so the term of a
composite $n$ is the product $f(d) f(n/d)$ with $d$ the smallest prime
factor of $n$, taken from a table of $f(k)$ for odd $k \le N/3$; only
the prime terms need logarithms, sines, cosines and exponentials. The
even $n$ are handled by Horner's rule in $2^{-s}$ over the dyadic ranges
of the odd part.

The implementation rearranges this so that all the bulk work goes
through vector operations of the generic rings:

- the odd numbers are processed in blocks of 4096 with a segmented
  sieve that carries, for each small prime, its next multiple and the
  cofactor of that multiple, so that no divisions are needed;
- the primes of a block are evaluated together: `_gr_vec_log` in the
  phase ring $P$, the phase $(t/2\pi)\log p$ reduced by subtracting a
  nearby integer exactly, rounded to the term ring $T$ with
  `_gr_vec_set_other`, then `_gr_vec_sin_cos`, `_gr_vec_exp` and
  `_gr_vec_mul`; on the critical line $\sigma = 1/2$ the magnitudes
  $p^{-1/2}$ come from `_gr_vec_rsqrt` instead of the exponential;
- the composite terms of a block are gathered from the table with
  `_gr_vec_gather` (real and imaginary parts interleaved, so that a
  lookup touches one cache line), multiplied with four `_gr_vec_mul`,
  and scattered back with `_gr_vec_scatter`;
- the terms of each dyadic range are summed pairwise with
  `_gr_vec_add`.

The phase ring $P$ needs about $\log_2(t \log N)$ more bits than the
term ring $T$, because the phase $t \log p$ is reduced modulo $2\pi$.

## Bounded memory: the two-phase version {#two-phase}

The table of $f(k)$ for $k \le N/3$ takes $8N/3$ bytes per double of
the table ring, which limits the plain algorithm. With a table bound
$N^{2/3} \le M < N/3$, the odd $k \le M$ are processed as above, and
among the odd $c > M$ only the *roots* are computed: the primes and the
composites with cofactor $c/d \le M$. Every $n \le N$ whose odd part
exceeds $M$ is uniquely $n = mc$ with $c$ a root and $m \le N/c < N/M$
an integer whose odd part has no prime factor above $d$ (the chain of
divisions of $n$ by its smallest prime factor reaches $c$ as its last
element above $M$). Hence
$$\sum_{n \le N} f(n) = \sum_j f(2)^j \sum_{\substack{k \text{ odd} \\ k \le \min(M, N/2^j)}} f(k)
+ \sum_{c \text{ root}} f(c)\, G(\lfloor N/c \rfloor, d_c),$$
where $G(z, p)$ is the sum of $f(m)$ over the $m \le z$ whose odd part
has no prime factor above $p$, and $d_c = \infty$ for a prime $c$. With
$Z = \lfloor N/(M+1) \rfloor \le N^{1/3}$, the values of $G$ form a table
of about $Z^2/(2\log Z)$ entries. The second phase gathers
$G(\lfloor N/c \rfloor, d_c)$ for each root and needs one complex product
per root; all other composites (all but a fraction of about
$1.12/\log(N/M)$) are never computed.

The default table bound (`M = 0`) is the full table for $N < 2048$ and
$\lceil N^{2/3} \rceil$ otherwise. Besides bounding the memory by
$O(N^{2/3})$, this is faster, since the table fits in the cache and
fewer products are needed: at $N = 10^8$, with $(T, P) = (1, 3)$, the
time per term drops from 24.6 to 15.8 ns and the memory from 262 MB to
11 MB.

The second phase, all but $O(N^{2/3})$ of the work, is split into
chunks of at least $2^{18}$ odd numbers, each with its own sieve state,
and run with `flint_parallel_do` on `flint_get_num_threads` threads.
The partial sums are added in the order of the chunks, so the result
does not depend on the number of threads.

## Error bounds {#errors}

With `F = NULL`, every term is a ball and the result is rigorous by
construction. With a plain ring $F$, only the prime terms are balls, and
the error of the composite terms and block sums is bounded a priori and
added per dyadic range. With $\varepsilon$ the largest relative radius
of a prime term (obtained with `_gr_vec_get_interval_mid_rad`) and
$\mu = \sqrt{2}(\delta_m + \delta_a + \delta_m \delta_a)$ the relative
error of a plain complex product, a term $f(n)$ with $\Omega(n)$ prime
factors has a relative error of at most
$(1 + \varepsilon)^{\Omega}(1 + \mu)^{\Omega - 1} - 1$, where
$\Omega \le \log_3 N$; the pairwise sums of at most 4096 terms add
$(1 + \delta_a)^{12} - 1$ times the sum of the magnitudes, which is
bounded by an integral. In the two-phase version, the sum over the roots
of $g(c) = c^{-\sigma}$ times the sum of $g(m)$ over the corresponding
$m$ is exactly the sum of $g(n)$ over the $n \le N$ with odd part above
$M$, so the bound is as tight as in the first phase.

The `dfloat` instantiation uses $\delta_m = \delta_a = u$ for `d1` and
$\delta_m = 8u^2$, $\delta_a = 3u^2$ for `d2`, $u = 2^{-53}$, from the
error analysis of the kernels
([](algorithms/floating-point/double-word)); the phase context is
strong (canonical), so that rounding the reduced phases to $T$ is a
truncation.

## Choosing the precisions {#precision}

`dfloat_powsum_sieved` chooses $T$ so that the terms have
$b = \mathit{prec} + \log_2(N)/2 + 8$ bits (the accumulated error grows
like $\sqrt{N} u^T$, so that the result has about `prec` bits after the
binary point), and $P$ so that the phase $t \log N$ leaves $b + 6$ bits
after the binary point. With $P = 4$, a phase short by up to 8 bits is
accepted. It returns 0 if no pair $(T, P)$ in $\{(1, 2), (1, 3), (2, 3),
(2, 4), (3, 4)\}$ fits, if $|\operatorname{Re}(s)| > 2^{10}$,
$|t| \ge 10^{60}$ or $N \ge 2^{50}$, if the tables would take more than
`DFLOAT_POWSUM_MAX_BYTES`, or if `dfloat_is_supported` returns 0. The
plain composites are used for $T \le 2$; `_dfloat_powsum_sieved_ball`
keeps every term a ball (and is the only variant for $T = 3$): its radii
are about 10 times smaller and it is 20 to 30% slower.

## Performance {#performance}

At the height of the $10^{15}$-th zero ($t \approx 2.1 \cdot 10^{14}$,
$N = 5\,760\,732$, one core), `acb_dirichlet_powsum_sieved` at 179 bits
takes 2.4 s, while the full-table `dfloat` version takes 0.13 s with
$(1, 2)$ (radius $2^{-30}$), 0.27 s with $(2, 3)$ (radius $2^{-79}$)
and 0.60 s with $(3, 4)$ (radius $2^{-133}$); the default two-phase
version is another 25 to 30% faster. About half of the time goes to the
prime terms, a quarter to the composites (dominated by cache misses of
the table lookups), and the rest to the sieve and the sums. At low
heights the fixed setup cost of a few microseconds makes the `arb_t`
version faster for $N < 16$, which is why `acb_dirichlet_zeta_rs_r`
uses the `dfloat` version only from $N = 16$ on.
