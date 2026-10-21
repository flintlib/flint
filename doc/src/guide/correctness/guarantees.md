---
title: What FLINT guarantees
summary: Provably correct algorithms as the rule; the documented exceptions (probable primes, Monte Carlo methods, conjectures, floating point); randomised algorithms with certified output; and answers that may be "unknown".
order: 1
---

The design rule of FLINT is that a function returns the mathematically
correct result for every valid input, using algorithms with a proof of
correctness. Heuristics are used everywhere, but only to choose between
algorithms or parameters: tuning cutoffs, `double` estimates of the
number of terms of a series or of the number of primes needed, random
evaluation points. They affect the running time, not the result. The
sections below describe the exceptions to the rule, which are always stated in
the documentation of the function concerned, and the cases where a
function is allowed to answer "unknown" instead.

Correctness is meant for valid input in the sense of
[](guide/correctness/input-validation), within the size limits of
[](guide/correctness/limits), and in the absence of bugs; see
[](guide/correctness/testing).

## Primes and probable primes {#primes}

Primality is the most visible exception. Many functions come in a
proved and a probable version, and the probable version is used when a
small risk is acceptable in exchange for speed:

| Proved | Not proved |
|---|---|
| `n_is_prime`: deterministic for every word-size $n$ | `n_is_probabprime`: obsolete, now the same as `n_is_prime` |
| `n_ll_is_prime`: proves or disproves primality of two-limb integers up to about 81 bits, and returns $-1$ (unknown) above | |
| `fmpz_is_prime`: a proof (Pocklington, Morrison or Brillhart–Lehmer–Selfridge tests, APRCL or ECPP) | `fmpz_is_probabprime`: exact up to 64 bits, and up to about 81 bits where `n_ll_is_prime` decides; above that the Baillie–PSW test |
| `fmpzi_is_prime` | `fmpzi_is_probabprime` |
| `fmpz_nextprime(res, n, 1)`, `fmpz_randprime(f, state, bits, 1)` | the same with `proved = 0`, which for large numbers use GMP's `mpz_nextprime` |
| `fmpz_factor`: every factor is proved prime | `fmpz_factor_smooth(..., proved = 0)`: factors are probable primes |

The Baillie–PSW test [@BaiWag1980] [@PSW1980] combines a strong
probable prime test to base 2 with a strong Lucas test. No composite
number is known to pass it, and none exists below $2^{64}$
[@FeiGal2013], but there is no proof that none exists, and heuristic
arguments suggest that infinitely many do. For word-size integers this
does not matter: `n_is_prime` combines strong probable prime tests
with bases chosen (through a hash table) so that the test is correct
for all $n < 2^{64}$, which is why `n_is_probabprime` became
redundant. See [](algorithms/number-theory/primality-word),
[](algorithms/number-theory/strong-probable-prime) and
[](algorithms/number-theory/primality-proving).

```c example name=guide_correctness_primes
#include <flint/flint.h>
#include <flint/ulong_extras.h>
#include <flint/fmpz.h>

int main(void)
{
    fmpz_t n;
    fmpz_init(n);

    /* a word-size prime: the answer of n_is_prime is a proof */
    flint_printf("n_is_prime(4294967291) = %d\n", n_is_prime(UWORD(4294967291)));

    /* the Mersenne prime 2^89 - 1 */
    fmpz_one(n);
    fmpz_mul_2exp(n, n, 89);
    fmpz_sub_ui(n, n, 1);
    flint_printf("fmpz_is_probabprime: %d (Baillie-PSW, not a proof)\n",
        fmpz_is_probabprime(n));
    flint_printf("fmpz_is_prime:       %d (proved)\n", fmpz_is_prime(n));

    /* 2^89 + 1 is divisible by 3 */
    fmpz_add_ui(n, n, 2);
    flint_printf("2^89 + 1: %d %d\n", fmpz_is_probabprime(n), fmpz_is_prime(n));

    fmpz_clear(n);
    flint_cleanup_master();
    return 0;
}
```

```output
n_is_prime(4294967291) = 1
fmpz_is_probabprime: 1 (Baillie-PSW, not a proof)
fmpz_is_prime:       1 (proved)
2^89 + 1: 0 0
```

A "composite" answer is always certain: a number that fails a probable
prime test is composite. Only the answer "prime" of the probable
versions lacks a proof.

## Other documented exceptions {#exceptions}

- *Monte Carlo methods.* Some multimodular algorithms can stop as soon
  as the reconstructed result has not changed for a while, instead of
  using enough primes to cover a proven bound. They take a flag
  `proved`: `fmpz_mat_det_modular`,
  `fmpz_mat_det_modular_accelerated` and
  `fmpz_mat_det_modular_given_divisor` with `proved = 0`, and
  `gr_poly_resultant_modular` with `proved = 0`, return a result that
  is correct with high probability. The main functions (`fmpz_mat_det`,
  `fmpz_poly_resultant`, ...) always pass `proved = 1`.
- *Heuristic results.* `qfb_exponent` returns the exponent of a class
  group with high probability, without a proof.
- *Conjectures.* `qfb_exponent_grh` is correct if the generalised
  Riemann hypothesis holds.
- *Bounds checked numerically.* A few auxiliary bounds are verified for
  all word-size arguments but not proved in general, for example
  `arith_bell_number_size`.
- *Floating point.* The floating-point types `arf_t`, `nfloat`, `acf_t`
  and `double` (the `d_vec` and `d_mat` modules) compute approximations
  without error bounds, as do the functions with `approx` in their
  name in the ball modules, such as `arb_mat_approx_solve` and
  `acb_mat_approx_eig_qr`. Each `arf_t` operation is correctly rounded,
  but errors accumulate over a computation as in any floating-point
  arithmetic. See [](algorithms/floating-point/overview).
- *Assumed properties.* Some functions take a flag that lets them assume
  a property of the input instead of checking it, for example
  `fmpz_mpoly_sqrt_heap` with `check = 0`, which assumes that the input
  is a perfect square. With an input that lacks the property the result
  is meaningless, not just unproved.

## Randomised algorithms {#randomised}

Many algorithms in FLINT make random choices. Almost all of them are
*Las Vegas* algorithms: the random choices influence the running time,
but the result is always correct, either because a wrong choice is
detected and the computation repeated, or because the result is
verified at the end. Examples are the factorisation of polynomials over
finite fields (Cantor–Zassenhaus, Berlekamp and Kaltofen–Shoup, all via
`gr_poly_factor_finite_field`), root finding over finite fields, the
modular and sparse interpolation algorithms for multivariate GCDs (bad
primes and bad evaluation points are detected), and integer factoring
with ECM and the quadratic sieve, whose factors are proved prime by
`fmpz_factor`. *Monte Carlo* algorithms, which may return a wrong result
with small probability, are only used where documented, as listed
above.

The random states these functions create internally start from the
fixed seed of `flint_rand_init`, so a computation is reproducible: the
same input gives the same output. Functions that take a `flint_rand_t` are
equally deterministic for a given state.

Heuristic methods with a certified result follow the same pattern.
`fmpz_lll` tries `double` arithmetic and heuristic dot products first
and accepts the result only when `fmpz_lll_is_reduced` proves that the
basis is reduced, falling back to multiprecision arithmetic otherwise.
`fmpz_poly_gcd_heuristic` returns a GCD only after checking that it
divides both inputs. The rational solvers of `fmpz_mat` guess the
denominator of the solution from a random linear combination of its
entries, and accept the solution only when a norm bound proves it.

## Rigorous numerics {#balls}

Real and complex numbers are computed with ball arithmetic
([](guide/types/ball-arithmetic)): every result is a ball that is guaranteed
to contain the exact value. This guarantee is the same kind of guarantee
as for exact computations, but the information content of the answer
varies: a ball may be too wide to be useful, or even contain all real
numbers. Predicates on balls return 1 only when the property is certain
and 0 when it is false *or* cannot be decided at the current precision;
see [](guide/types/ball-arithmetic#predicates). The `double` functions
of [](arb_fpwrap) build on this: a result they return is accurate to
within one unit in the last place (or correctly rounded on request),
and when that accuracy cannot be reached they return NaN with the
status `FPWRAP_UNABLE` instead of an inaccurate value.

## Hardware floating point {#floating-point}

Many parts of FLINT use `double` arithmetic, under different
assumptions.

- *Heuristics.* `double` arithmetic and `libm` functions choose
  parameters, such as the number of terms of a series, or compute
  approximations that are verified afterwards, such as the starting
  values of the root finder `acb_poly_find_roots` or the
  floating-point stage of `fmpz_lll`. A flawed floating-point
  implementation can then only make the computation slower or make it
  fail. Results computed this way may differ slightly between
  platforms, for example with and without fused multiply-adds; a
  rigorous enclosure may have a different radius, but it is still an
  enclosure.
- *Rigorous bounds in the ball modules.* `mag_t` arithmetic and some
  error bounds of `arb` (such as `arb_sin_cos_wide`) are computed with
  doubles under weak assumptions: `double` is IEEE 754 binary64; the
  compiler does not reorder floating-point operations; basic operations
  have an error of at most 1.1 units in the last place in any rounding
  mode (so x87 arithmetic and a changed rounding mode are allowed); only
  normal numbers and explicit infinities are used; and `libm`
  transcendental functions are never used for bounds (logarithms are
  bounded with `mag_d_log_upper_bound` and `mag_d_log_lower_bound`).
- *Exact or correctly rounded arithmetic.* Other code relies on IEEE
  arithmetic exactly: the floating-point FFTs of [](fft_small), the
  `_precomp` functions of [](ulong_extras), the floating-point matrix
  kernels, and the [`dfloat`](dfloat) numbers and balls (`d1b_t` to
  `d4b_t`), which also require round-to-nearest, a correctly rounded
  `fma`, evaluation of doubles in double precision and no contraction
  of expressions into fused multiply-adds. `dfloat_is_supported`
  checks these conditions where possible, and the `dfloat` functions
  must not be used where it returns 0.

FLINT must therefore not be compiled with unsafe floating-point
optimisations such as `-ffast-math`. These are off by default in GCC
and Clang, but some other compilers, such as the classic Intel
compiler, enable them by default and must be told not to.

## Unknown answers {#unknown}

Some questions cannot be decided, or not decided at reasonable cost, and
the functions answering them have a way to say so:

- The generic ring functions of [](gr) return a status code:
  `GR_SUCCESS`, `GR_DOMAIN` when the operation is mathematically
  undefined (division by zero, $1/2$ in $\Z$), and `GR_UNABLE` when the
  implementation cannot compute the result, for example the inverse of
  a real ball that contains zero. Predicates return a `truth_t`, one of
  `T_TRUE`, `T_FALSE` and `T_UNKNOWN`. See [](gr#status) and
  [](gr#truth).
- The exact real and complex numbers of [](ca) represent numbers
  symbolically, and equality of such numbers is not decidable in
  general. Predicates such as `ca_check_is_zero` return `T_UNKNOWN`
  when no proof is found; see [](ca/check-is-zero#unknown) and
  [](algorithms/exact/zero-testing#in-practice). The
  algebraic numbers of [](qqbar) are an exception: every question about
  them is decided.
- Other functions that may fail document which return value means
  failure, such as $-1$ for `n_ll_is_prime` above, or 0 for
  `fmpz_poly_gcd_heuristic` and for `fmpz_factor_smooth` when the
  factorisation is incomplete.

A program must treat "unknown" as a third possibility. A common mistake
is to test `if (gr_is_zero(x, ctx) == T_FALSE)` in one branch and to
assume `x` is zero in the other.
