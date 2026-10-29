---
title: ulong_extras.h
header: ulong_extras.h
description: Arithmetic and number theory for single-word unsigned integers.
category: Word-size integers
order: 1
---

This module works with `ulong`, an unsigned integer of `FLINT_BITS` bits
(64 on all common platforms). It provides division and modular arithmetic
with precomputed inverses, primality testing, factoring, integer roots and
elementary number-theoretic functions. Most higher-level modules
(`nmod`, `fmpz`, the polynomial and matrix modules) rely on it.

For background, [](algorithms/arithmetic/choosing) compares the ways of computing modulo a
word-size integer, and [](algorithms/number-theory/primality-word) and
[](algorithms/number-theory/factoring-word) describe the algorithms behind
`n_is_prime` and `n_factor`.

## Conventions {#conventions}

Throughout, $B = 2^{\bits}$.

- Functions ending in `_preinv` take a precomputed inverse of type `ulong`
  computed by `n_preinvert_limb`. Functions ending in `_precomp` take a
  `double` inverse computed by `n_precompute_inverse`.
- A `2` in the name (`n_mod2_preinv`, `n_mulmod2_preinv`) marks the variant
  with no size restrictions on its inputs. Variants without it are faster
  but may require reduced operands, a normalised modulus or a modulus
  below $2^{53}$.
- `n_ll_` and `n_lll_` functions take operands of two or three limbs.
- Precomputed inverses are assumed to match the modulus, and operands of
  modular functions are assumed to be reduced. None of this is checked.

::: summary algorithms/arithmetic/choosing :::

## Basic arithmetic {#basic}

::: function
ulong n_mulhi(ulong a, ulong b)

Returns the high word $\lfloor ab / B \rfloor$ of the double-word product
$ab$.
:::

::: function
int n_add_checked(ulong * a, ulong b, ulong c)
int n_sub_checked(ulong * a, ulong b, ulong c)
int n_mul_checked(ulong * a, ulong b, ulong c)

Sets `*a` to $b + c$, $b - c$ or $bc$ modulo $B$ and returns 1 if the
exact result does not fit in a `ulong`, otherwise 0.
:::

::: function
ulong n_pow(ulong n, ulong exp)

Returns $n^{\mathit{exp}}$, with $0^0 = 1$. Overflow is not detected.
:::

::: function
ulong n_flog(ulong n, ulong b)
ulong n_clog(ulong n, ulong b)

Returns $\lfloor \log_b n \rfloor$ and $\lceil \log_b n \rceil$
respectively. Requires $n \ge 1$ and $b \ge 2$.
:::

::: function
ulong n_clog_2exp(ulong n, ulong b)

Returns $\lceil \log_b 2^n \rceil$. Requires $b \ge 2$.
:::

::: function
ulong n_revbin(ulong n, ulong b)

Returns $n$ with its lowest $b$ bits in reverse order; for example the
6-bit number $010110_2$ becomes $011010_2$.
:::

::: function
int n_sizeinbase(ulong n, int base)
slong n_nonzero_sizeinbase10(ulong n)

Returns the number of digits of $n$ in the given base
($2 \le \mathit{base} \le 36$); `n_sizeinbase` returns 1 for $n = 0$.
`n_nonzero_sizeinbase10` is a branch-free version for base 10 and
requires $n \ne 0$.
:::

## Integer roots and perfect powers {#roots}

::: function
ulong n_sqrt(ulong a)
ulong n_sqrtrem(ulong * r, ulong a)

Returns $\lfloor \sqrt a \rfloor$. `n_sqrtrem` also sets `*r` to
$a - \lfloor \sqrt a \rfloor^2$.

::: details
`n_sqrt` computes the square root in IEEE double precision. If $a$ is one below a
square the rounding may be up, giving a result one too large, which is then
corrected; care is taken that squaring the too-large value does not
overflow.
:::
:::

::: function
int n_is_square(ulong x)

Returns 1 if $x$ is a perfect square, otherwise 0.

::: details
`n_is_square` first checks that $x$ is a square modulo 64, $63 = 3^2 \cdot 7$ and
$65 = 5 \cdot 13$ using lookup tables; only then is a square root taken.
:::
:::

::: function
ulong n_root(ulong n, ulong root)
ulong n_rootrem(ulong * remainder, ulong n, ulong root)

Returns $r = \lfloor n^{1/k} \rfloor$ for $k = \mathit{root}$, or 0 if $n$
or $k$ is 0. `n_rootrem` also sets the remainder $n - r^k$.
:::

::: function
ulong n_cbrt(ulong n)
ulong n_cbrtrem(ulong * remainder, ulong n)

Returns $\lfloor n^{1/3} \rfloor$; `n_cbrtrem` also sets the remainder.

::: details
`n_cbrt` takes a first approximation from `n_cbrt_estimate`. For $n > 2^{46}$ the
result is computed with `n_cbrt_chebyshev_approx`; otherwise the iteration
$x \leftarrow x - (x^3 - a) x / (2 x^3 + a)$ of Kahan [@Kahan1991] is used.
:::
:::

::: function
double n_cbrt_estimate(double a)
ulong n_cbrt_binary_search(ulong n)
ulong n_cbrt_chebyshev_approx(ulong n)

Alternative cube root algorithms used by `n_cbrt`: a bit-manipulation
estimate in the spirit of the fast inverse square root, a binary search,
and a Chebyshev approximation (requires $n \ne 0$).

::: details
`n_cbrt_chebyshev_approx` writes $n = x \cdot 2^e$ with
$x \in [1/2, 1]$, evaluates a Chebyshev approximation of $x^{1/3}$
and multiplies by
$2^{\lfloor e/3 \rfloor}$ and the cube root of $2^{e \bmod 3}$.
:::
:::

::: function
int n_is_perfect_power235(ulong n)
ulong n_factor_power235(ulong * exp, ulong n)

`n_is_perfect_power235` returns whether $n$ is a square, cube or fifth
power. `n_factor_power235` returns the root $r$ and sets `*exp` to
$k \in \{2, 3, 5\}$ if $n = r^k$, and otherwise returns 0.

::: details
Both functions first apply a sequence of modular tests, which rejects most
inputs: each test looks up a 3-bit mask saying whether $n$ can be a square,
cube or fifth power modulo some small modulus, and the masks are combined
with bitwise AND. For any bit that
survives, the corresponding root is computed and checked exactly.
:::
:::

::: function
int n_is_perfect_power(ulong * root, ulong n)

If $n = r^k$ with $k \ge 2$, returns $k$ and sets `*root` to $r$; otherwise
returns 0. Both 0 and 1 count as squares. The $k$ found is not necessarily
the largest possible.
:::

## Division and remainder {#division}

These functions compute $a \bmod n$ (and possibly the quotient) from a
precomputed approximation of $1/n$, which is much faster than a hardware
division when the same $n$ is used repeatedly. The approximation is either
a `ulong` (`_preinv` functions), a `double` (`_precomp` functions) or the
data for Barrett or Lemire reduction. For a comparison of the methods see
[](algorithms/arithmetic/choosing).

### With a ulong inverse {#division-preinv}

::: function
ulong n_preinvert_limb(ulong n)
ulong n_preinvert_limb_prenorm(ulong n)

Returns the precomputed inverse of $n \ne 0$ used by all `_preinv`
functions [@GraMol2010]. `n_preinvert_limb_prenorm` is a faster version for normalised
$n$ (top bit set).
:::

::: function
ulong n_mod2_preinv(ulong a, ulong n, ulong ninv)

Returns $a \bmod n$. No restrictions on $a$; $n \ne 0$.
:::

::: function
ulong n_divrem2_preinv(ulong * q, ulong a, ulong n, ulong ninv)
ulong n_divrem_preinv(ulong * q, ulong a, ulong n, ulong ninv, unsigned int norm)
ulong n_divrem_preinv_unnorm(ulong * q, ulong a, ulong n, ulong ninv, unsigned int norm)

Returns $a \bmod n$ and sets `*q` to $\lfloor a/n \rfloor$. The last two
variants take the number of leading zero bits of $n$ as `norm`;
`n_divrem_preinv_unnorm` requires $\mathit{norm} > 0$.
:::

::: function
ulong n_div2_preinv(ulong a, ulong n, ulong ninv)

Returns the quotient $\lfloor a/n \rfloor$. No restrictions on $a$;
$n \ne 0$.
:::

::: function
ulong n_ll_mod_preinv(ulong a_hi, ulong a_lo, ulong n, ulong ninv)

Returns $a \bmod n$ for the two-limb integer $a = a_{hi} B + a_{lo}$.
No restrictions on $a$; $n \ne 0$.
:::

::: function
ulong n_lll_mod_preinv(ulong a_hi, ulong a_mi, ulong a_lo, ulong n, ulong ninv)

Returns $a \bmod n$ for the three-limb integer
$a = a_{hi} B^2 + a_{mi} B + a_{lo}$. Requires $a_{hi} < n$.
:::

::: function
ulong n_divrem_norm(ulong * q, ulong a, ulong n)

Returns $a \bmod n$ and sets `*q` to the quotient (0 or 1), assuming $n$ is
normalised. No inverse is needed: the quotient is a single comparison.
:::

### With a double inverse {#division-precomp}

::: function
double n_precompute_inverse(ulong n)

Returns $1/n$ as a `double`, the precomputed inverse used by all `_precomp`
functions. Requires $n > 0$.
:::

::: function
ulong n_mod_precomp(ulong a, ulong n, double ninv)

Returns $a \bmod n$. Requires $n < 2^{53}$, $a < n^2$ and $a < 2^{\bits - 1}$.

::: details
`n_mod_precomp` takes the quotient to be `a * ninv` truncated to an
integer, and corrects the remainder by at most one multiple of $n$. This
suffices: the processor is assumed to round to nearest, so `ninv` is
correct to 53 bits (call the weight of the last bit a *place*) and is at most half a
place out. Since $a$ is exact and the mantissa of `ninv` is below 2, the
product `a * ninv` is at most one place out; truncating it to an integer
gives either the exact quotient or one less, because half a place is less
than 1 when the product is below $n < 2^{53}$. The only other case is that
the product was rounded up to an integer, which implies that $a/n$ lies
within $2^{-54}$ of an integer; this is handled by subtracting 1 from the
quotient. The computed quotient is therefore exact or one too small, and a
single correction step suffices.
:::
:::

::: function
ulong n_mod2_precomp(ulong a, ulong n, double ninv)

Returns $a \bmod n$ with no restrictions on $a$ or $n$.

::: details
For $n < 2^{53}$ and $a < n^2$, `n_mod2_precomp` computes the quotient as
`n_mod_precomp` does. Otherwise the top 52 bits of the quotient are correct; the
remainder is then multiplied by `ninv` again to obtain a correction, which
may be one out in either direction and is fixed up at the end.
:::
:::

::: function
ulong n_divrem2_precomp(ulong * q, ulong a, ulong n, double npre)

Returns $a \bmod n$ and sets `*q` to the quotient. No restrictions on $a$
or $n$.
:::

### Barrett and Lemire reduction {#division-other}

::: function
ulong n_barrett_precomp(ulong n)

Returns the precomputed inverse for `n_mod_barrett`. Requires $n \ge 2$.
:::

::: function
ulong n_mod_barrett(ulong x, ulong n, ulong npre)
ulong n_mod_barrett_lazy(ulong x, ulong n, ulong npre)

Returns $x \bmod n$ using Barrett reduction [@Bar1987]. The `lazy` version returns a
residue in $[0, 2n)$ instead of $[0, n)$. No restrictions on $x$. Both
support $n = 1$ if `npre` is set to `UWORD_MAX`.
:::

::: function
ulong n_lemire_precomp(ulong n)

Returns the precomputed inverse for `n_mod_lemire`. Requires $n \ge 1$.
:::

::: function
ulong n_mod_lemire(ulong x, ulong n, ulong npre)

Returns $x \bmod n$ using the method of Lemire, Kaser and Kurz
[@LKK2019]. Only
correct for $n, x < 2^{\bits / 2}$.
:::

## Divisibility {#divisibility}

::: function
int n_divides(ulong * q, ulong n, ulong p)

Returns 1 and sets `*q` to $n/p$ if $p$ divides $n$; otherwise returns 0 and
sets `*q` to 0.
:::

::: function
int n_divisible_odd_gm(ulong n, ulong inv1, ulong inv2)

Returns whether the odd number $d$ divides $n$, given
$\mathit{inv1} = d^{-1} \bmod B$ and
$\mathit{inv2} = \lfloor (B - 1)/d \rfloor$. Costs one multiplication.
See [](algorithms/arithmetic/divisibility-test).
:::

::: function
int n_remove(ulong * n, ulong p)
int n_remove2_precomp(ulong * n, ulong p, double ppre)

Divides `*n` by the largest possible power of $p > 1$ and returns the
exponent. Requires $n \ne 0$. `n_remove2_precomp` takes
`ppre = n_precompute_inverse(p)`.

::: details
For $p = 2$ trailing zeros are counted. Otherwise `n_remove` repeatedly
squares $p$, storing a table of powers and removing the highest power at
each step, then walks back down the table. `n_remove2_precomp` repeatedly
calls `n_divrem2_precomp`.
:::
:::

## Greatest common divisor {#gcd}

::: function
ulong n_gcd(ulong x, ulong y)

Returns $\gcd(x, y)$, with $\gcd(x, 0) = \gcd(0, x) = x$.

::: details
`n_gcd` removes common factors of 2 and calls GMP's `mpn_gcd_11` (or
`mpn_gcd_1` when GMP internals are not available), a binary algorithm
[@Knu1997].
:::
:::

::: function
ulong n_xgcd(ulong * a, ulong * b, ulong x, ulong y)

Returns $g = \gcd(x, y)$ and sets `*a`, `*b` to unsigned values with
$ax - by = g$. Requires $x \ge y$.

::: details
`n_xgcd` uses the extended Euclidean algorithm, which always produces
cofactors $\lvert a \rvert < x/2$ and $\lvert b \rvert < x/2$, apart from a
few degenerate cases. The proof is by induction. Suppose that at some step
$x_n = q y_n + r$ with $r \ge 1$ and $1 = s y_n - t r$ with $s < y_n/2$
and $t < y_n/2$. Then $1 = s y_n - t (x_n - q y_n) = (s + tq) y_n - t x_n$,
and it suffices to show $s + tq < x_n/2$, since $t < y_n/2 < x_n/2$. At the
previous step of the back-substitution we had $1 = s r - c d$ with
$s < r/2$ and $c < r/2$, so
$s + tq < r/2 + q y_n/2 = (r + q y_n)/2 = x_n/2$.
:::
:::

## Modular addition and inversion {#modular}

::: function
ulong n_addmod(ulong a, ulong b, ulong n)
ulong n_submod(ulong a, ulong b, ulong n)
ulong n_negmod(ulong a, ulong n)

Returns $(a + b) \bmod n$, $(a - b) \bmod n$ and $-a \bmod n$ for reduced
operands.
:::

::: function
ulong n_invmod(ulong x, ulong y)

Returns $x^{-1} \bmod y$. Throws an exception if the inverse does not exist.
:::

::: function
ulong n_gcdinv(ulong * a, ulong x, ulong y)

Returns $g = \gcd(x, y)$ and sets `*a` to the $0 \le a < y$ with
$ax \equiv g \pmod y$. Requires $x < y$. When $y = 1$, returns 1 and sets
`*a` to 0.
:::

::: function
ulong n_binvert(ulong n)

Returns the inverse of odd $n$ modulo $B$; the result is meaningless for
even $n$.

::: details
`n_binvert` is an inline function using Hurchalla's variant of Newton
iteration
[@Hur2022], starting from an 8-bit inverse read from a 128-byte
table.
:::
:::

## Modular multiplication {#mulmod}

::: function
ulong n_mulmod2_preinv(ulong a, ulong b, ulong n, ulong ninv)

Returns $ab \bmod n$. No restrictions on $a$, $b$ or $n$.
:::

::: function
ulong n_mulmod2(ulong a, ulong b, ulong n)

Returns $ab \bmod n$, computing the inverse of $n$ on the fly. Use
`n_mulmod2_preinv` when $n$ is reused.
:::

::: function
ulong n_mulmod_preinv(ulong a, ulong b, ulong n, ulong ninv, ulong norm)

Returns $\lfloor a / 2^{\mathit{norm}} \rfloor \, b \bmod n$ for
normalised $n$ and reduced $a, b$. If $a$, $b$ and $n$ are the operands and
modulus shifted left by `norm` bits, this is their product modulo the
unshifted modulus, shifted left by `norm` bits.
:::

::: function
ulong n_mulmod_precomp(ulong a, ulong b, ulong n, double ninv)

Returns $ab \bmod n$. Requires $n < 2^{53}$ and $a, b < n$.

::: details
`n_mulmod_precomp` computes the product $ab$ in floating point with an
error of at most half a place (as defined for `n_mod_precomp`), and multiplying by `ninv` gives an approximate quotient within two
places of the exact one. Since the quotient is below $n < 2^{53}$, after
truncation the computed integer quotient is at most two above and one below
the true quotient; an argument like the one for `n_mod_precomp` excludes
the case of two below. The remainder is corrected accordingly.
:::
:::

## Multiplication by a fixed operand {#shoup}

These functions speed up computing $a b_i \bmod n$ for fixed $a$ and $n$
and many $b_i$, as in scalar multiplication of vectors and matrices
[@Sho2021; @Har2014].
They require $n < 2^{\bits - 1}$.

::: function
ulong n_mulmod_precomp_shoup(ulong a, ulong n)

Returns $\lfloor a B / n \rfloor$, the precomputed data for
`n_mulmod_shoup`. Requires $a < n$.
:::

::: function
ulong n_mulmod_shoup(ulong a, ulong b, ulong a_precomp, ulong n)

Returns $ab \bmod n$ given `a_precomp` from `n_mulmod_precomp_shoup`.
Requires $a < n$; $b$ is arbitrary.
:::

::: function
void n_mulmod_precomp_shoup_quo_rem(ulong * a_pr_quo, ulong * a_pr_rem, ulong a, ulong n)
ulong n_mulmod_precomp_shoup_rem_from_quo(ulong a_pr_quo, ulong n)

Compute both the quotient $\lfloor aB/n \rfloor$ and the remainder
$aB \bmod n$, as needed by `n_mulmod_and_precomp_shoup`. The second
function recovers the remainder from a known quotient.
:::

::: function
ulong n_mulmod_precomp_shoup_preinv(ulong a, ulong n, ulong ninv)
void n_mulmod_precomp_shoup_quo_rem_preinv(ulong * a_pr_quo, ulong * a_pr_rem, ulong a, ulong n, ulong ninv)

Same as `n_mulmod_precomp_shoup` and `n_mulmod_precomp_shoup_quo_rem`,
given `ninv` from `n_preinvert_limb(n)` (the `ninv` of an `nmod_t`). They
use multiplications instead of a hardware division: much faster where the
division is slow (e.g. Skylake or Cascade Lake, and 64-bit ARM), about as fast
on recent x86-64 processors. Require $a < n$; any $n \ne 0$.
:::

::: function
void n_mulmod_and_precomp_shoup(ulong * ab, ulong * ab_precomp, ulong a, ulong b, ulong a_pr_quo, ulong a_pr_rem, ulong b_precomp, ulong n)

Sets `*ab` to $ab \bmod n$ and `*ab_precomp` to its precomputed data
$\lfloor (ab \bmod n) B / n \rfloor$, without a division. Requires
$a, b < n$. Useful for building tables of powers $a, a^2, a^3, \ldots$
together with their precomputed data.
:::

## Modular exponentiation {#powmod}

::: function
ulong n_powmod2_preinv(ulong a, slong exp, ulong n, ulong ninv)
ulong n_powmod2(ulong a, slong exp, ulong n)

Returns $a^{\mathit{exp}} \bmod n$ for $a < n$. Negative exponents are
allowed; an exception is thrown if $a$ is then not invertible.
:::

::: function
ulong n_powmod2_ui_preinv(ulong a, ulong exp, ulong n, ulong ninv)
ulong n_powmod_ui_preinv(ulong a, ulong exp, ulong n, ulong ninv, ulong norm)

Returns $a^{\mathit{exp}} \bmod n$ for an unsigned exponent and $a < n$.
`n_powmod_ui_preinv` requires normalised $n$, like `n_mulmod_preinv`.
:::

::: function
ulong n_powmod2_fmpz_preinv(ulong a, const fmpz_t exp, ulong n, ulong ninv)

Returns $a^{\mathit{exp}} \bmod n$ for a nonnegative multiprecision
exponent and $a < n$.
:::

::: function
ulong n_powmod_precomp(ulong a, slong exp, ulong n, double npre)
ulong n_powmod_ui_precomp(ulong a, ulong exp, ulong n, double npre)
ulong n_powmod(ulong a, slong exp, ulong n)

Returns $a^{\mathit{exp}} \bmod n$ using a `double` inverse. Require
$n < 2^{53}$ and $a < n$.
:::

## Chinese remaindering {#crt}

::: function
ulong n_CRT(ulong r1, ulong m1, ulong r2, ulong m2)

Returns the unique $0 \le x < m_1 m_2$ with $x \equiv r_1 \pmod{m_1}$ and
$x \equiv r_2 \pmod{m_2}$. Requires coprime $m_1, m_2 > 1$ with $m_1 m_2$
fitting in a `ulong`, and reduced $r_1, r_2$. Nothing is checked.
:::

## Jacobi symbols {#jacobi}

::: function
int n_jacobi(slong x, ulong y)
int n_jacobi_unsigned(ulong x, ulong y)

Returns the Jacobi symbol $\left(\frac{x}{y}\right)$ for odd $y$ (not
checked), with a binary algorithm [@ShaSor1993]. `n_jacobi_unsigned` accepts a full-limb $x$.
:::

## Primes {#primes}

Finding, counting and enumerating primes. To run through consecutive
primes, use a prime iterator (`n_primes_t`); the first primes and their
inverses are also available from a per-thread cache.

::: function
ulong n_nextprime(ulong n, int proved)

Returns the smallest prime greater than $n$. Throws an exception if it
does not fit in a `ulong`. The flag `proved` has no effect: the result is
always a proved prime.
:::

::: function
ulong n_prime_pi(ulong n)
void n_prime_pi_bounds(ulong * lo, ulong * hi, ulong n)

Return the number $\pi(n)$ of primes $p \le n$, or set
$\mathit{lo} \le \pi(n) \le \mathit{hi}$.

::: details
`n_prime_pi` looks up small $n$ in a table; otherwise it extends the
table of cached primes and performs a binary search.

`n_prime_pi_bounds` looks up small values in a table and otherwise uses the
bounds $n / \ln n < \pi(n) < 1.25506\, n / \ln n$ of Rosser and Schoenfeld
[@RosSch1962], valid for $n \ge 17$. The bit length of $n$ is used to bound
$\ln n$ from the appropriate side. If `lo` and `hi` point to the same
location, the upper bound is stored.
:::
:::

::: function
ulong n_nth_prime(ulong n)
void n_nth_prime_bounds(ulong * lo, ulong * hi, ulong n)

Return the $n$-th prime $p_n$ (with $p_1 = 2$; $n = 0$ throws an
exception), or set $\mathit{lo} \le p_n \le \mathit{hi}$.

::: details
`n_nth_prime` makes sure the table of cached primes is large enough and
looks up the entry. `n_nth_prime_bounds` does not check for overflow, so
its bounds are meaningless for large $n$. It uses the estimates, valid for
$n > 5$,

$$
\begin{aligned}
p_n &> n (\ln n + \ln \ln n - 1), \\
p_n &< n (\ln n + \ln \ln n), \\
p_n &< n (\ln n + \ln \ln n - 0.9427) \quad (n \ge 15985).
\end{aligned}
$$

The first inequality is due to Dusart [@Dus1999] and the last to Massias
and Robin [@MasRob1996].
:::
:::

::: function
void n_primes_init(n_primes_t iter)
void n_primes_clear(n_primes_t iter)

Initialise and clear a prime iterator.
:::

::: function
ulong n_primes_next(n_primes_t iter)

Returns the next prime and advances the iterator. The first call
returns 2.

::: details
`n_primes_next` reads small primes from the table `flint_primes_small`
(or from the iterator's own extended copy). When that table is exhausted,
it calls `n_primes_jump_after`, which sieves the next block of odd numbers
with `n_primes_sieve_range`.
:::
:::

::: function
void n_primes_jump_after(n_primes_t iter, ulong n)

Makes the iterator continue with the first prime greater than $n$.
:::

::: function
void n_compute_primes(ulong num_primes)
const ulong * n_primes_arr_readonly(ulong num_primes)
const double * n_prime_inverses_arr_readonly(ulong n)
void n_cleanup_primes(void)

Maintain a per-thread cache of the first primes and their `double`
inverses. The `_readonly` functions return pointers into the cache, which
stay valid until `n_cleanup_primes` is called in the same thread.
:::

::: function
void n_primes_extend_small(n_primes_t iter, ulong bound)
void n_primes_sieve_range(n_primes_t iter, ulong a, ulong b)

Internal helpers of the prime iterator: extend its table of small primes
until the second-largest entry is at least `bound`, or sieve the odd
numbers in $[a, b]$ and point the iterator at the start of that block.
:::

## Primality testing {#primality}

`n_is_prime`, and `n_ll_is_prime` for two-limb integers, are the functions
to use. The others are individual probable prime tests and two alternative
primality proofs, of interest mainly for testing and experiments.

::: summary algorithms/number-theory/primality-word :::

::: function
int n_is_prime(ulong n)

Returns 1 if $n$ is prime and 0 otherwise. Deterministic for all $n$.
:::

::: function
int n_is_prime_odd_no_trial(ulong n)

Like `n_is_prime` for odd $n$, but skips trial division. Faster when $n$ is
already known to have no small factors.
:::

::: function
int n_ll_is_prime(ulong nhi, ulong nlo)

Primality test for the two-limb integer $n = n_{hi} B + n_{lo}$ with
$n_{hi} \ne 0$. Returns 1 (prime), 0 (composite) or $-1$ (unknown). The
answer is always definite below about $2^{81}$ on 64-bit machines.

::: details
For $n$ up to about 81 bits, trial division is followed by a strong probable
prime test with the first 13 primes as bases, which proves primality in this
range [@SorWeb2016]. For larger $n$, trial division and a base-2 test are
done and the function returns 0 or $-1$; use `fmpz_is_prime` for a proof.
On 32-bit machines the function currently always returns $-1$.
:::
:::

::: function
int n_is_strong_probabprime2_preinv(ulong n, ulong ninv, ulong a, ulong d)
int n_is_strong_probabprime_precomp(ulong n, double npre, ulong a, ulong d)

Test whether odd $n$ is a strong probable prime to base $a$
[@Mil1976; @Rab1980], given the odd part $d$ of $n - 1$ and $0 < a < n$. The `_precomp` version requires
$n < 2^{53}$. See [](algorithms/number-theory/strong-probable-prime).
:::

::: function
int n_is_probabprime_fermat(ulong n, ulong i)

Returns whether $i^{n-1} \equiv 1 \pmod n$. Requires $1 < i < n$ and
$i \nmid n$.
:::

::: function
int n_is_probabprime_fibonacci(ulong n)

Fibonacci probable prime test: checks $F_{n - (n/5)} \equiv 0 \pmod n$,
where $(n/5)$ is a Jacobi symbol [@CraPom2005]. Requires $\gcd(n, 10) = 1$.
:::

::: function
int n_is_probabprime_lucas(ulong n)

Lucas probable prime test, similar to the one of Baillie and Wagstaff
[@BaiWag1980].
:::

::: function
int n_is_probabprime_BPSW(ulong n)

Baillie–PSW probable prime test [@PSW1980; @BaiWag1980]. No composite passing it is known, and it
has been checked against tables of pseudoprimes up to $2^{64}$, so it is a
primality test for word-size $n$.

::: details
This variant of BPSW, following [@Chen2003], uses a Fibonacci test instead
of a Lucas test when $n \equiv 2, 3 \pmod 5$ (see also [@CraPom2005],
p. 143).
:::
:::

::: function
int n_is_prime_pocklington(ulong n, ulong iterations)

Pocklington–Lehmer test [@Poc1914; @BLS1975]. Returns 1 if $n$ is proved
prime, 0 if it is composite, and $-1$ if nothing was proved because the
bases $2, \ldots, \mathit{iterations} - 1$ were not enough. If
$\mathit{iterations}$ is 0 or 1, the number of bases is not limited and
$-1$ is never returned for this reason.

::: details
Most of the time goes into factoring $n - 1$ with `n_factor_partial`
(trial division, then Hart's one-line algorithm [@Har2012] or SQUFOF),
which stops once the factored part $F$ of $n - 1$ is at least
$(n - 1)^{1/3}$. If $F \le \sqrt n$, the cube root criterion of
Brillhart, Lehmer and Selfridge is applied in addition. Then, for each
prime $q$ dividing $F$, `n_is_prime_pocklington` tries the bases
$a = 2, 3, \ldots, \mathit{iterations} - 1$ until $a^{(n-1)/q} \ne 1$,
returning 0 if some $a^{n-1} \not\equiv 1 \pmod n$ and $-1$ if no base
works. Finally it checks that $n$ is coprime to the product of the values
$a^{(n-1)/q} - 1$. In most cases a larger number of iterations hardly
affects the running time.
:::
:::

::: function
int n_is_prime_pseudosquare(ulong n)

Primality test based on pseudosquares, following Theorem 2.7 of
[@LukPatWil1996].

::: details
`n_is_prime_pseudosquare` first trial-divides $n$ by up to
`FLINT_PSEUDOSQUARES_CUTOFF` primes; if no factor is found and the square
of the last prime does not exceed $n$, let $B_0$ be the next prime. The
first pseudosquare $L_p > n/B_0$ is then looked up in a static table (OEIS
A002189), and step 3 of the test checks that
$p_j^{(n-1)/2} \equiv \pm 1 \pmod n$ for the primes $p_j$ up to $p$.
Every prime passes this, and most composites fail. If
$n \equiv 3, 7 \pmod 8$, $n$ is then prime. Otherwise step 4 proves
primality: for $n \equiv 5 \pmod 8$ by checking
$2^{(n-1)/2} \equiv -1 \pmod n$, and for $n \equiv 1 \pmod 8$ by finding
some $p_j$ with $p_j^{(n-1)/2} \equiv -1 \pmod n$, trying further primes
if necessary. Step 4 cannot fail for prime $n$; if it ever fails, an
exception is thrown, since a composite passing step 3 would be of
independent interest.
:::
:::

::: function
int n_is_probabprime(ulong n)

Obsolete; calls `n_is_prime`.
:::

## Factoring {#factoring}

::: summary algorithms/number-theory/factoring-word :::

::: type
n_factor_t

A factorisation $p_1^{e_1} \cdots p_k^{e_k}$ into distinct primes, with
fields `num` (the number $k$ of distinct primes), `p[]` (the primes) and
`exp[]` (the exponents).
:::

::: function
void n_factor_init(n_factor_t * factors)

Initialises `factors` to the empty factorisation.
:::

::: function
void n_factor(n_factor_t * factors, ulong n, int proved)

Sets `factors` to the complete factorisation of $n$ (empty for $n \le 1$);
`factors` need not be initialised. The flag `proved` has no effect, as all
primality tests used are deterministic for word-size integers.
:::

::: function
ulong n_factor_evaluate(const n_factor_t * factors)

Returns $p_1^{e_1} \cdots p_k^{e_k}$, or 0 if it does not fit in a `ulong`.
:::

::: function
void n_factor_insert(n_factor_t * factors, ulong p, ulong exp)

Multiplies the factorisation by $p^{\mathit{exp}}$, adding to the exponent
of $p$ if it is already present.
:::

::: function
ulong n_factor_trial(n_factor_t * factors, ulong n, ulong num_primes)
ulong n_factor_trial_range(n_factor_t * factors, ulong n, ulong start, ulong num_primes)

Trial division by the first `num_primes` primes (starting at the prime with
index `start`). Factors found are added to `factors`, which must be
initialised; the unfactored cofactor is returned.

::: details
The functions call `n_compute_primes` as needed and divide with
`n_remove2_precomp` using the cached `double` inverses. They stop as soon
as the square of the current prime exceeds the cofactor.
:::
:::

::: function
ulong n_factor_trial_partial(n_factor_t * factors, ulong n, ulong * prod, ulong num_primes, ulong limit)
ulong n_factor_partial(n_factor_t * factors, ulong n, ulong limit, int proved)

Partial factorisation: like `n_factor_trial` and `n_factor`, but stop as
soon as the product of the prime factors found exceeds `limit`. Return the
cofactor; `n_factor_trial_partial` also sets `*prod` to the product of the
factors found.
:::

::: function
ulong n_factor_one_line(ulong n, ulong iters)

Hart's one-line factoring algorithm [@Har2012], a variant of Fermat's
method. Returns a factor of $n$ or 0. Faster than SQUFOF below about
$2^{40}$.
:::

::: function
ulong n_factor_lehman(ulong n)

Lehman's algorithm [@Leh1974]. Always returns a factor of $n$, but works only up to
$10^{16}$ and is not used by `n_factor`.
:::

::: function
ulong n_factor_SQUFOF(ulong n, ulong iters)
ulong n_ll_factor_SQUFOF(ulong nhi, ulong nlo, ulong iters)

Shanks' square form factorisation as described by Gower and Wagstaff
[@GowWag2008]. Returns a factor of $n$ that fits in one limb, or 0 if none
was found within `iters` iterations. The `ll` version takes a two-limb
$n = n_{hi} B + n_{lo}$.

::: details
SQUFOF is first tried on $n$ itself and then on $kn$ for the small primes
$k$ in `flint_primes_small`. The product may need two limbs, which the
implementation supports; SQUFOF itself only works with numbers about the
size of $\sqrt n$.
:::
:::

::: function
ulong n_factor_pp1(ulong n, ulong B1, ulong c)
ulong n_factor_pp1_wrapper(ulong n)

Williams' $p + 1$ method [@Wil1982] with prime bound $B_1$ and a random starting value
$2 < c < n$. Returns a factor, or 0 or 1 on failure. The wrapper chooses
parameters tuned for 31 to 64-bit $n$.

::: details
Each choice of $c$ gives another chance of success, with roughly
exponentially decreasing probability. If neither $p + 1$ nor $p - 1$ is
smooth for any factor $p$ of $n$, the method cannot succeed.
:::
:::

::: function
int n_factor_pollard_brent_single(ulong * factor, ulong n, ulong ai, ulong xi, ulong max_iters)
int n_factor_pollard_brent(ulong * factor, flint_rand_t state, ulong n_in, ulong max_tries, ulong max_iters)

Pollard's rho method [@Pol1975] with Brent's cycle detection [@Bre1980],
for composite $n$.
Return 1 and set `*factor` to a nontrivial (not necessarily prime) factor
on success, otherwise 0. The first function uses the polynomial
$x^2 + a_i$ and starting value $x_i$; the second tries up to `max_tries`
random choices.
:::

::: function
int n_ll_factor_rho(nn_ptr factor, ulong nhi, ulong nlo, ulong max_tries, ulong max_iters)

Pollard–Brent rho for an odd two-limb $n$ with $n_{hi} \ne 0$. Returns 1
and writes a nontrivial two-limb factor on success, otherwise 0.

::: details
`max_iters` bounds Brent's outer doubling parameter, so each attempt
performs about $2\,\mathit{max\_iters}$ iterations and reaches factors up
to roughly $(2\,\mathit{max\_iters})^2$. The inner loop uses two-limb
Montgomery arithmetic and is several times faster than
`fmpz_factor_pollard_brent` at this size.
:::
:::

::: function
int n_ll_factor_ecm(nn_ptr factor, ulong nhi, ulong nlo, ulong curves, ulong B1, ulong B2, flint_rand_t state)

The elliptic curve method for an odd two-limb $n$ with $n_{hi} \ne 0$,
trying up to `curves` curves with bounds $B_1$, $B_2$ (if $B_2 < B_1$, then
$100 B_1$ is used). Returns 1 and writes a factor on success, otherwise 0.

::: details
Curves use Suyama's parametrisation. Both stages run in two-limb Montgomery
arithmetic, roughly twice as fast as `fmpz_factor_ecm` with the same
parameters for $n$ of this size.
:::
:::

::: function
int n_factor_ecm(ulong * f, ulong curves, ulong B1, ulong B2, flint_rand_t state, ulong n)

The elliptic curve method [@Len1987; @ZD2006] for a single-limb $n$,
trying up to `curves`
random curves with stage bounds $B_1$ and $B_2$. Returns 1 or 2 if a factor
was found in stage I or II, $-1$ if one was found while selecting a curve,
and 0 on failure. The factor is stored in `*f`.
:::

## Arithmetic functions {#arith}

::: function
int n_moebius_mu(ulong n)
void n_moebius_mu_vec(int * mu, ulong len)

Compute the Möbius function $\mu(n)$, with $\mu(0) = 0$, or the values
$\mu(0), \ldots, \mu(\mathit{len} - 1)$.

::: details
`n_moebius_mu` uses $\mu(4n) = 0$ and $\mu(2n) = -\mu(n)$ for even $n$,
looks up odd $n$ up to a cutoff in a table storing $\mu(n) + 1$ in two
bits, and otherwise checks for small odd square factors before calling
`n_factor`. `n_moebius_mu_vec` sieves: it flips the sign of $\mu$ at
every multiple of each prime $p$ and zeroes it at multiples of $p^2$.
:::
:::

::: function
int n_is_squarefree(ulong n)

Returns whether $n$ is squarefree, that is, $\mu(n) \ne 0$. By convention
1 is squarefree and 0 is not.
:::

::: function
ulong n_euler_phi(ulong n)

Returns Euler's totient $\varphi(n)$, the number of $1 \le k \le n$ coprime
to $n$.
:::

## Factorials {#factorials}

::: function
ulong n_factorial_mod2_preinv(ulong n, ulong p, ulong pinv)
ulong n_factorial_fast_mod2_preinv(ulong n, ulong p, ulong pinv)

Return $n! \bmod p$ for any modulus $p$ (no optimisations for composite
$p$). The `fast` version uses fast multipoint evaluation in about $O(n^{1/2})$
operations up to logarithmic factors [@Pol1974; @BGS2007]; the other chooses between table lookup, a direct product and the fast
algorithm depending on $n$.
:::

## Square roots modulo n {#sqrtmod}

::: function
ulong n_sqrtmod(ulong a, ulong p)

Returns the square root of $a$ modulo the prime $p$ that lies in
$[0, p/2]$, or 0 if $a$ is not a square. Requires $a < p$.

::: details
If $p$ is not prime, the result is either 0 or a genuine square root of $a$
modulo $p$. This is a wrapper of `nmod_sqrt`, which should be called
directly when an `nmod_t` for $p$ is available.
:::
:::

::: function
slong n_sqrtmod_2pow(ulong ** sqrt, ulong a, slong exp)
slong n_sqrtmod_primepow(ulong ** sqrt, ulong a, ulong p, slong exp)
slong n_sqrtmodn(ulong ** sqrt, ulong a, n_factor_t * fac)

Compute all square roots of $a$ modulo $2^{\mathit{exp}}$, modulo
$p^{\mathit{exp}}$, or modulo an integer $m$ given by its factorisation, and
return their number. The roots are written to a newly allocated array
stored in `*sqrt`, which the caller frees with `flint_free`. If $a$ is not a
square, 0 is returned and `*sqrt` is set to `NULL`.
:::

## Primitive roots and discrete logarithms {#roots-of-unity}

::: function
ulong n_quadratic_nonresidue(ulong n)

For odd non-square $n$, returns some $a$ with Jacobi symbol
$\left(\frac{a}{n}\right) = -1$; for prime $n$ this is a quadratic
nonresidue.
:::

::: function
ulong n_primitive_root_prime(ulong p)
ulong n_primitive_root_prime_prefactor(ulong p, n_factor_t * factors)

Return a primitive root modulo the prime $p$. The `prefactor` version takes
the prime factors $q$ of $p - 1$ (exponents are ignored) and returns an $a$
with $a^{(p-1)/q} \ne 1$ for every $q$ given.
:::

::: function
ulong n_discrete_log_bsgs(ulong b, ulong a, ulong n)

Returns $x$ with $a^x \equiv b \pmod n$ by baby-step giant-step
[@Sha1971], assuming
$(\Z/n\Z)^\times$ is cyclic, that is, $n$ is 2, 4, $p^k$ or $2p^k$ for an
odd prime $p$.
:::

## Random numbers {#random}

All functions draw from a `flint_rand_t` state. Those with `randtest` in
the name make corner cases such as 0, 1 and values with few set bits more
likely, and are meant for test code.

::: function
ulong n_randlimb(flint_rand_t state)
ulong _n_randlimb(flint_rand_t state)

Returns a uniformly random limb. `_n_randlimb` is an inline version for
tight loops.

::: details
`_n_randlimb` advances two linear congruential generators
$v \leftarrow (a v + c) \bmod B$ with fixed constants $a$, $c$, whose
states are stored in `state`. The high half of the first state gives the
low half of the result and the high half of the second state gives its
high half.
:::
:::

::: function
ulong n_randbits(flint_rand_t state, unsigned int bits)

Returns a uniformly random integer with exactly `bits` bits (the top bit
is set), or 0 if `bits` is 0.
:::

::: function
ulong n_randtest_bits(flint_rand_t state, int bits)

Like `n_randbits`, but values with sparse binary representations are more
likely.
:::

::: function
ulong n_randint(flint_rand_t state, ulong limit)
ulong n_urandint(flint_rand_t state, ulong limit)
ulong _n_randint(flint_rand_t state, ulong limit)

Returns a uniformly random integer in $[0, \mathit{limit})$, or a random
limb if `limit` is 0. The three functions are equivalent; `_n_randint` is
inline.
:::

::: function
ulong n_randtest(flint_rand_t state)

Returns a random integer of random bit length between 0 and `FLINT_BITS`,
with 0, 1, `COEFF_MAX`, `WORD_MAX`, `UWORD_MAX` and sparse values made
more likely.
:::

::: function
ulong n_randtest_not_zero(flint_rand_t state)

Like `n_randtest`, but never returns 0.
:::

::: function
ulong n_randprime(flint_rand_t state, ulong bits, int proved)

Returns a random prime with exactly `bits` bits, where
$2 \le \mathit{bits} \le \bits$. The result is always a proved prime;
the flag `proved` has no effect.
:::

::: function
ulong n_randtest_prime(flint_rand_t state, int proved)

Returns a random prime whose bit length is chosen at random between 2 and
`FLINT_BITS`.
:::

## Double-limb modular exponentiation {#double-limb}

Low-level helpers for the base-2 and small-base strong probable prime tests
in `n_ll_is_prime`. The modulus $m$ has two limbs ($m > B$) but must be
below roughly $2^{3 \bits / 2}$.

::: function
void n_ll_small_preinv(nn_ptr minv, nn_srcptr m)

Computes the two-limb inverse $\lfloor 2^{3 \bits} / m \rfloor$ of a
two-limb $m > B$, as needed by `n_ll_small_2_powmod` and
`n_ll_small_powmod_triple`.
:::

::: function
void n_ll_small_2_powmod(nn_ptr res, nn_srcptr exp, nn_srcptr m, nn_srcptr minv)
void n_ll_small_powmod_triple(nn_ptr res1, nn_ptr res2, nn_ptr res3, ulong b1, ulong b2, ulong b3, nn_srcptr exp, nn_srcptr m, nn_srcptr minv)

Compute $b^e \bmod m$ for small single-limb bases $b$ and two-limb $e$ and
$m > B$: the first function for $b = 2$, the second for three bases
$b_1 \le b_2 \le b_3$ at once. Helpers for `n_ll_is_prime`.

::: details
Both $m$ and $b$ must be small: these functions do not support the whole
two-limb range. Correctness is guaranteed if $(12 m)^2 < 2^{3 \bits}$ for the base-2
function and $(6 b_3 m)^2 < 2^{3 \bits}$ for the general one.
:::
:::

## ECM internals {#ecm}

Building blocks of `n_factor_ecm`. Points are in Montgomery projective
coordinates $(x : z)$ on curves in Montgomery form [@Mon1987], and arithmetic is
modulo $n$ using the parameters stored in `n_ecm_inf`.

::: function
void n_factor_ecm_double(ulong * x, ulong * z, ulong x0, ulong z0, ulong n, n_ecm_t n_ecm_inf)

Sets $(x : z) = 2 (x_0 : z_0)$, using
$x = (x_0 + z_0)^2 (x_0 - z_0)^2$ and
$z = 4 x_0 z_0 \left((x_0 - z_0)^2 + 4 a_{24} x_0 z_0\right)$.
:::

::: function
void n_factor_ecm_add(ulong * x, ulong * z, ulong x1, ulong z1, ulong x2, ulong z2, ulong x0, ulong z0, ulong n, n_ecm_t n_ecm_inf)

Sets $(x : z) = (x_1 : z_1) + (x_2 : z_2)$ given their difference
$(x_0 : z_0)$.
:::

::: function
void n_factor_ecm_mul_montgomery_ladder(ulong * x, ulong * z, ulong x0, ulong z0, ulong k, ulong n, n_ecm_t n_ecm_inf)

Sets $(x : z) = k (x_0 : z_0)$ with the Montgomery ladder.
:::

::: function
int n_factor_ecm_select_curve(ulong * f, ulong sigma, ulong n, n_ecm_t n_ecm_inf)

Selects a curve and starting point from `sigma` with Suyama's
parametrisation, storing $x_0$, $z_0$ and $(a + 2)/4$ in `n_ecm_inf`.
Returns 1 if a factor of $n$ turned up along the way (stored in `*f`),
otherwise 0.
:::

::: function
int n_factor_ecm_stage_I(ulong * f, const ulong * prime_array, ulong num, ulong B1, ulong n, n_ecm_t n_ecm_inf)
int n_factor_ecm_stage_II(ulong * f, ulong B1, ulong B2, ulong P, ulong n, n_ecm_t n_ecm_inf)

Stage I with the `num` primes up to $B_1$ in `prime_array`, and stage II
up to $B_2$ with primorial $P \approx \sqrt{B_2}$. Return 1 and set `*f`
if a factor is found, otherwise 0.
:::
