---
title: aprcl.h
header: aprcl.h
description: The APRCL primality test with Jacobi sums, and arithmetic in cyclotomic rings modulo n.
category: Integers
order: 9
---

This module implements the APRCL primality test of Adleman, Pomerance and
Rumely [@APR1983] in the Jacobi sum form of Cohen and Lenstra [@CL1984;
@CL1987], following Cohen's book [@Coh1996, §9.1]. It proves or disproves
the primality of integers of up to about 6000 digits, in time
$(\log n)^{O(\log \log \log n)}$, and its main step consists of many
independent tests that run in parallel when threads are available.
`fmpz_is_prime` uses it for its final proof below the size from which
[ECPP](ecpp) is faster (250 bits with one thread, more with several), and
as a fallback.

Most users should call `fmpz_is_prime`, which also handles the easy cases
by $n - 1$ and $n + 1$ tests. `aprcl_is_prime` is the entry point of this
module; the rest of it exposes the configuration of the test, the
individual checks, and the arithmetic in the rings
$\Z[\zeta_{p^k}]/(n)$ (`unity_zp`) and $\Z[\zeta_q, \zeta_p]/(n)$
(`unity_zpq`) on which the test is built. The module also contains the
older Gauss sum version of the test, which is slower and is kept for
reference. The module was written by Vladimir Glazachev (Google Summer of
Code 2015).

The test is described in [](algorithms/number-theory/primality-proving).

::: summary algorithms/number-theory/primality-proving :::

## Primality tests {#tests}

::: function
int aprcl_is_prime(const fmpz_t n)
int aprcl_is_prime_jacobi(const fmpz_t n)

Return 1 if $n$ is prime and 0 if it is not (in particular for $n < 2$),
using the Jacobi sum test with the parameters chosen by
`aprcl_config_jacobi_init`. `aprcl_is_prime` is the same function. If
the test cannot decide (which would require that no suitable auxiliary
prime exists for one of the additional tests, and has not been observed),
an exception is thrown; `_aprcl_is_prime_jacobi` returns `UNKNOWN`
instead. Above about 20000 bits (6000 digits) the configuration extends
the tabulated $R$ so that $s^2 > n$ still holds, but `ecpp_is_prime` is
then usually much faster.
:::

::: type
primality_test_status

The result of `_aprcl_is_prime_jacobi` and `_aprcl_is_prime_gauss`: one
of `PRIME`, `COMPOSITE`, `PROBABPRIME` (all checks passed but primality
could not be concluded) and `UNKNOWN` (the test could not decide).
:::

::: function
primality_test_status _aprcl_is_prime_jacobi(const fmpz_t n, const aprcl_config config)

The Jacobi sum test for $n$ with the given configuration. Returns
`PRIME`, `COMPOSITE` or `UNKNOWN`, the latter also when the configuration
does not satisfy $s^2 > n$.

::: details
The steps are those of Cohen's Algorithm 9.1.28: a check that $n$ is
coprime to $Rs$; for every prime $q \mid s$ and every prime power
$p^k \,\|\, q - 1$, the test that $J_0^u J_v$ is a root of unity in
$\Z[\zeta_{p^k}]/(n)$ (the pairs $(p, q)$ are independent tasks,
distributed over the available threads in decreasing order of cost); the
conditions $(\mathcal L_p)$, completed by additional tests with an extra
prime $q$ where needed; and the final division
(`aprcl_is_prime_final_division`). Moduli of 2 to 16 limbs use the
`unity_zp_mpn` arithmetic, larger ones the `unity_zp` arithmetic.
:::
:::

::: function
int aprcl_is_prime_gauss(const fmpz_t n)
int aprcl_is_prime_gauss_min_R(const fmpz_t n, ulong R)

Return 1 if $n$ is prime and 0 otherwise, using the Gauss sum test in the
form described by Schoof [@Sch2008], with the minimum value $R$
respectively 180. If the conditions $(\mathcal L_p)$ cannot be verified,
`aprcl_is_prime_gauss` retries with $R$ multiplied by 2, 3 and 5 before
throwing an exception. This version is much slower than the Jacobi sum
test.
:::

::: function
primality_test_status _aprcl_is_prime_gauss(const fmpz_t n, const aprcl_config config)

The Gauss sum test for $n$ with the given configuration. Returns `PRIME`,
`COMPOSITE` or `PROBABPRIME` (if primality could not be proved with this
$R$, or the configuration does not satisfy $s^2 > n$).
:::

::: function
int aprcl_is_prime_final_division(const fmpz_t n, const fmpz_t s, ulong r)

The final step of the test: given that every divisor of $n$ is congruent
to $n^i \bmod s$ for some $0 \le i < r$, and $\gcd(n, s) = 1$, returns 0
if $n$ is divisible by one of the residues $a = n^i \bmod s$ with
$1 < a \le \sqrt n$, $1 \le i < r$, and 1 otherwise. If $s^2 > n$, a
return value of 1 proves that $n$ is prime. The walk through the residues
is parallelised when threads are available.
:::

## Configuration {#config}

The test needs an integer $R$ and the product $s$ of prime powers $q^e$
with $q - 1 \mid R$, such that $a^R \equiv 1 \pmod s$ for all $a$ coprime
to $s$ and $s^2 > n$.

::: type
_aprcl_config
aprcl_config

The parameters of the test: `R`, `s`, the factorisation `rs` of $R$
(an `n_factor_t`), the factorisation `qs` of $s$, and flags `qs_used`
telling which prime powers of `qs` are used in the reduced $s$.
:::

::: function
ulong aprcl_R_value(const fmpz_t n)

Returns the value of $R$ from a precomputed table, chosen by the bit size
of $n$ so that the corresponding $s$ satisfies $s^2 > n$ and the test is
fast: from $R = 6$ for $n < 2^{17}$ to $R = 1396755360$ up to 12713 bits,
and $R = 6983776800$ (for which $s$ has about 10000 bits, enough up to
20004 bits) beyond. For larger $n$ this $s$ is too small, and
`aprcl_config_jacobi_init` extends $R$. On 32-bit machines inputs beyond
12713 bits raise an exception.
:::

::: function
void aprcl_config_jacobi_init(aprcl_config conf, const fmpz_t n)
void aprcl_config_jacobi_init_R(aprcl_config conf, const fmpz_t n, ulong R)
void aprcl_config_jacobi_clear(aprcl_config conf)

Initialise `conf` for the Jacobi sum test of $n$, with $R$ from
`aprcl_R_value` respectively the given $R$, and free it. If the value of
`aprcl_R_value` does not give $s^2 > n$, `aprcl_config_jacobi_init`
multiplies $R$ by the next primes after its largest prime factor until
it does (throwing an exception if $R$ would overflow a word). The full $s$
for $R$ is reduced by dropping the prime powers $q^e$ whose tests are the
most expensive per bit contributed to $s$ (estimated from the sizes of
$\varphi(p^k)$ for $p^k \| q - 1$, and very large $q$ first) as long as
$s^2 > n$ remains.
:::

::: function
void aprcl_config_gauss_init(aprcl_config conf, const fmpz_t n)
void aprcl_config_gauss_init_min_R(aprcl_config conf, const fmpz_t n, ulong R)
void aprcl_config_gauss_clear(aprcl_config conf)

Initialise `conf` for the Gauss sum test of $n$, choosing the least
suitable $R$ (respectively the least suitable $R$ that is a multiple of
the given one) with $s^2 > n$, and free it.
:::

::: function
int aprcl_is_mul_coprime_ui_ui(ulong x, ulong y, const fmpz_t n)
int aprcl_is_mul_coprime_ui_fmpz(ulong x, const fmpz_t y, const fmpz_t n)

Return 1 if $\gcd(xy, n) = 1$ and 0 otherwise.
:::

::: function
ulong aprcl_p_power_in_q(ulong q, ulong p)

Returns the exponent of the prime $p$ in $q > 0$.
:::

::: function
int _aprcl_p_ind(const aprcl_config conf, ulong p)

Returns the index of the prime $p$ in the factorisation of $R$ stored in
`conf`, or $-1$ if $p \nmid R$.
:::

## The individual checks {#checks}

Here $p^k \| q - 1$, $J(p, q)$ denotes a Jacobi sum (see
`unity_zp_jacobi_sum_pq`), and $u, v$ are the quantities of Cohen's
Algorithm 9.1.28, step 4: $n = p^k u + v$ with $0 \le v < p^k$. Each check computes $J_0^u J_v$ and returns $h$ if it
equals $\zeta_{p^k}^h$, and $-1$ if it is not a root of unity (which
proves $n$ composite).

::: function
slong _aprcl_is_prime_jacobi_check_pk(const unity_zp j, const fmpz_t u, ulong v)
slong _aprcl_is_prime_jacobi_check_21(ulong q, const fmpz_t n)
slong _aprcl_is_prime_jacobi_check_22(const unity_zp j, const fmpz_t u, ulong v, ulong q)
slong _aprcl_is_prime_jacobi_check_2k(const unity_zp j, const unity_zp j2_1, const unity_zp j2_2, const fmpz_t u, ulong v)

The check for odd $p$ (with $j = J(p, q)$), for $p^k = 2$ (which reduces
to $(-q)^{(n-1)/2} \equiv \pm 1 \pmod n$, returning 0 for $+1$, 1 for
$-1$ and $-1$ otherwise), for $p^k = 4$ (with $j = J(2, q)$; $h$ gives
$\zeta_4^h \in \{1, i, -1, -i\}$), and for $p = 2$, $k \ge 3$ (with
$j = J(2, q)$, $j_{2,1} = J_3(q)$ and $j_{2,2} = J_2(q)$).
:::

::: function
int _aprcl_is_prime_jacobi_check_pk_mpn(slong * h, const unity_zp j, const fmpz_t u, ulong v)
int _aprcl_is_prime_jacobi_check_22_mpn(slong * h, const unity_zp j, const fmpz_t u, ulong v, ulong q)
int _aprcl_is_prime_jacobi_check_2k_mpn(slong * h, const unity_zp j, const unity_zp j2_1, const unity_zp j2_2, const fmpz_t u, ulong v)

The same checks with `unity_zp_mpn` arithmetic. They return 0 without
computing anything if the modulus is not supported by `mpn_mod` (fewer
than 2 or more than 16 limbs), in which case the caller falls back to the
`unity_zp` versions; otherwise they return 1 and store the result in
`*h`.
:::

::: function
int _aprcl_is_prime_jacobi_additional_test(const fmpz_t n, ulong p)

The additional test of step 5 for a prime $p \mid R$ for which
$(\mathcal L_p)$ was not established: looks for a prime $q \equiv 1
\pmod{2p}$ (with $4 \| q - 1$ if $p = 2$) such that
$n^{(q-1)/p} \not\equiv 1 \pmod q$ and checks that the corresponding $h$ is
a primitive root of unity (and for $p = 2$ that
$q^{(n-1)/2} \equiv -1 \pmod n$). Returns 2 if $n$ is proved composite, 1
if the condition holds, and 0 if no suitable $q$ was found.
:::

::: function
int _aprcl_is_gausspower_2q_equal_first(ulong q, const fmpz_t n)
int _aprcl_is_gausspower_2q_equal_second(ulong q, const fmpz_t n)
slong _aprcl_is_gausspower_from_unity_p(ulong q, ulong r, const fmpz_t n)

Checks of the Gauss sum test, for the Gauss sum $\tau(\chi)$ of a
character $\chi$ of conductor $q$: whether
$(\chi(-1) q)^{(n-1)/2} \equiv -1 \pmod n$ and whether
$q^{(n-1)/2} \equiv -1 \pmod n$ (for the quadratic character, returning
1 or 0), and for a character of order $r$, returns $i$ such that
$\tau(\chi^n) = \zeta_r^i \tau(\chi)^n$, or $-1$ if there is none.
:::

## Arithmetic in $\Z[\zeta_{p^k}]/(n)$ {#unity-zp}

::: type
_unity_zp
unity_zp

An element of $\Z[\zeta_{p^k}]/(n)$, represented by an `fmpz_mod_poly_t`
(field `poly`) together with $p$ (`p`), $k$ (`exp`) and the modulus
context `ctx`. The polynomial is reduced modulo $x^{p^k} - 1$, and modulo
the cyclotomic polynomial $\Phi_{p^k}$ after `unity_zp_reduce_cyclotomic`.
Unless stated otherwise, the operands of a function must have been
initialised with the same $p$, $k$ and $n$.
:::

::: macro
SQUARING_SPACE
UNITY_ZP_MUL_BEGIN(t, glen, hlen)
UNITY_ZP_MUL_END(t)

The size (70) of the temporary arrays of `fmpz` used by the in-place
multiplications, and macros that set up and release such an array of
$\mathit{glen} + \mathit{hlen} - 1$ entries, on the stack if it fits.
:::

::: function
void unity_zp_init(unity_zp f, ulong p, ulong exp, const fmpz_t n)
void unity_zp_clear(unity_zp f)

Initialise $f$ as the zero element of $\Z[\zeta_{p^{exp}}]/(n)$, and free
it.
:::

::: function
void unity_zp_copy(unity_zp f, const unity_zp g)
void unity_zp_swap(unity_zp f, unity_zp g)
void unity_zp_set_zero(unity_zp f)

Copy, swap and set to zero.
:::

::: function
slong unity_zp_is_unity(unity_zp f)
int unity_zp_equal(unity_zp f, unity_zp g)

`unity_zp_is_unity` returns $h$ if $f = \zeta^h$ and $-1$ otherwise;
`unity_zp_equal` returns whether $f = g$ after reduction by
$\Phi_{p^{exp}}$.
:::

::: function
void unity_zp_coeff_set_fmpz(unity_zp f, ulong ind, const fmpz_t x)
void unity_zp_coeff_set_ui(unity_zp f, ulong ind, ulong x)
void unity_zp_coeff_add_fmpz(unity_zp f, ulong ind, const fmpz_t x)
void unity_zp_coeff_add_ui(unity_zp f, ulong ind, ulong x)
void unity_zp_coeff_inc(unity_zp f, ulong ind)
void unity_zp_coeff_dec(unity_zp f, ulong ind)

Set, add to, increment or decrement the coefficient of $\zeta^{ind}$, for
$\mathit{ind} < p^{exp}$. The added $x$ must be less than $n$.
:::

::: function
void unity_zp_mul_scalar_ui(unity_zp f, const unity_zp g, ulong s)
void unity_zp_add(unity_zp f, const unity_zp g, const unity_zp h)

Set $f$ to $sg$, respectively $g + h$.
:::

::: function
void unity_zp_mul(unity_zp f, const unity_zp g, const unity_zp h)
void unity_zp_sqr(unity_zp f, const unity_zp g)
void unity_zp_mul_inplace(unity_zp f, const unity_zp g, const unity_zp h, fmpz_t * t)
void unity_zp_sqr_inplace(unity_zp f, const unity_zp g, fmpz_t * t)

Set $f$ to $gh$, respectively $g^2$, reduced by $\Phi_{p^{exp}}$ and $n$.
The `inplace` versions use the preallocated array $t$ of
`SQUARING_SPACE` values as temporary space when it is large enough.
:::

::: function
void _unity_zp_mul_reduce(unity_zp f, const unity_zp g, const unity_zp h, fmpz * t)

Sets $f$ to $gh$ (or $g^2$ if $h$ is $g$) reduced by $\Phi_{p^{exp}}$ and
modulo $n$: the product over the integers is computed by
`_fmpz_poly_mul` or `_fmpz_poly_sqr` into $t$ (which must hold at least
`g->poly->length + h->poly->length - 1` initialised values), the
cyclotomic reduction is done there with additions and subtractions (first
by $x^{p^{exp}} = 1$, then by $\Phi_{p^{exp}}(x) = \Phi_p(x^{p^{exp-1}})$),
and each of the $\varphi(p^{exp})$ remaining coefficients is reduced
modulo $n$ once. The operands need not be reduced by $\Phi_{p^{exp}}$.
:::

::: function
int _unity_zp_mul_special(unity_zp f, const unity_zp g, const unity_zp h, fmpz * t)
int _unity_zp_sqr_special(unity_zp f, const unity_zp g, fmpz * t)

If a hand-optimised straight-line program exists for $p^{exp}$
(currently 3, 4, 5, 7, 8, 9 and 16), the operands have full length
$\varphi(p^{exp})$ and the modulus is small enough (a few thousand bits)
that the smaller number of multiplications pays off, compute the reduced
product (or square) using $t$ (at least `SQUARING_SPACE` values) as
temporary space and return 1; otherwise return 0.
:::

::: function
void unity_zp_pow_fmpz(unity_zp f, const unity_zp g, const fmpz_t pow)
void unity_zp_pow_ui(unity_zp f, const unity_zp g, ulong pow)
void unity_zp_pow_2k_fmpz(unity_zp f, const unity_zp g, const fmpz_t pow)
void unity_zp_pow_2k_ui(unity_zp f, const unity_zp g, ulong pow)
void unity_zp_pow_sliding_fmpz(unity_zp f, unity_zp g, const fmpz_t pow)

Set $f$ to $g^{\mathit{pow}}$, by binary exponentiation, by the $2^k$-ary
method with $k$ from `_unity_zp_pow_select_k`, and by sliding window
exponentiation.
:::

::: function
ulong _unity_zp_pow_select_k(const fmpz_t n)

Returns the least $k$ with
$\log n < k(k + 1) 2^{2k} / (2^{k + 1} - k - 2) + 1$, the window size for
an exponent of the size of $n$.
:::

::: function
void unity_zp_reduce_cyclotomic(unity_zp f, const unity_zp g)
void _unity_zp_reduce_cyclotomic(unity_zp f)
void _unity_zp_reduce_cyclotomic_divmod(unity_zp f)

Set $f = g \bmod \Phi_{p^{exp}}$, respectively reduce $f$ in place (by
additions, or by polynomial division); the in-place versions require $f$
to be reduced modulo $x^{p^{exp}} - 1$.
:::

::: function
void unity_zp_aut(unity_zp f, const unity_zp g, ulong x)
void unity_zp_aut_inv(unity_zp f, const unity_zp g, ulong x)

Set $f = \sigma_x(g)$ for the automorphism $\sigma_x(\zeta) = \zeta^x$,
respectively $f = \sigma_x^{-1}(g)$; for the latter $g$ must be reduced by
$\Phi_{p^{exp}}$.
:::

### Jacobi sums {#jacobi-sums}

Here $\chi_{p,q}$ is the character modulo the prime $q$ defined by
$\chi_{p,q}(g^x) = \zeta_{p^k}^x$, where $g$ is a primitive root modulo
$q$ and $p^k \| q - 1$, and $j(\chi_1, \chi_2) = \sum_x \chi_1(x)
\chi_2(1 - x)$.

::: function
void unity_zp_jacobi_sum_pq(unity_zp f, ulong q, ulong p)

Sets $f$ to $J(p, q) = j(\chi_{p,q}, \chi_{p,q})$.
:::

::: function
void unity_zp_jacobi_sum_2q_one(unity_zp f, ulong q)
void unity_zp_jacobi_sum_2q_two(unity_zp f, ulong q)

For $p = 2$, set $f$ to $J_2(q) = j(\chi_{2,q}^{2^{k-3}},
\chi_{2,q}^{3 \cdot 2^{k-3}})^2$, respectively $J_3(q) =
j(\chi_{2,q}, \chi_{2,q}, \chi_{2,q}) = J(2, q)\, j(\chi_{2,q}^2,
\chi_{2,q})$.
:::

::: function
nn_ptr aprcl_f_table(const ulong q)

Returns a table (allocated with `_nmod_vec_init`) of the discrete
logarithms $f(x)$ with $g^{f(x)} = 1 - g^x \pmod q$ for a primitive root
$g$ modulo the prime $q$, from which the Jacobi sums are computed.
:::

::: function
void _unity_zp_jacobi_sum_pq_general(unity_zp f, const nn_ptr table, ulong p, ulong q, ulong k, ulong a, ulong b)

Sets $f$ to $\sum_{x=1}^{q-2} \zeta_{p^k}^{ax + b f(x)}$, given the table
of `aprcl_f_table`.
:::

### Fixed-size coefficients {#unity-zp-mpn}

::: type
_unity_zp_mpn
unity_zp_mpn

A dense analogue of `unity_zp` with coefficients in an `mpn_mod` ring,
used by the Jacobi sum checks for moduli of 2 to 16 limbs. An element is
a vector of $d = \varphi(p^{exp})$ residues (field `coeffs`), always
reduced by $\Phi_{p^{exp}}$, together with $p$, $k$, $p^k$, $p^{k-1}$,
$d$, the borrowed `mpn_mod` context `ctx`, and data for the products,
which are computed with a single modular reduction per output coefficient
on top of `_flint_mpn_poly_mulmid` (classical or Karatsuba).
:::

::: function
void unity_zp_mpn_init(unity_zp_mpn f, ulong p, ulong exp, gr_ctx_t ctx)
void unity_zp_mpn_clear(unity_zp_mpn f)

Initialise $f$ as the zero element of $\Z[\zeta_{p^{exp}}]/(n)$, where
`ctx` is an `mpn_mod` context with modulus $n$ (borrowed: it must outlive
$f$), and free it.
:::

::: function
void unity_zp_mpn_set_zero(unity_zp_mpn f)
void unity_zp_mpn_swap(unity_zp_mpn f, unity_zp_mpn g)
void unity_zp_mpn_copy(unity_zp_mpn f, const unity_zp_mpn g)
void unity_zp_mpn_coeff_set_ui(unity_zp_mpn f, ulong ind, ulong x)
void unity_zp_mpn_set_unity_zp(unity_zp_mpn f, const unity_zp g)

Set to zero, swap, copy, set the coefficient of $\zeta^{ind}$ to $x$ (for
any $\mathit{ind} < p^{exp}$, exponents of at least $\varphi(p^{exp})$
being folded by the cyclotomic relation), and convert from a `unity_zp`
with the same $p$, $k$ and $n$.
:::

::: function
void unity_zp_mpn_mul(unity_zp_mpn f, const unity_zp_mpn g, const unity_zp_mpn h)
void unity_zp_mpn_sqr(unity_zp_mpn f, const unity_zp_mpn g)
void unity_zp_mpn_mul_scalar_ui(unity_zp_mpn f, const unity_zp_mpn g, ulong s)
void unity_zp_mpn_pow_ui(unity_zp_mpn f, const unity_zp_mpn g, ulong pow)
void unity_zp_mpn_pow_sliding_fmpz(unity_zp_mpn f, const unity_zp_mpn g, const fmpz_t pow)

Product, square, scalar multiple and powers (the last by sliding window
exponentiation). The output may alias the inputs.
:::

::: function
void unity_zp_mpn_aut_inv(unity_zp_mpn f, const unity_zp_mpn g, ulong x)
slong unity_zp_mpn_is_unity(const unity_zp_mpn f)

$f = \sigma_x^{-1}(g)$ ($f$ and $g$ must not be aliased), and the
exponent $h$ with $f = \zeta^h$, or $-1$.
:::

## Arithmetic in $\Z[\zeta_q, \zeta_p]/(n)$ {#unity-zpq}

These are used by the Gauss sum test only.

::: type
_unity_zpq
unity_zpq

An element of $\Z[\zeta_q, \zeta_p]/(n)$, represented by $p$
polynomials in $\zeta_q$ (field `polys`, of type `fmpz_mod_poly_t`), the
$j$-th one holding the coefficient of $\zeta_p^j$.
:::

::: function
void unity_zpq_init(unity_zpq f, ulong q, ulong p, const fmpz_t n)
void unity_zpq_clear(unity_zpq f)
void unity_zpq_copy(unity_zpq f, const unity_zpq g)
void unity_zpq_swap(unity_zpq f, unity_zpq g)
int unity_zpq_equal(const unity_zpq f, const unity_zpq g)

Initialise as zero, free, copy, swap and compare.
:::

::: function
void unity_zpq_coeff_set_fmpz(unity_zpq f, slong i, slong j, const fmpz_t x)
void unity_zpq_coeff_set_ui(unity_zpq f, slong i, slong j, ulong x)
void unity_zpq_coeff_add(unity_zpq f, slong i, slong j, const fmpz_t x)
void unity_zpq_coeff_add_ui(unity_zpq f, slong i, slong j, ulong x)

Set or add to the coefficient of $\zeta_q^i \zeta_p^j$, for $i < q$ and
$j < p$; the added $x$ must be less than $n$.
:::

::: function
void unity_zpq_add(unity_zpq f, const unity_zpq g, const unity_zpq h)
void unity_zpq_mul(unity_zpq f, const unity_zpq g, const unity_zpq h)
void _unity_zpq_mul_unity_p(unity_zpq f)
void unity_zpq_mul_unity_p_pow(unity_zpq f, const unity_zpq g, slong p)
void unity_zpq_pow(unity_zpq f, const unity_zpq g, const fmpz_t p)
void unity_zpq_pow_ui(unity_zpq f, const unity_zpq g, ulong pow)

Sum, product, multiplication of $f$ by $\zeta_p$ in place, $f = g
\zeta_p^k$ (the last argument being $k$), and powers.
:::

::: function
void unity_zpq_gauss_sum(unity_zpq f, ulong q, ulong p)
void unity_zpq_gauss_sum_character_pow(unity_zpq f, ulong q, ulong p, ulong pow)
void unity_zpq_gauss_sum_sigma_pow(unity_zpq f, ulong q, ulong p)

Set $f$ to the Gauss sum $\tau(\chi_{p,q}) = \sum_x \chi_{p,q}(x)
\zeta_q^x$, to $\tau(\chi_{p,q}^{\mathit{pow}})$, and to
$\sigma_n(\tau(\chi_{p,q}))$.
:::
