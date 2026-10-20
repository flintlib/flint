---
title: fmpz_factor.h
header: fmpz_factor.h
description: Factorisation of integers into proven primes, and the factoring methods behind it: trial division, Pollard rho, $p + 1$ and the elliptic curve method.
category: Integers
order: 4
---

This module computes the prime factorisation of an `fmpz` integer and
provides the factored-form type `fmpz_factor_t` together with the
individual factoring methods for multiprecision integers. `fmpz_factor`
combines trial division, perfect power detection, the elliptic curve
method, Pollard–Brent rho, SQUFOF and the [quadratic sieve](qsieve), and
proves the primality of every factor it returns. It is practical for
integers of up to about 90 digits, or larger ones whose second largest
prime factor has at most about 40 digits. `fmpz_factor_smooth` finds the
small factors only.

[](fmpz_factor/factor) explains the strategy of `fmpz_factor`
with its cutoffs; the topic page
[](algorithms/number-theory/integer-factoring) compares the methods.

::: summary algorithms/number-theory/integer-factoring :::

## The factored-form type {#type}

::: type
fmpz_factor_struct
fmpz_factor_t

An integer in factored form $n = s \prod_{i} p_i^{e_i}$. The fields are
`sign` ($-1$, 0 or 1), `p` (an array of `fmpz` bases), `exp` (an array of
`ulong` exponents), `num` (the number of factors) and `alloc` (the
allocated length of the arrays). Zero is represented by sign 0 and no
factors, $\pm 1$ by sign $\pm 1$ and no factors.

The factorisation functions produce distinct primes with positive
exponents, but *not* sorted in any particular order: small primes found by
trial division come first in increasing order, followed by the larger
factors in the order in which they were split off.
:::

::: function
void fmpz_factor_init(fmpz_factor_t factor)
void fmpz_factor_clear(fmpz_factor_t factor)

Initialise `factor` (representing 1) and free the memory it uses.
:::

::: function
void fmpz_factor_get_fmpz(fmpz_t z, const fmpz_factor_t factor, slong i)

Sets $z$ to the $i$-th base `factor->p + i`, for $0 \le i <$ `factor->num`.
:::

::: function
void fmpz_factor_expand(fmpz_t n, const fmpz_factor_t factor)
void fmpz_factor_expand_multiexp(fmpz_t n, const fmpz_factor_t factor)
void fmpz_factor_expand_iterative(fmpz_t n, const fmpz_factor_t factor)

Set $n$ to the integer represented by `factor`, including the sign.
`fmpz_factor_expand` calls `fmpz_factor_expand_multiexp`, which processes
the bits of all exponents simultaneously, from the top: at each step it
squares the result and multiplies it by the product of the bases whose
exponent has that bit set. This is much faster than
`fmpz_factor_expand_iterative`,
which raises each base to its power and multiplies the powers one at a
time.
:::

## Input and output {#io}

::: function
int fmpz_factor_print(const fmpz_factor_t factor)
int fmpz_factor_fprint(FILE * fs, const fmpz_factor_t factor)

Print `factor` to `stdout` or to `fs` as `p1^e1 * p2^e2 * ...` (omitting
exponents equal to 1), preceded by `-1 * ` if the sign is negative.
Zero is printed as `0`, and $1$ and $-1$ as `1` and `-1`. Return the
number of characters written.
:::

## Factoring {#factoring}

::: function
void fmpz_factor(fmpz_factor_t factor, const fmpz_t n)

Sets `factor` to the factorisation of $n$ into primes. The sign of $n$
goes into the sign field; for $n = 0$ the sign is 0 and there are no
factors. Every base is a prime whose primality has been proved.
:::

::: function
void fmpz_factor_si(fmpz_factor_t factor, slong n)

`fmpz_factor` for a machine integer $n$ (using `n_factor`).
:::

::: function
void fmpz_factor_no_trial(fmpz_factor_t factor, const fmpz_t n)

Appends the factorisation of $n > 1$ to `factor`, without the initial
trial division of `fmpz_factor`: it is meant for integers already known
to have no small prime factors. It does not reset `factor` or set its sign.
:::

::: function
int fmpz_factor_smooth(fmpz_factor_t factor, const fmpz_t n, slong bits, int proved)

Factors $n$ into primes of up to roughly `bits` bits and possibly one
remaining cofactor, which may or may not be prime. Returns 1 if the
factorisation is known to be complete; otherwise returns 0, and the last
entry (whose exponent may exceed 1) is a composite cofactor that needs
further work. Requires `bits` $\ge 1$.

::: details
The function trial-divides by the primes up to $2^{\mathit{bits}}$ (by
the first 3512 primes, up to 32749, if `bits` is at least 15), then tests
the cofactor for primality and for being a perfect power, and if
`bits` $\ge 16$ runs the elliptic curve method with parameters tuned to
find a factor of the given size with probability at least about $2/3$
(factors of about 18 to 22 bits first, then increasing by 6 bits at a
time, up to 100 bits). Smaller factors are found with much higher
probability. If a cofactor of at most one word remains (or $n$ fits in
one word), it is factored completely with `n_factor`.

The parameter `proved` controls the primality tests: with 1, all factors
are proved prime (except the last one if the return value is 0); with 0,
they are only checked to be probable primes; with $-1$, no primality test
is done after trial division (but the perfect power test is still done).

Choosing `bits` larger than about a third of the size of $n$ is usually
slower than finishing with the quadratic sieve; `fmpz_factor` uses
$\max(b/3 - 17, 2)$ for a $b$-bit number.
:::
:::

::: function
int _fmpz_factor_smooth(fmpz_factor_t factor, const fmpz_t n, slong bits, int proved, int adaptive)

`fmpz_factor_smooth` (which calls this with `adaptive` = 0). If
`adaptive` is nonzero, the bound `bits` is lowered whenever the
unfactored part gets smaller (after trial division and after each factor
found by the elliptic curve method) to the bound that
`fmpz_factor_no_trial` would choose for a number of that size, so that
the search stops once finishing the cofactor with the quadratic sieve is
expected to be faster. This is used by `fmpz_factor_no_trial`.
:::

::: function
int fmpz_factor_trial_range(fmpz_factor_t factor, const fmpz_t n, ulong start, ulong num_primes)

Sets `factor` to the sign of $n$ and the prime factors of $n$ found by
trial division: the power of 2 and the primes `flint_primes[i]` with
`start` $\le i <$ `start + num_primes`. Returns 1 if the remaining
cofactor is 1 and 0 otherwise; the cofactor itself is *not* stored. If
$n$ fits in a small `fmpz` (at most $\bits - 2$ bits), it is instead
factored completely, and 1 is returned.
:::

::: function
int fmpz_factor_trial(fmpz_factor_t factor, const fmpz_t n, slong num_primes)

Sets `factor` to the sign of $n$ and the prime factors of $n$ among the
first `num_primes` primes, found with the product tree method of
`flint_mpn_factor_trial_tree`. If the cofactor is not 1, it is appended
as the last entry (with exponent 1, and not necessarily prime) and 0 is
returned; otherwise 1 is returned. Requires $0 \le$ `num_primes`
$\le 3512$ (an exception is raised otherwise). If $n$ fits in a small
`fmpz`, it is instead factored completely.
:::

::: function
void fmpz_factor_refine(fmpz_factor_t res, const fmpz_factor_t f)

Given a partial factorisation `f` of an integer (bases need not be prime,
coprime or in any order, and may be 1), sets `res` to a factorisation of
the same integer whose bases are pairwise coprime and greater than 1, by
the factor refinement algorithm of Bach, Driscoll and Shallit
[@BDS1993], which uses repeated GCDs. The result need not consist of
primes and is not sorted.
:::

## Factoring methods {#methods}

These functions look for one factor of a composite number. They are
probabilistic: a return value of 0 means only that no factor was found
with the given effort.

::: function
int fmpz_factor_pollard_brent(fmpz_t factor, flint_rand_t state, fmpz_t n, ulong max_tries, ulong max_iters)

Looks for a factor of the composite odd $n$ by Pollard's rho method with
Brent's cycle detection [@Pol1975; @Bre1980], using up to `max_tries`
random polynomials $x^2 + a$ and starting points, each for about
`max_iters` iterations. If a nontrivial factor is found, sets `factor` to
it and returns its number of limbs; otherwise returns 0. The factor need
not be prime. For even $n$ it returns 1 and the factor 2. Single-word $n$
are passed to `n_factor_pollard_brent`. The expected number of
iterations to find a prime factor $p$ is $O(\sqrt p)$.
:::

::: function
int fmpz_factor_pollard_brent_single(fmpz_t p_factor, fmpz_t n_in, fmpz_t yi, fmpz_t ai, ulong max_iters)

One attempt of `fmpz_factor_pollard_brent` with the polynomial
$x^2 + a$ and the starting point $y$ given by `ai` and `yi` (with
$1 \le a \le n - 3$ and $1 \le y < n$). Returns the number of limbs of
the factor found, or 0.
:::

::: function
int fmpz_factor_pp1(fmpz_t factor, const fmpz_t n, ulong B1, ulong B2_sqrt, ulong c)

Williams' $p + 1$ method [@Wil1982]: looks for a factor of $n$ with a
stage 1 bound `B1` and a stage 2 bound of at least $\mathit{B2\_sqrt}^2$,
using the Lucas sequence with starting value $c > 2$. If a factor is
found, sets `factor` to it and returns 1; otherwise returns 0. For even
$n$ it returns the factor 2.

::: details
The method finds a prime factor $p$ when $p - \left(\frac{c^2 - 4}{p}\right)$
is `B1`-smooth apart from at most one prime below the stage 2 bound. Since
the Legendre symbol depends on $c$, which is unknown in advance, the method
succeeds either for $p + 1$ or, as a $p - 1$ method, for $p - 1$; repeated
calls with different random $c$ give additional chances with roughly
exponentially decaying probability of finding a missed factor. If neither
$p + 1$ nor $p - 1$ is smooth for any prime factor $p$ of $n$, the method
never succeeds.
:::
:::

::: function
int fmpz_factor_ecm(fmpz_t f, ulong curves, ulong B1, ulong B2, flint_rand_t state, const fmpz_t n_in)

Looks for a factor of $n$ = `n_in` by Lenstra's elliptic curve method
[@Len1987] with Montgomery's curves and arithmetic [@Mon1987], trying up
to `curves` random curves with stage 1 bound `B1` and stage 2 bound `B2`
($B_1 < B_2$). Returns 1 if a factor is found in stage 1, 2 if it is found
in stage 2, $-1$ if it is found while setting up a curve, and 0
otherwise; the factor (not necessarily prime) is written to `f`.
Single-word $n$ are passed to `n_factor_ecm`.

::: details
The curves are distributed over the available threads (the result is the
factor found by the curve of lowest index, so it does not depend on the
number of threads). See [](algorithms/number-theory/ecm) for the stages
and the parameter choice.
:::
:::

## Multiplicative functions {#multiplicative}

These evaluate arithmetic functions from a factorisation; `fmpz_euler_phi`,
`fmpz_moebius_mu` and `fmpz_divisor_sigma` in the `fmpz` module factor
their argument and call them.

::: function
void fmpz_factor_euler_phi(fmpz_t res, const fmpz_factor_t fac)

Sets `res` to Euler's totient $\varphi(n) = \prod p_i^{e_i - 1}(p_i - 1)$
of the integer whose factorisation into distinct primes is `fac` (the
sign is ignored).
:::

::: function
int fmpz_factor_moebius_mu(const fmpz_factor_t fac)

Returns the Möbius function $\mu(n)$: 0 if some exponent exceeds 1, and
$(-1)^k$ for $k$ distinct primes otherwise.
:::

::: function
void fmpz_factor_divisor_sigma(fmpz_t res, ulong k, const fmpz_factor_t fac)

Sets `res` to $\sigma_k(n) = \sum_{d \mid n} d^k = \prod_i
(p_i^{k(e_i + 1)} - 1)/(p_i^k - 1)$ (the number of divisors for $k = 0$).
:::

## Elliptic curve method: components {#ecm}

The stages of `fmpz_factor_ecm`, exposed for experimentation. They work
with limb arrays of `n_size` limbs in the *shifted* representation of
`flint_mpn_mulmod_preinvn`: the modulus $n$ is stored as the normalised
$n 2^{\mathrm{norm}}$ (`ecm_inf->normbits` is the shift) with its inverse
`ecm_inf->ninv`, and all residues are multiplied by $2^{\mathrm{norm}}$.
The argument `n` of these functions is the shifted modulus. The curves
are in Montgomery form $B y^2 = x^3 + A x^2 + x$, and points are
represented by projective $x$-coordinates $(X : Z)$.

::: type
ecm_s
ecm_t

Working data for the elliptic curve method: temporary arrays (`t`, `u`,
`v`, `w`), the current point (`x`, `z`), $a_{24} = (A + 2)/4$ (`a24`),
the inverse of the modulus (`ninv`), the shifted value of 1 (`one`), the
stage 2 tables (`GCD_table`, `prime_table`), and `n_size` and `normbits`.
:::

::: function
void fmpz_factor_ecm_init(ecm_t ecm_inf, ulong sz)
void fmpz_factor_ecm_clear(ecm_t ecm_inf)

Initialise `ecm_inf` for moduli of `sz` limbs, and free it. The stage 2
tables are not allocated by `fmpz_factor_ecm_init`.
:::

::: function
void fmpz_factor_ecm_double(nn_ptr x, nn_ptr z, nn_ptr x0, nn_ptr z0, nn_ptr n, ecm_t ecm_inf)

Sets $(x : z)$ to $2 (x_0 : z_0)$ by Montgomery's doubling formula
$x = (x_0 + z_0)^2 (x_0 - z_0)^2$,
$z = 4 x_0 z_0 \bigl((x_0 - z_0)^2 + a_{24} \cdot 4 x_0 z_0\bigr)$ modulo $n$.
:::

::: function
void fmpz_factor_ecm_add(nn_ptr x, nn_ptr z, nn_ptr x1, nn_ptr z1, nn_ptr x2, nn_ptr z2, nn_ptr x0, nn_ptr z0, nn_ptr n, ecm_t ecm_inf)

Sets $(x : z)$ to $(x_1 : z_1) + (x_2 : z_2)$ given their difference
$(x_0 : z_0)$, by Montgomery's differential addition
$x = z_0 (x_1 x_2 - z_1 z_2)^2$, $z = x_0 (x_2 z_1 - x_1 z_2)^2$ modulo
$n$ (up to a common factor 4).
:::

::: function
void fmpz_factor_ecm_mul_montgomery_ladder(nn_ptr x, nn_ptr z, nn_ptr x0, nn_ptr z0, ulong k, nn_ptr n, ecm_t ecm_inf)

Sets $(x : z)$ to $k (x_0 : z_0)$ by the Montgomery ladder (one doubling
and one differential addition per bit of $k$).
:::

::: function
int fmpz_factor_ecm_select_curve(nn_ptr f, nn_ptr sig, nn_ptr n, ecm_t ecm_inf)

Sets up a curve and a starting point by Suyama's parametrisation from the
parameter $\sigma$ = `sig`: with $u = \sigma^2 - 5$ and $v = 4\sigma$,
the point is $(u^3 : v^3)$ and $a_{24} = (v - u)^3 (3u + v) / (16 u^3 v)$;
the starting point is stored in `ecm_inf->x`, `ecm_inf->z` and $a_{24}$ in
`ecm_inf->a24`. Returns 0 on success; if the inversion needed for
$a_{24}$ reveals a factor of $n$, writes it to `f` and returns its number
of limbs; returns $-1$ if the curve is degenerate.
:::

::: function
int fmpz_factor_ecm_stage_I(nn_ptr f, const ulong * prime_array, ulong num, ulong B1, nn_ptr n, ecm_t ecm_inf)

Stage 1: multiplies the point by $p^{\lfloor \log_p B_1 \rfloor}$ for each
of the `num` primes $p \le B_1$ in `prime_array`, taking the GCD of
$z$ with $n$ after each prime. If a nontrivial factor is found, writes it
to `f` and returns its number of limbs; otherwise returns 0.
:::

::: function
int fmpz_factor_ecm_stage_II(nn_ptr f, ulong B1, ulong B2, ulong P, nn_ptr n, ecm_t ecm_inf)

Stage 2 (the standard continuation, by baby steps and giant steps): with
the primorial $P \approx \sqrt{B_2}$, accumulates the product of
$X_{mP} Z_j - X_j Z_{mP}$ over all $m$ and $j$ such that $mP \pm j$ is a
prime in $(B_1, B_2]$ (from the tables in `ecm_inf`), and takes one GCD
with $n$ at the end. Returns the number of limbs of the factor found, or
0.
:::

## Internal functions {#internal}

::: function
void _fmpz_factor_fit_length(fmpz_factor_t factor, slong len)
void _fmpz_factor_set_length(fmpz_factor_t factor, slong newlen)

Make room for `len` factors, and set the number of factors (freeing the
bases beyond the new length).
:::

::: function
void _fmpz_factor_append(fmpz_factor_t factor, const fmpz_t p, ulong exp)
void _fmpz_factor_append_ui(fmpz_factor_t factor, ulong p, ulong exp)

Append the factor $p^{\mathit{exp}}$. No check is made that $p$ is not
already present.
:::

::: function
void _fmpz_factor_concat(fmpz_factor_t factor1, fmpz_factor_t factor2, ulong exp)

Appends the factors of `factor2` to `factor1`, with their exponents
multiplied by `exp`. The bases must be distinct from those of `factor1`;
they are not merged.
:::

::: function
void _fmpz_factor_extend_factor_ui(fmpz_factor_t factor, ulong n)

Appends the factorisation of the word $n$ (computed by `n_factor`); for
$n = 0$ it instead clears `factor` and sets its sign to 0.
:::

::: function
void flint_mpn_sqr_and_add_a(nn_ptr y, nn_ptr a, nn_ptr n, ulong n_size, nn_ptr ninv, ulong normbits)

Sets $y \gets y^2 + a$ modulo the shifted modulus $n$ with inverse
`ninv`, the iteration of Pollard's rho method. Requires $y, a < n$; the
result is reduced.
:::

::: function
int flint_mpn_factor_pollard_brent_single(nn_ptr factor, nn_ptr n, nn_ptr ninv, nn_ptr a, nn_ptr y, ulong n_size, ulong normbits, ulong max_iters)

The limb-level implementation of `fmpz_factor_pollard_brent_single`, in
the shifted representation, for $\mathit{n\_size} \ge 2$ limbs. GCDs are
taken once every 100 iterations on the accumulated product of the
differences, with backtracking if the product becomes $0 \bmod n$.
Returns the number of limbs of the (shifted) factor, or 0.
:::
