---
title: Factoring integers
kind: topic
summary: The factoring methods available in FLINT (trial division, Pollard rho, $p \pm 1$, ECM, SQUFOF and the quadratic sieve), their running times, and how `fmpz_factor` combines them.
order: 30
---

No factoring algorithm is best for all inputs. The cost of the *special
purpose* methods (trial division, Pollard rho, $p - 1$, $p + 1$ and the
elliptic curve method) depends mainly on the size of the factor $p$ they
find, while that of the *general purpose* methods (SQUFOF, the quadratic
sieve) depends only on the size of the number $N$ itself. A factoring
routine therefore removes small factors first and calls the general
methods only on what remains.

## The methods {#methods}

We write $L_x[\alpha, c] = \exp\bigl((c + o(1)) (\log x)^\alpha
(\log \log x)^{1 - \alpha}\bigr)$. The running times below are in
arithmetic operations modulo $N$ (or, for the sieve, in operations on
small numbers), and those of the randomised methods are heuristic expected
values [@CraPom2005, Chapters 5–6].

| Method | Finds | Time | In FLINT |
|---|---|---|---|
| Trial division | primes $p \le b$ | $\pi(b)$ divisions, or one product tree | `n_factor_trial`, `fmpz_factor_trial`, `flint_mpn_factor_trial_tree` |
| Pollard–Brent rho | any $p$ | $O(\sqrt p)$ | `n_factor_pollard_brent`, `n_ll_factor_rho`, `fmpz_factor_pollard_brent` |
| Pollard $p - 1$ | $p$ with $p - 1$ smooth | $O(B_1)$ | `n_factor_pp1` (as a degenerate $p + 1$), `fmpz_factor_pp1` |
| Williams $p + 1$ | $p$ with $p \pm 1$ smooth | $O(B_1)$ | `n_factor_pp1`, `fmpz_factor_pp1` |
| Elliptic curves | any $p$ | $L_p[1/2, \sqrt 2]$ | `n_factor_ecm`, `n_ll_factor_ecm`, `fmpz_factor_ecm` |
| Fermat-type (Hart, Lehman) | factors near $\sqrt N$ | $O(N^{1/3})$ | `n_factor_lehman`, Hart's one-line method in `n_factor` |
| SQUFOF | any $N$ | $O(N^{1/4})$ | `n_factor_SQUFOF`, `n_ll_factor_SQUFOF` |
| Quadratic sieve | any $N$ | $L_N[1/2, 1]$ | `qsieve_factor` |

**Trial division** by all primes up to $b$ finds the factors below $b$.
For multiprecision $N$, FLINT divides by the primes one at a time
(`flint_mpn_factor_trial`, a division of $N$ by a word for each prime) or
computes $\gcd(N, \prod_{p \le b} p)$ with a product tree of the primes,
which is cached (`flint_mpn_factor_trial_tree`), and then identifies the
primes dividing the GCD.

**Pollard's rho method** [@Pol1975] iterates $x \mapsto x^2 + a \bmod N$.
Modulo a prime factor $p$ the sequence becomes periodic after about
$\sqrt{\pi p/2}$ steps, and a collision $x_i \equiv x_j \pmod p$ is detected
as $\gcd(x_i - x_j, N) > 1$. Brent's variant [@Bre1980] finds the cycle
with powers of two and multiplies the differences together, taking a GCD
only every 100 steps.

**The $p - 1$ and $p + 1$ methods** [@Pol1974; @Wil1982] compute
$a^E \bmod N$ (respectively a Lucas sequence) for $E$ the product of all
prime powers up to a bound $B_1$; if $p - 1$ (respectively $p + 1$)
divides $E$, then $p$ divides $\gcd(a^E - 1, N)$. A second stage allows
one larger prime up to $B_2$. They succeed only for special $p$.

**The elliptic curve method** [@Len1987] replaces the group
$(\Z/p\Z)^*$ of order $p - 1$ by the group of points of a random elliptic
curve modulo $p$, whose order varies in $[p + 1 - 2\sqrt p, p + 1 +
2\sqrt p]$: trying many curves gives many chances of a smooth order. See
[](algorithms/number-theory/ecm).

**SQUFOF** [@GowWag2008] walks through the continued fraction expansion
of $\sqrt{kN}$ looking for a square form, using only arithmetic on numbers
of half the size of $N$; it is the method of choice for $N$ of up to
about 70 bits.

**The quadratic sieve** [@Pom1985] finds many $x$ for which
$x^2 - kN$ factors over a base of small primes, combines them into a
congruence of squares $X^2 \equiv Y^2 \pmod N$ by linear algebra over
$\mathbb F_2$, and splits $N$ with $\gcd(X - Y, N)$. See
[](algorithms/number-theory/quadratic-sieve).

FLINT has no implementation of the number field sieve, whose running time
$L_N[1/3, (64/9)^{1/3}]$ makes it faster than the quadratic sieve above
about 100 digits.

## Words: `n_factor` {#words}

For $N < 2^{64}$, `n_factor` removes the primes up to $2^{16}$ by trial
division and then splits the remaining cofactor, which has at most three
prime factors, with rho, Hart's one-line method and SQUFOF. See
[](algorithms/number-theory/factoring-word). Two-word versions of rho,
ECM and SQUFOF (`n_ll_factor_rho`, `n_ll_factor_ecm`,
`n_ll_factor_SQUFOF`) are used by `fmpz_factor` for cofactors of up to 128
bits.

## Multiprecision integers: `fmpz_factor` {#fmpz}

`fmpz_factor` proceeds as follows ([](fmpz_factor/factor) has more
details).

1. Remove the powers of 2, then trial divide by the primes from 3 in
   blocks of 1000 primes, continuing with the next block as long as
   factors are found. A cofactor of one word is finished with `n_factor`.
2. Check whether the cofactor is a proven prime (`fmpz_is_prime`) or a
   perfect power (`fmpz_is_perfect_power`); in the latter case factor the
   root.
3. Run `_fmpz_factor_smooth` with the bound $\max(b/3 - 17, 2)$ bits for
   a $b$-bit cofactor: trial division by the first 3512 primes with a
   product tree, then the elliptic curve method with the parameters of
   the table below for factors of increasing size (from about 18 bits, in
   steps of 6 bits), lowering the bound as the cofactor shrinks.
4. Split what remains: on 64-bit machines, for cofactors of at most 128
   bits, by two-word rho (2048 or 8192 iterations), from 118 bits a few
   two-word ECM curves, and up to 72 bits by SQUFOF; otherwise, or if
   these fail, by the quadratic sieve.
5. Factor the pieces recursively, from step 2.

The choice of the bound $b/3 - 17$ in step 3 reflects the costs: the
elliptic curve method finds a factor of $f$ bits in time
$L_{2^f}[1/2, \sqrt 2]$, while the sieve takes $L_{2^b}[1/2, 1]$; looking
for factors of up to about a third of the size of the number before
switching to the sieve is a good compromise, and the constant was tuned
on random inputs. The table of elliptic curve parameters
(`fmpz_factor_smooth`, with $B_2 = 100 B_1$) is designed to find a
factor of the given size with probability about $1/3$ at each size, so
that over the overlapping sizes a factor of $f \le$ the bound bits is
found with probability at least about $2/3$:

| Factor bits | 20 | 30 | 40 | 50 | 60 | 70 | 80 | 90 | 100 |
|---|---|---|---|---|---|---|---|---|---|
| $B_1$ | 19 | 102 | 610 | 3850 | 14000 | 26000 | 70000 | 240000 | 470000 |
| Curves | 5 | 13 | 18 | 21 | 42 | 117 | 190 | 291 | 430 |

All prime factors returned by `fmpz_factor` are proved prime, by the
methods of [](algorithms/number-theory/primality-proving).
`fmpz_factor_smooth` with `proved` = 0 checks only probable primality,
which is much faster for large factors.

::: functions fmpz_factor fmpz_factor_smooth n_factor qsieve_factor fmpz_factor_ecm :::
