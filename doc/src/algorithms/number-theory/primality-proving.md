---
title: Primality proving: APRCL and ECPP
summary: How FLINT proves the primality of large integers, by $n \pm 1$ tests when the factorisations allow it, by the APRCL test with Jacobi sums, and by elliptic curve primality proving with verifiable certificates.
complexity: APRCL $(\log n)^{O(\log \log \log n)}$; ECPP heuristically $\tilde O(\log^4 n)$ (fastECPP), with certificates verifiable in $\tilde O(\log^3 n)$.
order: 25
---

A *probable prime* test such as BPSW
([](algorithms/number-theory/strong-probable-prime)) can show that a
number is composite, but not that it is prime. For word-size integers
the distinction disappears, because the tests have been verified for all
inputs ([](algorithms/number-theory/primality-word)). For larger integers
`fmpz_is_probabprime` runs BPSW, while `fmpz_is_prime` returns a proof:
it decides small inputs directly, tries the classical $n \pm 1$ tests,
which succeed when $n - 1$ or $n + 1$ is sufficiently factored, and
otherwise runs one of two general-purpose algorithms: the APRCL test,
based on Jacobi sums in cyclotomic rings, or elliptic curve primality
proving (ECPP), which also produces a certificate.

## Tests based on $n \pm 1$ {#n-pm-1}

::: theorem Pocklington {#thm-pocklington}
Let $n > 1$ and let $F \mid n - 1$ with $F > \sqrt n$ and $F$ completely
factored. If for every prime $q \mid F$ there is an $a$ with
$a^{n-1} \equiv 1 \pmod n$ and $\gcd(a^{(n-1)/q} - 1, n) = 1$, then $n$
is prime [@Poc1914].
:::

::: proof
Let $p$ be a prime divisor of $n$ and $q^e \| F$. With $a$ as in the
hypothesis, the order of $a$ modulo $p$ divides $n - 1$ but not
$(n-1)/q$, so $q^e$ divides it, and hence $q^e \mid p - 1$. Thus
$F \mid p - 1$ and $p > F > \sqrt n$. A composite $n$ would have a prime
divisor at most $\sqrt n$.
:::

The same idea with Lucas sequences gives Morrison's test for factored
parts of $n + 1$ [@Mor1975], and refinements allow smaller factored
parts: Brillhart, Lehmer and Selfridge [@BLS1975] need only
$F^3 > n$ (with a final check that a quadratic has no integer root), and
Lenstra's algorithm for divisors in residue classes [@Len1984] combines
the factored parts $F_1$ of $n - 1$ and $F_2$ of $n + 1$ when
$\mathrm{lcm}(F_1, F_2)^3 > n$.

`fmpz_is_prime` finds the small prime factors of $n - 1$ and $n + 1$ by
computing $n$ modulo all primes up to a bound (about
$(\log n)^3/100$, ten times more below 260 bits, increased up to three
times by factors of 10 if needed). Before running any of the tests, it
checks whether they can succeed at all: if neither cofactor is a probable
prime and the smooth parts do not satisfy $\mathrm{lcm}(F_1, F_2)^3 > n$,
it goes directly to APRCL or ECPP. For random inputs the $n \pm 1$ tests
succeed often below about 100 bits and essentially never above 260 bits.

## APRCL {#aprcl}

The test of Adleman, Pomerance and Rumely [@APR1983], in the form of
Cohen and Lenstra with Jacobi sums [@CL1984; @CL1987], generalises
Fermat's little theorem to cyclotomic rings. Choose an integer $R$ (a
product of small primes) and let $s$ be the product of the primes $q$
with $q - 1 \mid R$ (and suitable powers of 2 and 3), so that
$a^R \equiv 1 \pmod s$ for all $a$ coprime to $s$; $R$ is chosen so that
$s^2 > n$. Such an $s$ is huge compared with $R$: in FLINT, $R = 5040$
suffices for $n$ of up to 344 bits, and the largest tabulated value
$R = 6983776800$ gives $s$ of about 10000 bits.

For every pair $(p, q)$ with $q \mid s$ prime and $p^k \| q - 1$, let
$\chi$ be a character of order $p^k$ modulo $q$. The test checks, in the
ring $\Z[\zeta_{p^k}]/(n)$, that a certain product of Jacobi sums
$J_0^u J_v$, where $n = p^k u + v$, is a root of unity; if $n$ is prime
this is a consequence of the congruence $\tau(\chi)^{n} \equiv
\tau(\chi^n) \pmod n$ for Gauss sums. An additional condition
$(\mathcal L_p)$ must hold for every prime $p \mid R$; it usually follows
from the roots of unity found (when one of them is primitive), and is
otherwise established by an extra test with another prime $q$.

::: theorem {#thm-aprcl}
If $\gcd(n, Rs) = 1$, all the Jacobi sum conditions hold and
$(\mathcal L_p)$ holds for every prime $p \mid R$, then every divisor
$r$ of $n$ satisfies $r \equiv n^i \pmod s$ for some $0 \le i < R$
[@Coh1996, §9.1].
:::

::: proposition Final division {#prop-final}
Under the conclusion of [](#thm-aprcl), if $s > \sqrt n$ and no residue
$r_i = n^i \bmod s$ with $1 \le i < R$ satisfies $1 < r_i \le \sqrt n$ and
$r_i \mid n$, then $n$ is prime.
:::

::: proof
If $n$ were composite it would have a prime divisor $r$ with
$1 < r \le \sqrt n < s$. Then $r \equiv n^i \pmod s$ for some $i < R$,
and since $0 < r < s$, $r = n^i \bmod s$; $i \ne 0$ because $r > 1$. So
$r$ is one of the residues tested, a contradiction.
:::

The cost is dominated by the powerings $J^u$ in rings of degree
$\varphi(p^k)$ over $\Z/n\Z$, and is $(\log n)^{O(\log \log \log n)}$
[@APR1983]: in theory not polynomial, but in practice comparable to a few
hundred to a few thousand modular exponentiations. The tests for the
different pairs $(p, q)$ are independent and run in parallel.

In FLINT (`aprcl_is_prime`): $R$ comes from a table by the bit size of
$n$ (`aprcl_R_value`), and $s$ is then reduced by removing the prime
powers whose tests cost the most per bit of $s$, as long as $s^2 > n$
remains. The arithmetic in $\Z[\zeta_{p^k}]/(n)$ uses fixed-size
`mpn_mod` coefficients and delayed reduction for moduli of 2 to 16 limbs,
and hand-written formulas for the multiplications in degrees 2 to 8 for
$p^k \in \{3, 4, 5, 7, 8, 9, 16\}$. The final division walks through the
residues $n^i \bmod s$ in parallel. The table of $R$ goes up to inputs of
about 20000 bits; for larger $n$, `aprcl_config_jacobi_init` multiplies
$R$ by further primes until $s^2 > n$. The test functions check
$s^2 > n$ before they report a prime.

## Elliptic curve primality proving {#ecpp}

::: theorem Goldwasser–Kilian {#thm-gk}
Let $n > 1$ be coprime to 6, $E: y^2 = x^3 + ax + b$ with
$\gcd(4a^3 + 27b^2, n) = 1$, and $m, q$ integers with $q \mid m$ and
$q > (n^{1/4} + 1)^2$. Suppose there is a point $P$ on $E$ modulo $n$
such that $mP = O$ and $(m/q) P$ is a point whose $Z$-coordinate is a
unit modulo $n$ (that is, $(m/q)P \ne O$ modulo every prime divisor of
$n$). If $q$ is prime, then $n$ is prime [@GolKil1986].
:::

::: proof
Let $p$ be a prime divisor of $n$ with $p \le \sqrt n$ and consider the
reduction of $E$ and $P$ modulo $p$, an elliptic curve over $\mathbb F_p$.
The point $Q = (m/q) P$ satisfies $qQ = mP = O$ and $Q \ne O$ modulo $p$,
so $Q$ has order $q$, as $q$ is prime. Hence
$q \le \#E(\mathbb F_p) \le p + 1 + 2\sqrt p = (\sqrt p + 1)^2
\le (n^{1/4} + 1)^2$ by Hasse's theorem, contradicting the hypothesis on
$q$. So $n$ has no prime divisor at most $\sqrt n$.
:::

An ECPP *certificate* is a chain of such steps $n = n_0, n_1 = q_0,
n_2 = q_1, \ldots$ ending with a number small enough to be checked
directly (below $2^{64}$ in FLINT). Checking a step costs two scalar
multiplications on the curve, so the whole certificate is verified in
time $\tilde O(\log^3 n)$, much faster than it was found.

The difficulty is to find curves whose order $m$ is a small multiple of a
probable prime. Counting points on random curves is too slow; the method
of Atkin and Morain [@AtkMor1993] constructs curves with known order by
complex multiplication. If $4n = t^2 + |D| v^2$ for a negative
discriminant $D$ (found by Cornacchia's algorithm from $\sqrt D \bmod n$),
then the curves over $\mathbb F_n$ with complex multiplication by the
order of discriminant $D$ have $n + 1 \pm t$ points (and further twists
for $D = -3, -4$). Their $j$-invariants are the roots modulo $n$ of the
Hilbert class polynomial $H_D$, whose degree is the class number $h(D)$.
A step therefore consists of finding a $D$ for which one of the orders
$n + 1 \pm t$ is $q$ times a small number, with $q$ a probable prime
larger than $(n^{1/4} + 1)^2$, then computing a root of $H_D$ modulo $n$,
the curve, its correct twist and a point.

The fastECPP organisation of the search [@FKMW2004; @Morain2007] uses
many discriminants composed of a set of small primes, so that $\sqrt D$
modulo $n$ is a product of precomputed square roots, and removes the
small factors of all candidate orders at once with a remainder tree.
FLINT also uses the class field tower of Enge and Morain
[@EngeMorain2003] and Weber's functions as class invariants to find the
root of $H_D$ in steps of prime degree, and factors of $H_D$ over the
genus field. The heuristic running time of fastECPP is
$\tilde O(\log^4 n)$, against $\tilde O(\log^5 n)$ for the original
Atkin–Morain algorithm. The details are in [](ecpp/prove) and
[](ecpp).

## Choice of the method {#choice}

For $n$ beyond two words, `fmpz_is_prime` proceeds as follows:

1. Trial division by about as many primes as $n$ has bits, a test for
   being a square, and a base-2 strong probable prime test, which reject
   almost all composites cheaply.
2. The $n \pm 1$ tests, if the feasibility check says they can succeed.
3. Otherwise ECPP if $n$ has at least as many bits as a threshold that
   depends on the number of threads available (250 bits for one thread,
   550 for two, 1000 for three, 2000 for four, 8000 for up to eight and
   12000 beyond), and APRCL below it or if ECPP fails. APRCL parallelises
   better than the largely serial ECPP chain, which moves the crossover
   up with more threads.

For two-word $n$, `n_ll_is_prime` first does trial division and a base-2
strong probable prime test and can certify many inputs directly; the
others go through steps 2 and 3. Word-size $n$ are decided by
`n_is_prime`. To prove a number prime and keep the proof,
use `ecpp_prove`; a certificate can be checked by `ecpp_verify` or by
PARI/GP.

::: functions fmpz_is_prime aprcl_is_prime ecpp_prove ecpp_verify :::
