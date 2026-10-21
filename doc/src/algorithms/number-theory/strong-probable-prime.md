---
title: Strong probable prime test
summary: The Miller–Rabin test for a single base. Every odd prime passes it, at most a quarter of the bases let a composite pass, and with tabulated pseudoprimes it becomes a proof of primality for word-size integers.
complexity: One modular exponentiation with an exponent below $n$, about $\log_2 n$ modular squarings.
order: 40
---

The strong probable prime test, also called the Miller–Rabin test,
checks with one modular exponentiation a condition that every odd prime
satisfies. It goes back to Miller [@Mil1976], who showed that it proves
primality under the generalised Riemann hypothesis when enough bases are
tried, and Rabin [@Rab1980], who turned it into a randomised test.

## Definition {#definition}

::: definition {#def-sprp}
Let $n > 2$ be odd and write $n - 1 = 2^s d$ with $d$ odd. Then $n$ is a
*strong probable prime to base $a$*, or an $a$-SPRP, if

$$
a^d \equiv 1 \pmod n \quad\text{or}\quad a^{2^r d} \equiv -1 \pmod n
\text{ for some } 0 \le r < s.
$$

A composite $a$-SPRP is a *strong pseudoprime* to base $a$.
:::

::: theorem {#thm-prime-passes}
If $n$ is an odd prime and $n \nmid a$, then $n$ is an $a$-SPRP.
:::

::: proof
Let $x_r = a^{2^r d} \bmod n$ for $0 \le r \le s$, so that $x_{r+1} = x_r^2$.
By Fermat's little theorem $x_s = a^{n-1} \bmod n = 1$. Let $r$ be the least
index with $x_r = 1$. If $r = 0$ we are done. Otherwise $x_{r-1}^2 = 1$ and
$x_{r-1} \ne 1$. Since $\Z/n\Z$ is a field, the only square roots of 1 are
$\pm 1$, so $x_{r-1} = n - 1$.
:::

The proof shows what the test detects: a composite $n$ is exposed when the
sequence $x_0, x_1, \ldots$ reaches 1 through a square root of 1 other than
$\pm 1$, or never reaches 1 at all.

::: theorem Monier, Rabin {#thm-rabin}
If $n > 9$ is odd and composite, then $n$ is a strong pseudoprime to at most
$\varphi(n)/4$ bases $1 \le a < n$.
:::

See [@Mon1980; @Rab1980] for the original proofs and [@CraPom2005] for an
exposition. Hence $k$ random bases
detect a composite with probability at least $1 - 4^{-k}$, but FLINT never
relies on this bound for word-size integers.

## The test {#test}

::: algorithm Strong probable prime test {#alg-sprp}
**Input:** odd $n > 2$, a base $1 < a < n - 1$, and the odd part $d$ of $n - 1$.

**Output:** whether $n$ is an $a$-SPRP.

1. $y \gets a^d \bmod n$.
2. If $y = 1$ or $y = n - 1$, return true.
3. Repeat $s - 1$ times: set $y \gets y^2 \bmod n$; if $y = n - 1$, return
   true; if $y = 1$, return false.
4. Return false.
:::

Returning false early when $y = 1$ in step 3 is valid because 1 is then
reached through a square root of 1 other than $\pm 1$.

## Making it deterministic {#deterministic}

A single base is not a proof, but strong pseudoprimes are rare enough to be
listed completely in a range. Feitsma and Galway [@FeiGal2013] computed all
31894014 strong pseudoprimes to base 2 below $2^{64}$. FLINT uses this
in two ways.

- Below $2^{32}$, there are only 2314 base-2 strong pseudoprimes. After a
  base-2 test, `n_is_prime` looks $n$ up in a hash table of them.
- Above $2^{32}$, a second base is chosen from a table indexed by a hash of
  $n$, in the spirit of Forišek and Jančina [@ForJan2015]. Each table entry
  was chosen so that it is a witness for every base-2 pseudoprime that
  hashes to its slot. See [](algorithms/number-theory/primality-word) for the details.

The function `n_ll_is_prime` uses a different known result: for $n$ below
about $3.3 \cdot 10^{24}$, testing the first 13 primes as bases is a proof
[@SorWeb2016].

## Arithmetic {#arithmetic}

The test is one modular exponentiation followed by up to $s - 1$
squarings, so its speed is that of modular multiplication.
`n_is_strong_probabprime2_preinv` uses [Möller–Granlund
reduction](algorithms/arithmetic/moller-granlund) and
`n_is_strong_probabprime_precomp` floating-point reduction. Inside
`n_is_prime`, the tests for $n \ge 2^{32}$ run entirely in
[Montgomery form](algorithms/arithmetic/montgomery), which avoids the
normalisation shifts and lets the comparisons with $\pm 1$ be done on
Montgomery representatives directly.

::: functions n_is_strong_probabprime2_preinv n_is_strong_probabprime_precomp :::
