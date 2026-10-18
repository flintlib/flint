---
title: Primality testing for word-size integers
kind: topic
summary: How `n_is_prime` proves primality of any 64-bit integer with at most two strong probable prime tests, using table lookups, vectorised trial division and a hashed table of witnesses.
order: 2
---

Every primality test in `ulong_extras` that answers yes or no is
deterministic: for word-size integers, probable primality and primality
coincide because the relevant pseudoprimes have been enumerated. The
default function `n_is_prime` chooses its method by the size of $n$.

## The method by size {#by-size}

| Range of $n$ | What `n_is_prime` does |
|---|---|
| $n < 2^{15}$ | Looks up odd $n$ in a bit table. |
| all larger $n$ | Trial division by 2, 3 and 5. |
| $2^{15} \le n < 2^{20}$ | Looks up $n$ in a table indexed modulo 30. |
| $2^{20} \le n$ | Trial division by the 32 primes from 7 to 149. |
| $2^{20} \le n < 2^{32}$ | Strong probable prime test to base 2 in 32-bit arithmetic, then a hash lookup among the 2314 base-2 strong pseudoprimes below $2^{32}$. |
| $2^{32} \le n < 2^{64}$ | Base-2 test in [Montgomery arithmetic](algorithms/arithmetic/montgomery), a check for one special form of semiprime, and a second test with a base taken from a table of $2^{15}$ witnesses. |

All trial division uses the multiplication trick of
[](algorithms/arithmetic/divisibility-test), so it involves no division instructions.
For $n < 2^{32}$ the 32 divisibility tests have no branches and are
vectorised by the compiler.

::: summary algorithms/number-theory/strong-probable-prime :::

## Choosing the second base {#witness}

For $n > 2^{32}$ that pass the base-2 test, `n_is_prime` computes a 15-bit
hash of $n$ and reads a base $a < 2^{16}$ from a table of 32768 entries
(64 KiB). The table was built from the Feitsma–Galway list of base-2 strong
pseudoprimes [@FeiGal2013] so that, in each hash bucket, the stored base
exposes every pseudoprime in that bucket. This follows the approach of
Forišek and Jančina [@ForJan2015], whose table is eight times larger.

One family of pseudoprimes makes this hard: semiprimes $n = pq$ with
$q = k(p - 1) + 1$ for $k \in \{2, \ldots, 9, 12\}$, for which reliable
witnesses are rare. These are detected before the second test. Writing
$p = x + 1$ gives $n = (x + 1)(kx + 1)$, so for each $k$ the candidate
$x = \lfloor \sqrt{n/k} \rfloor$ is computed from a single floating-point
square root and checked. This costs about a twentieth of a Fermat test.

The correctness of the whole procedure up to $2^{64}$ can be verified with
the program `examples/check_n_is_prime`, which runs `n_is_prime` on every
entry of the Feitsma–Galway table.

## Other tests {#others}

The other tests in the module are kept for comparison, for other
applications, or as certificates:

- `n_is_probabprime_BPSW` (Baillie–PSW [@PSW1980; @BaiWag1980]) has no
  known counterexample at any size;
- `n_is_prime_pocklington` [@Poc1914; @BLS1975] and
  `n_is_prime_pseudosquare` [@LukPatWil1996] produce proofs that do not
  depend on tables;
- `n_ll_is_prime` handles two-word integers, deterministically below about
  $2^{81}$, with the double-word exponentiations `n_ll_small_2_powmod` and
  `n_ll_small_powmod_triple` reduced with the two-word Barrett
  inverse $\lfloor B^3/n \rfloor$ computed by `n_ll_small_preinv`.
