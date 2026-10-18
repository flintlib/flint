---
title: Factoring word-size integers
kind: topic
summary: How `n_factor` combines trial division, primality testing, perfect power detection, Pollard–Brent rho, Hart's one-line method and SQUFOF to factor any 64-bit integer.
order: 3
---

`n_factor` completely factors any word-size integer. Small factors are
removed by trial division; what remains has no prime factor below
$2^{16}$, so it has at most three prime factors and is split by one of
several special-purpose methods.

## The pipeline {#pipeline}

1. Remove the power of 2 by counting trailing zeros.
2. If $n < 2^{15}$ and $n$ is prime (by table lookup), stop.
3. Trial-divide by the primes from 3 to 311, using the multiplication test
   of [](algorithms/arithmetic/divisibility-test) and exact division. Stop as soon as
   the cofactor is at most the square of the last prime tried: it is then
   1 or prime.
4. Test the cofactor with `n_is_prime_odd_no_trial`; stop if it is prime.
5. Continue trial division with the remaining primes below $2^{16}$
   (6542 primes in total), again stopping early when possible. On 64-bit
   machines a cofactor below $65537^2$ is now prime.
6. Split the remaining composite cofactors, keeping a stack of factors that
   still need work. For each one:
   - replace it by its root if it is a square, cube or fifth power
     (`n_factor_power235`);
   - accept it if it is prime;
   - otherwise find a factor with the first of these that succeeds:
     Pollard–Brent rho [@Pol1975; @Bre1980], iterating in
     [Montgomery arithmetic](algorithms/arithmetic/montgomery), with one attempt of 32768 iterations, if
     the number exceeds $2^{38}$; Hart's one-line method [@Har2012] with
     40000 iterations, if the number is below $2^{50}$; SQUFOF
     [@GowWag2008] with 50000 iterations.

If all three methods fail, an exception is thrown. This has never been
observed.

## Why these methods {#choices}

After step 5, a composite cofactor $m < 2^{64}$ has all prime factors above
$2^{16}$, so it is a product of two or three primes. Rho finds a factor $p$
in about $\sqrt p$ steps, which is fast for the unbalanced splits typical of
large $m$. Hart's one-line method is a Fermat-type method that is very fast
for small $m$ and for factors close to $\sqrt m$, which is exactly the case
where rho is slowest. SQUFOF has a running time of about $m^{1/4}$ that
does not depend on the shape of the factorisation, which makes it a
dependable last resort.

The $p + 1$ method [@Wil1982] (`n_factor_pp1`) and the elliptic curve method
[@Len1987; @Mon1987] (`n_factor_ecm`) are available but not currently used
by `n_factor`. For background on all of these methods see [@CraPom2005,
Chapters 5–7].

For two-word integers, `n_ll_factor_rho` and `n_ll_factor_ecm` run rho and
ECM in two-word Montgomery arithmetic, and `n_ll_factor_SQUFOF` runs SQUFOF on
$kn$ for a sequence of small multipliers $k$. Their use by `fmpz_factor`,
and the methods for larger integers, are described in
[](algorithms/number-theory/integer-factoring).

::: functions n_factor n_factor_partial n_factor_trial :::
