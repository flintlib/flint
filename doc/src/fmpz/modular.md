---
title: Modular arithmetic with integers
functions: fmpz_mod fmpz_smod fmpz_invmod fmpz_powm_ui fmpz_sqrtmod
---

These functions reduce, invert, exponentiate and take square roots of
integers modulo an integer $m$ given as an `fmpz_t`. Each call is
self-contained: no data about $m$ is kept between calls. They are the
right choice for a few operations; for many operations with the same
modulus the module [](fmpz_mod) (which stores precomputed data in a
context) or, for a word-size modulus, [](nmod) is faster.

## Residues {#residues}

`fmpz_mod` returns the residue in $[0, |m|)$, and `fmpz_smod` the
symmetric residue $y$ with $-|m|/2 < y \le |m|/2$. Both require $m \ne 0$
without checking it. `fmpz_mod_ui`, `fmpz_fdiv_ui` and `fmpz_get_nmod`
reduce modulo a word. For the results of `fmpz_invmod`, `fmpz_powm` and
`fmpz_sqrtmod`, the inputs may be any integers; they are reduced first.

## Inversion {#invmod}

`fmpz_invmod` returns 1 and the inverse in $[0, |m|)$ if $\gcd(g, m) = 1$,
and 0 otherwise, so it doubles as a coprimality test. Modulo $\pm 1$
every $g$ is invertible with inverse 0. Word-size moduli use `n_gcdinv`;
larger ones GMP's `mpz_invert`, an extended GCD
([](algorithms/integers/gcd)). For a modulus whose factors are needed
when the inverse does not exist, use `fmpz_gcdinv` or the extended GCD
functions of [](fmpz/gcd).

## Exponentiation {#powm}

`fmpz_powm` and `fmpz_powm_ui` require $m \ge 1$ and return $g^e \bmod m$
in $[0, m)$; for $m = 1$ the result is 0 even for $e = 0$. `fmpz_powm`
also accepts $e < 0$ and then inverts $g$ first, raising an exception if
$g$ is not invertible. The algorithm depends on the sizes:

- a word-size modulus with a word-size base: `n_powmod2_ui_preinv` (or
  GMP for an exponent that is not a word);
- a large base with a word-size modulus, or a large modulus below 70000
  bits: GMP's `mpz_powm` and `mpz_powm_ui`, which use Montgomery
  reduction [@Mon1985] with a sliding window;
- $m$ of at least 70000 bits, in builds with `fft_small`: binary or
  sliding-window exponentiation in the `gr` context of $\Z/m\Z$
  (`gr_ctx_init_fmpz_mod`), where each modular multiplication is a
  product followed by a division with a precomputed inverse of $m$
  (`fmpz_fdiv_r_preinvn`), so that both steps use FLINT's FFT
  multiplication ([](fft_small/usage#callers));
- exponents $e \le 4$ with a large modulus are done by one to three
  multiplications and reductions.

## Square roots {#sqrtmod}

`fmpz_sqrtmod` computes a square root modulo a prime $p$. It returns 1
and a root if $a$ is a square modulo $p$ (0 and 1 are their own roots),
and 0 with $b = 0$ otherwise. It chooses, as `nmod_sqrt` does for
words ([](nmod/sqrt)), between $a^{(p+1)/4}$ for $p \equiv 3 \pmod 4$,
Atkin's formula for $p \equiv 5 \pmod 8$, Tonelli–Shanks when $2^s \|
p - 1$ for $s < 32$, and Müller's Lucas sequence variant of the
Cipolla–Lehmer method otherwise ([](algorithms/number-theory/sqrt-mod-p)).
Primality of $p$ is not checked: for even $p$ and for perfect squares $p$
the function returns 0, and for other composite $p$ a return value of 1
does not guarantee anything.

## Example {#example}

Computations modulo the Mersenne prime $p = 2^{127} - 1$, and modulo 12:

```c example name=fmpz_modular
#include <flint/flint.h>
#include <flint/fmpz.h>

int main(void)
{
    fmpz_t p, a, b, e, x;

    fmpz_init(p); fmpz_init(a); fmpz_init(b); fmpz_init(e); fmpz_init(x);

    /* p = 2^127 - 1 is prime */
    fmpz_one_2exp(p, 127);
    fmpz_sub_ui(p, p, 1);

    /* Fermat: 3^(p-1) = 1 mod p */
    fmpz_set_ui(a, 3);
    fmpz_sub_ui(e, p, 1);
    fmpz_powm(x, a, e, p);
    flint_printf("3^(p-1) mod p = "); fmpz_print(x); flint_printf("\n");

    /* inverses: 3^-1 mod p, and a negative exponent */
    fmpz_invmod(x, a, p);
    flint_printf("3^-1 mod p = "); fmpz_print(x); flint_printf("\n");
    fmpz_set_si(e, -2);
    fmpz_powm(b, a, e, p);
    fmpz_mul(x, x, x);
    fmpz_mod(x, x, p);
    flint_printf("3^-2 mod p agrees: %d\n", fmpz_equal(x, b));

    /* not invertible modulo 12 */
    fmpz_set_ui(b, 12);
    fmpz_set_si(a, -9);
    flint_printf("invmod(-9, 12) = %d\n", fmpz_invmod(x, a, b));
    fmpz_mod(x, a, b);
    flint_printf("-9 mod 12 = "); fmpz_print(x);
    fmpz_smod(x, a, b);
    flint_printf(", smod = "); fmpz_print(x); flint_printf("\n");

    /* square roots modulo p */
    fmpz_set_ui(a, 2);
    if (fmpz_sqrtmod(x, a, p))
    {
        fmpz_mul(b, x, x);
        fmpz_mod(b, b, p);
        flint_printf("sqrt(2) mod p = "); fmpz_print(x);
        flint_printf(", squared: "); fmpz_print(b); flint_printf("\n");
    }
    fmpz_set_si(a, -1);
    flint_printf("-1 is a square mod p: %d\n", fmpz_sqrtmod(x, a, p));

    fmpz_clear(p); fmpz_clear(a); fmpz_clear(b); fmpz_clear(e); fmpz_clear(x);
    return 0;
}
```
```output
3^(p-1) mod p = 1
3^-1 mod p = 113427455640312821154458202477256070485
3^-2 mod p agrees: 1
invmod(-9, 12) = 0
-9 mod 12 = 3, smod = 3
sqrt(2) mod p = 18446744073709551616, squared: 2
-1 is a square mod p: 0
```

Since $2^{128} = 2 \cdot 2^{127} \equiv 2 \pmod p$, the root found is
$2^{64}$. As $p \equiv 3 \pmod 4$, $-1$ is not a square.
