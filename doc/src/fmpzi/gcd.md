---
title: Greatest common divisors of Gaussian integers
functions: fmpzi_gcd
---

Since $\mathbb{Z}[i]$ is a Euclidean domain, any two Gaussian integers
$x, y$ have a greatest common divisor $g$, unique up to multiplication
by one of the units $\pm 1, \pm i$, and $(x, y) = g\mathbb{Z}[i]$. All
five functions return the same $g$: the associate with positive real
part $a$ and imaginary part in $(-a, a]$, so their results can be
compared directly. They differ only in the algorithm.

::: summary algorithms/number-theory/gaussian-integers :::

## Algorithms {#algorithms}

- `fmpzi_gcd_euclidean` runs Euclid's algorithm with the rounded
  quotients of `fmpzi_divrem`; each step at least halves the norm.
- `fmpzi_gcd_euclidean_improved` uses the cheaper approximate quotients
  of `fmpzi_divrem_approx`, and as soon as all parts are below $2^{50}$
  it finishes in double-precision arithmetic, where the remainders are
  still computed exactly.
- `fmpzi_gcd_binary` is Weilert's $(1+i)$-ary analogue of the binary GCD
  [@Wei2000], which replaces divisions by additions and divisions by
  $1 + i$. It is implemented for inputs of up to 500 bits and calls
  `fmpzi_gcd_euclidean` for larger ones.
- `fmpzi_gcd_shortest` computes the Hermite normal form of the lattice
  spanned by $x, ix, y, iy$ and reads off $g$ as a shortest vector,
  using a quasi-linear continued fraction algorithm.
- `fmpzi_gcd` uses the double-precision Euclidean algorithm when all
  parts are below $2^{50}$, `fmpzi_gcd_euclidean_improved` when an input
  fits in single words or has at most 30000 bits, and otherwise
  `fmpzi_gcd_shortest`, preceded by one division if the sizes are very
  different.

## Example {#example}

With $x = 3(1 + i)(2 + i)^3$ and $y = 3(2 + i)(2 - i)^2$, the GCD is
$3(2 + i) = 6 + 3i$, since $2 + i$ and $2 - i$ are non-associate primes
above 5. The program also shows a division with remainder, the removal
of the factor $1 + i$, and primality tests.

```c example name=fmpzi_gcd_basic
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpzi.h>

int main(void)
{
    fmpzi_t x, y, g, h, q, r, t;
    fmpz_t n;
    slong k;

    fmpzi_init(x); fmpzi_init(y); fmpzi_init(g); fmpzi_init(h);
    fmpzi_init(q); fmpzi_init(r); fmpzi_init(t);
    fmpz_init(n);

    /* x = 3 (1 + i) (2 + i)^3,  y = 3 (2 + i) (2 - i)^2 */
    fmpzi_set_si_si(t, 2, 1);
    fmpzi_pow_ui(x, t, 3);
    fmpzi_set_si_si(t, 3, 3);
    fmpzi_mul(x, x, t);
    fmpzi_set_si_si(t, 2, -1);
    fmpzi_sqr(y, t);
    fmpzi_set_si_si(t, 6, 3);
    fmpzi_mul(y, y, t);
    flint_printf("x = "); fmpzi_print(x);
    flint_printf(", y = "); fmpzi_print(y); flint_printf("\n");

    fmpzi_gcd(g, x, y);
    flint_printf("gcd = "); fmpzi_print(g); flint_printf("\n");
    fmpzi_gcd_binary(h, x, y);
    flint_printf("binary: %d", fmpzi_equal(g, h));
    fmpzi_gcd_shortest(h, x, y);
    flint_printf(", shortest: %d", fmpzi_equal(g, h));
    fmpzi_gcd_euclidean(h, x, y);
    flint_printf(", euclidean: %d\n", fmpzi_equal(g, h));

    fmpzi_divrem(q, r, x, y);
    fmpzi_norm(n, r);
    flint_printf("x = q*y + r with q = "); fmpzi_print(q);
    flint_printf(", r = "); fmpzi_print(r);
    flint_printf(", N(r) = "); fmpz_print(n);
    fmpzi_norm(n, y);
    flint_printf(", N(y) = "); fmpz_print(n); flint_printf("\n");

    k = fmpzi_remove_one_plus_i(t, x);
    flint_printf("x = (1+i)^%wd * (", k); fmpzi_print(t); flint_printf(")\n");

    fmpzi_set_si_si(t, 2, 1);
    flint_printf("prime: 2+i %d", fmpzi_is_prime(t));
    fmpzi_set_si_si(t, 3, 0);
    flint_printf(", 3 %d", fmpzi_is_prime(t));
    fmpzi_set_si_si(t, 5, 0);
    flint_printf(", 5 %d\n", fmpzi_is_prime(t));

    fmpzi_clear(x); fmpzi_clear(y); fmpzi_clear(g); fmpzi_clear(h);
    fmpzi_clear(q); fmpzi_clear(r); fmpzi_clear(t);
    fmpz_clear(n);
    return 0;
}
```
```output
x = -27+39*I, y = 30-15*I
gcd = 6+3*I
binary: 1, shortest: 1, euclidean: 1
x = q*y + r with q = -1+1*I, r = -12-6*I, N(r) = 180, N(y) = 1125
x = (1+i)^1 * (6+33*I)
prime: 2+i 1, 3 1, 5 0
```

The remainder satisfies $N(r) = 180 \le N(y)/2$. The rational prime 3
stays prime in $\mathbb{Z}[i]$ because $3 \equiv 3 \pmod 4$, while
$5 = (2 + i)(2 - i)$ splits.
