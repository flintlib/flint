---
title: Square roots, n-th roots and perfect powers
functions: fmpz_sqrt fmpz_perfect_sqrt fmpz_is_square fmpz_root fmpz_is_perfect_power
---

`fmpz_sqrt` and `fmpz_sqrtrem` compute the integer square root
$s = \lfloor \sqrt g \rfloor$ and the remainder $g - s^2 \in [0, 2s]$.
When only squares are of interest, `fmpz_perfect_sqrt` is the function
to use: it decides whether $a$ is a square and returns the root if it is,
rejecting almost all non-squares without computing a root.
`fmpz_is_square` only answers the question. `fmpz_root` computes
$n$-th roots rounded towards zero, and `fmpz_is_perfect_power` decides
whether an integer is a nontrivial power.

## Square roots {#sqrt}

Inline values use `n_sqrt` (a floating-point square root with a
correction). Large values use FLINT's own square root code
(`_flint_mpn_sqrtrem`) on 64-bit machines:

- one to four limbs: the hardware square root of a `double`, a Newton
  correction and one step of Zimmermann's divide-and-conquer square root,
  specialised to a two-limb root;
- up to `FLINT_MPN_SQRTREM_NEWTON_CUTOFF` limbs (7423 on recent x86-64):
  the divide-and-conquer ("Karatsuba") square root of Zimmermann
  [@Zim1999], as in GMP, which computes the root of the top half
  recursively and the remaining limbs by one division, at a cost of a
  small multiple of a multiplication;
- above: Newton iteration for $1/\sqrt{\alpha}$ in fixed-point
  arithmetic with the Karp–Markstein trick [@KM1997], followed by a
  check of the last limb that certifies the integer part.

The algorithms are described in [](algorithms/integers/division#sqrt).
`fmpz_sqrt` does not form the remainder, which saves a squaring.

## Recognising squares {#squares}

`fmpz_is_square` and `fmpz_perfect_sqrt` first reject a large input
unless it is a square modulo 256 (low limb), then modulo the factors of
$2^{48} - 1$ ($4095 = 9 \cdot 5 \cdot 7 \cdot 13$, $17$, $257$, $97$ and
$241$, from a single reduction modulo $2^{48} - 1$), which lets about one
random non-square in 900 through. A surviving candidate is decided by
computing the root: for inputs below $2^{100}$, rounding the
floating-point root is exact for squares; otherwise the divide-and-conquer
root is computed with an early exit as soon as its low bits show that the
input cannot be a square. `fmpz_perfect_sqrt` returns this root, so for a
square it costs one root computation, about half the cost of
`fmpz_is_square` followed by `fmpz_sqrt`.

## n-th roots and perfect powers {#nth}

`fmpz_root` uses `n_sqrtrem`, `n_cbrtrem` and `n_rootrem` for inline
values and GMP's `mpz_root` (Newton iteration) otherwise. The root of a
negative number is negative for odd $n$ and an error for even $n$.

`fmpz_is_perfect_power` returns an exponent $k \ge 2$ with $f = r^k$ if
there is one. For one-limb values it uses `n_is_perfect_power`
repeatedly, which finds the largest exponent. For larger values it
follows the algorithm of GMP's `mpz_perfect_power_p`: the multiplicities
of 2 and of small primes in $f$ restrict the possible exponents to the
divisors of their GCD, and the remaining candidates are tested with
`mpz_root`. The exponent found is then not necessarily the largest:
$2^{120}$ is reported as the square of $2^{60}$. Apply the function again to the root to obtain a larger
exponent.

## Example {#example}

Fifty digits of $\sqrt 2$, checked square roots, a cube root of a
negative number and a perfect power:

```c example name=fmpz_roots
#include <flint/flint.h>
#include <flint/fmpz.h>

int main(void)
{
    fmpz_t a, s, r;
    int k;

    fmpz_init(a); fmpz_init(s); fmpz_init(r);

    /* 50 digits of sqrt(2): floor(sqrt(2 * 10^98)) */
    fmpz_ui_pow_ui(a, 10, 98);
    fmpz_mul_ui(a, a, 2);
    fmpz_sqrtrem(s, r, a);
    flint_printf("isqrt(2*10^98) = "); fmpz_print(s); flint_printf("\n");
    flint_printf("remainder has %wu bits\n", fmpz_bits(r));

    /* checked square roots */
    fmpz_mul(a, s, s);
    flint_printf("perfect_sqrt(s^2) = %d, ", fmpz_perfect_sqrt(r, a));
    flint_printf("equal to s: %d\n", fmpz_equal(r, s));
    fmpz_add_ui(a, a, 1);
    flint_printf("is_square(s^2 + 1) = %d\n", fmpz_is_square(a));

    /* n-th roots round towards zero */
    fmpz_set_si(a, -1000001);
    k = fmpz_root(r, a, 3);
    flint_printf("root(-1000001, 3) = "); fmpz_print(r);
    flint_printf(", exact: %d\n", k);

    /* perfect powers */
    fmpz_ui_pow_ui(a, 7, 45);
    k = fmpz_is_perfect_power(r, a);
    flint_printf("7^45 = r^%d with r = ", k); fmpz_print(r); flint_printf("\n");

    fmpz_clear(a); fmpz_clear(s); fmpz_clear(r);
    return 0;
}
```
```output
isqrt(2*10^98) = 14142135623730950488016887242096980785696718753769
remainder has 164 bits
perfect_sqrt(s^2) = 1, equal to s: 1
is_square(s^2 + 1) = 0
root(-1000001, 3) = -100, exact: 0
7^45 = r^45 with r = 7
```

The remainder $2 \cdot 10^{98} - s^2$ is at most $2s$; here both $s$ and
the remainder have 164 bits.
