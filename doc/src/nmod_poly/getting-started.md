---
title: Getting started with nmod_poly
summary: How to create, fill, print and free polynomials over Z/nZ, the rules for normalisation and aliasing, how the underscore functions work, and what changes when the modulus is not prime.
order: 1
---

An `nmod_poly_t` holds a polynomial over $\Z/n\Z$ for a modulus $1 \le n
< 2^{\bits}$. It owns its coefficient array and its modulus, grows as
needed, and must be initialised before use and cleared afterwards.

## A first program {#first}

The program below works modulo 7. It builds $f = 5x^3 + 6$ coefficient by
coefficient, squares it, divides the square by $x^2 + 1$ and evaluates
$f$ at 2.

```c example name=nmod_poly_first
#include <flint/flint.h>
#include <flint/nmod_poly.h>

int main(void)
{
    nmod_poly_t f, g, h, q, r;

    nmod_poly_init(f, 7);
    nmod_poly_init(g, 7);
    nmod_poly_init_mod(h, f->mod);   /* the same modulus, no recomputation */
    nmod_poly_init_mod(q, f->mod);
    nmod_poly_init_mod(r, f->mod);

    nmod_poly_set_coeff_ui(f, 3, 5);
    nmod_poly_set_coeff_ui(f, 0, 13);   /* reduced: 13 = 6 mod 7 */

    nmod_poly_mul(g, f, f);

    nmod_poly_set_coeff_ui(h, 2, 1);
    nmod_poly_set_coeff_ui(h, 0, 1);
    nmod_poly_divrem(q, r, g, h);

    flint_printf("f     = "); nmod_poly_print_pretty(f, "x"); flint_printf("\n");
    flint_printf("f^2   = "); nmod_poly_print_pretty(g, "x"); flint_printf("\n");
    flint_printf("plain = "); nmod_poly_print(g); flint_printf("\n");
    flint_printf("q     = "); nmod_poly_print_pretty(q, "x"); flint_printf("\n");
    flint_printf("r     = "); nmod_poly_print_pretty(r, "x"); flint_printf("\n");
    flint_printf("f(2)  = %wu\n", nmod_poly_evaluate_nmod(f, 2));
    flint_printf("deg f = %wd, length = %wd\n",
                 nmod_poly_degree(f), nmod_poly_length(f));

    nmod_poly_clear(f);
    nmod_poly_clear(g);
    nmod_poly_clear(h);
    nmod_poly_clear(q);
    nmod_poly_clear(r);
    return 0;
}
```
```output
f     = 5*x^3+6
f^2   = 4*x^6+4*x^3+1
plain = 7 7  1 0 0 4 0 0 4
q     = 4*x^4+3*x^2+4*x+4
r     = 3*x+4
f(2)  = 4
deg f = 3, length = 4
```

The plain format printed by `nmod_poly_print` is the length, the modulus
and the coefficients from the constant term up. `nmod_poly_init_mod`
copies the `nmod_t` of an existing polynomial; `nmod_poly_init` computes
its precomputed inverse, which costs about one division.

## Storage and normalisation {#storage}

The coefficients are stored in `f->coeffs[0]`, ..., `f->coeffs[f->length
- 1]`, with $x^i$ at index $i$, in an array with room for `f->alloc`
entries. Every function expects and produces polynomials that are

- *normalised*: the leading coefficient `coeffs[length - 1]` is nonzero,
  or the length is 0 (the zero polynomial, of degree $-1$);
- *reduced*: every coefficient lies in $[0, n)$.

Neither condition is checked. `nmod_poly_set_coeff_ui` maintains both, so
it is the safe way to build a polynomial. To fill many coefficients
quickly, write them directly and normalise at the end:

```c
nmod_poly_fit_length(f, len);           /* room for len coefficients */
for (i = 0; i < len; i++)
    f->coeffs[i] = nmod_set_ui(data[i], f->mod);
_nmod_poly_set_length(f, len);
_nmod_poly_normalise(f);                /* strip leading zeros */
```

Forgetting `_nmod_poly_normalise` leaves a polynomial with a zero leading
coefficient, which makes `nmod_poly_degree` wrong and can break division
and GCD functions, which divide by the leading coefficient.
`nmod_poly_fit_length` only grows the allocation, by at least a factor of
two, so filling a polynomial term by term costs amortised linear time.

## Aliasing {#aliasing}

Functions without a leading underscore allow any output to be the same
object as any input: `nmod_poly_mul(f, f, f)` squares $f$ in place, and
`nmod_poly_divrem(q, f, f, g)` replaces $f$ by its remainder. They
allocate a temporary when the algorithm cannot work in place. Several
*outputs* of one function, such as the quotient and remainder, or the
GCD and the cofactors of `nmod_poly_xgcd`, must be distinct objects.

All polynomials passed to one function must have the same modulus. The
modulus is not checked; the one used is usually that of the first input
or of the output.

## The underscore functions {#underscore}

Most functions come in two versions. `nmod_poly_mul(res, a, b)` handles
zero polynomials, allocates memory, deals with aliasing and normalises the
result. `_nmod_poly_mul(res, a, alen, b, blen, mod)` does only the
arithmetic, on plain arrays: it requires $\mathit{alen} \ge \mathit{blen}
\ge 1$, writes exactly $\mathit{alen} + \mathit{blen} - 1$ coefficients
(leading zeros included) to an array that the caller has allocated, and
does not allow `res` to overlap `a` or `b`. The conditions differ from
function to function and are stated in each description; typical ones are
nonzero lengths, an ordering of the lengths, zero padding to a fixed
length (for example, residues modulo $f$ padded to length
$\deg f$), and no aliasing.

The underscore functions are useful inside loops, where the checks and
allocations of the wrappers would dominate, and for working with
coefficient arrays that are not `nmod_poly_t` objects, such as rows of a
matrix or parts of another polynomial. The example multiplies two
polynomials stored in plain arrays and computes the middle coefficients
$3, \ldots, 6$ of the same product with `_nmod_poly_mulmid`:

```c example name=nmod_poly_underscore
#include <flint/flint.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_poly.h>

int main(void)
{
    nmod_t mod;
    ulong a[5] = { 1, 2, 3, 4, 5 };     /* 1 + 2x + 3x^2 + 4x^3 + 5x^4 */
    ulong b[3] = { 6, 0, 1 };           /* 6 + x^2 */
    ulong c[7], m[4];
    slong i;

    nmod_init(&mod, 11);

    _nmod_poly_mul(c, a, 5, b, 3, mod);           /* requires 5 >= 3 >= 1 */
    _nmod_poly_mulmid(m, a, 5, b, 3, 3, 7, mod);  /* coefficients 3..6 */

    flint_printf("product:");
    for (i = 0; i < 7; i++)
        flint_printf(" %wu", c[i]);
    flint_printf("\nmiddle: ");
    for (i = 0; i < 4; i++)
        flint_printf(" %wu", m[i]);
    flint_printf("\n");
    return 0;
}
```
```output
product: 6 1 8 4 0 4 5
middle:  4 0 4 5
```

## Prime and composite moduli {#composite}

Ring operations work for every modulus: addition, multiplication,
powering, evaluation, composition, Taylor shifts by Horner's rule, and
division by a polynomial whose leading coefficient is a unit, for example
a monic one. Many other functions invert elements along the way and
assume that $n$ is prime: GCDs and resultants, inverses modulo a
polynomial, square roots, interpolation, power series functions such as
`nmod_poly_exp_series` (which also need $n$ larger than the length) and
the root-finding and factoring functions. The assumption is not checked.
When an inverse does not exist, these functions usually throw an
exception (`FLINT_IMPINV` for an impossible inverse), which by default
aborts the program, but some can return a wrong result without warning.

The next program works modulo 12. Division by the monic $x^2 + 5$ is
exact arithmetic in $(\Z/12\Z)[x]$ and works, and so does the
multiplication modulo this polynomial; a GCD modulo 12 should not be
attempted.

```c example name=nmod_poly_composite
#include <flint/flint.h>
#include <flint/nmod_poly.h>

int main(void)
{
    nmod_poly_t a, b, q, r;

    nmod_poly_init(a, 12);
    nmod_poly_init(b, 12);
    nmod_poly_init(q, 12);
    nmod_poly_init(r, 12);

    nmod_poly_set_str(a, "5 12  7 0 3 9 4");   /* 4x^4 + 9x^3 + 3x^2 + 7 */
    nmod_poly_set_str(b, "3 12  5 0 1");       /* x^2 + 5 */

    nmod_poly_divrem(q, r, a, b);
    flint_printf("q = "); nmod_poly_print_pretty(q, "x"); flint_printf("\n");
    flint_printf("r = "); nmod_poly_print_pretty(r, "x"); flint_printf("\n");

    nmod_poly_mulmod(q, a, a, b);
    flint_printf("a^2 mod b = "); nmod_poly_print_pretty(q, "x"); flint_printf("\n");

    /* 6 * (2x + 2) = 0 modulo 12: the product is normalised to zero */
    nmod_poly_set_str(r, "2 12  2 2");
    nmod_poly_scalar_mul_nmod(q, r, 6);
    flint_printf("6 (2x + 2) = "); nmod_poly_print_pretty(q, "x");
    flint_printf(", length %wd\n", nmod_poly_length(q));

    nmod_poly_clear(a);
    nmod_poly_clear(b);
    nmod_poly_clear(q);
    nmod_poly_clear(r);
    return 0;
}
```
```output
q = 4*x^2+9*x+7
r = 3*x+8
a^2 mod b = 7
6 (2x + 2) = 0, length 0
```

Over a composite modulus the degree of a product can be smaller than the
sum of the degrees, since the leading coefficients can multiply to zero;
the functions normalise their results accordingly.

## Where to go next {#next}

- [](nmod_poly/mul) explains how products are computed and how to compute
  only part of a product.
- [](nmod_poly/divrem) and [](nmod_poly/inv_series) cover division with
  remainder and power series division.
- [](nmod_poly/modular-arithmetic) shows how to compute efficiently modulo
  a fixed polynomial, for example in a finite field $\mathbb{F}_p[x]/(f)$.
- [](nmod_poly/gcd), [](nmod_poly/resultant), [](nmod_poly/multipoint) and
  [](nmod_poly/compose_mod) cover the main advanced operations.
- The tutorial [](tutorials/algebra/reed-solomon) uses many of these functions in a
  complete program.
