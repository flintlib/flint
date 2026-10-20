---
title: Getting started with fmpz_poly
summary: How to create, fill, print and free polynomials over the integers, the rules for normalisation and aliasing, how the underscore functions work, how division behaves over Z, and how to move between Z[x], (Z/nZ)[x] and Q[x].
order: 1
---

An `fmpz_poly_t` holds a polynomial with integer coefficients of any
size. It owns its array of `fmpz` coefficients, grows as needed, and
must be initialised before use and cleared afterwards.

## A first program {#first}

The program builds $f = 5x^3 - 1$ coefficient by coefficient, squares
it, divides the square by $x^2 + 1$, factors $f^2 - 1$ and evaluates
$f$ at 2.

```c example name=fmpz_poly_first
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_factor.h>

int main(void)
{
    fmpz_poly_t f, g, h, q, r;
    fmpz_poly_factor_t fac;
    fmpz_t x, y;

    fmpz_poly_init(f);
    fmpz_poly_init(g);
    fmpz_poly_init(h);
    fmpz_poly_init(q);
    fmpz_poly_init(r);
    fmpz_poly_factor_init(fac);
    fmpz_init(x);
    fmpz_init(y);

    fmpz_poly_set_coeff_ui(f, 3, 5);
    fmpz_poly_set_coeff_si(f, 0, -1);

    fmpz_poly_mul(g, f, f);
    fmpz_poly_set_str(h, "3  1 0 1");          /* x^2 + 1 */
    fmpz_poly_divrem(q, r, g, h);

    flint_printf("f     = "); fmpz_poly_print_pretty(f, "x"); flint_printf("\n");
    flint_printf("f^2   = "); fmpz_poly_print_pretty(g, "x"); flint_printf("\n");
    flint_printf("plain = "); fmpz_poly_print(g); flint_printf("\n");
    flint_printf("q     = "); fmpz_poly_print_pretty(q, "x"); flint_printf("\n");
    flint_printf("r     = "); fmpz_poly_print_pretty(r, "x"); flint_printf("\n");

    fmpz_poly_sub_si(g, g, 1);                 /* f^2 - 1 = (f - 1)(f + 1) */
    fmpz_poly_factor(fac, g);
    flint_printf("f^2 - 1 = "); fmpz_poly_factor_print(fac);

    fmpz_set_ui(x, 2);
    fmpz_poly_evaluate_fmpz(y, f, x);
    flint_printf("f(2)  = "); fmpz_print(y); flint_printf("\n");
    flint_printf("deg f = %wd, length = %wd\n",
                 fmpz_poly_degree(f), fmpz_poly_length(f));

    fmpz_poly_clear(f);
    fmpz_poly_clear(g);
    fmpz_poly_clear(h);
    fmpz_poly_clear(q);
    fmpz_poly_clear(r);
    fmpz_poly_factor_clear(fac);
    fmpz_clear(x);
    fmpz_clear(y);
    return 0;
}
```
```output
f     = 5*x^3-1
f^2   = 25*x^6-10*x^3+1
plain = 7  1 0 0 -10 0 0 25
q     = 25*x^4-25*x^2-10*x+25
r     = 10*x-24
f^2 - 1 = 5
2  0 1 ^ 3
4  -2 0 0 5 ^ 1
f(2)  = 39
deg f = 3, length = 4
```

The plain format printed by `fmpz_poly_print` is the length, two
spaces, and the coefficients from the constant term up; it is also
the format read by `fmpz_poly_set_str`. `fmpz_poly_factor_print`
prints the content followed by one line per irreducible factor in
plain format with its multiplicity: here $f^2 - 1 = 25x^6 - 10x^3 = 5
\cdot x^3 \cdot (5x^3 - 2)$.

## Storage and normalisation {#storage}

The coefficients are stored in `f->coeffs[0]`, ...,
`f->coeffs[f->length - 1]`, with $x^i$ at index $i$, in an array with
room for `f->alloc` entries. Each entry is an `fmpz`, which holds small
values inline and larger values in a separately allocated GMP integer,
so every coefficient has its own size. Every function expects and
produces *normalised* polynomials: the leading coefficient
`coeffs[length - 1]` is nonzero, or the length is 0 (the zero
polynomial, of degree $-1$). In addition, the allocated entries beyond
the length are always zero.

`fmpz_poly_set_coeff_ui`, `fmpz_poly_set_coeff_si` and
`fmpz_poly_set_coeff_fmpz` maintain these rules. To fill many
coefficients quickly, allocate, write the entries with the `fmpz`
functions and normalise at the end:

```c
fmpz_poly_fit_length(f, len);           /* room for len coefficients */
for (i = 0; i < len; i++)
    fmpz_set_si(f->coeffs + i, data[i]);
_fmpz_poly_set_length(f, len);
_fmpz_poly_normalise(f);                /* strip leading zeros */
```

Forgetting `_fmpz_poly_normalise` leaves a zero leading coefficient,
which makes `fmpz_poly_degree` wrong and breaks the functions that
divide by the leading coefficient. A coefficient can also be read or
modified in place through `fmpz_poly_get_coeff_ptr`.

## Aliasing {#aliasing}

Functions without a leading underscore allow any output to be the same
object as any input: `fmpz_poly_mul(f, f, f)` squares $f$ in place and
`fmpz_poly_gcd(f, f, g)` replaces $f$ by a GCD. Several *outputs* of
one function, such as the quotient and remainder of
`fmpz_poly_divrem`, must be distinct objects. The few exceptions are
stated in the descriptions; for example, the Hensel lifting functions
allow only some aliasing.

## The underscore functions {#underscore}

Most functions come in two versions. `fmpz_poly_mul(res, a, b)`
handles zero polynomials, allocates memory, deals with aliasing and
normalises. `_fmpz_poly_mul(res, a, alen, b, blen)` only does the
arithmetic on arrays: it requires positive lengths, writes exactly
$\mathit{alen} + \mathit{blen} - 1$ coefficients to an array of
initialised `fmpz` that the caller provides, and does not allow `res`
to overlap the inputs. The arrays are usually managed with
`_fmpz_vec_init` and `_fmpz_vec_clear`. The conditions vary and are
stated with each function.

The example computes the product of two arrays and its middle
coefficients 2 to 4 with `_fmpz_poly_mulmid`, and the GCD of two
polynomials stored in one array:

```c example name=fmpz_poly_underscore
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_vec.h>
#include <flint/fmpz_poly.h>

int main(void)
{
    fmpz * a, * b, * c;
    slong i;

    a = _fmpz_vec_init(4);
    b = _fmpz_vec_init(3);
    c = _fmpz_vec_init(6);

    for (i = 0; i < 4; i++)
        fmpz_set_si(a + i, i + 1);              /* 1 + 2x + 3x^2 + 4x^3 */
    fmpz_set_si(b + 0, -1);
    fmpz_set_si(b + 2, 1);                      /* x^2 - 1 */

    _fmpz_poly_mul(c, a, 4, b, 3);              /* requires lengths > 0 */
    flint_printf("product:");
    for (i = 0; i < 6; i++)
        flint_printf(" %wd", fmpz_get_si(c + i));

    _fmpz_poly_mulmid(c, a, 4, b, 3, 2, 5);     /* coefficients 2, 3, 4 */
    flint_printf("\nmiddle: ");
    for (i = 0; i < 3; i++)
        flint_printf(" %wd", fmpz_get_si(c + i));

    /* gcd(x^2 - 1, x^2 + 2x + 1) = x + 1, written to (c, 3) */
    fmpz_set_si(a + 0, 1);
    fmpz_set_si(a + 1, 2);
    fmpz_set_si(a + 2, 1);
    _fmpz_poly_gcd(c, a, 3, b, 3);              /* requires len1 >= len2 > 0 */
    flint_printf("\ngcd:    ");
    for (i = 0; i < 3; i++)
        flint_printf(" %wd", fmpz_get_si(c + i));
    flint_printf("\n");

    _fmpz_vec_clear(a, 4);
    _fmpz_vec_clear(b, 3);
    _fmpz_vec_clear(c, 6);
    return 0;
}
```
```output
product: -1 -2 -2 -2 3 4
middle:  -2 -2 3
gcd:     1 1 0
```

The GCD is written with zero padding to the length of the shorter
input.

## Division over the integers {#division}

$\Z$ is not a field, so a polynomial can only be divided by another in
the sense of $\mathbb{Q}[x]$ if the leading coefficient of the divisor
is $\pm 1$ or the division is exact. `fmpz_poly_divrem` always returns
integer polynomials with $A = BQ + R$; when the leading coefficient
$\ell$ of $B$ is not $\pm 1$, the remainder can have degree at least
$\deg B$, with coefficients smaller than $|\ell|$ in those degrees.
`fmpz_poly_divexact` and `fmpz_poly_divides` are the right tools when
$B$ divides $A$, and `fmpz_poly_pseudo_divrem` computes $\ell^d A = BQ
+ R$ with $\deg R < \deg B$. [](fmpz_poly/divrem) explains the details.

The same applies to the other operations of the Euclidean algorithm.
`fmpz_poly_gcd` returns the GCD in $\Z[x]$ with positive leading
coefficient, which includes the GCD of the contents, and
`fmpz_poly_xgcd` returns the resultant $r$ and $s, t$ with $sf + tg =
r$ instead of a Bézout identity for the GCD.

## Other polynomial rings {#conversions}

- `fmpz_poly_get_nmod_poly` reduces a polynomial modulo the modulus of
  an `nmod_poly_t`, and `fmpz_poly_set_nmod_poly` lifts back with
  coefficients in the symmetric range ($-n/2 \le c < n/2$), which is the
  right choice when the integer coefficients are known to be small in
  absolute value. `fmpz_poly_CRT_ui` combines images modulo coprime
  moduli. These are the building blocks of multimodular algorithms.
- For polynomials over $\mathbb{Q}$, `fmpq_poly_set_fmpz_poly` and
  `fmpq_poly_get_numerator` convert; an `fmpq_poly_t` is an
  `fmpz_poly` numerator with a common denominator.
- For multiprecision moduli, use `fmpz_mod_poly_set_fmpz_poly`.

```c example name=fmpz_poly_modular
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/nmod_poly.h>
#include <flint/fmpz_poly.h>

int main(void)
{
    fmpz_poly_t f, g;
    nmod_poly_t a, b;
    fmpz_t m;

    fmpz_poly_init(f);
    fmpz_poly_init(g);
    nmod_poly_init(a, 101);
    nmod_poly_init(b, 103);
    fmpz_init(m);

    fmpz_poly_set_str(f, "4  -1000 3 0 -2500");   /* -2500x^3 + 3x - 1000 */

    fmpz_poly_get_nmod_poly(a, f);
    fmpz_poly_get_nmod_poly(b, f);
    flint_printf("mod 101: "); nmod_poly_print_pretty(a, "x"); flint_printf("\n");

    fmpz_poly_set_nmod_poly(g, a);                 /* symmetric lift */
    flint_printf("lift:    "); fmpz_poly_print_pretty(g, "x"); flint_printf("\n");

    fmpz_set_ui(m, 101);
    fmpz_poly_CRT_ui(g, g, m, b, 1);               /* modulo 101 * 103 */
    flint_printf("CRT:     "); fmpz_poly_print_pretty(g, "x"); flint_printf("\n");

    fmpz_poly_clear(f);
    fmpz_poly_clear(g);
    nmod_poly_clear(a);
    nmod_poly_clear(b);
    fmpz_clear(m);
    return 0;
}
```
```output
mod 101: 25*x^3+3*x+10
lift:    25*x^3+3*x+10
CRT:     -2500*x^3+3*x-1000
```

One prime is not enough to recover coefficients of size 2500; the
product $101 \cdot 103 = 10403$ is, and the symmetric lift returns the
original polynomial.

## Where to go next {#next}

- [](fmpz_poly/mul) explains how products are computed and how to
  compute only part of a product.
- [](fmpz_poly/divrem), [](fmpz_poly/gcd) and [](fmpz_poly/resultant)
  cover division, GCDs and resultants over $\Z$.
- [](fmpz_poly/evaluate), [](fmpz_poly/compose) and [](fmpz_poly/series)
  cover evaluation, composition and power series.
- [](fmpz_poly_factor/factor) explains factorisation into irreducible
  polynomials.
- The algorithms are described in
  [](algorithms/polynomials/multiplication-integers),
  [](algorithms/polynomials/gcd-integers),
  [](algorithms/polynomials/resultant-integers),
  [](algorithms/polynomials/hensel-lifting),
  [](algorithms/polynomials/factoring-integers) and
  [](algorithms/polynomials/real-roots-integers).
