---
title: Symbolic expressions and closed forms
functions: qqbar_get_fexpr_formula qqbar_set_fexpr qqbar_get_fexpr_repr qqbar_get_fexpr_root_nearest qqbar_get_quadratic QQBAR_FORMULA_GAUSSIANS QQBAR_FORMULA_EXP_FORM
---

A `qqbar_t` can be converted to a symbolic expression (see [](fexpr)) in
three ways, and converted back with
`qqbar_set_fexpr`.

- `qqbar_get_fexpr_repr` writes the internal representation exactly,
  as `AlgebraicNumberSerialized(List(c_0, ..., c_d), enclosure)`. Reading
  it back costs no root isolation, so this is the format for storage and
  for transferring numbers between programs or contexts (`ca_transfer`
  uses it).
- `qqbar_get_fexpr_root_nearest` and `qqbar_get_fexpr_root_indexed`
  describe the number as a root of its minimal polynomial, selected by a
  short decimal approximation or by its index in the canonical order.
  These are readable and unambiguous; reading them back requires
  isolating the roots.
- `qqbar_get_fexpr_formula` tries to find a closed form with radicals,
  trigonometric functions or exponentials. It succeeds for all numbers
  of degree at most 2, for elements of cyclotomic fields that it
  recognises, and for some numbers built from these by roots and by
  combining real and imaginary parts; otherwise it fails and returns 0.

## Finding closed forms {#formulas}

Every algebraic number of degree at most 4 can be written with
radicals, and every element of an abelian extension of $\mathbb{Q}$ is,
by the Kronecker–Weber theorem, an element of a cyclotomic field
$\mathbb{Q}(\zeta_n)$, which gives closed forms with $e^{2 \pi i/n}$ or
with $\cos$ and $\sin$ of rational multiples of $\pi$. Most algebraic
numbers of degree 5 and higher have no expression by radicals at all.
`qqbar_get_fexpr_formula` uses the following methods, selected by the
`QQBAR_FORMULA_...` flags.

- Gaussian rationals $a + bi$ and quadratic numbers $(a + b\sqrt{c})/q$
  are written directly, using `qqbar_get_quadratic` with smooth
  factorisation to make $c$ small.
- Roots of unity are recognised from their cyclotomic minimal
  polynomials. For other numbers whose degree divides $\varphi(n)$ for a
  plausible $n$, an expression as a polynomial in $\zeta_n$ (or in $2
  \cos(2\pi/n)$ for real numbers) is sought with LLL and verified
  exactly; the result is written in exponential or trigonometric form.
- Deflation: if the minimal polynomial is a polynomial in $x^k$, a closed
  form for $x^k$ is sought recursively and the $k$-th root taken, with
  the correct branch.
- Separation: as a last resort, a formula is sought for the real and
  imaginary parts, or for the modulus and the sign, separately.

The cubic, quartic and quintic formulas and depression are not yet
implemented, so a cubic irrationality that is not cyclotomic and not a
pure root (such as a root of $x^3 - x - 1$) gets no formula.

## Example {#example}

The program prints closed forms for $\cos(\pi/5)$, $\cos(2\pi/7)$
(found in a cyclotomic field and written with trigonometric functions)
and $\sqrt[3]{2} + i$ (by separation into real and imaginary parts). A
root of $x^5 - x - 1$ has no closed form and is written as the root of
its minimal polynomial nearest to a decimal approximation. Finally the
exact representation is written and read back.

```c example name=qqbar_formulas
#include <flint/fmpz_poly.h>
#include <flint/fexpr.h>
#include <flint/qqbar.h>

static void show(const qqbar_t x)
{
    fexpr_t e;
    fexpr_init(e);
    if (!qqbar_get_fexpr_formula(e, x, QQBAR_FORMULA_ALL))
        qqbar_get_fexpr_root_nearest(e, x);
    fexpr_print(e);
    flint_printf("\n");
    fexpr_clear(e);
}

int main(void)
{
    qqbar_t x, y;
    qqbar_ptr r;
    fmpz_poly_t f;
    fexpr_t e;
    char * s;

    qqbar_init(x);
    qqbar_init(y);
    fexpr_init(e);
    fmpz_poly_init(f);
    r = _qqbar_vec_init(5);

    qqbar_cos_pi(x, 1, 5);            /* cos(pi/5) = (1 + sqrt(5))/4 */
    show(x);

    qqbar_cos_pi(x, 2, 7);            /* cos(2 pi/7), degree 3 */
    show(x);

    qqbar_set_ui(x, 2);               /* 2^(1/3) + i, degree 6 */
    qqbar_root_ui(x, x, 3);
    qqbar_i(y);
    qqbar_add(x, x, y);
    show(x);

    fmpz_poly_set_str(f, "6  -1 -1 0 0 0 1");   /* x^5 - x - 1 */
    qqbar_roots_fmpz_poly(r, f, 0);
    show(r + 0);

    qqbar_get_fexpr_repr(e, r + 0);   /* exact serialization */
    s = fexpr_get_str(e);
    flint_printf("%.48s...\n", s);
    flint_free(s);
    qqbar_set_fexpr(y, e);
    flint_printf("read back: %d\n", qqbar_equal(r + 0, y));

    _qqbar_vec_clear(r, 5);
    fmpz_poly_clear(f);
    qqbar_clear(x);
    qqbar_clear(y);
    fexpr_clear(e);
    return 0;
}
```
```output
Div(Add(1, Sqrt(5)), 4)
Div(Add(Mul(2, Cos(Div(Pi, 7))), Mul(2, Sin(Div(Pi, 14))), -1), 2)
Add(Pow(2, Div(1, 3)), Mul(1, NumberI))
PolynomialRootNearest(List(-1, -1, 0, 0, 0, 1), Decimal("1.16730"))
AlgebraicNumberSerialized(List(-1, -1, 0, 0, 0, ...
read back: 1
```

The formula for $\cos(2\pi/7)$ is correct but not the simplest one: the
LLL search finds a short relation in the real cyclotomic field, which is
then written term by term. Formulas are meant for presentation; they
are not canonical, and two equal numbers can get different formulas.
Use `qqbar_get_fexpr_repr` when the expression will be read by a
program.
