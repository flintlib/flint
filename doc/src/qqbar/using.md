---
title: Using qqbar
summary: How to compute with algebraic numbers in canonical form: creating them from polynomials, radicals and roots of unity, arithmetic and its cost, printing and closed forms, predicates, and the qqbar contexts of the generic rings interface.
order: 1
---

The module [](qqbar) represents an algebraic number by its minimal
polynomial over
$\Z$ and a ball that isolates it from the other roots (see
[](algorithms/exact/qqbar-representation)). Unlike [](ca), it needs
no context, its representation is canonical, and every predicate is
decided. The price is that operations compute minimal polynomials,
whose degrees can grow quickly.

## Creating algebraic numbers {#creating}

A `qqbar_t` is initialised with `qqbar_init` and cleared with
`qqbar_clear`. Integers and rational numbers are set with
`qqbar_set_si`, `qqbar_set_fmpq` and similar functions, and
`qqbar_i` gives the imaginary unit. Irrational numbers are obtained:

- as roots of polynomials with `qqbar_roots_fmpz_poly`, which returns all
  roots in a canonical order (see [](qqbar/roots));
- as radicals with `qqbar_sqrt`, `qqbar_root_ui` and `qqbar_pow`, which
  take principal branches;
- as values of trigonometric functions and roots of unity at rational
  multiples of $\pi$: `qqbar_cos_pi`, `qqbar_exp_pi_i`,
  `qqbar_root_of_unity`, ...;
- from a symbolic expression with `qqbar_set_fexpr`, or as eigenvalues
  of an integer matrix with `qqbar_eigenvalues_fmpz_mat`.

## Arithmetic and degrees {#arithmetic}

`qqbar_add`, `qqbar_mul`, `qqbar_div` and the other arithmetic
functions return the exact result, with its minimal polynomial. If $x$
and $y$ have degrees $m$ and $n$, the result has degree at most $mn$,
and in general an operation costs a resultant computation and a
factorisation of a polynomial of degree $mn$ (see
[](algorithms/exact/qqbar-arithmetic)). Long computations therefore
become expensive: the sum of the square roots of the first $k$ primes
has degree $2^k$. Many structured cases are recognised and cheaper, for
example operations on quadratic numbers and on numbers that lie in a
common small number field or cyclotomic field. Functions such as
`qqbar_binop_within_limits` estimate the size of a result before
computing it.

## Printing and closed forms {#printing}

`qqbar_print` prints the minimal polynomial (as its coefficient list,
constant term first) and the enclosure; `qqbar_printnd` and
`qqbar_get_str_nd` print a decimal approximation and the degree, and
`qqbar_get_acb` gives an enclosure of any precision.
`qqbar_get_fexpr_formula` tries to express a number with radicals and
trigonometric functions and returns 0 if it does not find a formula (see
[](qqbar/formulas)); `qqbar_get_fexpr_root_nearest` and
`qqbar_get_fexpr_repr` always succeed, describing the number by its
minimal polynomial and an approximation or its enclosure. These
functions are declared only when `fexpr.h` is included before `qqbar.h`.

```c example name=qqbar_using
#include <flint/fexpr.h>
#include <flint/qqbar.h>

int main(void)
{
    qqbar_t x, y, z;
    fexpr_t e;
    ulong p[] = {2, 3, 5, 7, 11};
    int i;

    qqbar_init(x);
    qqbar_init(y);
    qqbar_init(z);
    fexpr_init(e);

    /* x = 2^(1/3) + i */
    qqbar_set_ui(x, 2);
    qqbar_root_ui(x, x, 3);
    qqbar_i(y);
    qqbar_add(x, x, y);
    qqbar_print(x); flint_printf("\n");
    qqbar_printnd(x, 20); flint_printf("\n");

    /* the minimal polynomial determines the number up to conjugates */
    flint_printf("degree %wd, real: %d\n", qqbar_degree(x), qqbar_is_real(x));
    qqbar_pow_ui(z, x, 3);
    qqbar_print(z); flint_printf("\n");

    /* cos(pi/5) is (1 + sqrt(5))/4 */
    qqbar_cos_pi(y, 1, 5);
    qqbar_get_fexpr_formula(e, y, QQBAR_FORMULA_ALL);
    fexpr_print(e); flint_printf("\n");

    /* degrees multiply: sqrt(2) + sqrt(3) + ... + sqrt(11) */
    qqbar_zero(z);
    for (i = 0; i < 5; i++)
    {
        qqbar_sqrt_ui(y, p[i]);
        qqbar_add(z, z, y);
        flint_printf("%d terms: degree %wd\n", i + 1, qqbar_degree(z));
    }

    fexpr_clear(e);
    qqbar_clear(x);
    qqbar_clear(y);
    qqbar_clear(z);
    return 0;
}
```

```output
deg 6 [5, 12, 3, -4, 3, 0, 1] [1.2599210498948731647672106072782283506 +/- 4.95e-38] + [1.0000000000000000000000000000000000000 +/- 1.79e-38]*I
1.2599210498948731648 + 1.0000000000000000000*I (deg 6)
degree 6, real: 0
deg 6 [125, 888, 1935, 140, 63, -12, 1] [-1.77976314968461949430163182183468505 +/- 2.39e-36] + [3.762203155904598424255116917816924781 +/- 9.47e-37]*I
Div(Add(1, Sqrt(5)), 4)
1 terms: degree 2
2 terms: degree 4
3 terms: degree 8
4 terms: degree 16
5 terms: degree 32
```

## Predicates {#predicates}

Equality (`qqbar_equal`), comparisons of real and imaginary parts and
absolute values (`qqbar_cmp_re`, `qqbar_cmpabs`, ...), signs and tests
such as `qqbar_is_real` and `qqbar_is_rational` return ordinary
integers: they are always decided. Most of them only compare minimal
polynomials and enclosures; when two enclosures cannot be separated, an
exact computation decides (see
[](algorithms/exact/qqbar-representation#predicates)).

## Errors {#errors}

Division by zero and other operations whose result is not defined
throw a FLINT exception (see [](flint)). Functions whose result may fail
to be an algebraic number or may be too expensive to compute return an
`int` success flag instead: `qqbar_pow` returns 0 when $x^y$ is not
algebraic (for irrational $y$, by the Gelfond–Schneider theorem), the
trigonometric functions such as `qqbar_tan_pi` return 0 at poles, and
the functions with explicit degree and bit limits, such as
`qqbar_evaluate_fmpz_mpoly`, return 0 when a limit would be exceeded.

## The generic rings interface {#gr}

The contexts `gr_ctx_init_real_qqbar` and `gr_ctx_init_complex_qqbar`
make the real algebraic numbers $\overline{\mathbb{Q}} \cap \mathbb{R}$
and the algebraic numbers $\overline{\mathbb{Q}}$ available to the
generic algorithms of [](gr) (see [](gr_domains#algebraic)).
`_gr_ctx_qqbar_set_limits` sets limits on degrees and bit sizes, beyond
which operations return `GR_UNABLE` instead of attempting expensive
computations. With these contexts, functions such as `gr_poly_roots`
and `gr_mat_eigenvalues` return exact algebraic roots and eigenvalues.
For long computations in a fixed number field, [](ca) or [](nf_elem)
are usually faster; numbers can be converted to and from `ca_t` with
`ca_set_qqbar` and `ca_get_qqbar`.
