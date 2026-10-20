---
title: Canonical form and the fmpq interface
summary: How an fmpq_t stores a fraction, the canonical form that all functions rely on, how arithmetic keeps fractions canonical without large GCDs, and how to use the underscore functions on pairs of fmpz.
order: 1
---

An `fmpq_t` holds a rational number as a numerator and a denominator, two
`fmpz` integers.

## The canonical form {#canonical}

A fraction $p/q$ is *canonical* if $q > 0$ and $\gcd(p, q) = 1$, with zero
written as $0/1$. Every rational number has exactly one canonical
representation, so that equality is a comparison of components
(`fmpq_equal`), the sign is the sign of the numerator (`fmpq_sgn`) and the
height $\max(|p|, q)$ measures the size of the number (`fmpq_height`).

All functions of `fmpq` assume canonical inputs and produce canonical
outputs, unless their documentation says otherwise. The assumption is not
checked. A non-canonical input does not crash, but the results are not
specified: `fmpq_equal` can report $1/2 \ne 2/4$, and `fmpq_sgn` of
$1/(-2)$ is positive.

The components can be accessed through `fmpq_numref` and `fmpq_denref`,
which give `fmpz *` pointers. Writing to them is allowed, after which the
caller is responsible for restoring the invariant, typically by calling
`fmpq_canonicalise`:

```
fmpz_set_ui(fmpq_numref(x), 6);
fmpz_set_si(fmpq_denref(x), -4);
fmpq_canonicalise(x);              /* x = -3/2 */
```

Functions that build a fraction from two integers, such as
`fmpq_set_si`, `fmpq_set_ui` and `fmpq_set_fmpz_frac`, canonicalise. The
exceptions are the conversions that copy components verbatim:
`fmpq_set`, `fmpq_set_mpq`, `fmpq_get_mpq` and `fmpq_set_str` (which, like
`mpq_set_str`, accepts `"2/4"` as it is).

## Arithmetic on canonical fractions {#arithmetic}

Keeping fractions in lowest terms costs GCD computations. Arithmetic on
canonical inputs can arrange them so that the GCDs involve the smaller
numbers, a method due to Henrici [@Hen1956]; see also
[@Knu1997, §4.5.1].

::: proposition Addition {#prop-add}
Let $p/q$ and $r/s$ be canonical, $g = \gcd(q, s)$, $q = ga$ and $s = gb$.
Put $t = pb + ra$ and $h = \gcd(t, g)$. Then
$$
\frac{p}{q} + \frac{r}{s} = \frac{t/h}{ab \cdot g/h}
$$
and the right-hand side is canonical.
:::

::: proof
The identity is clear since $qs/g = gab$. As $\gcd(p, q) = 1$, $p$ is
coprime to $a$, and $b$ is coprime to $a$ by the choice of $g$, so
$\gcd(t, a) = \gcd(pb, a) = 1$; similarly $\gcd(t, b) = 1$. Hence $t/h$ is
coprime to $ab$, and to $g/h$ by the definition of $h$.
:::

So `fmpq_add` computes $\gcd(q, s)$, of the size of the denominators, and
$\gcd(t, g)$, where $g$ is often small, instead of a GCD of the full
numerator and denominator of the result. When $g = 1$ the sum $(ps +
qr)/(qs)$ is canonical without any further GCD, and an integer
summand needs none at all: $(p + qr)/q$ is canonical.

::: proposition Multiplication {#prop-mul}
Let $p/q$ and $r/s$ be canonical, $g_1 = \gcd(p, s)$ and $g_2 = \gcd(q,
r)$. Then
$$
\frac{p}{q} \cdot \frac{r}{s} = \frac{(p/g_1)(r/g_2)}{(q/g_2)(s/g_1)}
$$
and the right-hand side is canonical.
:::

::: proof
Any prime dividing both the numerator and the denominator on the right
divides one of $p/g_1$, $r/g_2$ and one of $q/g_2$, $s/g_1$. It cannot
divide both $p$ and $q$, or $r$ and $s$, which are coprime, and it cannot
divide both $p/g_1$ and $s/g_1$, or $r/g_2$ and $q/g_2$, by the
definition of $g_1$ and $g_2$.
:::

The two GCDs are between a numerator and a denominator of the inputs,
again smaller than the product. Squaring, equal denominators and integer
factors are special cases that need fewer or no GCDs.

When all four components are small `fmpz` values (at most 62 bits in
absolute value on a 64-bit machine), `fmpq_add` and `fmpq_mul` call
`_fmpq_add_small` and `_fmpq_mul_small`, which compute in single and
double words with `n_gcd` and allocate no memory unless the result needs
more than 62 bits.

## The underscore interface {#underscore}

Most functions have a version with a leading underscore that takes the
numerator and denominator as separate `fmpz_t` arguments. These are
useful when rationals are stored in other ways, for example as a vector of
numerators with a common denominator, or the content of a polynomial.
The underscore functions still expect canonical inputs (with a few
documented exceptions, such as `_fmpq_gcd` and the comparison functions
`_fmpq_cmp`, which only need positive denominators). Outputs may be
aliased with inputs, but a numerator must never be aliased with a
denominator.

## Example {#example}

The program computes $\sum_{k=1}^{10} 1/(k(k+1)) = 10/11$ with
`fmpq_add`, builds a fraction from components and canonicalises it, and
adds two fractions given as pairs of `fmpz` with `_fmpq_add`.

```c example name=fmpq_canonical_basic
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>

int main(void)
{
    fmpq_t s, t;
    fmpz_t a, b, c, d, rn, rd;
    ulong k;

    fmpq_init(s);
    fmpq_init(t);

    /* telescoping sum: 1/(1*2) + 1/(2*3) + ... + 1/(10*11) */
    for (k = 1; k <= 10; k++)
    {
        fmpq_set_ui(t, 1, k * (k + 1));
        fmpq_add(s, s, t);
    }
    flint_printf("sum = ");
    fmpq_print(s);
    flint_printf("\n");

    /* writing the components directly */
    fmpz_set_ui(fmpq_numref(t), 6);
    fmpz_set_si(fmpq_denref(t), -4);
    flint_printf("canonical before: %d\n", fmpq_is_canonical(t));
    fmpq_canonicalise(t);
    flint_printf("after: ");
    fmpq_print(t);
    flint_printf("\n");

    /* 5/6 + 7/10 on pairs of fmpz: g = 2, t = 5*5 + 7*3 = 46, h = 2 */
    fmpz_init_set_ui(a, 5);
    fmpz_init_set_ui(b, 6);
    fmpz_init_set_ui(c, 7);
    fmpz_init_set_ui(d, 10);
    fmpz_init(rn);
    fmpz_init(rd);
    _fmpq_add(rn, rd, a, b, c, d);
    _fmpq_print(rn, rd);
    flint_printf("\n");

    fmpz_clear(a); fmpz_clear(b); fmpz_clear(c); fmpz_clear(d);
    fmpz_clear(rn); fmpz_clear(rd);
    fmpq_clear(s);
    fmpq_clear(t);
    return 0;
}
```
```output
sum = 10/11
canonical before: 0
after: -3/2
23/15
```
