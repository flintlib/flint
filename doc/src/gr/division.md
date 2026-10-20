---
title: Division in generic rings
functions: gr_is_invertible gr_div gr_div_nonunique gr_divexact gr_euclidean_div
---

A generic ring offers several division functions because "$x$ divided
by $y$" has different meanings in different rings. In a field every
nonzero $y$ is invertible and $x / y$ is unique. In an integral domain
such as $\Z$ the quotient is unique when it exists, but it may not exist.
In a ring with zero divisors such as $\Z/10\Z$ the equation $q y = x$
may have several solutions. And in a Euclidean ring there is also
division with remainder. The functions below cover these cases:

| Function | Computes | Fails with `GR_DOMAIN` when |
|---|---|---|
| `gr_inv` | $x^{-1}$ | $x$ is not a unit |
| `gr_div` | the unique $q$ with $q y = x$ | no solution, or several |
| `gr_div_nonunique` | some $q$ with $q y = x$ | no solution |
| `gr_divexact` | $q$ with $q y = x$, assuming one exists | (may not check) |
| `gr_euclidean_div`, `gr_euclidean_rem` | $q$, $r$ with $x = q y + r$ | $y = 0$ |
| `gr_divides` (predicate) | whether $y \mid x$ | |

Each function returns `GR_UNABLE` when the answer cannot be determined:
the ring has no algorithm for it, or the input is an inexact ball that
contains zero. Over rings that pretend to be fields
(`gr_ctx_set_is_pretend_field`), an inversion that meets a zero divisor
also returns `GR_UNABLE` and records the zero divisor.

## Uniqueness {#uniqueness}

`gr_div` returns `GR_DOMAIN` for $0/0$ in every ring except the zero
ring, since every $q$ solves $q \cdot 0 = 0$; in the zero ring $0/0 = 0$
is the unique quotient. In $\Z/10\Z$, both $3$ and $8$ solve
$2 q = 6$, so `gr_div` fails there too, while `gr_div_nonunique` returns
one of the solutions. `gr_divides(y, x)` answers the question that
`gr_div_nonunique` answers: whether some $q$ exists, so $0 \mid 0$ is
true in every ring. (Note the order of the arguments:
`gr_divides(d, x, ctx)` tests whether $d$ divides $x$.)

For ball arithmetic, a quotient by a ball that contains zero but is not
exactly zero is `GR_UNABLE`, since the ball might represent a nonzero
number, while a quotient by an exact zero is `GR_DOMAIN`.

```c example name=gr_division
#include <flint/flint.h>
#include <flint/gr.h>

const char * st(int status)
{
    return status == GR_SUCCESS ? "GR_SUCCESS" : status == GR_DOMAIN ?
        "GR_DOMAIN" : status == GR_UNABLE ? "GR_UNABLE" : "both";
}

const char * tr(truth_t t)
{
    return t == T_TRUE ? "T_TRUE" : t == T_FALSE ? "T_FALSE" : "T_UNKNOWN";
}

void divide(const char * xs, const char * ys, gr_ctx_t ctx)
{
    gr_ptr x, y, q;
    int status;

    GR_TMP_INIT3(x, y, q, ctx);
    GR_MUST_SUCCEED(gr_set_str(x, xs, ctx));
    GR_MUST_SUCCEED(gr_set_str(y, ys, ctx));

    flint_printf("%s / %s in %{gr_ctx}\n", xs, ys, ctx);
    status = gr_div(q, x, y, ctx);
    flint_printf("  div:           %s", st(status));
    if (status == GR_SUCCESS)
        flint_printf(", %{gr}", q, ctx);
    status = gr_div_nonunique(q, x, y, ctx);
    flint_printf("\n  div_nonunique: %s", st(status));
    if (status == GR_SUCCESS)
        flint_printf(", %{gr}", q, ctx);
    flint_printf("\n  divides:       %s\n", tr(gr_divides(y, x, ctx)));

    GR_TMP_CLEAR3(x, y, q, ctx);
}

int main(void)
{
    gr_ctx_t ZZ, ZZ10, RR;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_nmod(ZZ10, 10);
    gr_ctx_init_real_arb(RR, 53);

    divide("7", "3", ZZ);
    divide("0", "0", ZZ);
    divide("6", "2", ZZ10);
    divide("7", "2", ZZ10);
    divide("1", "[0 +/- 0.5]", RR);
    divide("1", "0", RR);

    gr_ctx_clear(ZZ);
    gr_ctx_clear(ZZ10);
    gr_ctx_clear(RR);
    return 0;
}
```
```output
7 / 3 in Integer ring (fmpz)
  div:           GR_DOMAIN
  div_nonunique: GR_DOMAIN
  divides:       T_FALSE
0 / 0 in Integer ring (fmpz)
  div:           GR_DOMAIN
  div_nonunique: GR_SUCCESS, 0
  divides:       T_TRUE
6 / 2 in Integers mod 10 (_gr_nmod)
  div:           GR_DOMAIN
  div_nonunique: GR_SUCCESS, 3
  divides:       T_TRUE
7 / 2 in Integers mod 10 (_gr_nmod)
  div:           GR_DOMAIN
  div_nonunique: GR_DOMAIN
  divides:       T_FALSE
1 / [0 +/- 0.5] in Real numbers (arb, prec = 53)
  div:           GR_UNABLE
  div_nonunique: GR_UNABLE
  divides:       T_UNKNOWN
1 / 0 in Real numbers (arb, prec = 53)
  div:           GR_DOMAIN
  div_nonunique: GR_DOMAIN
  divides:       T_FALSE
```

Which of the solutions `gr_div_nonunique` returns is not specified;
here it is 3 rather than 8.

## Exact division {#divexact}

`gr_divexact` is for the common situation in which the caller knows that
$y$ divides $x$, for example when removing a known content from a
polynomial or in fraction-free elimination. The ring may then use a
faster algorithm and is allowed to skip the check: `fmpz` uses
`fmpz_divexact`, so in $\Z$ the call `gr_divexact(q, 7, 3)` returns
`GR_SUCCESS` with $q = 2$, a meaningless value. Use `gr_div` when the
quotient might not exist. The default implementation simply calls
`gr_div`, and the `_ui`, `_si`, `_fmpz`, `_fmpq` and `_other` variants
default to the corresponding `gr_div` variants.

## Inverses {#inverses}

`gr_is_invertible` and `gr_inv` decide and compute units: in
$\Z/10\Z$, $3^{-1} = 7$ and `gr_inv` of $4$ returns `GR_DOMAIN`. The
default implementations, used by rings that do not provide an inverse,
handle $x = \pm 1$ and $x = 0$ (`GR_DOMAIN`) and return `GR_UNABLE` for
any other element; [](gr/implementing) shows a ring that relies on this
default. `gr_div` has no default implementation (it returns `GR_UNABLE`
in rings that do not implement it), while the `_ui`, `_si` and `_fmpz`
variants convert the divisor to the ring and call `gr_div`, and
`gr_div_fmpq` multiplies by the inverse of the rational number.

Negative powers (`gr_pow_si`, `gr_pow_fmpz`) are computed by default by
inverting and then raising to the positive power, so they fail exactly
when $x$ is not invertible.

## Division with remainder {#euclidean}

`gr_euclidean_div`, `gr_euclidean_rem` and `gr_euclidean_divrem` compute
$q$ and $r$ with $x = q y + r$ in rings with a Euclidean function, with
a choice of quotient defined by the ring: in $\Z$ the floor quotient
($-7 = (-4) \cdot 2 + 1$), in $\Z[i]$ the nearest Gaussian integer to
$x / y$ ($7 + 3i = 3 (2 + i) + 1$), and the usual division with remainder
in $F[x]$. There is no default implementation; rings without one return
`GR_UNABLE`.
