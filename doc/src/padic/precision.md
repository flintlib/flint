---
title: Precision in the padic module
summary: How a padic_t stores p^v u modulo p^N, what the precision of a variable means, why results are reduced to the precision of the output and what this implies for precision loss, with a first program.
order: 1
---

A `padic_t` holds an integer $u$, a valuation $v$ and a precision $N$, and
represents the class of $x = p^v u$ modulo $p^N$. The precision is a
property of the variable, like the precision of an MPFR variable: it is
chosen when the variable is initialised and stays fixed. Every function
treats the stored values of its inputs as exact elements of
$\Z[1/p] \subset \mathbf{Q}_p$, computes the exact result, and reduces it
modulo $p^N$ for the precision $N$ of the output. The precisions of the
inputs play no role.

## A first program {#first}

The program works in $\mathbf{Q}_7$. It stores $1/3$ and $12/49$ to
precision 5, multiplies them into a variable of precision 3, inverts
$1/3$ and converts back to a rational number.

```c example name=padic_precision_first
#include <flint/flint.h>
#include <flint/fmpq.h>
#include <flint/padic.h>

int main(void)
{
    fmpz_t p;
    fmpq_t q;
    padic_ctx_t ctx;
    padic_t x, y, z;

    fmpz_init_set_ui(p, 7);
    padic_ctx_init(ctx, p, 0, 10, PADIC_SERIES);
    fmpq_init(q);

    padic_init2(x, 5);              /* x is stored modulo 7^5 */
    padic_init2(y, 5);
    padic_init2(z, 3);              /* z is stored modulo 7^3 */

    fmpq_set_si(q, 1, 3);
    padic_set_fmpq(x, q, ctx);      /* 1/3 mod 7^5 */
    flint_printf("x = "); padic_print(x, ctx); flint_printf("\n");

    fmpq_set_si(q, 12, 49);
    padic_set_fmpq(y, q, ctx);      /* 12/49 = 12 * 7^-2 */
    flint_printf("y = "); padic_print(y, ctx);
    flint_printf("   (u = "); fmpz_print(padic_unit(y));
    flint_printf(", v = %wd)\n", padic_val(y));

    padic_mul(z, x, y, ctx);        /* exact product, reduced mod 7^3 */
    flint_printf("x*y mod 7^3 = "); padic_print(z, ctx); flint_printf("\n");

    padic_inv(z, x, ctx);
    flint_printf("1/x mod 7^3 = "); padic_print(z, ctx); flint_printf("\n");

    ctx->mode = PADIC_TERSE;
    padic_get_fmpq(q, y, ctx);
    flint_printf("y as a rational: "); fmpq_print(q); flint_printf("\n");
    flint_printf("x in terse mode: "); padic_print(x, ctx); flint_printf("\n");

    padic_clear(x); padic_clear(y); padic_clear(z);
    fmpq_clear(q);
    padic_ctx_clear(ctx);
    fmpz_clear(p);
    return 0;
}
```
```output
x = 5 + 4*7^1 + 4*7^2 + 4*7^3 + 4*7^4
y = 5*7^-2 + 1*7^-1   (u = 12, v = -2)
x*y mod 7^3 = 4*7^-2
1/x mod 7^3 = 3
y as a rational: 12/49
x in terse mode: 11205
```

The context stores the powers $p^0, \ldots, p^9$, which covers every
reduction in this program; other powers would be computed on demand.
The value of $y$ is exact, since $12/49 = 12 \cdot 7^{-2}$ terminates.
In terse mode $x$ prints as the integer $11205 \equiv 1/3 \pmod{7^5}$,
the representative actually stored, not as $1/3$.

## Canonical and reduced form {#forms}

A value $p^v u$ is *canonical* if $p \nmid u$ or $u = v = 0$; the
valuation of a nonzero canonical element is then $v$. It is *reduced* if
in addition $0 < u < p^{N - v}$, so that it is the smallest nonnegative
representative of its class modulo $p^N$. The functions of the module
return reduced values and accept any canonical values. Code that writes
the fields directly, through `padic_unit`, `padic_val` and `padic_prec`,
must restore this with `padic_reduce` (or `_padic_canonicalise` and
`_padic_reduce`).

Since the unit is reduced modulo $p^{N - v}$, the number of stored
digits depends on the valuation: an element of valuation $-10$ in a
variable of precision 5 stores 15 digits. Conversely, any element of
valuation at least $N$ is stored as zero, and `padic_one` gives zero in
a variable of precision $N \le 0$.

## Precision is not tracked {#not-tracked}

Because the inputs are treated as exact, the module never reports a
loss of precision. If $x$ and $y$ agree modulo $p^{20}$ and are stored to
precision 20, `padic_sub` returns $x - y$ reduced to the precision of the
output, which may be 20 as well, although in a computation with genuine
$p$-adic numbers only the class of $x - y$ modulo $p^{20}$ is known, and
its leading digits are meaningless. In the same way an inverse
$1/(p^v u)$ with $v > 0$ is only determined modulo $p^{N - 2v}$ by its
input modulo $p^N$, but `padic_inv` returns it to the full precision of
the output. The rules for how much precision each operation preserves
are given in [](algorithms/p-adic/precision); following them is the
responsibility of the caller. The module `padic_radix.h` does track
error bounds and should be preferred when this matters.

A related consequence is that `padic_equal` compares stored values.
Two variables of different precisions holding the same $p$-adic number
compare equal only if their stored representatives coincide; reducing
both to the smaller precision first gives the expected comparison.

## Choosing the context {#context}

The context stores a table of $p^e$ for $\mathtt{min} \le e <
\mathtt{max}$. Almost every function reduces modulo $p^{N - v}$, and a
power found in the table saves an exponentiation, so the table should
cover the precisions in use, for example $\mathtt{min} = 0$ and
$\mathtt{max}$ a little larger than the largest $N$. The printing mode
is a field of the context that may be changed at any time, as in the
program above. A context may be shared by any number of variables of
different precisions, and is read-only after initialisation.
