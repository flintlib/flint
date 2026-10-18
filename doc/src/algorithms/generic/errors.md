---
title: Domain errors and unable errors
summary: The meaning of the status codes GR_DOMAIN and GR_UNABLE, how statuses compose, why there are no junk values, how algorithms use the distinction for fallback strategies, and dynamic evaluation in rings that pretend to be fields.
order: 30
---

Almost every function of the generic interface returns a status code:
`GR_SUCCESS`, or a combination of the flags `GR_DOMAIN` and
`GR_UNABLE` [@Joh2025gr, §9]. The two flags answer different questions,
and generic algorithms depend on the difference: a domain error is a
mathematical fact about the input, while an unable error is a statement
about the implementation.

## The two kinds of failure {#kinds}

`GR_DOMAIN` means that the requested result has no value in the target
structure:

| Situation | Example |
|---|---|
| division by zero | $1/0$ in any nonzero ring (but $0/0 = 0$ in the zero ring) |
| no quotient | $1/2$ in $\Z$; $3/2$ in $\Z/6\Z$ |
| no root | $\sqrt{2}$ in $\mathbb{Q}$; $\sqrt{-1}$ in $\mathbb{R}$; $\sqrt{3}$ in $\Z/7\Z$ |
| no solution | an inconsistent linear system; the inverse of a singular matrix |
| pole | $\Gamma(0)$, $\zeta(1)$ |
| conversion | $1/7$ in $\Z/7\Z$; $-1$ to `ulong`; $i$ to $\mathbb{R}$ |
| shapes | adding matrices of different shapes; `gr_one` among matrices of all shapes |

`GR_UNABLE` means that the result may exist but was not computed:

| Situation | Example |
|---|---|
| not implemented | $\exp(2)$ in $\Z/7\Z$; `gr_div` in a ring without division |
| undecidable from the data | $1 / [0 \pm 0.5]$; the rank of a matrix of balls with an undecided pivot |
| too large | $2^{2^{100}}$ in $\Z$ |
| too expensive | a `qqbar` or `ca` operation beyond the evaluation limits of the context |
| a zero divisor in a ring pretending to be a field | see [below](#dynamic-evaluation) |

Conversions to C types are treated as conversions to a structure:
`gr_get_ui` and `gr_get_si` return `GR_DOMAIN` when the value is out of
range of the type. For inexact inputs, `GR_DOMAIN` is returned only if
the result is
undefined for every exact value in the enclosures of the inputs: $1/0$
in ball arithmetic is a domain error, but $1/[0 \pm 0.5]$ is unable,
since the ball also contains invertible numbers. When the flag
`GR_UNABLE` is set, a domain error may also be present without being
known; `GR_DOMAIN | GR_UNABLE` therefore means "unable".

The framework has no *junk values*. Proof assistants often define
$1/0 = 0$ so that division is total, and floating-point arithmetic
returns infinities and NaNs; a `gr` function instead reports the
failure, and the output variable holds an unspecified (valid) element
that must not be used. Ball arithmetic rings do not have infinite or NaN
elements, but they do have balls of infinite radius, which are legitimate
enclosures.

## Composition of statuses {#composition}

Statuses are combined with bitwise or. Code that needs every
intermediate result accumulates the statuses of a sequence of operations
and checks once at the end, which keeps generic code free of error
handling in the common case:

```c
status = GR_SUCCESS;
status |= gr_mul(t, x, y, ctx);
status |= gr_add(t, t, z, ctx);
status |= gr_sqrt(res, t, ctx);
if (status != GR_SUCCESS)
    ...
```

::: proposition {#prop-or}
Let a straight-line program compute $F(x)$ by operations each of which
is needed for the result and returns `GR_DOMAIN` only if its result is
undefined for all exact values of its inputs. If the OR of the statuses
is `GR_SUCCESS`, the result encloses $F(x)$. If it is exactly
`GR_DOMAIN`, then $F$ is undefined at every exact input in the
enclosures.
:::

::: proof
The first statement is [](algorithms/generic/inexact#thm-composition).
For the second, no step returned `GR_UNABLE`, so all steps before the
first failing one succeeded and enclose the exact intermediate values;
the failing step reported that its result is undefined for all of those
values, so the expression $F$, which needs that result, is undefined.
:::

Once a status is nonzero the later operations work on unspecified
values; they still run (and may return any status), which is harmless
because the combined status already reports the failure. Code that
wants to stop early tests the status after the critical step.

## Fallback strategies {#fallbacks}

The distinction between the two flags tells an algorithm what to do
next. A domain error is final: no other method will produce the
result, so it is propagated to the caller. An unable error invites
another attempt: higher precision, an exact representation, or a
different algorithm. A typical pattern is to try a fast inexact method
first and to fall back to a slower exact one only when it is unable to
decide. The program below inverts $2 \times 2$ matrices with real ball
arithmetic and retries with exact real algebraic numbers (`qqbar`) only
on `GR_UNABLE`. The first matrix is invertible, which the balls prove;
the second is singular, which the balls cannot decide and the exact
computation proves; the third is exactly zero, which the balls already
prove singular.

```c example name=gr_fallback
#include <flint/flint.h>
#include <flint/gr.h>
#include <flint/gr_mat.h>

/* invert a 2 x 2 matrix given as a string: try balls, then exact numbers */
int inverse(const char * s)
{
    gr_ctx_t RR, QQbar;
    gr_mat_t A, B;
    int status;

    gr_ctx_init_real_arb(RR, 64);
    gr_mat_init(A, 2, 2, RR);
    gr_mat_init(B, 2, 2, RR);
    GR_MUST_SUCCEED(gr_mat_set_str(A, s, 0, RR));
    status = gr_mat_inv(B, A, RR);
    flint_printf("%s\n  balls: status %d\n", s, status);
    gr_mat_clear(A, RR);
    gr_mat_clear(B, RR);
    gr_ctx_clear(RR);

    if (status == GR_UNABLE)     /* not decided: fall back to exact */
    {
        gr_ctx_init_real_qqbar(QQbar);
        gr_mat_init(A, 2, 2, QQbar);
        gr_mat_init(B, 2, 2, QQbar);
        GR_MUST_SUCCEED(gr_mat_set_str(A, s, 0, QQbar));
        status = gr_mat_inv(B, A, QQbar);
        flint_printf("  exact: status %d\n", status);
        gr_mat_clear(A, QQbar);
        gr_mat_clear(B, QQbar);
        gr_ctx_clear(QQbar);
    }

    if (status == GR_DOMAIN)
        flint_printf("  singular\n");
    return status;
}

int main(void)
{
    inverse("[[1, sqrt(2)], [sqrt(2), 3]]");
    inverse("[[1, sqrt(2)], [sqrt(2), 2]]");
    inverse("[[0, 0], [0, 0]]");
    return 0;
}
```
```output
[[1, sqrt(2)], [sqrt(2), 3]]
  balls: status 0
[[1, sqrt(2)], [sqrt(2), 2]]
  balls: status 2
  exact: status 1
  singular
[[0, 0], [0, 0]]
  balls: status 1
  singular
```

Limits serve the same purpose inside a computation: the evaluation
limits of `qqbar` contexts (`_gr_ctx_qqbar_set_limits`) and the options
of `ca` contexts turn operations that would take too long into
`GR_UNABLE`, so that an algorithm can recover and try another strategy
instead of running indefinitely.

## Dynamic evaluation {#dynamic-evaluation}

Many algorithms for fields (Euclid's algorithm, Gaussian elimination)
work over a commutative ring as long as every element they need to
invert is invertible. Dynamic evaluation [@DDDD1985] runs such an
algorithm over a ring that is only *conjectured* to be a field, such as
$\Z/n\Z$ for an $n$ believed prime or $K[x]/(m)$ for an $m$ believed
irreducible; when an inversion fails, the non-invertible element splits
the ring.

The generic interface supports this through *pretend fields*
(`gr_ctx_set_is_pretend_field`). In a ring that pretends to be a field,
the inversion of a nonzero element that is not invertible returns
`GR_UNABLE`, not `GR_DOMAIN`, and records the element in the context,
from which `gr_ctx_recover_zero_divisor` retrieves it. Algorithms that
only need inversions to fail cleanly check `gr_ctx_is_pretend_field`
instead of `gr_ctx_is_field`, and propagate the unable error.

::: lemma {#lem-split}
Let $n \ge 2$ and let $a \in \Z/n\Z$ be nonzero and not invertible.
Then $g = \gcd(a, n)$ satisfies $1 < g < n$. Similarly, if $K$ is a field,
$m \in K[x]$ is monic of positive degree and $a \in K[x]/(m)$ is nonzero
and not invertible, then the monic $g = \gcd(a, m)$ is a factor of $m$
with $0 < \deg g < \deg m$.
:::

::: proof
If $g = 1$ then $a$ is invertible modulo $n$ by Bézout's identity, a
contradiction; if $g = n$ then $n \mid a$ and $a = 0$. The polynomial
case is identical, with degrees in place of sizes.
:::

::: proposition {#prop-dynamic}
Let an algorithm for fields be run over $\Z/n\Z$ pretending to be a
field, where every element that the algorithm finds to be nonzero is
subsequently inverted (as with pivots in Gaussian elimination and
leading coefficients in Euclid's algorithm). If it returns
`GR_SUCCESS`, then for every prime $p \mid n$ the reduction of its
output modulo $p$ is the output of the same algorithm over the field
$\Z/p\Z$ on the reduced input.
:::

::: proof
Every element found nonzero was inverted successfully, so it is a unit
modulo $n$ and hence nonzero modulo every $p \mid n$; every element found
zero is zero modulo $n$ and hence modulo $p$. So the zero tests, and
therefore all branches, agree with those of the computation over
$\Z/p\Z$, and all ring operations commute with reduction modulo $p$.
:::

If the algorithm fails instead, the lemma gives a factorisation of $n$
(or of $m$), and the computation can be continued on the factors, using
the Chinese remainder theorem to combine results. The program below
computes a determinant by LU decomposition modulo $1003 = 17 \cdot 59$.
`gr_ctx_init_nmod` knows that 1003 is composite, so LU decomposition,
which needs a field, fails. Pretending that the ring is a field, the
first pivot 17 cannot be inverted, and the recovered zero divisor gives
the factor 17.

```c example name=gr_dynamic_evaluation
#include <flint/flint.h>
#include <flint/ulong_extras.h>
#include <flint/gr.h>
#include <flint/gr_mat.h>

int main(void)
{
    gr_ctx_t R;
    gr_mat_t A;
    gr_ptr d, z;
    ulong n = 1003;    /* = 17 * 59, but we pretend that it is prime */
    int status;

    gr_ctx_init_nmod(R, n);
    gr_mat_init(A, 3, 3, R);
    d = gr_heap_init(R);
    z = gr_heap_init(R);
    GR_MUST_SUCCEED(gr_mat_set_str(A, "[[17, 1, 2], [3, 4, 5], [6, 7, 9]]", 0, R));

    status = gr_mat_det_lu(d, A, R);
    flint_printf("det_lu without pretence: status %d\n", status);

    GR_MUST_SUCCEED(gr_ctx_set_is_pretend_field(R, T_TRUE));
    status = gr_mat_det_lu(d, A, R);
    flint_printf("det_lu as a pretend field: status %d\n", status);
    if (status == GR_UNABLE && gr_ctx_recover_zero_divisor(z, R) == GR_SUCCESS)
    {
        ulong a = *(ulong *) z;
        flint_printf("zero divisor %wu, factor %wu of %wu\n", a, n_gcd(a, n), n);
    }

    gr_heap_clear(d, R);
    gr_heap_clear(z, R);
    gr_mat_clear(A, R);
    gr_ctx_clear(R);
    return 0;
}
```
```output
det_lu without pretence: status 3
det_lu as a pretend field: status 2
zero divisor 17, factor 17 of 1003
```

The contexts `gr_ctx_init_nmod`, `gr_ctx_init_fmpz_mod`,
`gr_ctx_init_mpn_mod` and `gr_ctx_init_gr_poly_quotient` support
pretending; the quotient rings record nontrivial monic factors of the
modulus and can be refined in place with
`gr_poly_quotient_ctx_refine`. Recording is serialised by a global lock,
so pretend fields can be used from several threads.
