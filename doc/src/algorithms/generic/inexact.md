---
title: Inexact elements and truth values
summary: The enclosure semantics that make ball arithmetic, power series with error terms and exact real numbers rigorous in generic code: inclusion-preserving operations, three-valued predicates, and weak normalisation of polynomials.
order: 20
---

Generic algorithms are written for exact rings, but many of FLINT's
rings represent elements inexactly: a real ball $[m \pm r]$ stands for
an unknown real number in an interval, and a power series
$a_0 + a_1 x + O(x^n)$ for any series with those first terms. The
framework makes such rings usable by generic code without losing rigour
by fixing a semantics for inexact elements and a three-valued logic for
predicates [@Joh2025gr, §8].

## Enclosure semantics {#enclosures}

::: definition Enclosure {#def-enclosure}
Let $R$ be a ring and $\tilde R$ the set of representations of a
context for $R$. Each $\tilde x \in \tilde R$ denotes a nonempty set
$S(\tilde x) \subseteq R$, its *enclosure*. An exact representation
denotes a single element. An implementation $\tilde f$ of an operation
$f\colon R^k \to R$ is *inclusion-preserving* if, whenever it succeeds,
$$
f(x_1, \ldots, x_k) \in S(\tilde f(\tilde x_1, \ldots, \tilde x_k))
\quad \text{for all } x_i \in S(\tilde x_i).
$$
:::

For real and complex balls (`arb`, `acb`), $S$ is the ball. For power
series with an error term (`gr_series`), an element
$a_0 + \cdots + a_{n-1} x^{n-1} + O(x^n)$ with ball coefficients denotes
all series whose first $n$ coefficients lie in the balls; the error term
encloses the unknown tail. Exact rings, including the exact real numbers
of `ca`, have singleton enclosures. Every method of a ring context must
be inclusion-preserving; this is the meaning of "a ring may have inexact
representations, but must handle inexactness rigorously" in [`gr`](gr).

::: theorem Soundness of composition {#thm-composition}
Let $P$ be a program that computes $y = F(x_1, \ldots, x_k)$ from its
inputs by a finite sequence of operations, each implemented by an
inclusion-preserving method, with branches only on predicates that are
sound in the sense of [](#def-sound). If $P$ returns `GR_SUCCESS` on
inputs $\tilde x_1, \ldots, \tilde x_k$, then its output $\tilde y$
satisfies $F(x_1, \ldots, x_k) \in S(\tilde y)$ for all
$x_i \in S(\tilde x_i)$.
:::

::: proof
Fix exact values $x_i \in S(\tilde x_i)$ and follow the exact computation
of $F$ alongside $P$. By induction on the number of operations executed,
every intermediate exact value lies in the enclosure of the
corresponding representation: this holds for the inputs, and an
inclusion-preserving operation applied to enclosures of its exact
arguments returns an enclosure of the exact result. A branch is taken on
a predicate value `T_TRUE` or `T_FALSE`, which by soundness is the truth
value of the predicate at the exact values, so the exact computation
takes the same branch; on `T_UNKNOWN` the program does not branch on
the predicate but fails (or takes a path valid in both cases). Hence the
two computations stay in step and the final exact value lies in
$S(\tilde y)$.
:::

The theorem is why the generic algorithms of `gr_poly`, `gr_mat` and
`gr_series` are rigorous over balls without any code specific to
balls. It also explains the rule that predicates may never guess.

## Three-valued predicates {#truth}

::: definition Sound predicate {#def-sound}
An implementation $\tilde P$ of a predicate $P$ on $R^k$ is *sound* if
it returns `T_TRUE` only when $P(x_1, \ldots, x_k)$ holds for all
$x_i \in S(\tilde x_i)$, and `T_FALSE` only when it fails for all such
$x_i$. Otherwise it returns `T_UNKNOWN`.
:::

`T_UNKNOWN` arises for three reasons: the enclosures do not decide the
question (is $[0 \pm 0.1]$ zero?), the question is undecidable for the
representation (equality of exact real numbers involving
transcendental constants, where `ca` uses heuristics that may fail
[@Joh2021calcium; @Ric1997]), or deciding it is too expensive (is a
2048-bit modulus prime?). Sound predicates are closed under Kleene's
connectives, which `truth_and`, `truth_or` and `truth_not` implement:
if $\tilde P$ and $\tilde Q$ are sound, so are `truth_and`$(\tilde P,
\tilde Q)$ and the others, since for example `T_FALSE` from
`truth_and` requires one operand to be `T_FALSE`, so that the
conjunction fails at all points.

A generic algorithm that meets `T_UNKNOWN` must either continue in a way
that is correct in both cases or stop with `GR_UNABLE`. Gaussian
elimination is an example: `gr_mat_find_nonzero_pivot` returns the first
row whose entry is known to be nonzero (`T_FALSE` from `gr_is_zero`),
returns `GR_DOMAIN` if all entries are known to be zero (the column has
no pivot), and returns `GR_UNABLE` if no entry is known to be nonzero but
some are not known to be zero. A rank or a determinant over balls is
therefore either correct or not computed, never wrong.

## Weak normalisation {#weak-normalisation}

Dense polynomials are normally stored with a nonzero leading coefficient,
so that the length is the degree plus one. Over an inexact ring this
cannot be guaranteed: the leading coefficient may be a ball that
contains zero. `gr_poly` therefore normalises *weakly*: trailing
coefficients are removed only while `gr_is_zero` returns `T_TRUE`
(`_gr_vec_normalise_weak`). A leading coefficient $[0 \pm 0.001]$ is kept.

::: proposition {#prop-weak}
If every coefficient operation is inclusion-preserving and polynomials
are weakly normalised, then for every polynomial $f$ in the enclosure of
a stored polynomial of length $\ell$, $\deg f \le \ell - 1$.
:::

::: proof
The stored coefficients of index at least $\ell$ are exactly zero, so
every $f$ in the enclosure has zero coefficients there; and weak
normalisation removes a stored coefficient only if it is known to be
exactly zero, which does not change the enclosure.
:::

So the length is an upper bound for the degree, and operations that need
the exact degree must check the leading coefficient. Division with
remainder by $B$ requires the leading coefficient of $B$ to be invertible:
`gr_poly_divrem` (also called by `gr_euclidean_divrem` in a polynomial
ring) returns `GR_DOMAIN` when $B = 0$ and `GR_UNABLE` when the
leading coefficient of $B$ is not known to be nonzero, instead of
computing a quotient for a degree that may be wrong. Polynomial GCDs,
which depend on the degrees of all remainders, fail in the same way. The
program below shows this, together with the cancellation $f - f$, which
is weakly normalised to a polynomial of length 3 that is not known to be
zero, and a power series computation in which the error term carries the
unknown tail.

```c example name=gr_inexact
#include <flint/flint.h>
#include <flint/gr.h>
#include <flint/gr_poly.h>
#include <flint/gr_series.h>
#include <flint/gr_special.h>

int main(void)
{
    gr_ctx_t RR, RRx, S;
    gr_ptr f, g, q, r, s;
    int status;

    gr_ctx_init_real_arb(RR, 53);
    gr_ctx_init_gr_poly(RRx, RR);
    f = gr_heap_init(RRx);
    g = gr_heap_init(RRx);
    q = gr_heap_init(RRx);
    r = gr_heap_init(RRx);

    /* a polynomial whose leading coefficient may be zero */
    GR_MUST_SUCCEED(gr_set_str(f, "[0 +/- 0.001]*x^2 + x + 1", RRx));
    GR_MUST_SUCCEED(gr_set_str(g, "x^3 + 2", RRx));
    flint_printf("f = %{gr}, length %wd\n", f, RRx, gr_poly_length(f, RR));
    status = gr_euclidean_divrem(q, r, g, f, RRx);
    flint_printf("g divrem f: status %d\n", status);
    GR_MUST_SUCCEED(gr_sub(q, f, f, RRx));
    flint_printf("f - f = %{gr}, length %wd\n", q, RRx, gr_poly_length(q, RR));
    flint_printf("f - f == 0: ");
    truth_println(gr_is_zero(q, RRx));

    /* a power series: the error term encloses the unknown tail */
    gr_series_ctx_init(S, RR, 4);
    s = gr_heap_init(S);
    GR_MUST_SUCCEED(gr_set_str(s, "[1 +/- 0.001] + x + O(x^3)", S));
    status = gr_exp(s, s, S);
    flint_printf("exp(s) = %{gr}, status %d\n", s, S, status);

    gr_heap_clear(s, S);
    gr_ctx_clear(S);
    gr_heap_clear(f, RRx);
    gr_heap_clear(g, RRx);
    gr_heap_clear(q, RRx);
    gr_heap_clear(r, RRx);
    gr_ctx_clear(RRx);
    gr_ctx_clear(RR);
    return 0;
}
```
```output
f = 1 + x + [+/- 1.01e-3]*x^2, length 3
g divrem f: status 2
f - f = [+/- 2.01e-3]*x^2, length 3
f - f == 0: T_UNKNOWN
exp(s) = [2.72 +/- 4.44e-3] + [2.72 +/- 4.44e-3]*x + [1.36 +/- 2.22e-3]*x^2 + O(x^3), status 0
```

The constant and linear coefficients of $f - f$ are exactly zero because
$1 - 1$ and $x - x$ are computed exactly with balls of radius zero; the
quadratic coefficient $[0 \pm 0.001] - [0 \pm 0.001]$ is $[0 \pm 0.002]$,
since the two operands are not known to be the same number. Precision
of the series is limited by the error term $O(x^3)$ of the input, not by
the precision 4 of the ring.

## Exact rings with undecidable predicates {#undecidable}

The exact real and complex numbers of [`ca`](ca) are exact in the sense
of [](#def-enclosure), but equality is undecidable in general for
expressions involving transcendental functions [@Ric1997]. `ca`
decides it when it can (by exact algebraic methods or by numerical
separation) and otherwise returns `T_UNKNOWN` [@Joh2021calcium]. For
generic algorithms the situation is the same as with balls: a pivot
whose vanishing cannot be decided makes elimination return
`GR_UNABLE`, and the result is never wrong.

## Floating-point arithmetic {#floating-point}

Floating-point numbers (`arf`, `acf`, `nfloat`, `dfloat`) do not satisfy
enclosure semantics: rounding errors are not tracked, so a result can be
arbitrarily wrong after cancellation. Their contexts report
`gr_ctx_is_ring` as `T_FALSE` and `gr_ctx_is_approx_commutative_ring` as
`T_TRUE`, and their predicates compare the floating-point values. They
can be used with generic algorithms to obtain fast approximate results,
which is how `nfloat` serves approximate root-finding and linear
solving, but the theorem above does not apply to them.
