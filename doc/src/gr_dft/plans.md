---
title: Constructing plans
functions: gr_dft_default_root gr_dft_precomp_init_root gr_dft_precomp_init gr_dft_precomp_init_karatsuba gr_dft_precomp_clear
---

Every transform in this module is applied through a plan, which is
constructed once for a length $n$, a ring and a root of unity $w$, and can
then be used for any number of forward and inverse transforms. Plan
construction chooses the algorithm, splits the length into the sub-plans
that the algorithm needs, and computes the table of powers
$w^0, \ldots, w^{n-1}$ (or the parts of it that the algorithm reads).
Its cost is comparable to that of a few transforms, so a plan should be
reused whenever several vectors of the same length are transformed.

## Which constructor {#which}

| Ring | Constructor | Root |
|---|---|---|
| $\Z/m\Z$, finite fields, other rings without $e^{-2\pi i/n}$ | `gr_dft_precomp_init_root` | supplied by the caller |
| `acb`, exact complex numbers, other rings with `gr_exp_pi_i` | `gr_dft_precomp_init` | $e^{-2\pi i/n}$ |
| `acb` over `arb`, fixed-point contexts, other complex rings stored as two real parts | `gr_dft_precomp_init_karatsuba` | $e^{-2\pi i/n}$ |

`gr_dft_default_root` computes the root used by the last two
constructors. It evaluates $e^{\pi i q}$ with $q = -2/n$ by `gr_exp_pi_i`,
or, if the ring lacks that method, $\exp(-2\pi i/n)$ from `gr_pi`, `gr_i`
and `gr_exp`. Over a ring without these operations, or without the
value (real numbers for $n \ge 3$), the constructors return its status,
`GR_UNABLE` or `GR_DOMAIN`. For $n \le 2$ the default root $\pm 1$
exists in every ring.

With the default root, the constructors do not compute the table from
powers of $w$: over `acb` each entry is rounded from a running product at
a precision raised by about $\log_2 n + 16$ bits, completed by the
symmetries of the unit circle, and over other rings it is computed entry
by entry with `gr_exp_pi_i` (over the fixed-point contexts, see
[](gr_dft/nfixed#plans)). The entries therefore have about the accuracy
of single evaluations, and over `acb` those equal to $\pm 1$ and $\pm i$
are exact.

## Roots of unity {#roots}

`gr_dft_precomp_init_root` requires a principal $n$-th root of unity:
$w^n = 1$, and $w^{n/p} - 1$ is not a zero divisor for every prime
$p \mid n$. Over a field this means that $w$ has order $n$. Over
$\Z/m\Z$ with composite $m$ it means that $w$ has order $n$ modulo every
prime factor of $m$. [](algorithms/fft/dft-algorithms#principal) explains
why the fast algorithms need this. The constructor checks two necessary
conditions, using `gr_is_neg_one` and `gr_is_one` on the table entries:

- if $n$ is even, $w^{n/2} = -1$;
- $w^{n/p} \ne 1$ for every odd prime $p \mid n$.

A definite failure gives `GR_DOMAIN`. When the ring cannot decide (an
inexact ring returns `T_UNKNOWN`), the root is accepted. The condition
$w^n = 1$ itself is not checked; the table is simply
$w^0, \ldots, w^{n-1}$. For $n$ a power of two the checks are complete
(and only $w^{n/2} = -1$ is needed); in a ring with zero divisors they
are not, as the [second example](#example-zero-divisor) shows: the
algorithms then return wrong results without any error. The
constructors with default roots perform no checks.

Over $\Z/p\Z$ with $p$ prime, a root of order $n$ exists if and only if
$n \mid p - 1$, and $g^{(p-1)/n}$ is one for any primitive root $g$
(`n_primitive_root_prime` computes one). Primes $p = c 2^k + 1$ with
large $k$, such as $998244353 = 119 \cdot 2^{23} + 1$, allow long
power-of-two transforms.

Over the complex fixed-point contexts (`gr_dft_ctx_init_nfixed_complex`),
the root table is always that of $e^{-2\pi i/n}$, computed to the full
precision of the context, and `w` only selects this root: it is accepted
if it lies within $1/n$ of $e^{-2\pi i/n}$, so a value rounded to a lower
precision is fine. `gr_dft_precomp_init_root` returns `GR_UNABLE` if `w`
approximates another primitive $n$-th root of unity, and `GR_DOMAIN` if
it approximates none.

## Complex mode {#complex}

`gr_dft_precomp_init_karatsuba` with a real ring `real_ctx` builds a plan
in complex mode. The ring `ctx` must store each element as a real part
followed by an imaginary part, both elements of `real_ctx`; only the
element sizes are checked: `GR_UNABLE` is returned unless an element of
`ctx` is exactly twice as large as an element of `real_ctx`. With `real_ctx` equal to `NULL` the function is
`gr_dft_precomp_init`. In complex mode:

- multiplications by $\pm 1$ and $\pm i$ (the powers $w^{jn/4}$) are
  done by swapping and negating parts, with no multiplication;
- multiplications by $w^{jn/8}$ for odd $j$, which have the form
  $c(1 \pm i)$, take two real multiplications;
- multiplications by other roots take three real multiplications with
  precomputed tables (the complex Karatsuba product), but only over
  exact rings and over rings whose `gr_ctx_get_real_prec` is at least
  1024 bits, and never over the fixed-point contexts. Otherwise they take
  one multiplication in `ctx`;
- `GR_DFT_ALG_AUTO` chooses split-radix rather than radix-2 for powers of
  two, unless `GR_DFT_SCRAMBLED` is set.

The special values are recognised from their exponents, which is valid
because complex mode always uses $w = e^{-2\pi i/n}$; over balls a
comparison could not prove that a computed value is exactly $-i$. Since
the formulas for $c(1 \pm i)$ use only the enclosure of $c$, they still
give enclosures of the exact products. The gain is largest for small
powers of two: a split-radix transform of length 16 in complex mode
needs 20 real multiplications.

## Failure and lifetime {#lifetime}

On any return value other than `GR_SUCCESS`, the plan is left in a
cleared state: it must not be used for transforms, and calling
`gr_dft_precomp_clear` on it is harmless. The plan stores pointers to
`ctx` and `real_ctx` (sub-plans included), so the contexts must outlive
the plan, and transforms must be given the same `ctx`. A plan is not
modified by transforms.

The status values of the constructors are:

| Status | Cause |
|---|---|
| `GR_DOMAIN` | $n = 0$; the algorithm does not support $n$ (see `GR_DFT_ALG_AUTO`); a failed root check; no default root (for example over `arb`) |
| `GR_UNABLE` | no default root (for example over `nmod`); `GR_DFT_ALG_BLUESTEIN` over a ring without default roots; complex mode with mismatched element sizes; an operation failed while building the tables |

## Example {#example}

The number 2 is a primitive root modulo 13, so it is a principal 12th
root of unity, while 4 has order 6. The automatic plan of length
$12 = 4 \cdot 3$ uses the prime factor algorithm.

```c example name=gr_dft_plan_nmod
#include <flint/flint.h>
#include <flint/gr.h>
#include <flint/gr_vec.h>
#include <flint/gr_dft.h>

int main(void)
{
    gr_ctx_t ctx;
    gr_dft_pre_t P;
    gr_ptr x, y, w;
    slong j, n = 12;
    int status;

    gr_ctx_init_nmod(ctx, 13);
    x = gr_heap_init_vec(n, ctx);
    y = gr_heap_init_vec(n, ctx);
    w = gr_heap_init(ctx);

    /* 4 has order 6 modulo 13: not a principal 12th root of unity */
    GR_MUST_SUCCEED(gr_set_ui(w, 4, ctx));
    status = gr_dft_precomp_init_root(P, w, n, GR_DFT_ALG_AUTO, 0, ctx);
    flint_printf("w = 4: %s\n", status == GR_DOMAIN ? "GR_DOMAIN" : "?");
    gr_dft_precomp_clear(P);   /* harmless after a failure */

    /* 2 is a primitive root modulo 13 */
    GR_MUST_SUCCEED(gr_set_ui(w, 2, ctx));
    GR_MUST_SUCCEED(gr_dft_precomp_init_root(P, w, n, GR_DFT_ALG_AUTO, 0, ctx));

    for (j = 0; j < n; j++)
        GR_MUST_SUCCEED(gr_set_si(GR_ENTRY(x, j, ctx->sizeof_elem), j, ctx));

    GR_MUST_SUCCEED(gr_dft_precomp(y, x, P, ctx));
    _gr_vec_print(y, n, ctx);
    flint_printf("\n");

    GR_MUST_SUCCEED(gr_dft_inverse_precomp(y, y, P, ctx));
    _gr_vec_print(y, n, ctx);
    flint_printf("\n");

    gr_dft_precomp_clear(P);
    gr_heap_clear(w, ctx);
    gr_heap_clear_vec(x, n, ctx);
    gr_heap_clear_vec(y, n, ctx);
    gr_ctx_clear(ctx);
    return 0;
}
```
```output
w = 4: GR_DOMAIN
[1, 12, 4, 11, 6, 5, 7, 9, 8, 3, 10, 2]
[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
```

## Example: a root that is not principal {#example-zero-divisor}

Modulo $91 = 7 \cdot 13$, the number $w = 79$ satisfies $w^3 = 1$ and
$w \ne 1$, so it passes the checks for $n = 3$. But $w \equiv 1 \pmod{13}$,
so $w - 1$ is a zero divisor and $1 + w + w^2 = 42 \ne 0$. The naive
algorithm evaluates the definition and is correct for any $w$ with
$w^n = 1$, while the mixed-radix kernel chosen by `GR_DFT_ALG_AUTO` uses
$1 + w + w^2 = 0$ and returns a different vector.

```c example name=gr_dft_plan_not_principal
#include <flint/flint.h>
#include <flint/gr.h>
#include <flint/gr_vec.h>
#include <flint/gr_dft.h>

/* w = 79 satisfies w^3 = 1 and w != 1 modulo 91 = 7 * 13,
   but w - 1 = 78 is a zero divisor */
int main(void)
{
    gr_ctx_t ctx;
    gr_dft_pre_t P1, P2;
    gr_ptr x, y, w;
    slong j, n = 3;

    gr_ctx_init_nmod(ctx, 91);
    x = gr_heap_init_vec(n, ctx);
    y = gr_heap_init_vec(n, ctx);
    w = gr_heap_init(ctx);
    GR_MUST_SUCCEED(gr_set_ui(w, 79, ctx));
    for (j = 0; j < n; j++)
        GR_MUST_SUCCEED(gr_set_si(GR_ENTRY(x, j, ctx->sizeof_elem), j + 1, ctx));

    /* both constructors accept w */
    GR_MUST_SUCCEED(gr_dft_precomp_init_root(P1, w, n, GR_DFT_ALG_NAIVE, 0, ctx));
    GR_MUST_SUCCEED(gr_dft_precomp_init_root(P2, w, n, GR_DFT_ALG_AUTO, 0, ctx));

    GR_MUST_SUCCEED(gr_dft_precomp(y, x, P1, ctx));
    flint_printf("naive: ");
    _gr_vec_print(y, n, ctx);
    GR_MUST_SUCCEED(gr_dft_precomp(y, x, P2, ctx));
    flint_printf("\nauto:  ");
    _gr_vec_print(y, n, ctx);
    flint_printf("\n");

    gr_dft_precomp_clear(P1);
    gr_dft_precomp_clear(P2);
    gr_heap_clear(w, ctx);
    gr_heap_clear_vec(x, n, ctx);
    gr_heap_clear_vec(y, n, ctx);
    gr_ctx_clear(ctx);
    return 0;
}
```
```output
naive: [6, 45, 71]
auto:  [6, 10, 78]
```

The naive result is the transform as defined, and the automatic one
differs from it by $-(1 + w + w^2)(0, x_2, x_1) = (0, -126, -84)$ modulo
91. Neither is useful for convolutions, because for this $w$ the inverse
transform does not undo the forward one.
