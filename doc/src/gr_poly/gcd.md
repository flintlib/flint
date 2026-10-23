---
title: Greatest common divisors over generic rings
functions: gr_poly_gcd gr_poly_xgcd
---

`gr_poly_gcd` and `gr_poly_xgcd` compute a greatest common divisor of
two polynomials, and with it the Bézout cofactors, over any ring whose
context says that the GCD makes sense: a field (or a ring that pretends
to be one, see `gr_ctx_set_is_pretend_field`) or a unique factorisation
domain. The algorithms are compared in [](algorithms/polynomials/gcd).

::: summary algorithms/polynomials/gcd :::

## Normalisation of the result {#normalisation}

A GCD is unique only up to a unit factor. The non-underscore functions
return the canonical associate (`gr_poly_canonical_associate`): the
monic GCD over a field, and over $\Z$ the GCD with positive leading
coefficient. $\gcd(0, 0) = 0$, and $\gcd(a, 0)$ is the canonical
associate of $a$. For the extended GCD, $G = SA + TB$ holds exactly;
over a field $G$ is the monic GCD, $\deg S < \deg B - \deg G$ and $\deg T
< \deg A - \deg G$. Over a ring that is not a field, there may be no
cofactors for the GCD itself (in $\Z[x]$, $1 = \gcd(2, x)$ is not a
combination $2S + xT$), and `gr_poly_xgcd` then returns a constant
multiple of the GCD: the subresultant pseudo-remainder sequence of the
primitive parts of $A$ and $B$, which carries the cofactors along,
ends with a multiple of their GCD, and this is multiplied by the least
common multiple of the contents of $A$ and $B$.

Both functions return `GR_UNABLE` if a leading coefficient of an input
is not provably nonzero, which can happen over inexact rings (see
[](gr_poly/getting-started#inexact)), and if the ring is neither a field
nor a unique factorisation domain. An algorithm that meets a
non-invertible leading coefficient returns `GR_DOMAIN`, or `GR_UNABLE`
when the ring pretends to be a field.

## Algorithm selection {#selection}

`gr_poly_gcd` and `gr_poly_xgcd` call the ring methods
`GR_METHOD_POLY_GCD` and `GR_METHOD_POLY_XGCD`. The defaults,
`_gr_poly_gcd_generic` and `_gr_poly_xgcd_generic`, choose:

- over a field (or a pretended field) that is finite, the half-gcd
  algorithm ([](algorithms/polynomials/half-gcd)) when both inputs have
  length at least 200, with the Euclidean algorithm for remainders
  shorter than 200 and inside the half-gcd recursion below length 100;
- over other fields, such as $\mathbb{Q}$, number fields, real and complex balls,
  the Euclidean algorithm, since the half-gcd algorithm suffers from
  coefficient growth over exact infinite fields and from numerical
  instability over approximate fields;
- over unique factorisation domains that are not fields, such as $\Z$
  and $\Z[y]$, the subresultant algorithm: the contents are removed,
  the primitive parts are run through the subresultant pseudo-remainder
  sequence, and the GCD of the contents is multiplied back;
- otherwise `GR_UNABLE`.

Several rings override the methods with tuned algorithms:

| Ring | GCD | Extended GCD |
|---|---|---|
| `nmod` (prime modulus) | `_nmod_poly_gcd`: Euclidean algorithm (in Montgomery form for some moduli) below a cutoff tabulated by the bit size of $n$ (470 to 2199), half-gcd above | default |
| `fmpz_mod` | Euclidean below length 256, half-gcd above (inner cutoff 128) | default |
| `mpn_mod` | Euclidean below length 240, half-gcd above (inner cutoff 80) | the same cutoffs |
| `fq`, `fq_nmod`, `fq_zech` | `_fq_poly_gcd` and its analogues, which call the generic Euclidean or half-gcd algorithm with their own cutoffs | `_fq_poly_xgcd` and analogues |
| `fmpz` | `_fmpz_poly_gcd`: heuristic and multimodular algorithms ([](algorithms/polynomials/gcd-integers)) | default |

There is no override for $\mathbb{Q}$: over `fmpq` the generic Euclidean
algorithm is used, whose intermediate coefficients can grow large; for
large inputs `fmpq_poly_gcd` is much faster.

The underscore functions `_gr_poly_gcd_euclidean`, `_gr_poly_gcd_hgcd`
and `_gr_poly_gcd_subresultant`, and the corresponding *xgcd*
functions, run one algorithm unconditionally; `_gr_poly_gcd_hgcd` and
`_gr_poly_xgcd_hgcd` take the half-gcd recursion cutoff `inner_cutoff`
(respectively `hgcd_cutoff`) and the length `cutoff` below which the
outer loop switches to the Euclidean algorithm.

## Conventions of the underscore functions {#underscore}

`_gr_poly_gcd` and the algorithm-specific underscore GCD functions
require $\mathit{lenA} \ge \mathit{lenB} \ge 1$ and nonzero leading
coefficients, write at most $\mathit{lenB}$ coefficients to $G$ and its
length to `*lenG`, and need not canonicalise the result.

`_gr_poly_xgcd` and its variants require $\mathit{lenA} \ge \mathit{lenB}
\ge 2$ and do not allow aliasing. The output arrays must have room for
$\mathit{lenB}$ coefficients ($G$), $\mathit{lenB} - 1$ ($S$) and
$\mathit{lenA} - 1$ ($T$), even if the results are shorter. The length
of $G$ is written to `*lenG`, and coefficients of $G$ beyond it may
contain garbage; on success the unused high coefficients of $S$ and $T$
are zeroed. A length $\mathit{lenG} = 0$ is possible if the algorithm
fails, for example at a non-invertible leading coefficient.

To wrap a different underscore implementation with the memory
management, aliasing rules and normalisation of the non-underscore
functions, use `gr_poly_gcd_wrapper` and `gr_poly_xgcd_wrapper`.

## Example {#example}

The program computes an extended GCD of $A = (x + 1)(x - 2)$ and $B =
(x + 1)(x + 3)$ over $\Z$ and over $\mathbb{F}_7$. The GCD over $\Z$ is
$x + 1$, but it is not a $\Z[x]$-combination of $A$ and $B$, since $B -
A = 5(x + 1)$; the subresultant algorithm returns $G = 5(x + 1)$ with
$S = -1$, $T = 1$. The factor $5$ is the resultant of the cofactors $x -
2$ and $x + 3$. Over $\mathbb{F}_7$ the monic GCD is returned.

```c example name=gr_poly_xgcd_z_f7
#include <flint/flint.h>
#include <flint/gr.h>
#include <flint/gr_poly.h>

int main(void)
{
    gr_ctx_t ZZ, F;
    gr_poly_t a, b, g, s, t;

    gr_ctx_init_fmpz(ZZ);
    gr_poly_init(a, ZZ);
    gr_poly_init(b, ZZ);
    gr_poly_init(g, ZZ);
    gr_poly_init(s, ZZ);
    gr_poly_init(t, ZZ);

    /* a = (x + 1)(x - 2), b = (x + 1)(x + 3) */
    GR_MUST_SUCCEED(gr_poly_set_str(a, "x^2 - x - 2", "x", ZZ));
    GR_MUST_SUCCEED(gr_poly_set_str(b, "x^2 + 4*x + 3", "x", ZZ));

    GR_MUST_SUCCEED(gr_poly_xgcd(g, s, t, a, b, ZZ));
    flint_printf("over Z:    g = "); gr_poly_print(g, ZZ);
    flint_printf(", s = "); gr_poly_print(s, ZZ);
    flint_printf(", t = "); gr_poly_print(t, ZZ); flint_printf("\n");

    /* the same inputs over GF(7) */
    gr_ctx_init_nmod(F, 7);
    {
        gr_poly_t a7, b7, g7, s7, t7;
        gr_poly_init(a7, F); gr_poly_init(b7, F); gr_poly_init(g7, F);
        gr_poly_init(s7, F); gr_poly_init(t7, F);
        GR_MUST_SUCCEED(gr_poly_set_gr_poly_other(a7, a, ZZ, F));
        GR_MUST_SUCCEED(gr_poly_set_gr_poly_other(b7, b, ZZ, F));
        GR_MUST_SUCCEED(gr_poly_xgcd(g7, s7, t7, a7, b7, F));
        flint_printf("over GF(7): g = "); gr_poly_print(g7, F);
        flint_printf(", s = "); gr_poly_print(s7, F);
        flint_printf(", t = "); gr_poly_print(t7, F); flint_printf("\n");
        gr_poly_clear(a7, F); gr_poly_clear(b7, F); gr_poly_clear(g7, F);
        gr_poly_clear(s7, F); gr_poly_clear(t7, F);
    }
    gr_ctx_clear(F);

    gr_poly_clear(a, ZZ);
    gr_poly_clear(b, ZZ);
    gr_poly_clear(g, ZZ);
    gr_poly_clear(s, ZZ);
    gr_poly_clear(t, ZZ);
    gr_ctx_clear(ZZ);
    return 0;
}
```
```output
over Z:    g = 5 + 5*x, s = -1, t = 1
over GF(7): g = 1 + x, s = 4, t = 3
```

The GCD over $\mathbb{F}_7$ satisfies $4A + 3B = 7x^2 + 8x + 1 \equiv x +
1$.
