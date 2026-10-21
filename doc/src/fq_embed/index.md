---
title: fq_embed.h
header: fq_embed.h
template-headers: fq_embed_templates
description: Embeddings of a finite field $\mathbb{F}_{p^m}$ into $\mathbb{F}_{p^n}$, $m \mid n$, given by different defining polynomials, as matrices over $\mathbb{F}_p$, with projections by the trace.
category: Finite fields
order: 10
---

If $m$ divides $n$, the field $\mathbb{F}_{p^n}$ contains a unique
subfield with $p^m$ elements, but two FLINT contexts for
$\mathbb{F}_{p^m}$ and $\mathbb{F}_{p^n}$ are defined by unrelated
polynomials, so the embedding has to be computed. This module finds an
embedding, by computing a root in $\mathbb{F}_{p^n}$ of the defining
polynomial of $\mathbb{F}_{p^m}$, and represents it, and a projection
back onto the subfield, by matrices over $\mathbb{F}_p$ acting on
coordinate vectors with respect to the bases $1, x, \ldots, x^{d-1}$. It
also provides the linear algebra on which this rests: multiplication
and composition matrices, and the change of basis between the monomial
basis and its dual basis for the trace form.

The interface is generated from `fq_embed_templates.h` and is the same
for `fq_nmod_embed` and `fq_zech_embed`, with matrices and polynomials
over `nmod` instead of `fmpz_mod`. Coordinate vectors can be converted
with `fq_get_fmpz_mod_mat` and `fq_set_fmpz_mod_mat`.

::: summary algorithms/finite-fields/embeddings :::

All functions require both contexts to have the same characteristic $p$,
and matrix outputs to be initialised with the stated dimensions over the
modulus $p$; neither is checked.

## Embeddings {#embeddings}

::: function
void fq_embed_gens(fq_t gen_sub, fq_t gen_sup, fmpz_mod_poly_t minpoly, const fq_ctx_t sub_ctx, const fq_ctx_t sup_ctx)

Given contexts for fields $k$ = `sub_ctx` of degree $m$ and $K$ =
`sup_ctx` of degree $n$ with $m \mid n$, computes a generator `gen_sub`
of $k$ over $\mathbb{F}_p$, its minimal polynomial `minpoly`, and a root
`gen_sup` of `minpoly` in $K$. These data define the embedding
$k \to K$ that maps `gen_sub` to `gen_sup`.

::: details
For $m = 1$ both generators are set to the generator of the prime field
context ($-c_0/c_1$ for the defining polynomial $c_1 x + c_0$) and
`minpoly` is not modified (it should be $x - \mathit{gen\_sub}$, which
the caller must set before calling `fq_embed_matrices`); otherwise
`_fq_embed_gens_naive` is called. Different roots give different
embeddings, which differ by a power of the Frobenius map; the one
returned is not in general compatible with the Conway polynomials of
the two fields.
:::
:::

::: function
void _fq_embed_gens_naive(fq_t gen_sub, fq_t gen_sup, fmpz_mod_poly_t minpoly, const fq_ctx_t sub_ctx, const fq_ctx_t sup_ctx)

Computes the data of `fq_embed_gens` as follows: `gen_sub` is the
generator $x$ of `sub_ctx`, `minpoly` is the defining polynomial of
`sub_ctx`, and `gen_sup` is a root of it in `sup_ctx`, found with
`fq_poly_factor_split_single`, that is, by root finding for a
polynomial of degree $m$ over $K$. Requires $m > 1$ to divide $n$.
:::

::: function
void fq_embed_matrices(fmpz_mod_mat_t embed, fmpz_mod_mat_t project, const fq_t gen_sub, const fq_ctx_t sub_ctx, const fq_t gen_sup, const fq_ctx_t sup_ctx, const fmpz_mod_poly_t gen_minpoly)

Given the output of `fq_embed_gens` for fields $k$ = `sub_ctx` and $K$ =
`sup_ctx` of degrees $m \mid n$, computes the $n \times m$ matrix
`embed` of the embedding $\varphi: k \to K$ that maps `gen_sub` to
`gen_sup`, and an $m \times n$ matrix `project` of a $k$-linear
projection $\pi: K \to k$ with $\pi \circ \varphi = \mathrm{id}$, so that
$\mathit{project} \cdot \mathit{embed}$ is the $m \times m$ identity.
The matrices act on column vectors of coordinates in the monomial bases.

::: details
The projection is $\pi(y) = \varphi^{-1}(\operatorname{Tr}_{K/\varphi(k)}(\theta y))$ for an element $\theta$
with $\operatorname{Tr}_{K/\varphi(k)}(\theta) = 1$: $\theta = 1/(n/m)$
if $p \nmid n/m$, and otherwise $\theta = x^i / \operatorname{Tr}(x^i)$
for a suitable power of the generator. The trace is computed with
`fq_embed_trace_matrix`.
:::
:::

::: function
void fq_embed_trace_matrix(fmpz_mod_mat_t res, const fmpz_mod_mat_t basis, const fq_ctx_t sub_ctx, const fq_ctx_t sup_ctx)

Given contexts for $k$ = `sub_ctx` and $K$ = `sup_ctx` of degrees
$m \mid n$ and the $n \times m$ matrix `basis` of an embedding
$\varphi: k \to K$, sets `res` to the $m \times n$ matrix of the map
$y \mapsto \varphi^{-1}(\operatorname{Tr}_{K/\varphi(k)}(y))$. It is
computed as $D_k \cdot \mathit{basis}^T \cdot M_K$, with $D_k$ =
`fq_embed_dual_to_mono_matrix` for $k$ and $M_K$ =
`fq_embed_mono_to_dual_matrix` for $K$. For $m = n$ the result is the
inverse of `basis`.
:::

## Matrices of elements {#matrices}

::: function
void fq_embed_composition_matrix(fmpz_mod_mat_t matrix, const fq_t gen, const fq_ctx_t ctx)
void fq_embed_composition_matrix_sub(fmpz_mod_mat_t matrix, const fq_t gen, const fq_ctx_t ctx, slong trunc)

Set `matrix` to the composition matrix of the element $a$ = `gen` of a
field of degree $d$: the $d \times d$ matrix whose columns are the
coordinates of $a^0, a^1, \ldots, a^{d-1}$, which maps the coordinates of
a polynomial $g$ of degree less than $d$ to those of $g(a)$. The `_sub`
version computes only the first `trunc` $\le d$ columns, into a
$d \times \mathit{trunc}$ matrix.
:::

::: function
void fq_embed_mul_matrix(fmpz_mod_mat_t matrix, const fq_t gen, const fq_ctx_t ctx)

Sets the $d \times d$ matrix `matrix` to the matrix of multiplication by
$a$ = `gen`, whose columns are the coordinates of
$a, ax, \ldots, ax^{d-1}$. The defining polynomial need not be monic.
:::

::: function
void fq_embed_mono_to_dual_matrix(fmpz_mod_mat_t res, const fq_ctx_t ctx)

Sets the $d \times d$ matrix `res` to the change of basis from the
monomial basis $1, x, \ldots, x^{d-1}$ to its dual basis with respect to
the trace form $(u, v) \mapsto \operatorname{Tr}(uv)$: the Hankel matrix
with entries $\operatorname{Tr}(x^{i+j})$, computed from the power series
expansion of $\operatorname{rev}(f')/\operatorname{rev}(f)$ (Newton's
identities). Mapping a coordinate vector of $y$ by this matrix gives the
traces $\operatorname{Tr}(y x^i)$.
:::

::: function
void fq_embed_dual_to_mono_matrix(fmpz_mod_mat_t res, const fq_ctx_t ctx)

Sets the $d \times d$ matrix `res` to the inverse of
`fq_embed_mono_to_dual_matrix`. By Euler's formula the dual basis of
$1, x, \ldots, x^{d-1}$ is $b_0 / f'(x), \ldots, b_{d-1} / f'(x)$, where
$f(X)/(X - x) = \sum b_j X^j$, so the matrix is the product of the
multiplication matrix of $1/f'(x)$ and the Hankel matrix of the
coefficients of $f$.
:::

::: function
void fq_modulus_pow_series_inv(fmpz_mod_poly_t res, const fq_ctx_t ctx, slong trunc)

Sets `res` to the inverse of the reversal of the defining polynomial, as
a power series modulo $x^{\mathit{trunc}}$.
:::

::: function
void fq_modulus_derivative_inv(fq_t m_prime, fq_t m_prime_inv, const fq_ctx_t ctx)

Sets `m_prime` to $f'(x)$, the derivative of the defining polynomial
evaluated at the generator, and `m_prime_inv` to its inverse. $f'(x) \ne
0$ because $f$ is irreducible, hence separable.
:::

## Example {#example}

The program embeds $\mathbb{F}_{5^2}$ into $\mathbb{F}_{5^6}$, both
defined by Conway polynomials, checks that the projection is a left
inverse of the embedding, and checks on one element that the embedding
is a ring homomorphism.

```c example name=fq_embed_gf25
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mod_poly.h>
#include <flint/fmpz_mod_mat.h>
#include <flint/fq.h>
#include <flint/fq_embed.h>

/* the image of a in the larger field, as the matrix product M * a */
void embed_elem(fq_t b, const fq_t a, const fmpz_mod_mat_t M,
                const fq_ctx_t k, const fq_ctx_t K)
{
    fmpz_mod_mat_t u, v;
    fmpz_mod_mat_init(u, fq_ctx_degree(k), 1, k->ctxp);
    fmpz_mod_mat_init(v, fq_ctx_degree(K), 1, k->ctxp);
    fq_get_fmpz_mod_mat(u, a, k);
    fmpz_mod_mat_mul(v, M, u, k->ctxp);
    fq_set_fmpz_mod_mat(b, v, K);
    fmpz_mod_mat_clear(u, k->ctxp);
    fmpz_mod_mat_clear(v, k->ctxp);
}

int main(void)
{
    fmpz_t p;
    fq_ctx_t k, K;
    fq_t g, G, a, b, c, d;
    fmpz_mod_poly_t minpoly;
    fmpz_mod_mat_t emb, proj, prod;

    fmpz_init_set_ui(p, 5);
    fq_ctx_init(k, p, 2, "a");                  /* GF(5^2) */
    fq_ctx_init(K, p, 6, "b");                  /* GF(5^6) */
    fq_init(g, k); fq_init(a, k); fq_init(c, k);
    fq_init(G, K); fq_init(b, K); fq_init(d, K);
    fmpz_mod_poly_init(minpoly, k->ctxp);
    fmpz_mod_mat_init(emb, 6, 2, k->ctxp);
    fmpz_mod_mat_init(proj, 2, 6, k->ctxp);
    fmpz_mod_mat_init(prod, 2, 2, k->ctxp);

    fq_embed_gens(g, G, minpoly, k, K);
    flint_printf("minpoly: ");
    fmpz_mod_poly_print_pretty(minpoly, "x", k->ctxp);
    flint_printf("\nimage of a: "); fq_print_pretty(G, K); flint_printf("\n");

    fq_embed_matrices(emb, proj, g, k, G, K, minpoly);
    fmpz_mod_mat_mul(prod, proj, emb, k->ctxp);
    flint_printf("project * embed = 1: %d\n", fmpz_mod_mat_is_one(prod, k->ctxp));

    /* phi(a^2 + 1) = phi(a)^2 + 1 */
    fq_gen(a, k);
    fq_sqr(c, a, k);
    fq_one(a, k);
    fq_add(c, c, a, k);
    embed_elem(b, c, emb, k, K);
    fq_sqr(d, G, K);
    fq_one(G, K);
    fq_add(d, d, G, K);
    flint_printf("homomorphism: %d\n", fq_equal(b, d, K));

    fmpz_mod_mat_clear(emb, k->ctxp);
    fmpz_mod_mat_clear(proj, k->ctxp);
    fmpz_mod_mat_clear(prod, k->ctxp);
    fmpz_mod_poly_clear(minpoly, k->ctxp);
    fq_clear(g, k); fq_clear(a, k); fq_clear(c, k);
    fq_clear(G, K); fq_clear(b, K); fq_clear(d, K);
    fq_ctx_clear(k);
    fq_ctx_clear(K);
    fmpz_clear(p);
    return 0;
}
```
```output
minpoly: x^2+4*x+2
image of a: b^5+3*b^4+b^2+4*b+1
project * embed = 1: 1
homomorphism: 1
```
