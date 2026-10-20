---
title: Multipoint evaluation and interpolation
functions: nmod_poly_evaluate_nmod_vec nmod_poly_evaluate_nmod_vec_fast nmod_poly_interpolate_nmod_vec nmod_poly_interpolate_nmod_vec_fast _nmod_poly_tree_build
---

Evaluating a polynomial of length $m$ at $m$ points with Horner's rule
costs $m^2$ multiplications, and so does Lagrange interpolation through
$m$ points. With a *subproduct tree*, both cost $O(M(m) \log m)$
operations [@GG2003, §10]. For points in geometric progression,
$O(M(m))$ is possible.

::: summary algorithms/polynomials/subproduct-tree :::

## Subproduct trees {#tree}

For points $x_0, \ldots, x_{m-1}$, the leaves of the tree are the linear
polynomials $x - x_i$ and each node holds the product of its two
children. `_nmod_poly_tree_build` computes all levels except the root,
the full product $P = \prod (x - x_i)$, which is not needed.

- *Evaluation* reduces $f$ modulo the products at a high level of the
  tree and then the remainders down the tree; the remainder modulo $x -
  x_i$ is $f(x_i)$. If $f$ is longer than the number of points, it is
  first reduced modulo the largest products.
- *Interpolation* computes the weights $w_i = 1 / P'(x_i)$ by evaluating
  $P'$ with the same tree (`_nmod_poly_interpolation_weights`), and then
  combines the scaled values $y_i w_i$ up the tree: a node with children
  $L, R$ and partial results $a, b$ gets $a R + b L$. The root gives the
  Lagrange interpolant $\sum_i y_i w_i P / (x - x_i)$.

Building the tree is a substantial part of the cost, so when many polynomials are evaluated at the same points, build it once with
`_nmod_poly_tree_alloc` and `_nmod_poly_tree_build` and call
`_nmod_poly_evaluate_nmod_vec_fast_precomp` and
`_nmod_poly_interpolate_nmod_vec_fast_precomp`.

`nmod_poly_evaluate_nmod_vec` evaluates point by point for polynomials of
length less than 32 and uses the tree otherwise. Point-by-point evaluation
uses `_nmod_poly_evaluate_nmod`, which for long polynomials uses
rectangular splitting with dot products and is fast in practice, so the
tree pays off only from a few hundred points on. `nmod_poly_interpolate_nmod_vec`
uses the Newton form (divided differences) for fewer than 6 points, the
barycentric Lagrange formula for fewer than 16 and the tree otherwise.
Interpolation requires distinct points with invertible differences, which
for prime $n$ just means distinct points.

## Geometric progressions {#geometric}

If the points are $1, q, q^2, \ldots, q^{m-1}$, evaluation is a single
product: from $ij = \binom{i+j}{2} - \binom{i}{2} - \binom{j}{2}$,
$$
f(q^j) = q^{-\binom{j}{2}} \sum_i \bigl(f_i q^{-\binom{i}{2}}\bigr)
q^{\binom{i+j}{2}},
$$
a middle product of the scaled coefficients with the fixed polynomial
$\sum_k q^{\binom{k}{2}} x^k$ (Bluestein's chirp transform [@Blu1970]
[@RSR1969]). FLINT uses the equivalent form with $q^{i^2/2} = r^{i^2}$,
which is why the functions take a value $r$ and evaluate at the powers of
$q = r^2$. Interpolation at these points is two products [@BS2005]: the
values are converted to the Newton basis of the progression and then to
the monomial basis. *Extrapolation* computes the values of a polynomial
at the next points of the progression from its values at consecutive
points, with one middle product and without computing the polynomial.

The data that depend only on $r$ and the number of points can be stored
in an `nmod_geometric_progression_t`, which also accepts the
`fft_small` transforms of the fixed polynomials, so that each product
costs two transforms instead of three (for evaluation and extrapolation,
transposed transforms of the length of the fixed polynomial; see
[](algorithms/polynomials/subproduct-tree#geometric)). For prime $n$, the points are
distinct if the multiplicative order of $q$ is at least $m$; a primitive
root modulo $n$, or a root of unity of large order for FFT-friendly
primes, is a safe choice.

## Examples {#examples}

Evaluating a polynomial of length 1000 at 1000 points with and without a
subproduct tree, and interpolating back:

```c example name=nmod_poly_multipoint
#include <flint/flint.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_poly.h>

int main(void)
{
    flint_rand_t state;
    nmod_poly_t f, g;
    nn_ptr xs, ys, zs;
    slong i, m = 1000;
    ulong p = 1000003;

    flint_rand_init(state);
    nmod_poly_init(f, p);
    nmod_poly_init(g, p);
    xs = _nmod_vec_init(m);
    ys = _nmod_vec_init(m);
    zs = _nmod_vec_init(m);

    nmod_poly_rand(f, state, m);
    for (i = 0; i < m; i++)
        xs[i] = 3 * i + 1;                 /* distinct points */

    nmod_poly_evaluate_nmod_vec_fast(ys, f, xs, m);
    nmod_poly_evaluate_nmod_vec_iter(zs, f, xs, m);
    flint_printf("fast == iter: %d\n", _nmod_vec_equal(ys, zs, m));

    nmod_poly_interpolate_nmod_vec(g, xs, ys, m);
    flint_printf("interpolant == f: %d\n", nmod_poly_equal(f, g));

    /* the product of the x - x_i vanishes at every point */
    nmod_poly_product_roots_nmod_vec(g, xs, m);
    nmod_poly_evaluate_nmod_vec(zs, g, xs, m);
    flint_printf("P(x_i) all zero: %d\n", _nmod_vec_is_zero(zs, m));

    _nmod_vec_clear(xs);
    _nmod_vec_clear(ys);
    _nmod_vec_clear(zs);
    nmod_poly_clear(f);
    nmod_poly_clear(g);
    flint_rand_clear(state);
    return 0;
}
```
```output
fast == iter: 1
interpolant == f: 1
P(x_i) all zero: 1
```

The same at the points $3^{2i}$ modulo the prime $p = 7 \cdot 2^{26} + 1$,
for which 3 is a primitive root, so that the first $p - 1$ powers of $3$
are distinct. One precomputation serves evaluation, interpolation and
extrapolation:

```c example name=nmod_poly_geometric
#include <flint/flint.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_poly.h>

int main(void)
{
    flint_rand_t state;
    nmod_t mod;
    nmod_poly_t f, g;
    nmod_geometric_progression_t G;
    nn_ptr ys, zs;
    slong m = 1000;
    ulong r = 3;

    flint_rand_init(state);
    nmod_init(&mod, 469762049);
    nmod_poly_init_mod(f, mod);
    nmod_poly_init_mod(g, mod);
    ys = _nmod_vec_init(2 * m);
    zs = _nmod_vec_init(2 * m);

    nmod_poly_rand(f, state, m);
    nmod_geometric_progression_init(G, r, 2 * m, mod);

    /* values at the 2m points 1, r^2, r^4, ... */
    _nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(ys, f->coeffs, f->length, G, 2 * m, mod);
    nmod_poly_evaluate_geometric_nmod_vec_iter(zs, f, r, 2 * m);
    flint_printf("fast == iter: %d\n", _nmod_vec_equal(ys, zs, 2 * m));

    /* interpolate from the first m values */
    nmod_poly_interpolate_geometric_nmod_vec_fast_precomp(g, ys, G, m);
    flint_printf("interpolant == f: %d\n", nmod_poly_equal(f, g));

    /* the next m values from the first m, without f */
    nmod_poly_extrapolate_geometric_precomp(zs, m, ys, m, m, G);
    flint_printf("extrapolated == values: %d\n", _nmod_vec_equal(zs, ys + m, m));

    nmod_geometric_progression_clear(G);
    _nmod_vec_clear(ys);
    _nmod_vec_clear(zs);
    nmod_poly_clear(f);
    nmod_poly_clear(g);
    flint_rand_clear(state);
    return 0;
}
```
```output
fast == iter: 1
interpolant == f: 1
extrapolated == values: 1
```
