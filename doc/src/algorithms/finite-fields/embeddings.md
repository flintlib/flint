---
title: Embeddings of finite fields
summary: To embed $\mathbb{F}_{p^m}$ into $\mathbb{F}_{p^n}$ given by unrelated defining polynomials, FLINT finds a root of the smaller field's polynomial in the larger field and expresses the embedding, and a left inverse built from the trace and the dual basis, as matrices over $\mathbb{F}_p$.
complexity: One root finding of a degree-$m$ polynomial over $\mathbb{F}_{p^n}$, plus $O(n^3)$ operations in $\mathbb{F}_p$ for the matrices.
order: 40
---

Let $k = \mathbb{F}_p[y]/(f_m)$ and $K = \mathbb{F}_p[x]/(f_n)$ be fields
of degrees $m \mid n$. Then $K$ has a unique subfield of order $p^m$,
isomorphic to $k$, but an isomorphism has to be computed: it is
determined by the image of the generator $y$, which must be a root of
$f_m$ in $K$, and there are $m$ such roots, giving the $m$ embeddings
$k \to K$, which differ by powers of the Frobenius map. Lenstra showed
that isomorphisms of finite fields can be found in deterministic
polynomial time [@Len1991]; practical algorithms either find a root of
$f_m$ in $K$, as FLINT does, or use Kummer theory and Gauss periods to
construct compatible generators without factoring [@All2002]
[@DFRR2019]. Systems such as Magma maintain a lattice of compatible
embeddings between all fields created [@BCS1997b]; FLINT computes
embeddings on demand and leaves compatibility to the user (Conway
polynomials, [](algorithms/finite-fields/defining-polynomials), provide
compatible generators by construction, but `fq_embed` does not use this).

## Finding the embedding {#gens}

`fq_embed_gens` returns the generator $y$ of $k$, its minimal polynomial
$f_m$, and a root $\beta$ of $f_m$ in $K$. The root is found by
`fq_poly_factor_split_single`, that is, by the root finding of
`gr_poly_roots_finite_field` applied to $f_m$ viewed as a polynomial over
$K$ ([](algorithms/polynomials/root-finding)), which splits $f_m$ (a
product of $m$ distinct linear factors over $K$) by random equal-degree
splitting; its cost is dominated by modular exponentiations of
polynomials of degree less than $m$ over $K$ to exponents of size about
$p^n$. The embedding is
$\varphi(\sum c_i y^i) = \sum c_i \beta^i$.

## The embedding matrix {#matrix}

With respect to the bases $1, y, \ldots, y^{m-1}$ of $k$ and
$1, x, \ldots, x^{n-1}$ of $K$, $\varphi$ is the $n \times m$ matrix
whose columns are the coordinates of $\beta^0, \ldots, \beta^{m-1}$, the
first $m$ columns of the *composition matrix* of $\beta$
(`fq_embed_composition_matrix_sub`). More generally, `fq_embed_matrices`
accepts any generator $\gamma$ of $k$ with minimal polynomial $g$ and a
root $\gamma'$ of $g$ in $K$; it then composes the matrix of
$\mathbb{F}_p[z]/(g) \to K$, $z \mapsto \gamma'$, with the inverse of the
matrix of $\mathbb{F}_p[z]/(g) \to k$, $z \mapsto \gamma$.

## Projecting back with the trace {#projection}

A map in the other direction, $\pi: K \to k$ with $\pi \circ \varphi =
\mathrm{id}$, is needed to test whether an element of $K$ lies in the
image and to compute its preimage. FLINT obtains it from the relative
trace $\operatorname{Tr}_{K/k'}$, where $k' = \varphi(k)$, which is
$k'$-linear and surjective.

::: lemma {#lem-dual}
Let $L = \mathbb{F}_p[x]/(f)$ with $\deg f = d$, and let $T$ be the
$d \times d$ matrix with entries $T_{ij} = \operatorname{Tr}(x^{i+j})$.
Then $T$ is invertible, and for $u \in L$ with coordinate vector
$\mathbf{u}$, $T\mathbf{u} = (\operatorname{Tr}(u x^i))_{i < d}$.
Moreover, if $f(X)/(X - x) = \sum_{j < d} b_j X^j$ in $L[X]$, then the
basis dual to $1, x, \ldots, x^{d-1}$ for the trace form is
$b_0/f'(x), \ldots, b_{d-1}/f'(x)$.
:::

::: proof
The trace form $(u, v) \mapsto \operatorname{Tr}(uv)$ is nondegenerate
because $L/\mathbb{F}_p$ is separable, so its Gram matrix $T$ in the
monomial basis is invertible, and $(T\mathbf{u})_i = \sum_j
\operatorname{Tr}(x^{i+j}) u_j = \operatorname{Tr}(x^i u)$. The second
statement is Euler's formula for the dual basis [@LN1997]
(a consequence of Lagrange interpolation at the conjugates of $x$).
:::

The matrix $T$ is `fq_embed_mono_to_dual_matrix`, computed from the power
sums $\operatorname{Tr}(x^i)$ of
[](algorithms/finite-fields/frobenius-norm-trace); its inverse,
`fq_embed_dual_to_mono_matrix`, is the product of the multiplication
matrix of $1/f'(x)$ (`fq_embed_mul_matrix`) and the Hankel matrix of the
coefficients of $f$, by the second part of the lemma.

::: proposition {#prop-trace}
Let $B$ be the $n \times m$ matrix of an embedding $\varphi: k \to K$,
$T_K$ the matrix of the trace form of $K$ and $D_k$ the inverse of that
of $k$. Then $D_k B^T T_K$ is the matrix of
$\varphi^{-1} \circ \operatorname{Tr}_{K/k'}$.
:::

::: proof
For $u \in K$, $T_K \mathbf{u} = (\operatorname{Tr}_{K/\mathbb{F}_p}(u
x^i))_i$, so $B^T T_K \mathbf{u} = (\operatorname{Tr}_{K/\mathbb{F}_p}(u
\varphi(y^j)))_{j < m}$, because column $j$ of $B$ is the coordinate
vector of $\varphi(y^j)$. By transitivity of the trace and
$k'$-linearity of $\operatorname{Tr}_{K/k'}$,
$\operatorname{Tr}_{K/\mathbb{F}_p}(u \varphi(y^j)) =
\operatorname{Tr}_{k'/\mathbb{F}_p}(\varphi(y^j)
\operatorname{Tr}_{K/k'}(u)) = \operatorname{Tr}_{k/\mathbb{F}_p}(y^j
w)$ with $w = \varphi^{-1}(\operatorname{Tr}_{K/k'}(u))$. These are the
dual coordinates of $w$ in $k$, and $D_k$ converts them to monomial
coordinates by [](#lem-dual).
:::

This is `fq_embed_trace_matrix`. Since $\operatorname{Tr}_{K/k'}$ is the
multiplication by $n/m$ on $k'$, the projection is
$\pi(u) = \varphi^{-1}(\operatorname{Tr}_{K/k'}(\theta u))$ for any
$\theta \in K$ with $\operatorname{Tr}_{K/k'}(\theta) = 1$:
$\theta = (n/m)^{-1}$ if $p \nmid n/m$, and otherwise
$\theta = x^i / \operatorname{Tr}_{K/k'}(x^i)$ for a power of the
generator with nonzero relative trace, which `fq_embed_matrices` finds by
looking for a nonzero entry in the first row of the trace matrix. Then
$\pi(\varphi(c)) = \varphi^{-1}(\varphi(c) \operatorname{Tr}_{K/k'}(\theta))
= c$. Using duality in this way, rather than solving a linear system for
each element, follows Shoup's approach to computations with linear forms
in extensions of finite fields [@Sho1994] [@Sho1999].

## Use in FLINT {#use}

- `fq_embed_gens`, `_fq_embed_gens_naive`: the generators.
- `fq_embed_matrices`: the embedding and projection matrices, from
  `fq_embed_composition_matrix_sub`, `fq_embed_trace_matrix`,
  `fq_embed_mono_to_dual_matrix`, `fq_embed_dual_to_mono_matrix`,
  `fq_embed_mul_matrix`, `fq_modulus_pow_series_inv` and
  `fq_modulus_derivative_inv`.
- The same functions exist for `fq_nmod` and `fq_zech` (with `nmod_mat`
  matrices); the `fq_zech` versions convert elements to the `fq_nmod`
  representation.
