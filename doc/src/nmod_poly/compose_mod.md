---
title: Modular composition
functions: nmod_poly_compose_mod nmod_poly_compose_mod_horner nmod_poly_compose_mod_brent_kung nmod_poly_compose_mod_brent_kung_preinv nmod_poly_precompute_matrix nmod_poly_compose_mod_brent_kung_precomp_preinv
---

Modular composition computes $f(g) \bmod h$ for polynomials $f, g, h$
with $\deg f, \deg g < \deg h = n$. It is the key operation of polynomial
factorisation over finite fields, where it computes the iterates
$x^{p^i} \bmod h$ of the Frobenius map, and of the computation of minimal
polynomials and of field embeddings. Horner's rule needs $n - 1$
multiplications modulo $h$, that is, $O(n M(n))$ operations; the
algorithm of Brent and Kung [@BrentKung1978] reduces this to $O(\sqrt n
M(n) + n^{(\omega + 1)/2})$, where $\omega$ is the exponent of matrix
multiplication.

::: summary algorithms/polynomials/composition :::

## The Brent–Kung algorithm {#brent-kung}

Let $m = \lfloor \sqrt n \rfloor + 1$ and write $f = \sum_{j} f_j(x)\,
x^{jm}$ with blocks $f_j$ of length $m$. Then
$$
f(g) = \sum_j f_j(g)\, (g^m)^j \bmod h, \qquad f_j(g) = \sum_{i < m}
f_{j,i}\, g^i.
$$
The baby steps $g^0, \ldots, g^{m-1} \bmod h$ form the rows of an $m
\times n$ matrix $A$, and all $f_j(g)$ together are the rows of the
product $C = BA$, where row $j$ of the $m \times m$ matrix $B$ holds the
coefficients of $f_j$. The giant steps evaluate $\sum_j C_j H^j$ for $H
= g^m \bmod h$ by Horner's rule or rectangular splitting
(`_nmod_poly_mod_matrix_rows_evaluate`), with about $\sqrt n$
multiplications modulo $h$. The matrix product is done by `nmod_mat_mul`
(see [](algorithms/linear-algebra/matrix-multiplication)).

The matrix $A$ depends only on $g$ and $h$. `nmod_poly_precompute_matrix`
computes it once (it must be initialised with $m$ rows and $n$ columns),
and `nmod_poly_compose_mod_brent_kung_precomp_preinv` then composes any
$f$ with $g$ at the cost of the matrix product and the giant steps. The
`preinv` versions take the inverse of the reversal of $h$, as explained in
[](nmod_poly/modular-arithmetic).

`nmod_poly_compose_mod` uses Horner's rule for $\deg h < 7$ and when
$\deg f \ge \deg h$ (Brent–Kung requires $\deg f < \deg h$), and the
Brent–Kung algorithm otherwise.

## Example: Rabin's irreducibility test {#example}

A polynomial $h$ of degree $n$ over $\mathbb{F}_p$ is irreducible if and
only if $x^{p^n} \equiv x \pmod h$ and $\gcd(x^{p^{n/q}} - x, h) = 1$ for
every prime $q$ dividing $n$ [@Rab1980b]. The powers $x^{p^i}$ are
computed by one exponentiation, $x^p \bmod h$, followed by repeated
composition with $x^p$, using one precomputed matrix:

```c example name=nmod_poly_rabin_irreducible
#include <flint/flint.h>
#include <flint/ulong_extras.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_factor.h>

/* x^(p^i) mod h for i = 0, ..., n, by composition with x^p */
int is_irreducible_rabin(const nmod_poly_t h)
{
    slong n = nmod_poly_degree(h), i, m = n_sqrt(n) + 1;
    ulong p = h->mod.n;
    nmod_poly_t hinv, xp, t, g;
    nmod_mat_t A;
    int result = 1;

    nmod_poly_init_mod(hinv, h->mod);
    nmod_poly_init_mod(xp, h->mod);
    nmod_poly_init_mod(t, h->mod);
    nmod_poly_init_mod(g, h->mod);
    nmod_mat_init(A, m, n, p);

    nmod_poly_reverse(hinv, h, h->length);
    nmod_poly_inv_series(hinv, hinv, h->length);

    nmod_poly_powmod_x_ui_preinv(xp, p, h, hinv);       /* x^p mod h */
    nmod_poly_precompute_matrix(A, xp, h, hinv);

    nmod_poly_set(t, xp);                               /* t = x^(p^1) */
    for (i = 1; i <= n; i++)
    {
        /* gcd test for the maximal proper divisors i = n/q of n */
        if (i < n && n % i == 0 && n_is_prime(n / i))
        {
            nmod_poly_set_coeff_ui(g, 1, 1);
            nmod_poly_sub(g, t, g);                     /* x^(p^i) - x */
            nmod_poly_gcd(g, g, h);
            if (!nmod_poly_is_one(g))
                result = 0;
            nmod_poly_zero(g);
        }
        if (i < n)                                      /* t = t(x^p) */
            nmod_poly_compose_mod_brent_kung_precomp_preinv(t, t, A, h, hinv);
    }
    if (!nmod_poly_is_gen(t))                           /* x^(p^n) == x */
        result = 0;

    nmod_mat_clear(A);
    nmod_poly_clear(hinv);
    nmod_poly_clear(xp);
    nmod_poly_clear(t);
    nmod_poly_clear(g);
    return result;
}

int main(void)
{
    flint_rand_t state;
    nmod_poly_t a, b, h;
    ulong p = 1000003;

    flint_rand_init(state);
    nmod_poly_init(a, p);
    nmod_poly_init(b, p);
    nmod_poly_init(h, p);

    nmod_poly_randtest_monic_irreducible(a, state, 61);   /* degree 60 */
    nmod_poly_randtest_monic_irreducible(b, state, 41);   /* degree 40 */

    flint_printf("a irreducible: %d (library: %d)\n",
                 is_irreducible_rabin(a), nmod_poly_is_irreducible(a));
    nmod_poly_mul(h, a, b);
    flint_printf("a*b irreducible: %d\n", is_irreducible_rabin(h));
    nmod_poly_mul(h, a, a);
    flint_printf("a^2 irreducible: %d\n", is_irreducible_rabin(h));

    nmod_poly_clear(a);
    nmod_poly_clear(b);
    nmod_poly_clear(h);
    flint_rand_clear(state);
    return 0;
}
```
```output
a irreducible: 1 (library: 1)
a*b irreducible: 0
a^2 irreducible: 0
```

`nmod_poly_is_irreducible` is declared in `nmod_poly_factor.h`. Each composition costs one $m \times m$ by $m
\times n$ matrix product and about $m$ modular multiplications, against
$\log_2 p \approx 20$ squarings modulo $h$ for an exponentiation, which is
where modular composition pays off for large $p$ [@GatSho1992]
[@KalSho1998].
