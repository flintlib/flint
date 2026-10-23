---
title: Computing modulo a fixed polynomial
summary: How to compute in (Z/nZ)[x]/(f) efficiently, with precomputed inverses of the reversed modulus, preconditioned multiplication, fast exponentiation and modular composition; with finite field arithmetic as the example.
order: 2
---

Many algorithms compute in the quotient ring $(\Z/n\Z)[x]/(f)$: arithmetic
in a finite field $\mathbb{F}_{p^d} = \mathbb{F}_p[x]/(f)$ for an
irreducible $f$ of degree $d$, distinct-degree factorisation, primality
tests for polynomials, and the computation of $x^{p^i} \bmod f$ in
general. An element of the quotient ring is represented by its remainder
modulo $f$, a polynomial of degree less than $d = \deg f$, and the basic
operation is a product followed by a reduction modulo $f$. The leading coefficient of $f$ must be invertible; for
a field, $n$ must be prime and $f$ irreducible.

## Reduction with a precomputed inverse {#preinv}

A remainder modulo $f$ can be computed with `nmod_poly_rem`, which for
large $d$ uses [Newton iteration](algorithms/polynomials/newton) to find
the quotient. That iteration computes an inverse of the reversed divisor
each time. When the same $f$ is used many times, this inverse should be
computed once. Write $\operatorname{rev}_k(a) = x^{k-1} a(1/x)$ for the
reversal of a polynomial viewed with length $k$. If $A = Qf + R$ with
$\deg R < d$ and $A$ of length $\ell$, then reversing gives
$$
\operatorname{rev}_\ell(A) = \operatorname{rev}_{\ell - d}(Q) \cdot
\operatorname{rev}_{d+1}(f) + x^{\ell - d} \operatorname{rev}_d(R),
$$
so the quotient is
$$
\operatorname{rev}_{\ell - d}(Q) = \operatorname{rev}_\ell(A) \cdot
\operatorname{rev}_{d+1}(f)^{-1} \bmod x^{\ell - d},
$$
a truncated product with the power series inverse of the reversal of $f$,
and the remainder is $R = A - Qf$, of which only the low $d$ coefficients
need to be computed. The reversal has constant term $\operatorname{lc}(f)$,
which is why the leading coefficient must be invertible.

FLINT calls this inverse `finv` (or `Binv`, `poly3inv`, ...) and computes
it with

```c
nmod_poly_reverse(finv, f, f->length);
nmod_poly_inv_series(finv, finv, f->length);
```

A product of two reduced elements has length at most $2d - 1$, so the
quotient has length at most $d$ and `finv` must be correct modulo $x^d$;
computing it modulo $x^{d+1}$ as above is the convention. Every function
with `preinv` in its name expects this `finv`: `nmod_poly_mulmod_preinv`,
`nmod_poly_divrem_newton_n_preinv`, `nmod_poly_div_newton_n_preinv`, the
`powmod` functions with `preinv`, `nmod_poly_powers_mod_bsgs` (which
computes it itself) and the Brent–Kung modular composition functions. A
reduction then costs two products of length about $d$, against the
inversion and two products of `nmod_poly_rem`.

For moduli with at most five nonzero coefficients below the leading one,
such as the trinomials and pentanomials returned by
`nmod_poly_minimal_irreducible`, `_nmod_poly_divrem_newton_n_preinv`
detects the sparsity (for $d > 20$) and reduces directly, with a few
multiply-adds per coefficient of the dividend and no products.

## Which function to use {#choice}

| Task | Function |
|---|---|
| One product, small $d$ | `nmod_poly_mulmod` |
| Many products modulo $f$ | `nmod_poly_mulmod_preinv` |
| Many products by the same $a$ | `nmod_poly_mulmod_precond` |
| $a^e \bmod f$ | `nmod_poly_powmod_ui_binexp_preinv`, `nmod_poly_powmod_fmpz_binexp_preinv` |
| $x^e \bmod f$ | `nmod_poly_powmod_x_ui_preinv`, `nmod_poly_powmod_x_fmpz_preinv` |
| $a^0, \ldots, a^{k-1} \bmod f$ | `nmod_poly_powers_mod_bsgs` |
| $a^{-1} \bmod f$ | `nmod_poly_invmod` |
| $g(a) \bmod f$ | `nmod_poly_compose_mod`, or with a precomputed matrix `nmod_poly_compose_mod_brent_kung_precomp_preinv` |

The `preinv` and `precond` functions require their inputs to be reduced
modulo $f$; the plain `mulmod` and `powmod` functions reduce them first.
`nmod_poly_mulmod_precond` precomputes either the matrix of multiplication
by $a$ (fastest for $d$ up to a few hundred) or the quotient $\lfloor a
x^d / f \rfloor$ (Shoup's method), and
`nmod_poly_mulmod_precond_init_num` chooses between them from the expected
number of multiplications. The precomputed object keeps pointers to $a$,
$f$ and `finv`, which must not change while it is used.

The function $x \mapsto x^p$ is the Frobenius automorphism of
$\mathbb{F}_{p^d}$, and $x^{p^i} \bmod f$ for successive $i$ is the core
of distinct-degree factorisation. Once $h = x^p \bmod f$ is known,
$x^{p^{i+1}} = (x^{p^i})^p = h(x^{p^i}) \bmod f$ is a modular composition,
which [Brent and Kung's algorithm](algorithms/polynomials/composition)
computes in $O(d^{1/2} M(d) + d^{(\omega+1)/2})$ operations instead of the
$O(M(d) \log p)$ of an exponentiation [@GatSho1992] [@KalSho1998]. The
matrix of powers of $h$ can be precomputed with
`nmod_poly_precompute_matrix` and reused for every step; see
[](nmod_poly/compose_mod).

## Example: arithmetic in a finite field {#example}

The program constructs $\mathbb{F}_{101^8}$ with the sparsest irreducible
polynomial of degree 8 and checks some identities: the Frobenius
automorphism has order 8, the multiplicative group has order $p^8 - 1$,
inverses multiply to 1, and 1000 multiplications by a fixed element with a
preconditioned object agree with one exponentiation.

```c example name=nmod_poly_finite_field
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/nmod_poly.h>

int main(void)
{
    nmod_poly_t f, finv, a, b, t;
    nmod_poly_mulmod_precond_t pre;
    fmpz_t e;
    slong i;

    nmod_poly_init(f, 101);
    nmod_poly_init_mod(finv, f->mod);
    nmod_poly_init_mod(a, f->mod);
    nmod_poly_init_mod(b, f->mod);
    nmod_poly_init_mod(t, f->mod);
    fmpz_init(e);

    nmod_poly_minimal_irreducible(f, 8);
    flint_printf("f = "); nmod_poly_print_pretty(f, "x"); flint_printf("\n");

    /* finv = 1 / rev(f) mod x^(d+1) */
    nmod_poly_reverse(finv, f, f->length);
    nmod_poly_inv_series(finv, finv, f->length);

    /* x^(p^8) = x in GF(p^8), but x^(p^4) != x */
    fmpz_ui_pow_ui(e, 101, 8);
    nmod_poly_powmod_x_fmpz_preinv(t, e, f, finv);
    flint_printf("x^(p^8) == x: %d\n", nmod_poly_is_gen(t));
    fmpz_ui_pow_ui(e, 101, 4);
    nmod_poly_powmod_x_fmpz_preinv(t, e, f, finv);
    flint_printf("x^(p^4) == x: %d\n", nmod_poly_is_gen(t));

    /* a^(p^8 - 1) = 1 for a = x + 2 */
    nmod_poly_set_coeff_ui(a, 1, 1);
    nmod_poly_set_coeff_ui(a, 0, 2);
    fmpz_ui_pow_ui(e, 101, 8);
    fmpz_sub_ui(e, e, 1);
    nmod_poly_powmod_fmpz_binexp_preinv(t, a, e, f, finv);
    flint_printf("a^(p^8 - 1) == 1: %d\n", nmod_poly_is_one(t));

    /* inverse */
    nmod_poly_invmod(b, a, f);
    nmod_poly_mulmod_preinv(t, a, b, f, finv);
    flint_printf("a * a^-1 == 1: %d\n", nmod_poly_is_one(t));

    /* 1000 multiplications by a with a preconditioned object */
    nmod_poly_mulmod_precond_init_num(pre, a, f, finv, 1000);
    nmod_poly_one(b);
    for (i = 0; i < 1000; i++)
        nmod_poly_mulmod_precond(b, pre, b);
    nmod_poly_mulmod_precond_clear(pre);
    nmod_poly_powmod_ui_binexp_preinv(t, a, 1000, f, finv);
    flint_printf("a^1000 agrees: %d\n", nmod_poly_equal(b, t));

    fmpz_clear(e);
    nmod_poly_clear(f);
    nmod_poly_clear(finv);
    nmod_poly_clear(a);
    nmod_poly_clear(b);
    nmod_poly_clear(t);
    return 0;
}
```
```output
f = x^8+2
x^(p^8) == x: 1
x^(p^4) == x: 0
a^(p^8 - 1) == 1: 1
a * a^-1 == 1: 1
a^1000 agrees: 1
```

The modulus is a binomial: an irreducible $x^8 + c$ exists because the
only prime factor 2 of 8 divides $p - 1 = 100$ and, as $4 \mid 8$, also
$4 \mid p - 1$. Reduction modulo $x^8 + 2$ costs one multiplication per
coefficient, and products modulo $f$ are then dominated by the
multiplication itself.
