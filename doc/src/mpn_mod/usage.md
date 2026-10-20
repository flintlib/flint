---
title: Using mpn_mod
summary: How to create an `mpn_mod` ring, store and operate on its elements, use generic polynomials and matrices over it, flag a prime modulus, and recover factors of a composite modulus with the pretend-field mechanism.
order: 1
---

`mpn_mod` provides $\Z/n\Z$ for moduli of 2 to 16 limbs (from $2^{64}$
up to $2^{1024}$ on 64-bit machines) as a generic ring. There is no
`mpn_mod_poly` or `mpn_mod_mat` type: polynomials and matrices are the
generic `gr_poly_t` and `gr_mat_t` over an `mpn_mod` context, and the
generic functions dispatch to the specialised code of this module through
the ring's method table.

## Creating the ring {#ring}

`gr_ctx_init_mpn_mod` creates the ring from an `fmpz` modulus. It returns
a status code that must be checked: `GR_UNABLE` means that the modulus
has one limb or more than 16 limbs, and then another ring must be used.
A program that accepts any modulus typically chooses as follows.

| Size of $n$ | Ring | Element type |
|---|---|---|
| one limb | `gr_ctx_init_nmod` | `ulong` |
| 2 to 16 limbs | `gr_ctx_init_mpn_mod` | array of $\ell$ limbs |
| larger | `gr_ctx_init_fmpz_mod` | `fmpz` |

The context must be cleared with `gr_ctx_clear`. Creating one is cheap
(it computes a precomputed inverse of $n$), so multimodular algorithms
can create many.

The ring does not know whether $n$ is prime. If it is, say so with
`gr_ctx_set_is_field(ctx, T_TRUE)`: algorithms that are only valid over
fields, such as polynomial factorisation, root finding and square
roots, otherwise return `GR_UNABLE`. Do not set the flag for a composite
modulus.

## Elements {#elements}

An element is an array of $\ell$ = `MPN_MOD_CTX_NLIMBS(ctx)` limbs
holding a residue in $[0, n)$. It needs no initialisation or clearing,
so a local variable can simply be an array of `MPN_MOD_MAX_LIMBS` limbs,
or `GR_TMP_INIT` can be used in code that should work for every ring.
The functions `mpn_mod_set_si`, `mpn_mod_mul` and so on can be called
directly, or through `gr_set_si`, `gr_mul` and so on, which cost one
indirect call more. Conversions to and from integers are
`mpn_mod_set_fmpz` (which reduces any integer) and `mpn_mod_get_fmpz`
(which returns the residue in $[0, n)$).

The program below computes in $\Z/p\Z$ for the Mersenne prime
$p = 2^{127} - 1$: a division of scalars, and a gcd of two polynomials
with the generic polynomial functions. The gcd is monic and its constant
term is printed as a residue in $[0, p)$, here $-2 \bmod p$.

```c example name=mpn_mod_usage
#include <flint/fmpz.h>
#include <flint/gr.h>
#include <flint/gr_poly.h>
#include <flint/mpn_mod.h>

int main(void)
{
    gr_ctx_t ctx;
    fmpz_t p;
    gr_poly_t f, g, h;
    ulong x[MPN_MOD_MAX_LIMBS], y[MPN_MOD_MAX_LIMBS];
    int status = GR_SUCCESS;

    /* Z/pZ with the Mersenne prime p = 2^127 - 1 */
    fmpz_init(p);
    fmpz_ui_pow_ui(p, 2, 127);
    fmpz_sub_ui(p, p, 1);
    status |= gr_ctx_init_mpn_mod(ctx, p);
    status |= gr_ctx_set_is_field(ctx, T_TRUE);

    /* elements are plain limb arrays */
    status |= mpn_mod_set_si(x, -3, ctx);
    status |= mpn_mod_set_ui(y, 10, ctx);
    status |= mpn_mod_div(x, x, y, ctx);        /* x = -3/10 */
    status |= mpn_mod_mul(y, x, y, ctx);        /* y = -3 */
    status |= gr_println(y, ctx);

    /* generic polynomials over the same ring */
    gr_poly_init(f, ctx);
    gr_poly_init(g, ctx);
    gr_poly_init(h, ctx);
    status |= gr_poly_set_str(f, "(x-1)*(x-2)^2", "x", ctx);
    status |= gr_poly_set_str(g, "(x-2)*(x+5)", "x", ctx);
    status |= gr_poly_gcd(h, f, g, ctx);
    status |= gr_poly_print(h, ctx);
    printf("\n");

    gr_poly_clear(f, ctx);
    gr_poly_clear(g, ctx);
    gr_poly_clear(h, ctx);
    gr_ctx_clear(ctx);
    fmpz_clear(p);
    printf("status: %d\n", status);
    return 0;
}
```

```output
170141183460469231731687303715884105724
170141183460469231731687303715884105725 + x
status: 0
```

## Polynomials and matrices {#generic}

Every `gr_poly` and `gr_mat` function works over the ring. The ring
provides its own multiplication of polynomials and matrices, division
with remainder, power series inversion, polynomial gcd, LU
decomposition, triangular solving, determinants and characteristic
polynomials; everything else (Newton iteration, half-gcd, composition,
factorisation, linear solving, ranks, inverses, ...) is generic code
that calls these. The specialised code keeps residues in fixed-size limb
arrays and delays modular reductions: sums of products are accumulated
exactly in a few extra limbs and reduced once
([](algorithms/arithmetic/multiword-delayed-reduction)). See
[](mpn_mod/poly-mul) and [](mpn_mod/mat-mul) for the multiplication
algorithms.

`fmpz_mod_poly` and `fmpz_mod_mat` do not use this module; code that
wants its speed for moduli of 2 to 16 limbs works with `gr_poly_t` and
`gr_mat_t` directly. Conversion between the two element types is
`gr_set_other`, which copies directly between an `fmpz_mod` ring and an
`mpn_mod` ring with the same modulus.

## Composite moduli and zero divisors {#zero-divisors}

Over a composite modulus, algorithms that divide (Euclid's algorithm, LU
decomposition, series inversion) may meet a nonzero element that is not
invertible. By default, `mpn_mod_inv` then returns `GR_DOMAIN`, which
the algorithms pass on. A ring can instead *pretend* to be a field
(`gr_ctx_set_is_pretend_field`): algorithms for fields then run, an
inversion of a nonzero non-unit $x$ returns `GR_UNABLE`, and the factor
$\gcd(x, n)$ of $n$ is recorded in the context, where
`gr_ctx_recover_zero_divisor` finds it. This is the dynamic evaluation
technique: compute as if $n$ were prime, and split $n$ when this fails.

```c example name=mpn_mod_pretend
#include <flint/fmpz.h>
#include <flint/gr.h>
#include <flint/mpn_mod.h>

int main(void)
{
    gr_ctx_t ctx;
    fmpz_t n, t;
    ulong x[MPN_MOD_MAX_LIMBS], y[MPN_MOD_MAX_LIMBS];
    int status;

    /* n = (2^64 - 59) (2^61 - 1), a product of two primes */
    fmpz_init(n);
    fmpz_init(t);
    fmpz_set_str(n, "42535295865117307778430344311653531707", 10);
    GR_MUST_SUCCEED(gr_ctx_init_mpn_mod(ctx, n));

    /* without the pretense, a non-unit gives GR_DOMAIN (1) */
    GR_MUST_SUCCEED(mpn_mod_set_ui(x, UWORD(2305843009213693951), ctx));
    status = mpn_mod_inv(y, x, ctx);
    printf("status %d\n", status);

    /* pretending to be a field: GR_UNABLE (2), and a factor is recorded */
    GR_MUST_SUCCEED(gr_ctx_set_is_pretend_field(ctx, T_TRUE));
    status = mpn_mod_inv(y, x, ctx);
    printf("status %d\n", status);
    if (gr_ctx_recover_zero_divisor(y, ctx) == GR_SUCCESS)
    {
        GR_MUST_SUCCEED(mpn_mod_get_fmpz(t, y, ctx));
        fmpz_print(t);
        printf("\n");
    }

    gr_ctx_clear(ctx);
    fmpz_clear(n);
    fmpz_clear(t);
    return 0;
}
```

```output
status 1
status 2
2305843009213693951
```

The determinant (`gr_mat_det`) works for every modulus: when the LU
decomposition meets a non-invertible pivot, `mpn_mod_mat_det` falls back
to the division-free Berkowitz algorithm. The characteristic polynomial
(`gr_mat_charpoly`) falls back in the same way when Danilevsky's
algorithm meets a non-invertible pivot.
