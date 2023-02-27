/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_poly.h"
#include "fmpz_mod.h"

/* split f assuming that f has degree(f) distinct nonzero roots in Fp */
void _fmpz_mod_poly_split_rabin(
    fmpz_mod_poly_t a,
    fmpz_mod_poly_t b,
    const fmpz_mod_poly_t f,
    const fmpz_t halfp, /* floor((p-1)/2) */
    fmpz_mod_poly_t t,
    fmpz_mod_poly_t t2,
    flint_rand_t randstate,
    const fmpz_mod_ctx_t ctx)
{
    FLINT_ASSERT(fmpz_mod_poly_degree(f, ctx) > 1);

    fmpz_mod_poly_fit_length(a, 2, ctx);
    fmpz_mod_poly_fit_length(b, 2, ctx);
    fmpz_mod_poly_fit_length(t, 3, ctx);

    if (fmpz_mod_poly_degree(f, ctx) == 2)
    {
        fmpz * A, * T;
        const fmpz * B;

        B = f->coeffs + 0;
        A = t->coeffs + 1;
        T = t->coeffs + 2;

        if (!fmpz_is_zero(halfp))
        {
            fmpz_mod_mul(A, f->coeffs + 1, halfp, ctx);
            fmpz_mod_neg(A, A, ctx);
            if (!fmpz_is_one(f->coeffs + 2))
            {
                fmpz_mod_inv(T, f->coeffs + 2, ctx);
                fmpz_mod_mul(A, A, T, ctx);
                fmpz_mod_mul(t->coeffs + 0, B, T, ctx);
                B = t->coeffs + 0;
            }

            fmpz_mod_mul(T, A, A, ctx);
            fmpz_mod_sub(T, T, B, ctx);
            if (!fmpz_sqrtmod(b->coeffs + 0, T, fmpz_mod_ctx_modulus(ctx)))
                flint_throw(FLINT_ERROR, "_fmpz_mod_poly_split_rabin: f is irreducible");

            fmpz_mod_add(a->coeffs + 0, A, b->coeffs + 0, ctx);
            fmpz_mod_sub(b->coeffs + 0, A, b->coeffs + 0, ctx);
        }
        else
        {
            fmpz_one(a->coeffs + 0);
            fmpz_zero(b->coeffs + 0);
        }

        fmpz_one(a->coeffs + 1);
        fmpz_one(b->coeffs + 1);
        _fmpz_mod_poly_set_length(a, 2);
        _fmpz_mod_poly_set_length(b, 2);

    #if FLINT_WANT_ASSERT
        fmpz_mod_add(T, a->coeffs + 0, b->coeffs + 0, ctx);
        fmpz_mod_mul(T, T, f->coeffs + 2, ctx);
        FLINT_ASSERT(fmpz_equal(T, f->coeffs + 1));
        fmpz_mod_mul(T, a->coeffs + 0, b->coeffs + 0, ctx);
        fmpz_mod_mul(T, T, f->coeffs + 2, ctx);
        FLINT_ASSERT(fmpz_equal(T, f->coeffs + 0));
    #endif

        return;
    }

    fmpz_mod_poly_reverse(t, f, f->length, ctx);
    fmpz_mod_poly_inv_series_newton(t2, t, t->length, ctx);

try_again:

    fmpz_randm(a->coeffs, randstate, fmpz_mod_ctx_modulus(ctx));
    fmpz_mod_poly_powmod_linear_fmpz_preinv(t, a->coeffs, halfp, f, t2, ctx);
    fmpz_mod_poly_sub_si(t, t, 1, ctx);
    fmpz_mod_poly_gcd(a, t, f, ctx);

    if (a->length <= 1 || a->length >= f->length)
        goto try_again;

    fmpz_mod_poly_divrem(b, t, f, a, ctx);
    FLINT_ASSERT(fmpz_mod_poly_is_zero(t, ctx));

    /* ensure deg a >= deg b */
    if (a->length < b->length)
        fmpz_mod_poly_swap(a, b, ctx);

    return;
}


/*
    If P has deg(P) distinct nonzero roots of P, fill them in and return 1.
    Otherwise return 0. Function is undefined for zero P.
    The modulus of P is assumed to be prime.
*/
int fmpz_mod_poly_find_distinct_nonzero_roots(
    fmpz * roots,
    const fmpz_mod_poly_t P,
    const fmpz_mod_ctx_t ctx)
{
    fmpz_t a0, a1;
    int success;
    slong i, roots_idx, sp;
    fmpz_t halfp;
    fmpz_mod_poly_struct * a , * b;
    fmpz_mod_poly_t f, t, t2;
    fmpz_mod_poly_struct stack[FLINT_BITS + 1];
    flint_rand_t randstate;
    slong d = fmpz_mod_poly_degree(P, ctx);

    FLINT_ASSERT(d >= 0);

    fmpz_init(a0);
    fmpz_init(a1);
    fmpz_init(halfp);

    if (d < 2)
    {
        if (d == 1)
        {
            fmpz_mod_poly_get_coeff_fmpz(a0, P, 0, ctx);
            fmpz_mod_poly_get_coeff_fmpz(a1, P, 1, ctx);
            if (fmpz_is_zero(a0))
            {
                success = 0;
                goto cleanup1;
            }
            fmpz_mod_inv(a1, a1, ctx);
            fmpz_mod_neg(a1, a1, ctx);
            fmpz_mod_mul(roots + 0, a0, a1, ctx);
        }
        success = 1;
        goto cleanup1;
    }

    if (fmpz_cmp_ui(fmpz_mod_ctx_modulus(ctx), 2) <= 0)
    {
        success = 0;
        goto cleanup1;
    }

    if (fmpz_is_zero(P->coeffs + 0))
    {
        success = 0;
        goto cleanup1;
    }

    flint_randinit(randstate);
    fmpz_mod_poly_init(t, ctx);
    fmpz_mod_poly_init(t2, ctx);
    fmpz_mod_poly_init(f, ctx);
    for (i = 0; i <= FLINT_BITS; i++)
        fmpz_mod_poly_init(stack + i, ctx);

    roots_idx = 0;

    fmpz_mod_poly_make_monic(f, P, ctx);
    fmpz_mod_poly_reverse(t, f, f->length, ctx);
    fmpz_mod_poly_inv_series_newton(t2, t, f->length, ctx);

    a = stack + 0;
    fmpz_sub_ui(halfp, fmpz_mod_ctx_modulus(ctx), 1);
    fmpz_divexact_ui(halfp, halfp, 2);
    fmpz_mod_poly_powmod_x_fmpz_preinv(t, halfp, f, t2, ctx);
    fmpz_mod_poly_sub_si(t, t, 1, ctx);
    fmpz_mod_poly_gcd(a, t, f, ctx);

    b = stack + 1;
    fmpz_mod_poly_add_si(t, t, 2, ctx);
    fmpz_mod_poly_gcd(b, t, f, ctx);

    if (fmpz_mod_poly_degree(b, ctx) + fmpz_mod_poly_degree(a, ctx) != d)
    {
        success = 0;
        goto cleanup;
    }
    /* deg a >= deg b */
    if (fmpz_mod_poly_degree(a, ctx) < fmpz_mod_poly_degree(b, ctx))
    {
        fmpz_mod_poly_swap(a, b, ctx);
    }

    sp = fmpz_mod_poly_degree(b, ctx) > 0 ? 2 : 1;
    while (sp > 0)
    {
        FLINT_ASSERT(sp < FLINT_BITS);
        sp--;
        fmpz_mod_poly_swap(f, stack + sp, ctx);

        FLINT_ASSERT(fmpz_mod_poly_degree(f, ctx) > 0);
        if (fmpz_mod_poly_degree(f, ctx) == 1)
        {
            fmpz_mod_poly_get_coeff_fmpz(a0, f, 0, ctx);
            fmpz_mod_poly_get_coeff_fmpz(a1, f, 1, ctx);
            FLINT_ASSERT(!fmpz_is_zero(a0));
            FLINT_ASSERT(fmpz_is_one(a1));
            fmpz_mod_neg(roots + roots_idx, a0, ctx);
            roots_idx++;
        }
        else
        {
            _fmpz_mod_poly_split_rabin(stack + sp + 0, stack + sp + 1, f,
                                                 halfp, t, t2, randstate, ctx);

            FLINT_ASSERT(
                FLINT_BIT_COUNT(fmpz_mod_poly_degree(stack + sp + 1, ctx))
                                                       <= FLINT_BITS - sp - 1);
            sp += 2;
        }
    }

    success = 1;

cleanup:

    flint_randclear(randstate);
    fmpz_mod_poly_clear(t, ctx);
    fmpz_mod_poly_clear(t2, ctx);
    fmpz_mod_poly_clear(f, ctx);
    for (i = 0; i <= FLINT_BITS; i++)
        fmpz_mod_poly_clear(stack + i, ctx);

    FLINT_ASSERT((!success) || roots_idx == d);

cleanup1:

    fmpz_clear(a0);
    fmpz_clear(a1);
    fmpz_clear(halfp);

    return success;
}
