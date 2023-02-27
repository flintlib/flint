/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_poly.h"
#include "fmpz_mod_poly_factor.h"
#include "fmpz.h"

/*
    Helper function for finding roots. The roots of a monic f are written
    with exponent given in mult to r. Uses Rabin's Las Vegas algorithm
    via gcd computations with (x + delta)^halfp - 1.
*/
static void _fmpz_mod_poly_push_roots(
    fmpz_mod_poly_factor_t r,
    fmpz_mod_poly_t f,              /* clobbered */
    slong mult,                     /* expoenent to write on the roots */
    const fmpz_t halfp,             /* (f->p - 1)/2 */
    fmpz_mod_poly_t t,              /* temp */
    fmpz_mod_poly_t t2,             /* more temp */
    fmpz_mod_poly_struct * stack,   /* temp of size FLINT_BITS */
    flint_rand_t randstate,
    const fmpz_mod_ctx_t ctx)
{
    slong i, sp;
    fmpz_mod_poly_struct * a, * b;
    const fmpz * p = fmpz_mod_ctx_modulus(ctx);

    FLINT_ASSERT(fmpz_mod_poly_degree(f, ctx) >= 1);
    FLINT_ASSERT(fmpz_is_one(f->coeffs + f->length - 1));
    FLINT_ASSERT(fmpz_is_probabprime(p));

    /* handle at least p = 2 */
    if (fmpz_cmp_ui(p, 10) < 0)
    {
        fmpz_t x, e;
        fmpz_init(x);
        fmpz_init(e);
        for (fmpz_zero(x); fmpz_cmp(x, p) < 0; fmpz_add_ui(x, x, 1))
        {
            fmpz_mod_poly_evaluate_fmpz(e, f, x, ctx);
            if (!fmpz_is_zero(e))
                continue;
            fmpz_mod_poly_factor_fit_length(r, r->num + 1, ctx);
            fmpz_mod_poly_fit_length(r->poly + r->num, 2, ctx);
            fmpz_mod_neg(r->poly[r->num].coeffs + 0, x, ctx);
            fmpz_one(r->poly[r->num].coeffs + 1);
            r->poly[r->num].length = 2;
            r->exp[r->num] = mult;
            r->num++;
        }
        fmpz_clear(e);
        fmpz_clear(x);
        return;
    }

    /* handle zero roots */
    if (fmpz_is_zero(f->coeffs + 0))
    {
        fmpz_mod_poly_factor_fit_length(r, r->num + 1, ctx);
        fmpz_mod_poly_fit_length(r->poly + r->num, 2, ctx);
        fmpz_zero(r->poly[r->num].coeffs + 0);
        fmpz_one(r->poly[r->num].coeffs + 1);
        r->poly[r->num].length = 2;
        r->exp[r->num] = mult;
        r->num++;

        i = 1;
        while (i < f->length && fmpz_is_zero(f->coeffs + i))
            i++;

        fmpz_mod_poly_shift_right(f, f, i, ctx);
    }

    if (fmpz_mod_poly_degree(f, ctx) <= 1)
    {
        if (fmpz_mod_poly_degree(f, ctx) == 1)
        {
            fmpz_mod_poly_factor_fit_length(r, r->num + 1, ctx);
            fmpz_mod_poly_swap(r->poly + r->num, f, ctx);
            r->exp[r->num] = mult;
            r->num++;
        }
        return;
    }

    FLINT_ASSERT(!fmpz_is_zero(f->coeffs + 0));
    fmpz_mod_poly_reverse(t, f, f->length, ctx);
    fmpz_mod_poly_inv_series_newton(t2, t, t->length, ctx);

    a = stack + 0;
    fmpz_mod_poly_powmod_x_fmpz_preinv(t, halfp, f, t2, ctx);
    fmpz_mod_poly_sub_si(t, t, 1, ctx);
    fmpz_mod_poly_gcd(a, t, f, ctx);

    b = stack + 1;
    fmpz_mod_poly_add_si(t, t, 2, ctx);
    fmpz_mod_poly_gcd(b, t, f, ctx);

    /* ensure deg a >= deg b */
    if (fmpz_mod_poly_degree(a, ctx) < fmpz_mod_poly_degree(b, ctx))
        fmpz_mod_poly_swap(a, b, ctx);

    fmpz_mod_poly_factor_fit_length(r, r->num +
             fmpz_mod_poly_degree(a, ctx) + fmpz_mod_poly_degree(b, ctx), ctx);

    /* initial split failed if b = 1 */
    sp = (fmpz_mod_poly_degree(b, ctx) > 0) ? 2 : 1;
    while (sp > 0)
    {
        sp--;
        FLINT_ASSERT(sp < FLINT_BITS);
        fmpz_mod_poly_swap(f, stack + sp, ctx);
        FLINT_ASSERT(fmpz_mod_poly_degree(f, ctx) >= 0);
        FLINT_ASSERT(FLINT_BIT_COUNT(fmpz_mod_poly_degree(f, ctx)) <= FLINT_BITS - sp);

        if (fmpz_mod_poly_degree(f, ctx) <= 1)
        {
            if (fmpz_mod_poly_degree(f, ctx) == 1)
            {
                FLINT_ASSERT(r->num < r->alloc);
                fmpz_mod_poly_set(r->poly + r->num, f, ctx);
                r->exp[r->num] = mult;
                r->num++;
            }
        }
        else
        {
            FLINT_ASSERT(sp + 1 < FLINT_BITS);

            _fmpz_mod_poly_split_rabin(stack + sp + 0, stack + sp + 1, f,
                                                 halfp, t, t2, randstate, ctx);

            FLINT_ASSERT(FLINT_BIT_COUNT(fmpz_mod_poly_degree(stack + sp + 1, ctx))
                                                     <= FLINT_BITS - (sp + 1));
            FLINT_ASSERT(FLINT_BIT_COUNT(fmpz_mod_poly_degree(stack + sp + 0, ctx))
                                                     <= FLINT_BITS - (sp + 0));
            sp += 2;
        }
    }
}

void fmpz_mod_poly_roots(fmpz_mod_poly_factor_t r, const fmpz_mod_poly_t f,
                                       int with_mult, const fmpz_mod_ctx_t ctx)
{
    slong i;
    fmpz_t p2;
    flint_rand_t randstate;
    fmpz_mod_poly_struct t[FLINT_BITS + 3];

    FLINT_ASSERT(fmpz_is_probabprime(fmpz_mod_ctx_modulus(ctx)));

    r->num = 0;

    if (fmpz_mod_poly_degree(f, ctx) < 2)
    {
        if (fmpz_mod_poly_degree(f, ctx) == 1)
        {
            fmpz_mod_poly_factor_fit_length(r, 1, ctx);
            fmpz_mod_poly_make_monic(r->poly + 0, f, ctx);
            r->exp[0] = 1;
            r->num = 1;
        }
        else if (fmpz_mod_poly_degree(f, ctx) < 0)
        {
            flint_throw(FLINT_ERROR,
                "Exception in fmpz_mod_poly_roots: input polynomial is zero.");
        }
        return;
    }

    fmpz_init_set(p2, fmpz_mod_ctx_modulus(ctx));
    fmpz_sub_ui(p2, p2, 1);
    fmpz_fdiv_q_2exp(p2, p2, 1);

    flint_randinit(randstate);

    for (i = 0; i < FLINT_BITS + 3; i++)
        fmpz_mod_poly_init(t + i, ctx);

    fmpz_mod_poly_make_monic(t + 0, f, ctx);

    if (with_mult)
    {
        fmpz_mod_poly_factor_t sqf;
        fmpz_mod_poly_factor_init(sqf, ctx);
        fmpz_mod_poly_factor_squarefree(sqf, t + 0, ctx);
        for (i = 0; i < sqf->num; i++)
        {
            _fmpz_mod_poly_push_roots(r, sqf->poly + i, sqf->exp[i],
                                      p2, t + 1, t + 2, t + 3, randstate, ctx);
        }
        fmpz_mod_poly_factor_clear(sqf, ctx);
    }
    else
    {
        _fmpz_mod_poly_push_roots(r, t + 0, 1,
                                      p2, t + 1, t + 2, t + 3, randstate, ctx);
    }

    fmpz_clear(p2);

    flint_randclear(randstate);

    for (i = 0; i < FLINT_BITS + 3; i++)
        fmpz_mod_poly_clear(t + i, ctx);
}

