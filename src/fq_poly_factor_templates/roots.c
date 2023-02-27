/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

#include "ulong_extras.h"

/* split f assuming that f has degree(f) distinct nonzero roots in Fq */
void _TEMPLATE(T, poly_split_rabin)(
    TEMPLATE(T, poly_t) a,
    TEMPLATE(T, poly_t) b,
    const TEMPLATE(T, poly_t) f,
    const fmpz_t halfq,         /* (q-1)/2  or 0 in characteristic 2 */
    TEMPLATE(T, poly_t) t,      /* temp space */
    TEMPLATE(T, poly_t) t2,     /* temp space */
    flint_rand_t randstate,
    const TEMPLATE(T, ctx_t) ctx)
{
    slong i;

    FLINT_ASSERT(TEMPLATE(T, poly_degree)(f, ctx) > 1);

    TEMPLATE(T, poly_reverse)(t, f, f->length, ctx);
    TEMPLATE(T, poly_inv_series_newton)(t2, t, t->length, ctx);

try_again:

    /* a = random linear */
    TEMPLATE(T, poly_fit_length)(a, 2, ctx);
    TEMPLATE(T, rand)(a->coeffs + 0, randstate, ctx);
    TEMPLATE(T, rand)(a->coeffs + 1, randstate, ctx);
    if (TEMPLATE(T, is_zero)(a->coeffs + 1, ctx))
        TEMPLATE(T, one)(a->coeffs + 1, ctx);
    a->length = 2;

    if (!fmpz_is_zero(halfq))
    {
        FLINT_ASSERT(fmpz_cmp_ui(TEMPLATE(T, ctx_prime)(ctx), 2) > 0);
        TEMPLATE(T, poly_powmod_fmpz_sliding_preinv)(t, a, halfq, 0, f, t2, ctx);
        TEMPLATE(T, poly_add_si)(t, t, -1, ctx);
    }
    else
    {
        /* it is important that coeff(a, x^1) is random */
        FLINT_ASSERT(fmpz_cmp_ui(TEMPLATE(T, ctx_prime)(ctx), 2) == 0);
        TEMPLATE(T, poly_set)(t, a, ctx);
        for (i = TEMPLATE(T, ctx_degree)(ctx); i > 1; i--)
        {
            TEMPLATE(T, poly_powmod_ui_binexp_preinv)(a, a, 2, f, t2, ctx);
            TEMPLATE(T, poly_add)(t, t, a, ctx);
        }
    }

    TEMPLATE(T, poly_gcd)(a, t, f, ctx);

    FLINT_ASSERT(!TEMPLATE(T, poly_is_zero)(a, ctx));

    if (0 >= TEMPLATE(T, poly_degree)(a, ctx) ||
        TEMPLATE(T, poly_degree)(a, ctx) >= TEMPLATE(T, poly_degree)(f, ctx))
    {
        goto try_again;
    }

    TEMPLATE(T, poly_div_basecase)(b, f, a, ctx);

    /* ensure deg a >= deg b */
    if (TEMPLATE(T, poly_degree)(a, ctx) < TEMPLATE(T, poly_degree)(b, ctx))
        TEMPLATE(T, poly_swap)(a, b, ctx);

    return;
}


/*
    Helper function for finding roots. The roots of a monic f are written
    with exponent given in mult to r. Uses Rabin's Las Vegas algorithm
    via gcd computations with (x + delta)^halfq - 1.
*/
static void _TEMPLATE(T, poly_push_roots)(
    TEMPLATE(T, poly_factor_t) r,
    TEMPLATE(T, poly_t) f,              /* clobbered */
    slong mult,                         /* expoenent to write on the roots */
    const fmpz_t halfq,                 /* (q-1)/2  or 0 in characteristic 2 */
    TEMPLATE(T, poly_t) t,              /* temp */
    TEMPLATE(T, poly_t) t2,             /* more temp */
    TEMPLATE(T, poly_struct) * stack,   /* temp of size FLINT_BITS */
    flint_rand_t randstate,
    const TEMPLATE(T, ctx_t) ctx)
{
    slong i, sp;
    TEMPLATE(T, poly_struct) * a, * b;

    FLINT_ASSERT(TEMPLATE(T, poly_degree)(f, ctx) >= 1);
    FLINT_ASSERT(TEMPLATE(T, is_one)(f->coeffs + TEMPLATE(T, poly_degree)(f, ctx), ctx));

    /* handle zero roots */
    if (TEMPLATE(T, is_zero)(f->coeffs + 0, ctx))
    {
        TEMPLATE(T, poly_factor_fit_length)(r, r->num + 1, ctx);
        TEMPLATE(T, poly_fit_length)(r->poly + r->num, 2, ctx);
        TEMPLATE(T, zero)(r->poly[r->num].coeffs + 0, ctx);
        TEMPLATE(T, one)(r->poly[r->num].coeffs + 1, ctx);
        r->poly[r->num].length = 2;
        r->exp[r->num] = mult;
        r->num++;

        i = 1;
        while (i < f->length && TEMPLATE(T, is_zero)(f->coeffs + i, ctx))
            i++;

        TEMPLATE(T, poly_shift_right)(f, f, i, ctx);
    }

    if (TEMPLATE(T, poly_degree)(f, ctx) <= 1)
    {
        if (TEMPLATE(T, poly_degree)(f, ctx) == 1)
        {
            TEMPLATE(T, poly_factor_fit_length)(r, r->num + 1, ctx);
            TEMPLATE(T, poly_swap)(r->poly + r->num, f, ctx);
            r->exp[r->num] = mult;
            r->num++;
        }
        return;
    }

    FLINT_ASSERT(!(TEMPLATE(T, is_zero)(f->coeffs + 0, ctx)));
    TEMPLATE(T, poly_reverse)(t, f, f->length, ctx);
    TEMPLATE(T, poly_inv_series_newton)(t2, t, t->length, ctx);

    a = stack + 0;
    b = stack + 1;

    TEMPLATE(T, poly_gen)(a, ctx);

    if (!fmpz_is_zero(halfq))
    {
        FLINT_ASSERT(fmpz_cmp_ui(TEMPLATE(T, ctx_prime)(ctx), 2) > 0);
        TEMPLATE(T, poly_powmod_fmpz_sliding_preinv)(t, a, halfq, 0, f, t2, ctx);
        TEMPLATE(T, poly_add_si)(t, t, -1, ctx);
        TEMPLATE(T, poly_gcd)(a, t, f, ctx);
        TEMPLATE(T, poly_add_si)(t, t, 1, ctx);
    }
    else
    {
        FLINT_ASSERT(fmpz_cmp_ui(TEMPLATE(T, ctx_prime)(ctx), 2) == 0);
        TEMPLATE(T, poly_set)(t, a, ctx);
        for (i = TEMPLATE(T, ctx_degree)(ctx); i > 1; i--)
        {
            TEMPLATE(T, poly_powmod_ui_binexp_preinv)(a, a, 2, f, t2, ctx);
            TEMPLATE(T, poly_add)(t, t, a, ctx);
        }
        TEMPLATE(T, poly_gcd)(a, t, f, ctx);
    }

    TEMPLATE(T, poly_add_si)(t, t, 1, ctx);
    TEMPLATE(T, poly_gcd)(b, t, f, ctx);

    /* ensure deg a >= deg b */
    if (TEMPLATE(T, poly_degree)(a, ctx) < TEMPLATE(T, poly_degree)(b, ctx))
        TEMPLATE(T, poly_swap)(a, b, ctx);

    TEMPLATE(T, poly_factor_fit_length)(r, r->num +
                                        TEMPLATE(T, poly_degree)(a, ctx) +
                                        TEMPLATE(T, poly_degree)(b, ctx), ctx);

    /* initial split failed if b = 1 */
    sp = (TEMPLATE(T, poly_degree)(b, ctx) > 0) ? 2 : 1;
    while (sp > 0)
    {
        sp--;
        FLINT_ASSERT(sp < FLINT_BITS);

        TEMPLATE(T, poly_swap)(f, stack + sp, ctx);

        FLINT_ASSERT(TEMPLATE(T, poly_degree)(f, ctx) >= 0);
        FLINT_ASSERT(
            FLINT_BIT_COUNT(TEMPLATE(T, poly_degree)(f, ctx))
            <= FLINT_BITS - sp);

        if (TEMPLATE(T, poly_degree)(f, ctx) <= 1)
        {
            if (TEMPLATE(T, poly_degree)(f, ctx) == 1)
            {
                FLINT_ASSERT(r->num < r->alloc);
                TEMPLATE(T, poly_set)(r->poly + r->num, f, ctx);
                r->exp[r->num] = mult;
                r->num++;
            }
        }
        else
        {
            FLINT_ASSERT(sp + 1 < FLINT_BITS);

            _TEMPLATE(T, poly_split_rabin)(stack + sp + 0, stack + sp + 1,
                                              f, halfq, t, t2, randstate, ctx);

            FLINT_ASSERT(
                FLINT_BIT_COUNT(TEMPLATE(T, poly_degree)(stack + sp + 1, ctx))
                <= FLINT_BITS - (sp + 1));

            FLINT_ASSERT(
                FLINT_BIT_COUNT(TEMPLATE(T, poly_degree)(stack + sp + 0, ctx))
                <= FLINT_BITS - (sp + 0));

            sp += 2;
        }
    }
}

void TEMPLATE(T, poly_roots)(
    TEMPLATE(T, poly_factor_t) r,
    const TEMPLATE(T, poly_t) f,
    int with_multiplicity,
    const TEMPLATE(T, ctx_t) ctx)
{
    slong i;
    fmpz_t q2;
    flint_rand_t randstate;
    TEMPLATE(T, poly_struct) t[FLINT_BITS + 3];

    r->num = 0;

    if (TEMPLATE(T, poly_degree)(f, ctx) < 2)
    {
        if (TEMPLATE(T, poly_degree)(f, ctx) == 1)
        {
            TEMPLATE(T, poly_factor_fit_length)(r, 1, ctx);
            TEMPLATE(T, poly_make_monic)(r->poly + 0, f, ctx);
            r->exp[0] = 1;
            r->num = 1;
        }
        else if (TEMPLATE(T, poly_degree)(f, ctx) < 0)
        {
            flint_throw(FLINT_ERROR, "Exception in fq_poly_roots: "
                                                  "input polynomial is zero.");
        }

        return;
    }

    fmpz_init(q2);
    TEMPLATE(T, ctx_order(q2, ctx));
    fmpz_sub_ui(q2, q2, 1);
    if (fmpz_is_even(q2))
        fmpz_fdiv_q_2exp(q2, q2, 1);
    else
        fmpz_zero(q2);

    flint_randinit(randstate);

    for (i = 0; i < FLINT_BITS + 3; i++)
        TEMPLATE(T, poly_init)(t + i, ctx);

    if (with_multiplicity)
    {
        TEMPLATE(T, poly_factor_t) sqf;
        TEMPLATE(T, poly_factor_init)(sqf, ctx);
        TEMPLATE(T, poly_factor_squarefree)(sqf, f, ctx);
        for (i = 0; i < sqf->num; i++)
        {
            _TEMPLATE(T, poly_push_roots)(r, sqf->poly + i, sqf->exp[i],
                                      q2, t + 1, t + 2, t + 3, randstate, ctx);
        }
        TEMPLATE(T, poly_factor_clear)(sqf, ctx);
    }
    else
    {
        TEMPLATE(T, poly_make_monic)(t + 0, f, ctx);
        _TEMPLATE(T, poly_push_roots)(r, t + 0, 1,
                                      q2, t + 1, t + 2, t + 3, randstate, ctx);
    }

    fmpz_clear(q2);

    flint_randclear(randstate);

    for (i = 0; i < FLINT_BITS + 3; i++)
        TEMPLATE(T, poly_clear)(t + i, ctx);
}

#endif
