/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_poly.h"
#include "fmpz_mod_poly_factor.h"
#include "fmpz.h"


void _fmpz_mod_poly_push_roots(
    fmpz_mod_poly_factor_t r,
    fmpz_mod_poly_t f,              /* clobbered */
    slong mult,                     /* expoenent to write on the roots */
    const fmpz_t halfp,             /* (f->p - 1)/2 */
    fmpz_mod_poly_t T,              /* temp */
    fmpz_mod_poly_struct * stack,   /* temp of size FLINT_BITS */
    flint_rand_t randstate)
{
    slong sp;
    fmpz_mod_poly_struct * a, * b;

    FLINT_ASSERT(fmpz_mod_poly_degree(f) >= 1);
    FLINT_ASSERT(fmpz_is_one(f->coeffs + fmpz_mod_poly_degree(f)));
    FLINT_ASSERT(fmpz_is_prime(&f->p));

    if (fmpz_mod_poly_degree(f) <= 1)
    {
        fmpz_mod_poly_factor_fit_length(r, r->num + 1);
        fmpz_set(&r->poly[r->num].p, &f->p);        /* bummer */
        fmpz_mod_poly_swap(r->poly + r->num, f);
        r->exp[r->num] = mult;
        r->num++;
        return;
    }

    /* handle at least p = 2 */
    if (fmpz_cmp_ui(&f->p, 10) < 0)
    {
        fmpz_t x, e;
        fmpz_init(x);
        fmpz_init(e);
        for (fmpz_zero(x); fmpz_cmp(x, &f->p) < 0; fmpz_add_ui(x, x, 1))
        {
            fmpz_mod_poly_evaluate_fmpz(e, f, x);
            if (!fmpz_is_zero(e))
                continue;
            fmpz_mod_poly_factor_fit_length(r, r->num + 1);
            fmpz_mod_poly_fit_length(r->poly + r->num, 2);
            fmpz_negmod(r->poly[r->num].coeffs + 0, x, &f->p);
            fmpz_one(r->poly[r->num].coeffs + 1);
            fmpz_set(&r->poly[r->num].p, &f->p);    /* bummer */
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
        fmpz_mod_poly_factor_fit_length(r, r->num + 1);
        fmpz_mod_poly_fit_length(r->poly + r->num, 2);
        fmpz_zero(r->poly[r->num].coeffs + 0);
        fmpz_one(r->poly[r->num].coeffs + 1);
        fmpz_set(&r->poly[r->num].p, &f->p);    /* bummer */
        r->poly[r->num].length = 2;
        r->exp[r->num] = mult;
        r->num++;
    }

    a = stack + 0;
    fmpz_mod_poly_zero(a);
    fmpz_mod_poly_set_coeff_ui(a, 1, 1);
    fmpz_mod_poly_powmod_fmpz_binexp(T, a, halfp, f);
    fmpz_mod_poly_zero(a);
    fmpz_mod_poly_set_coeff_ui(a, 0, 1);
    fmpz_mod_poly_sub(T, T, a);
    fmpz_mod_poly_gcd(a, T, f);

    b = stack + 1;
    fmpz_mod_poly_zero(b);
    fmpz_mod_poly_set_coeff_ui(b, 0, 2);
    fmpz_mod_poly_add(T, T, b);
    fmpz_mod_poly_gcd(b, T, f);

    /* ensure deg a >= deg b */
    if (fmpz_mod_poly_degree(a) < fmpz_mod_poly_degree(b))
        fmpz_mod_poly_swap(a, b);

    fmpz_mod_poly_factor_fit_length(r, r->num +
                 fmpz_mod_poly_degree(a) + fmpz_mod_poly_degree(b));

    /* initial split failed if b = 1 */
    sp = (fmpz_mod_poly_degree(b) > 0) ? 2 : 1;
    while (sp > 0)
    {
        sp--;
        FLINT_ASSERT(sp < FLINT_BITS);
        fmpz_mod_poly_swap(f, stack + sp);
        FLINT_ASSERT(fmpz_mod_poly_degree(f) >= 0);
        FLINT_ASSERT(FLINT_BIT_COUNT(fmpz_mod_poly_degree(f)) <= FLINT_BITS - sp);

        if (fmpz_mod_poly_degree(f) <= 1)
        {
            if (fmpz_mod_poly_degree(f) == 1)
            {
                FLINT_ASSERT(r->num < r->alloc);
                fmpz_set(&r->poly[r->num].p, &f->p);        /* bummer */
                fmpz_mod_poly_set(r->poly + r->num, f);
                r->exp[r->num] = mult;
                r->num++;
            }
        }
        else
        {
            FLINT_ASSERT(sp + 1 < FLINT_BITS);

            _fmpz_mod_poly_split_rabin(stack + sp + 0, stack + sp + 1, T, f,
                                                             halfp, randstate);

            FLINT_ASSERT(FLINT_BIT_COUNT(fmpz_mod_poly_degree(stack + sp + 1))
                                                     <= FLINT_BITS - (sp + 1));
            FLINT_ASSERT(FLINT_BIT_COUNT(fmpz_mod_poly_degree(stack + sp + 0))
                                                     <= FLINT_BITS - (sp + 0));
            sp += 2;
        }
    }
}

void fmpz_mod_poly_roots(fmpz_mod_poly_factor_t r, const fmpz_mod_poly_t f,
                                                                 int want_mult)
{
    slong i;
    fmpz_t p2;
    flint_rand_t randstate;
    fmpz_mod_poly_struct t[FLINT_BITS + 2];

    FLINT_ASSERT(fmpz_is_probabprime(&f->p));

    r->num = 0;

    if (fmpz_mod_poly_degree(f) < 2)
    {
        if (fmpz_mod_poly_degree(f) == 1)
        {
            fmpz_mod_poly_factor_fit_length(r, 1);
            fmpz_set(&r->poly[0].p, &f->p);             /* bummer */
            fmpz_mod_poly_make_monic(r->poly + 0, f);
            r->exp[0] = 1;
            r->num = 1;
        }
        else if (fmpz_mod_poly_degree(f) < 0)
        {
            flint_throw(FLINT_ERROR,
                "Exception in fmpz_mod_poly_roots: input polynomial is zero.");
        }
        return;
    }

    fmpz_init_set(p2, &f->p);
    fmpz_sub_ui(p2, p2, 1);
    fmpz_fdiv_q_2exp(p2, p2, 1);

    flint_randinit(randstate);

    for (i = 0; i < FLINT_BITS + 2; i++)
        fmpz_mod_poly_init(t + i, &f->p);

    if (want_mult)
    {
        fmpz_mod_poly_factor_t sqf;
        fmpz_mod_poly_factor_init(sqf);
        fmpz_mod_poly_factor_squarefree(sqf, f);
        for (i = 0; i < sqf->num; i++)
        {
            _fmpz_mod_poly_push_roots(r, sqf->poly + i, sqf->exp[i],
                                               p2, t + 1, t + 2, randstate);
        }
        fmpz_mod_poly_factor_clear(sqf);
    }
    else
    {
        fmpz_mod_poly_make_monic(t + 0, f);
        _fmpz_mod_poly_push_roots(r, t + 0, 1, p2, t + 1, t + 2, randstate);
    }

    fmpz_clear(p2);

    flint_randclear(randstate);

    for (i = 0; i < FLINT_BITS + 2; i++)
        fmpz_mod_poly_clear(t + i);
}

