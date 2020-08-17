/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"
#include "nmod_poly_factor.h"

/*
    Helper function for finding roots. The roots of a monic f are written
    with exponent given in mult to r. Uses Rabin's Las Vegas algorithm
    via gcd computations with (x + delta)^((p-1)/2) - 1.
*/
static void _nmod_poly_push_roots(
    nmod_poly_factor_t r,
    nmod_poly_t f,              /* clobbered */
    slong mult,                 /* expoenent to write on the roots */
    nmod_poly_t t,              /* temp */
    nmod_poly_t t2,             /* more temp */
    nmod_poly_struct * stack,   /* temp of size FLINT_BITS */
    flint_rand_t randstate)
{
    slong i, sp;
    nmod_poly_struct * a, * b;

    FLINT_ASSERT(nmod_poly_degree(f) >= 1);
    FLINT_ASSERT(f->coeffs[nmod_poly_degree(f)] == 1);
    FLINT_ASSERT(n_is_probabprime(f->mod.n));

    /* handle at least p = 2 */
    if (f->mod.n < 10)
    {
        ulong x;
        for (x = 0; x < f->mod.n; x++)
        {
            if (0 != nmod_poly_evaluate_nmod(f, x))
                continue;
            nmod_poly_factor_fit_length(r, r->num + 1);
            nmod_poly_fit_length(r->p + r->num, 2);
            r->p[r->num].mod = f->mod;        /* bummer */
            r->p[r->num].coeffs[0] = nmod_neg(x, f->mod);
            r->p[r->num].coeffs[1] = 1;
            r->p[r->num].length = 2;
            r->exp[r->num] = mult;
            r->num++;
        }
        return;
    }

    /* handle zero roots */
    if (f->coeffs[0] == 0)
    {
        nmod_poly_factor_fit_length(r, r->num + 1);
        nmod_poly_fit_length(r->p + r->num, 2);
        r->p[r->num].mod = f->mod;        /* bummer */
        r->p[r->num].coeffs[0] = 0;
        r->p[r->num].coeffs[1] = 1;
        r->p[r->num].length = 2;
        r->exp[r->num] = mult;
        r->num++;

        i = 1;
        while (i < f->length && f->coeffs[i] == 0)
            i++;

        nmod_poly_shift_right(f, f, i);
    }

    if (nmod_poly_degree(f) <= 1)
    {
        if (nmod_poly_degree(f) == 1)
        {
            nmod_poly_factor_fit_length(r, r->num + 1);
            r->p[r->num].mod = f->mod;        /* bummer */
            nmod_poly_swap(r->p + r->num, f);
            r->exp[r->num] = mult;
            r->num++;
        }
        return;
    }

    FLINT_ASSERT(f->coeffs[0] != 0);
    nmod_poly_reverse(t, f, f->length);
    nmod_poly_inv_series_newton(t2, t, t->length);

    a = stack + 0;
    b = stack + 1;

    nmod_poly_zero(a);
    nmod_poly_set_coeff_ui(a, 1, 1);
    nmod_poly_powmod_ui_binexp_preinv(t, a, (f->mod.n - 1)/2, f, t2);
    nmod_poly_sub_ui(t, t, 1);
    nmod_poly_gcd(a, t, f);
    nmod_poly_add_ui(t, t, 2);
    nmod_poly_gcd(b, t, f);

    /* ensure deg a >= deg b */
    if (nmod_poly_degree(a) < nmod_poly_degree(b))
        nmod_poly_swap(a, b);

    nmod_poly_factor_fit_length(r, r->num + nmod_poly_degree(a) +
                                            nmod_poly_degree(b));

    /* initial split failed if b = 1 */
    sp = (nmod_poly_degree(b) > 0) ? 2 : 1;
    while (sp > 0)
    {
        sp--;
        FLINT_ASSERT(sp < FLINT_BITS);
        nmod_poly_swap(f, stack + sp);
        FLINT_ASSERT(nmod_poly_degree(f) >= 0);
        FLINT_ASSERT(FLINT_BIT_COUNT(nmod_poly_degree(f)) <= FLINT_BITS - sp);

        if (nmod_poly_degree(f) <= 1)
        {
            if (nmod_poly_degree(f) == 1)
            {
                FLINT_ASSERT(r->num < r->alloc);
                r->p[r->num].mod = f->mod;       /* bummer */
                nmod_poly_set(r->p + r->num, f);
                r->exp[r->num] = mult;
                r->num++;
            }
        }
        else
        {
            FLINT_ASSERT(sp + 1 < FLINT_BITS);

            _nmod_poly_split_rabin(stack + sp + 0, stack + sp + 1, f, t, t2,
                                                                    randstate);

            FLINT_ASSERT(FLINT_BIT_COUNT(nmod_poly_degree(stack + sp + 1))
                                                 <= FLINT_BITS - (sp + 1));
            FLINT_ASSERT(FLINT_BIT_COUNT(nmod_poly_degree(stack + sp + 0))
                                                 <= FLINT_BITS - (sp + 0));
            sp += 2;
        }
    }
}

void nmod_poly_roots(nmod_poly_factor_t r, const nmod_poly_t f,
                                                         int with_multiplicity)
{
    slong i;
    flint_rand_t randstate;
    nmod_poly_struct t[FLINT_BITS + 3];

    FLINT_ASSERT(n_is_probabprime(f->mod.n));

    r->num = 0;

    if (nmod_poly_degree(f) < 2)
    {
        if (nmod_poly_degree(f) == 1)
        {
            nmod_poly_factor_fit_length(r, 1);
            r->p[0].mod = f->mod;            /* bummer */
            nmod_poly_make_monic(r->p + 0, f);
            r->exp[0] = 1;
            r->num = 1;
        }
        else if (nmod_poly_degree(f) < 0)
        {
            flint_throw(FLINT_ERROR, "Exception in nmod_poly_roots: "
                                                  "input polynomial is zero.");
        }
        return;
    }

    flint_randinit(randstate);

    for (i = 0; i < FLINT_BITS + 3; i++)
        nmod_poly_init_mod(t + i, f->mod);

    if (with_multiplicity)
    {
        nmod_poly_factor_t sqf;
        nmod_poly_factor_init(sqf);
        nmod_poly_factor_squarefree(sqf, f);
        for (i = 0; i < sqf->num; i++)
        {
            _nmod_poly_push_roots(r, sqf->p + i, sqf->exp[i],
                                               t + 1, t + 2, t + 3, randstate);
        }
        nmod_poly_factor_clear(sqf);
    }
    else
    {
        nmod_poly_make_monic(t + 0, f);
        _nmod_poly_push_roots(r, t + 0, 1, t + 1, t + 2, t + 3, randstate);
    }

    flint_randclear(randstate);

    for (i = 0; i < FLINT_BITS + 3; i++)
        nmod_poly_clear(t + i);
}

