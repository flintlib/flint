/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"
#include "ulong_extras.h"

/* split f assuming that f has degree(f) distinct nonzero roots in Fp */
static void _nmod_poly_rabinsplit(
    nmod_poly_t a,
    nmod_poly_t b,
    nmod_poly_t T,
    const nmod_poly_t f,
    flint_rand_t randstate)
{
    mp_limb_t delta;

    FLINT_ASSERT(nmod_poly_degree(f) > 1);

try_again:

    delta = n_randint(randstate, f->mod.n);

    nmod_poly_zero(a);
    nmod_poly_set_coeff_ui(a, 1, 1);
    nmod_poly_set_coeff_ui(a, 0, delta);
    nmod_poly_powmod_ui_binexp(T, a, (f->mod.n - 1)/2, f);
    nmod_poly_zero(a);
    nmod_poly_set_coeff_ui(a, 0, 1);
    nmod_poly_sub(T, T, a);
    nmod_poly_gcd(a, T, f);
    FLINT_ASSERT(!nmod_poly_is_zero(a));
    if (0 >= nmod_poly_degree(a) || nmod_poly_degree(a) >= nmod_poly_degree(f))
    {
        goto try_again;
    }
    nmod_poly_div(b, f, a);
    /* deg a >= deg b */
    if (nmod_poly_degree(a) < nmod_poly_degree(b))
    {
        nmod_poly_swap(a, b);
    }
    return;
}

/*
    If P has deg(P) distinct nonzero roots of P, fill them in and return 1.
    Otherwise return 0. Function is undefined for zero P.
    The modulus of P is assumed to be prime.
*/
int nmod_poly_find_distinct_nonzero_roots(
    mp_limb_t * roots,
    const nmod_poly_t P)
{
    mp_limb_t a0, a1;
    int success;
    slong i, roots_idx, sp;
    nmod_poly_struct * a , * b;
    nmod_poly_t f, T;
    nmod_poly_struct stack[FLINT_BITS + 1];
    flint_rand_t randstate;
    slong t = nmod_poly_degree(P);

    if (t < 2)
    {
        if (t == 1)
        {
            a0 = nmod_poly_get_coeff_ui(P, 0);
            a1 = nmod_poly_get_coeff_ui(P, 1);
            if (a0 == 0)
            {
                return 0;
            }
            roots[0] = nmod_mul(a0, nmod_inv(P->mod.n - a1, P->mod), P->mod);
        }
        return 1;
    }

    if (P->mod.n == 2)
    {
        return 0;
    }

    flint_randinit(randstate);
    nmod_poly_init_mod(T, P->mod);
    nmod_poly_init_mod(f, P->mod);
    for (i = 0; i <= FLINT_BITS; i++)
    {
        nmod_poly_init_mod(stack + i, P->mod);
    }

    roots_idx = 0;

    nmod_poly_make_monic(f, P);

    a = stack + 0;
    nmod_poly_zero(a);
    nmod_poly_set_coeff_ui(a, 1, 1);
    nmod_poly_powmod_ui_binexp(T, a, (P->mod.n - 1)/2, f);
    nmod_poly_zero(a);
    nmod_poly_set_coeff_ui(a, 0, 1);
    nmod_poly_sub(T, T, a);
    nmod_poly_gcd(a, T, f);

    b = stack + 1;
    nmod_poly_zero(b);
    nmod_poly_set_coeff_ui(b, 0, 2);
    nmod_poly_add(T, T, b);
    nmod_poly_gcd(b, T, f);

    if (nmod_poly_degree(b) + nmod_poly_degree(a) != t)
    {
        success = 0;
        goto cleanup;
    }
    /* deg a >= deg b */
    if (nmod_poly_degree(a) < nmod_poly_degree(b))
    {
        nmod_poly_swap(a, b);
    }

    sp = nmod_poly_degree(b) > 0 ? 2 : 1;
    while (sp > 0)
    {
        FLINT_ASSERT(sp < FLINT_BITS);
        sp--;
        nmod_poly_swap(f, stack + sp);

        FLINT_ASSERT(nmod_poly_degree(f) > 0);
        if (nmod_poly_degree(f) == 1)
        {
            a0 = nmod_poly_get_coeff_ui(f, 0);
            a1 = nmod_poly_get_coeff_ui(f, 1);
            FLINT_ASSERT(a0 != 0);
            FLINT_ASSERT(a1 == 1);
            roots[roots_idx] = P->mod.n - a0;
            roots_idx++;
        }
        else
        {
            _nmod_poly_rabinsplit(stack + sp + 0, stack + sp + 1, T, f, randstate);
            FLINT_ASSERT(FLINT_BIT_COUNT(nmod_poly_degree(stack + sp + 1)) <= FLINT_BITS - sp - 1);
            sp += 2;
        }
    }

    success = 1;

cleanup:

    flint_randclear(randstate);
    nmod_poly_clear(T);
    nmod_poly_clear(f);
    for (i = 0; i <= FLINT_BITS; i++)
    {
        nmod_poly_clear(stack + i);
    }

    if (success)
    {
        FLINT_ASSERT(roots_idx == t);
    }

    return success;
}
