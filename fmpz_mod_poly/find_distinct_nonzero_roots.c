/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_poly.h"
#include "fmpz_mod.h"

/* split f assuming that f has degree(f) distinct nonzero roots in Fp */
static void _fmpz_mod_poly_rabinsplit(
    fmpz_mod_poly_t a,
    fmpz_mod_poly_t b,
    fmpz_mod_poly_t T,
    const fmpz_mod_poly_t f,
    flint_rand_t randstate)
{
    fmpz_t delta;

    fmpz_init(delta);

    FLINT_ASSERT(fmpz_mod_poly_degree(f) > 1);

try_again:

    fmpz_randm(delta, randstate, &f->p);

    fmpz_mod_poly_zero(a);
    fmpz_mod_poly_set_coeff_ui(a, 1, 1);
    fmpz_mod_poly_set_coeff_fmpz(a, 0, delta);
    fmpz_sub_ui(delta, &f->p, 1);
    fmpz_divexact_ui(delta, delta, 2);
    fmpz_mod_poly_powmod_fmpz_binexp(T, a, delta, f);
    fmpz_mod_poly_zero(a);
    fmpz_mod_poly_set_coeff_ui(a, 0, 1);
    fmpz_mod_poly_sub(T, T, a);
    fmpz_mod_poly_gcd(a, T, f);
    FLINT_ASSERT(!fmpz_mod_poly_is_zero(a));
    if (0 >= fmpz_mod_poly_degree(a) || fmpz_mod_poly_degree(a) >= fmpz_mod_poly_degree(f))
    {
        goto try_again;
    }
    fmpz_mod_poly_div_basecase(b, f, a);
    /* deg a >= deg b */
    if (fmpz_mod_poly_degree(a) < fmpz_mod_poly_degree(b))
    {
        fmpz_mod_poly_swap(a, b);
    }

    fmpz_clear(delta);
    return;
}

/*
    If P has deg(P) distinct nonzero roots of P, fill them in and return 1.
    Otherwise return 0. Function is undefined for zero P.
    The modulus of P is assumed to be prime.
*/
int fmpz_mod_poly_find_distinct_nonzero_roots(
    fmpz * roots,
    const fmpz_mod_poly_t P)
{
    fmpz_t a0, a1;
    int success;
    slong i, roots_idx, sp;
    fmpz_t halfp;
    fmpz_mod_poly_struct * a , * b;
    fmpz_mod_poly_t f, T;
    fmpz_mod_poly_struct stack[FLINT_BITS + 1];
    flint_rand_t randstate;
    slong t = fmpz_mod_poly_degree(P);
    fmpz_mod_ctx_t fpctx;

    FLINT_ASSERT(t >= 0);

    fmpz_mod_ctx_init(fpctx, &P->p);
    fmpz_init(a0);
    fmpz_init(a1);
    fmpz_init(halfp);

    if (t < 2)
    {
        if (t == 1)
        {
            fmpz_mod_poly_get_coeff_fmpz(a0, P, 0);
            fmpz_mod_poly_get_coeff_fmpz(a1, P, 1);
            if (fmpz_is_zero(a0))
            {
                success = 0;
                goto cleanup1;
            }
            fmpz_mod_inv(a1, a1, fpctx);
            fmpz_mod_neg(a1, a1, fpctx);
            fmpz_mod_mul(roots + 0, a0, a1, fpctx);
        }
        success = 1;
        goto cleanup1;
    }

    if (fmpz_equal_ui(&P->p, 2))
    {
        success = 0;
        goto cleanup1;
    }

    flint_randinit(randstate);
    fmpz_mod_poly_init(T, &P->p);
    fmpz_mod_poly_init(f, &P->p);
    for (i = 0; i <= FLINT_BITS; i++)
    {
        fmpz_mod_poly_init(stack + i, &P->p);
    }

    roots_idx = 0;

    fmpz_mod_poly_make_monic(f, P);

    a = stack + 0;
    fmpz_mod_poly_zero(a);
    fmpz_mod_poly_set_coeff_ui(a, 1, 1);
    fmpz_sub_ui(halfp, &P->p, 1);
    fmpz_divexact_ui(halfp, halfp, 2);
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

    if (fmpz_mod_poly_degree(b) + fmpz_mod_poly_degree(a) != t)
    {
        success = 0;
        goto cleanup;
    }
    /* deg a >= deg b */
    if (fmpz_mod_poly_degree(a) < fmpz_mod_poly_degree(b))
    {
        fmpz_mod_poly_swap(a, b);
    }

    sp = fmpz_mod_poly_degree(b) > 0 ? 2 : 1;
    while (sp > 0)
    {
        FLINT_ASSERT(sp < FLINT_BITS);
        sp--;
        fmpz_mod_poly_swap(f, stack + sp);

        FLINT_ASSERT(fmpz_mod_poly_degree(f) > 0);
        if (fmpz_mod_poly_degree(f) == 1)
        {
            fmpz_mod_poly_get_coeff_fmpz(a0, f, 0);
            fmpz_mod_poly_get_coeff_fmpz(a1, f, 1);
            FLINT_ASSERT(!fmpz_is_zero(a0));
            FLINT_ASSERT(fmpz_is_one(a1));
            fmpz_mod_neg(roots + roots_idx, a0, fpctx);
            roots_idx++;
        }
        else
        {
            _fmpz_mod_poly_rabinsplit(stack + sp + 0, stack + sp + 1, T, f, randstate);
            FLINT_ASSERT(FLINT_BIT_COUNT(fmpz_mod_poly_degree(stack + sp + 1)) <= FLINT_BITS - sp - 1);
            sp += 2;
        }
    }

    success = 1;

cleanup:

    flint_randclear(randstate);
    fmpz_mod_poly_clear(T);
    fmpz_mod_poly_clear(f);
    for (i = 0; i <= FLINT_BITS; i++)
    {
        fmpz_mod_poly_clear(stack + i);
    }

    if (success)
    {
        FLINT_ASSERT(roots_idx == t);
    }

cleanup1:

    fmpz_mod_ctx_clear(fpctx);
    fmpz_clear(a0);
    fmpz_clear(a1);
    fmpz_clear(halfp);

    return success;
}
