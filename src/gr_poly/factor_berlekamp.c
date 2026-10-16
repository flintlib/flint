/*
    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "gr_vec.h"
#include "gr_mat.h"
#include "gr_poly.h"
#include "factor_impl.h"

/* Factor a monic squarefree polynomial f of degree >= 1; append the
   irreducible factors to fac. */
static int
_gr_poly_factor_berlekamp_squarefree(gr_poly_vec_t fac, const gr_poly_t f,
    flint_rand_t state, const fmpz_t q, const fmpz_t halfq, slong k, gr_ctx_t ctx)
{
    const slong n = f->length - 1;
    slong sz = ctx->sizeof_elem;
    gr_poly_t xq, xqi, g, factor, power, t;
    gr_poly_preinv_t P;
    gr_mat_t M, X;
    gr_ptr r;
    slong i, j, nullity;
    int status = GR_SUCCESS;

    if (n <= 1)
        return (n == 1) ? gr_poly_vec_append(fac, f, ctx) : GR_SUCCESS;

    gr_poly_preinv_init(P, ctx);
    gr_poly_init(xq, ctx);
    gr_poly_init(xqi, ctx);
    gr_poly_init(g, ctx);
    gr_poly_init(factor, ctx);
    gr_poly_init(power, ctx);
    gr_poly_init(t, ctx);
    gr_mat_init(M, n, n, ctx);
    gr_mat_init(X, n, 0, ctx);
    GR_TMP_INIT(r, ctx);

    /* Step 1: x^q mod f */
    status |= gr_poly_preinv_set(P, f, ctx);
    status |= gr_poly_preinv_powmod_x_fmpz(xq, q, P, ctx);
    GR_POLY_FACTOR_CHECK_STATUS()

    /* Step 2: matrix of the Berlekamp map, Q - I, where column i holds
       x^(q i) mod f */
    status |= gr_poly_one(xqi, ctx);
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < xqi->length; j++)
            status |= gr_set(gr_mat_entry_ptr(M, j, i, ctx), GR_ENTRY(xqi->coeffs, j, sz), ctx);
        for (; j < n; j++)
            status |= gr_zero(gr_mat_entry_ptr(M, j, i, ctx), ctx);
        status |= gr_sub_ui(gr_mat_entry_ptr(M, i, i, ctx), gr_mat_entry_ptr(M, i, i, ctx), 1, ctx);
        if (i + 1 < n)
            status |= gr_poly_preinv_mulmod(xqi, xqi, xq, P, ctx);
    }
    GR_POLY_FACTOR_CHECK_STATUS()

    /* Step 3: basis of the Berlekamp subalgebra */
    status |= gr_mat_nullspace(X, M, ctx);
    GR_POLY_FACTOR_CHECK_STATUS()
    nullity = X->c;

    if (nullity <= 1)
    {
        /* irreducible (nullity == 1); nullity == 0 cannot happen */
        status |= gr_poly_vec_append(fac, f, ctx);
        goto cleanup;
    }

    /* Step 4: find a proper factor using random elements of the subalgebra */
    while (status == GR_SUCCESS)
    {
        /* random linear combination of the basis vectors */
        do
        {
            gr_poly_fit_length(factor, n, ctx);
            status |= _gr_vec_zero(factor->coeffs, n, ctx);

            for (i = 0; i < nullity; i++)
            {
                status |= gr_randtest(r, state, ctx);
                for (j = 0; j < n; j++)
                    status |= gr_addmul(GR_ENTRY(factor->coeffs, j, sz),
                        gr_mat_entry_srcptr(X, j, i, ctx), r, ctx);
            }

            _gr_poly_set_length_normalise(factor, n, ctx);

            if (status != GR_SUCCESS)
                goto cleanup;
        }
        while (factor->length <= 1);

        status |= gr_poly_gcd(g, f, factor, ctx);
        GR_POLY_FACTOR_CHECK(g)

        if (g->length > 1 && g->length < f->length)
            break;

        if (fmpz_is_odd(q))
        {
            status |= gr_poly_preinv_powmod_fmpz_sliding(power, factor, halfq, 0, P, ctx);
            status |= gr_poly_sub_ui(power, power, 1, ctx);
        }
        else
        {
            /* absolute trace to F_2 */
            status |= gr_poly_set(power, factor, ctx);
            status |= gr_poly_set(t, factor, ctx);
            for (i = 1; i < k; i++)
            {
                status |= gr_poly_preinv_mulmod(t, t, t, P, ctx);
                status |= gr_poly_add(power, power, t, ctx);
            }
        }

        status |= gr_poly_gcd(g, power, f, ctx);
        GR_POLY_FACTOR_CHECK(g)

        if (g->length > 1 && g->length < f->length)
            break;
    }

    GR_POLY_FACTOR_CHECK_STATUS()

    /* Step 5: recurse on g and f/g */
    status |= gr_poly_divexact(t, f, g, ctx);
    GR_POLY_FACTOR_CHECK(t)

    status |= _gr_poly_factor_berlekamp_squarefree(fac, g, state, q, halfq, k, ctx);
    status |= _gr_poly_factor_berlekamp_squarefree(fac, t, state, q, halfq, k, ctx);

cleanup:
    gr_poly_preinv_clear(P, ctx);
    gr_poly_clear(xq, ctx);
    gr_poly_clear(xqi, ctx);
    gr_poly_clear(g, ctx);
    gr_poly_clear(factor, ctx);
    gr_poly_clear(power, ctx);
    gr_poly_clear(t, ctx);
    gr_mat_clear(M, ctx);
    gr_mat_clear(X, ctx);
    GR_TMP_CLEAR(r, ctx);

    return status;
}

int
gr_poly_factor_berlekamp(gr_ptr c, gr_poly_vec_t fac, fmpz_vec_t exp,
    const gr_poly_t F, gr_ctx_t ctx)
{
    gr_poly_vec_t sqf;
    fmpz_vec_t sqf_exp;
    fmpz_t q, halfq;
    flint_rand_t state;
    slong i, j, k, num;
    int status = GR_SUCCESS;

    gr_poly_vec_set_length(fac, 0, ctx);
    fmpz_vec_set_length(exp, 0);

    if (F->length == 0)
        return gr_zero(c, ctx);

    fmpz_init(q);
    fmpz_init(halfq);

    status |= _gr_poly_factor_ff_info(q, NULL, &k, ctx);
    if (status != GR_SUCCESS)
    {
        fmpz_clear(q);
        fmpz_clear(halfq);
        return status;
    }

    if (fmpz_is_odd(q))
    {
        fmpz_sub_ui(halfq, q, 1);
        fmpz_fdiv_q_2exp(halfq, halfq, 1);
    }

    gr_poly_vec_init(sqf, 0, ctx);
    fmpz_vec_init(sqf_exp, 0);
    flint_rand_init(state);

    status |= gr_poly_factor_squarefree(c, sqf, sqf_exp, F, ctx);

    for (i = 0; i < sqf->length && status == GR_SUCCESS; i++)
    {
        num = fac->length;
        status |= _gr_poly_factor_berlekamp_squarefree(fac, sqf->entries + i, state, q, halfq, k, ctx);
        for (j = num; j < fac->length; j++)
            fmpz_vec_append(exp, sqf_exp->entries + i);
    }

    fmpz_clear(q);
    fmpz_clear(halfq);
    gr_poly_vec_clear(sqf, ctx);
    fmpz_vec_clear(sqf_exp);
    flint_rand_clear(state);

    return status;
}
