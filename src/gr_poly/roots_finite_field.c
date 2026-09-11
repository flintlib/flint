/*
    Copyright (C) 2020 Daniel Schultz
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
#include "gr_poly.h"
#include "factor_impl.h"

/* Append the roots of the monic squarefree polynomial f (with multiplicity
   e) to roots. */
static int
_gr_poly_push_roots(gr_vec_t roots, fmpz_vec_t mult, const gr_poly_t f,
    const fmpz_t e, flint_rand_t state, const fmpz_t q, gr_ctx_t ctx)
{
    gr_poly_t g, t, x, h;
    gr_poly_preinv_t P;
    gr_poly_vec_t lin;
    slong i, n;
    int status = GR_SUCCESS;

    n = f->length - 1;

    if (n < 1)
        return GR_SUCCESS;

    gr_poly_preinv_init(P, ctx);
    gr_poly_init(g, ctx);
    gr_poly_init(t, ctx);
    gr_poly_init(x, ctx);
    gr_poly_init(h, ctx);
    gr_poly_vec_init(lin, 0, ctx);

    if (n == 1)
    {
        status |= gr_poly_vec_append(lin, f, ctx);
    }
    else
    {
        /* Zero root */
        if (gr_is_zero(f->coeffs, ctx) == T_TRUE)
        {
            gr_vec_fit_length(roots, roots->length + 1, ctx);
            status |= gr_zero(gr_vec_entry_ptr(roots, roots->length, ctx), ctx);
            roots->length++;
            fmpz_vec_append(mult, e);

            i = 1;
            while (i < f->length && gr_is_zero(GR_ENTRY(f->coeffs, i, ctx->sizeof_elem), ctx) == T_TRUE)
                i++;
            status |= gr_poly_shift_right(g, f, i, ctx);
        }
        else
        {
            status |= gr_poly_set(g, f, ctx);
        }

        n = g->length - 1;

        if (n == 1)
        {
            status |= gr_poly_vec_append(lin, g, ctx);
        }
        else if (n > 1)
        {
            /* The product of the distinct nonzero linear factors is
               gcd(g, x^(q-1) - 1) = gcd(g, t - 1) gcd(g, t + 1) where
               t = x^((q-1)/2) for odd q (Rabin's initial split), or
               gcd(g, T) gcd(g, T + 1) where T = x + x^2 + ... + x^(2^(k-1))
               is the absolute trace for q = 2^k. */
            status |= gr_poly_preinv_set(P, g, ctx);
            status |= gr_poly_gen(x, ctx);

            if (fmpz_is_odd(q))
            {
                fmpz_t halfq;
                fmpz_init(halfq);
                fmpz_sub_ui(halfq, q, 1);
                fmpz_fdiv_q_2exp(halfq, halfq, 1);
                status |= gr_poly_preinv_powmod_x_fmpz(t, halfq, P, ctx);
                fmpz_clear(halfq);
                status |= gr_poly_sub_ui(t, t, 1, ctx);
            }
            else
            {
                slong k = fmpz_val2(q);
                status |= gr_poly_set(t, x, ctx);
                status |= gr_poly_set(h, x, ctx);
                for (i = 1; i < k && status == GR_SUCCESS; i++)
                {
                    status |= gr_poly_preinv_mulmod(h, h, h, P, ctx);
                    status |= gr_poly_add(t, t, h, ctx);
                }
            }

            status |= gr_poly_gcd(h, t, g, ctx);
            GR_POLY_FACTOR_CHECK(h)
            if (h->length > 1)
                status |= _gr_poly_factor_equal_deg_with_frob(lin, h, 1, t, state, ctx);

            if (fmpz_is_odd(q))
                status |= gr_poly_add_ui(t, t, 2, ctx);
            else
                status |= gr_poly_add_ui(t, t, 1, ctx);

            status |= gr_poly_gcd(h, t, g, ctx);
            GR_POLY_FACTOR_CHECK(h)
            if (h->length > 1)
                status |= _gr_poly_factor_equal_deg_with_frob(lin, h, 1, t, state, ctx);
        }
    }

    GR_POLY_FACTOR_CHECK_STATUS()

    for (i = 0; i < lin->length; i++)
    {
        gr_poly_struct * l = lin->entries + i;
        gr_vec_fit_length(roots, roots->length + 1, ctx);
        status |= gr_neg(gr_vec_entry_ptr(roots, roots->length, ctx), l->coeffs, ctx);
        roots->length++;
        fmpz_vec_append(mult, e);
    }

cleanup:
    gr_poly_preinv_clear(P, ctx);
    gr_poly_clear(g, ctx);
    gr_poly_clear(t, ctx);
    gr_poly_clear(x, ctx);
    gr_poly_clear(h, ctx);
    gr_poly_vec_clear(lin, ctx);

    return status;
}

int
gr_poly_roots_finite_field(gr_vec_t roots, fmpz_vec_t mult, const gr_poly_t poly,
    int FLINT_UNUSED(flags), gr_ctx_t ctx)
{
    gr_poly_vec_t sqf;
    fmpz_vec_t sqf_exp;
    fmpz_t q;
    gr_ptr c;
    flint_rand_t state;
    slong i;
    int status = GR_SUCCESS;

    gr_vec_set_length(roots, 0, ctx);
    fmpz_vec_set_length(mult, 0);

    if (poly->length == 0)
        return GR_DOMAIN;

    fmpz_init(q);
    status |= _gr_poly_factor_ff_info(q, NULL, NULL, ctx);
    if (status != GR_SUCCESS)
    {
        fmpz_clear(q);
        return status;
    }

    if (poly->length == 1)
    {
        fmpz_clear(q);
        return GR_SUCCESS;
    }

    gr_poly_vec_init(sqf, 0, ctx);
    fmpz_vec_init(sqf_exp, 0);
    c = gr_heap_init(ctx);
    flint_rand_init(state);

    status |= gr_poly_factor_squarefree(c, sqf, sqf_exp, poly, ctx);

    for (i = 0; i < sqf->length && status == GR_SUCCESS; i++)
        status |= _gr_poly_push_roots(roots, mult, sqf->entries + i, sqf_exp->entries + i, state, q, ctx);

    fmpz_clear(q);
    gr_poly_vec_clear(sqf, ctx);
    fmpz_vec_clear(sqf_exp);
    gr_heap_clear(c, ctx);
    flint_rand_clear(state);

    return status;
}
