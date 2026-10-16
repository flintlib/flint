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
#include "gr_poly.h"
#include "factor_impl.h"

/*
    Factor a monic squarefree polynomial f of degree n >= 1; append the
    factors to fac with exponent e appended to exp.

    This is written at the level of coefficient arrays with a single
    scratch block rather than with gr_poly_t temporaries, which roughly
    halves the overhead for very small degrees (where the whole
    factorization takes a few hundred nanoseconds over nmod).
*/
static int
_gr_poly_factor_cantor_zassenhaus_squarefree(gr_poly_vec_t fac, fmpz_vec_t exp,
    gr_srcptr f, slong n, const fmpz_t e, flint_rand_t state, const fmpz_t q, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    gr_ptr W, v, h, t, g, xq, tmp;
    slong lenv, lenh, lent, leng, lenxq, i, j;
    gr_poly_preinv_t P;
    int have_P;
    gr_poly_vec_t ed;
    gr_poly_t gg, frob;
    int status = GR_SUCCESS;

    if (n == 1)
    {
        gr_poly_vec_fit_length(fac, fac->length + 1, ctx);
        gr_poly_fit_length(fac->entries + fac->length, 2, ctx);
        status |= _gr_vec_set(fac->entries[fac->length].coeffs, f, 2, ctx);
        _gr_poly_set_length(fac->entries + fac->length, 2, ctx);
        fac->length++;
        fmpz_vec_append(exp, e);
        return status;
    }

    gr_poly_preinv_init(P, ctx);
    have_P = 0;

    /* v (n+1), h (n), t (n), g (n+1), xq (n), tmp (n+1) */
    GR_TMP_INIT_VEC(W, 6 * n + 3, ctx);
    v = W;
    h = GR_ENTRY(v, n + 1, sz);
    t = GR_ENTRY(h, n, sz);
    g = GR_ENTRY(t, n, sz);
    xq = GR_ENTRY(g, n + 1, sz);
    tmp = GR_ENTRY(xq, n, sz);

    status |= _gr_vec_set(v, f, n + 1, ctx);
    lenv = n + 1;

    /* h = x (n >= 2 so this is reduced mod v) */
    status |= _gr_vec_zero(h, n, ctx);
    status |= gr_one(GR_ENTRY(h, 1, sz), ctx);
    lenh = 2;
    lenxq = 0;

    i = 0;
    do
    {
        i++;

        /* h = h^q mod v */
        if (!have_P)
        {
            status |= _gr_poly_preinv_set(P, v, lenv, ctx);
            have_P = 1;
        }

        status |= _gr_poly_preinv_powmod_fmpz_binexp(t, h, lenh, q, P, ctx);
        { gr_ptr s = h; h = t; t = s; }
        lenh = lenv - 1;
        GR_IGNORE(_gr_vec_normalise(&lenh, h, lenh, ctx));

        if (status != GR_SUCCESS)
            break;

        if (i == 1)
        {
            /* keep x^q mod (the original) v for equal degree factorization */
            status |= _gr_vec_set(xq, h, lenh, ctx);
            lenxq = lenh;
        }

        /* t = h - x */
        lent = FLINT_MAX(lenh, 2);
        status |= _gr_vec_set(t, h, lenh, ctx);
        status |= _gr_vec_zero(GR_ENTRY(t, lenh, sz), lent - lenh, ctx);
        status |= gr_sub_ui(GR_ENTRY(t, 1, sz), GR_ENTRY(t, 1, sz), 1, ctx);
        GR_IGNORE(_gr_vec_normalise(&lent, t, lent, ctx));

        /* g = gcd(v, t), made monic; g is the product of the irreducible
           factors of v of degree dividing i (and hence of degree i) */
        if (lent == 0)
        {
            status |= _gr_vec_set(g, v, lenv, ctx);
            leng = lenv;
        }
        else
        {
            status |= _gr_poly_gcd(g, &leng, v, lenv, t, lent, ctx);
            if (leng > 1)
                status |= _gr_poly_make_monic(g, g, leng, ctx);
        }

        if (status != GR_SUCCESS)
            break;

        if (leng > 1)
        {
            gr_poly_init(gg, ctx);
            gr_poly_fit_length(gg, leng, ctx);
            status |= _gr_vec_set(gg->coeffs, g, leng, ctx);
            _gr_poly_set_length(gg, leng, ctx);

            gr_poly_vec_init(ed, 0, ctx);

            if (leng - 1 == i)
            {
                gr_poly_vec_append_swap(ed, gg, ctx);
            }
            else
            {
                /* frob = x^q mod g */
                gr_poly_init(frob, ctx);
                gr_poly_fit_length(frob, leng - 1, ctx);
                if (lenxq >= leng)
                    status |= _gr_poly_divrem(tmp, frob->coeffs, xq, lenxq, g, leng, ctx);
                else
                    status |= _gr_vec_set(frob->coeffs, xq, lenxq, ctx);
                _gr_poly_set_length_normalise(frob, FLINT_MIN(lenxq, leng - 1), ctx);
                status |= _gr_poly_factor_equal_deg_with_frob(ed, gg, i, frob, state, ctx);
                gr_poly_clear(frob, ctx);
            }

            for (j = 0; j < ed->length; j++)
            {
                gr_poly_vec_append_swap(fac, ed->entries + j, ctx);
                fmpz_vec_append(exp, e);
            }

            gr_poly_vec_clear(ed, ctx);
            gr_poly_clear(gg, ctx);

            if (status != GR_SUCCESS)
                break;

            /* v = v / g */
            status |= _gr_poly_divexact(tmp, v, lenv, g, leng, ctx);
            lenv = lenv - leng + 1;
            status |= _gr_vec_set(v, tmp, lenv, ctx);
            have_P = 0;

            /* h = h mod v */
            if (lenv > 1 && lenh >= lenv)
            {
                status |= _gr_poly_divrem(tmp, t, h, lenh, v, lenv, ctx);
                status |= _gr_vec_set(h, t, lenv - 1, ctx);
                lenh = lenv - 1;
                GR_IGNORE(_gr_vec_normalise(&lenh, h, lenh, ctx));
            }
        }
    }
    while (lenv >= 2 * i + 3 && status == GR_SUCCESS);

    if (status == GR_SUCCESS && lenv > 1)
    {
        gr_poly_vec_fit_length(fac, fac->length + 1, ctx);
        gr_poly_fit_length(fac->entries + fac->length, lenv, ctx);
        status |= _gr_vec_set(fac->entries[fac->length].coeffs, v, lenv, ctx);
        _gr_poly_set_length(fac->entries + fac->length, lenv, ctx);
        fac->length++;
        fmpz_vec_append(exp, e);
    }

    GR_TMP_CLEAR_VEC(W, 6 * n + 3, ctx);
    gr_poly_preinv_clear(P, ctx);

    return status;
}

int
gr_poly_factor_cantor_zassenhaus(gr_ptr c, gr_poly_vec_t fac, fmpz_vec_t exp,
    const gr_poly_t F, gr_ctx_t ctx)
{
    gr_poly_vec_t sqf;
    fmpz_vec_t sqf_exp;
    fmpz_t q;
    flint_rand_t state;
    slong i;
    int status = GR_SUCCESS;

    gr_poly_vec_set_length(fac, 0, ctx);
    fmpz_vec_set_length(exp, 0);

    status |= _gr_poly_factor_ff_info(NULL, NULL, NULL, ctx);
    if (status != GR_SUCCESS)
        return status;

    if (F->length == 0)
        return gr_zero(c, ctx);

    fmpz_init(q);
    gr_poly_vec_init(sqf, 0, ctx);
    fmpz_vec_init(sqf_exp, 0);
    flint_rand_init(state);

    status |= gr_ctx_fq_order(q, ctx);
    status |= gr_poly_factor_squarefree(c, sqf, sqf_exp, F, ctx);

    for (i = 0; i < sqf->length && status == GR_SUCCESS; i++)
        status |= _gr_poly_factor_cantor_zassenhaus_squarefree(fac, exp,
            sqf->entries[i].coeffs, sqf->entries[i].length - 1,
            sqf_exp->entries + i, state, q, ctx);

    fmpz_clear(q);
    gr_poly_vec_clear(sqf, ctx);
    fmpz_vec_clear(sqf_exp);
    flint_rand_clear(state);

    return status;
}
