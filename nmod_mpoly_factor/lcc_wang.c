/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"


int nmod_mpoly_evaluate_all_n_poly(
    n_poly_t A,
    const nmod_mpoly_t B,
    const n_poly_struct * C,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, nvars = ctx->minfo->nvars;
    nmod_poly_t a;
    nmod_poly_struct * t1, ** t2;

    t1 = (nmod_poly_struct *) flint_malloc(nvars*sizeof(nmod_poly_struct));
    t2 = (nmod_poly_struct **) flint_malloc(nvars*sizeof(nmod_poly_struct *));

    nmod_poly_mock(a, A, ctx->ffinfo->mod);
    for (i = 0; i < nvars; i++)
    {
        nmod_poly_mock(t1 + i, C + i, ctx->ffinfo->mod);
        t2[i] = t1 + i;
    }

    nmod_mpoly_compose_nmod_poly(a, B, t2, ctx);

    n_poly_mock(A, a);

    flint_free(t1);
    flint_free(t2);

    return 1;
}

int nmod_mpoly_factor_lcc_wang(
    nmod_mpoly_struct * lc_divs,
    const nmod_mpoly_factor_t lcAfac,
    const n_poly_t Auc,
    const n_bpoly_struct * Auf,
    slong r,
    const n_poly_struct * alpha,
    const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong i, j, k;
    const slong n = ctx->minfo->nvars - 1;
    n_poly_struct * lcAfaceval;
    n_poly_struct * salpha;
    n_poly_struct * d;
    n_poly_t Q, R;
    nmod_mpoly_t t;

    salpha = (n_poly_struct *) flint_malloc((n + 1)*sizeof(n_poly_struct));
    n_poly_init(salpha + 0);
    for (i = 0; i < n; i++)
        salpha[i + 1] = alpha[i];

    lcAfaceval = (n_poly_struct *) flint_malloc(lcAfac->num*sizeof(n_poly_struct));
    for (i = 0; i < lcAfac->num; i++)
        n_poly_init(lcAfaceval + i);

    d = (n_poly_struct *) flint_malloc((lcAfac->num + 1)*sizeof(n_poly_struct));
    for (i = 0; i < lcAfac->num + 1; i++)
        n_poly_init(d + i);

    n_poly_init(Q);
    n_poly_init(R);

    nmod_mpoly_init(t, ctx);

    for (j = 0; j < lcAfac->num; j++)
        nmod_mpoly_evaluate_all_n_poly(lcAfaceval + j, lcAfac->poly + j, salpha, ctx);

    n_poly_set(d + 0, Auc);
    for (i = 0; i < lcAfac->num; i++)
    {
        n_poly_mod_make_monic(Q, lcAfaceval + i, ctx->ffinfo->mod);
        if (n_poly_degree(Q) < 1)
        {
            success = 0;
            goto cleanup;
        }
        for (j = i; j >= 0; j--)
        {
            n_poly_set(R, d + j);
            while (n_poly_degree(R) > 0)
            {
                n_poly_mod_gcd(R, R, Q, ctx->ffinfo->mod);
                n_poly_mod_div(Q, Q, R, ctx->ffinfo->mod);
                if (n_poly_degree(Q) < 1)
                {
                    success = 0;
                    goto cleanup;
                }
            }
        }
        n_poly_set(d + i + 1, Q);
    }

    for (j = 0; j < r; j++)
    {
        nmod_mpoly_one(lc_divs + j, ctx);
        n_poly_mod_mul(R, Auf[j].coeffs + Auf[j].length - 1, Auc, ctx->ffinfo->mod);
        for (i = lcAfac->num - 1; i >= 0; i--)
        {
            n_poly_mod_make_monic(Q, lcAfaceval + i, ctx->ffinfo->mod);
            if (n_poly_degree(Q) < 1)
                continue;
            k = n_poly_mod_remove(R, Q, ctx->ffinfo->mod);
            nmod_mpoly_pow_ui(t, lcAfac->poly + i, k, ctx);
            nmod_mpoly_mul(lc_divs + j, lc_divs + j, t, ctx);
        }
    }

    success = 1;

cleanup:

    n_poly_clear(Q);
    n_poly_clear(R);
    nmod_mpoly_clear(t, ctx);

    for (i = 0; i < lcAfac->num; i++)
        n_poly_clear(lcAfaceval + i);
    flint_free(lcAfaceval);

    for (i = 0; i < lcAfac->num + 1; i++)
        n_poly_clear(d + i);
    flint_free(d);

    n_poly_clear(salpha + 0);
    flint_free(salpha);

    return success;
}
