/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2020 Claus Fieker
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "nmod_poly.h"
#include "fmpz_vec.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"
#include "padic.h"
#include "qadic.h"

int _qadic_ctx_init_conway_ui(qadic_ctx_t ctx, ulong p, slong d,
        slong min, slong max, const char * var, enum padic_print_mode mode)
{
    int result;
    ulong tmp[410]; /* Largest degree is 409 */
    slong num_nzcoeffs;
    slong * idx_nzcoeffs;
    nn_ptr nzcoeffs;
    slong ix, jx;

    result = _nmod_poly_conway(tmp, p, d);

    if (!result)
        return 0;

    /* Find number of non-zero coefficients */
    num_nzcoeffs = 2;

    for (ix = 1; ix < d; ix++)
        if (tmp[ix])
            num_nzcoeffs++;

    idx_nzcoeffs = flint_malloc(sizeof(slong) * num_nzcoeffs);
    nzcoeffs = flint_malloc(sizeof(fmpz) * num_nzcoeffs);

    /* Copy the polynomial */
    for (ix = 0, jx = 0; ix < d; ix++)
    {
        if (tmp[ix])
        {
            nzcoeffs[jx] = tmp[ix];
            idx_nzcoeffs[jx] = ix;
            jx++;
        }
    }

    nzcoeffs[jx] = 1;
    idx_nzcoeffs[jx] = d;

    ctx->len = num_nzcoeffs;
    ctx->j = idx_nzcoeffs;
    ctx->a = (fmpz *) nzcoeffs;

    /* Complete the initialisation of the context */
    padic_ctx_init(&ctx->pctx, (fmpz *) &p, min, max, mode);

    ctx->var = flint_malloc(strlen(var) + 1);
    strcpy(ctx->var, var);

    return 1;
}

void qadic_ctx_init_conway(qadic_ctx_t ctx, const fmpz_t p, slong d,
        slong min, slong max, const char * var, enum padic_print_mode mode)
{
    int result;

    if (*p < 2 || *p > 109987) /* NOTE: This works. */
        flint_throw(FLINT_ERROR, "Exception (qadic_ctx_init_conway).  "
                "Conway polynomials are only available for primes up to 109987.\n");

    result = _qadic_ctx_init_conway_ui(ctx, *p, d, min, max, var, mode);

    if (!result)
        flint_throw(FLINT_ERROR, "Exception (qadic_ctx_init_conway).  "
                "The polynomial for (p, d) = (%{fmpz}, %wd) is not present "
                "in the database.\n", p, d);
}

void qadic_ctx_init(qadic_ctx_t ctx, const fmpz_t p, slong d,
        slong min, slong max, const char * var, enum padic_print_mode mode)
{
    flint_rand_t state;
    fmpz_mod_poly_t poly;
    slong i, j;
    fmpz_mod_ctx_t ctxp;

    if (*p >= 2 && *p <= 109987)
        if (_qadic_ctx_init_conway_ui(ctx, *p, d, min, max, var, mode))
            return;

    flint_rand_init(state);

    fmpz_mod_ctx_init(ctxp, p);
    fmpz_mod_poly_init2(poly, d + 1, ctxp);

    fmpz_mod_poly_randtest_sparse_irreducible(poly, state, d + 1, ctxp);

    flint_rand_clear(state);

    /* Find number of non-zero coefficients */
    ctx->len = 1;

    for (i = 0; i < d; i++)
    {
        if (!fmpz_is_zero(poly->coeffs + i))
            ctx->len ++;
    }

    ctx->a = _fmpz_vec_init(ctx->len);
    ctx->j = flint_malloc(ctx->len*sizeof(slong));

    /* Copy the polynomial */
    j = 0;

    for (i = 0; i < d; i++)
    {
        if (!fmpz_is_zero(poly->coeffs+i))
        {
            fmpz_set(ctx->a + j, poly->coeffs + i);
            ctx->j[j] = i;
            j++;
        }
    }

    fmpz_set_ui(ctx->a + j, 1);
    ctx->j[j] = d;

    /* Complete the initialisation of the context */
    padic_ctx_init(&ctx->pctx, p, min, max, mode);

    ctx->var = flint_malloc(strlen(var) + 1);
    strcpy(ctx->var, var);

    fmpz_mod_poly_clear(poly, ctxp);
    fmpz_mod_ctx_clear(ctxp);
}
