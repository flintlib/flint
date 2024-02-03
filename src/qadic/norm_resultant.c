/*
    Copyright (C) 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "qadic.h"
#include "fmpz_mod.h"
#include "fmpz_mod_mat.h"

void _qadic_norm_resultant(fmpz_t rop, const fmpz *op, slong len,
                           const fmpz *a, const slong *j, slong lena,
                           const fmpz_t p, slong N)
{
    const slong d = j[lena - 1];

    fmpz_t pN;

    fmpz_init(pN);
    fmpz_pow_ui(pN, p, N);

    if (len == 1)
    {
        fmpz_powm_ui(rop, op + 0, d, pN);
    }
    else  /* len >= 2 */
    {
        fmpz_mod_ctx_t ctx;
        fmpz_mod_mat_t M;
        slong n = d + len - 1, i, k;

        /* todo: should use fmpz_mod_poly_resultant when that exists */
        fmpz_mod_ctx_init(ctx, pN);
        fmpz_mod_mat_init(M, n, n, ctx);

        for (k = 0; k < len-1; k++)
            for (i = 0; i < lena; i++)
                fmpz_mod_set_fmpz(fmpz_mod_mat_entry(M, k, k + d - j[i]), a + i, ctx);
        for (k = 0; k < d; k++)
            for (i = 0; i < len; i++)
                fmpz_mod_set_fmpz(fmpz_mod_mat_entry(M, len - 1 + k, k + len -1 - i), op + i, ctx);

        fmpz_mod_mat_det(rop, M, ctx);

        fmpz_mod_mat_clear(M, ctx);
        fmpz_mod_ctx_clear(ctx);

        /*
            XXX:  This part of the code is currently untested as the Conway
            polynomials used for the extension Qq/Qp are monic.
         */
        if (!fmpz_is_one(a + (lena - 1)))
        {
            fmpz_t f;

            fmpz_init(f);
            fmpz_powm_ui(f, a + (lena - 1), len - 1, pN);
            _padic_inv(f, f, p, N);
            fmpz_mul(rop, f, rop);
            fmpz_mod(rop, rop, pN);
            fmpz_clear(f);
        }
    }
    fmpz_clear(pN);
}

void qadic_norm_resultant(padic_t rop, const qadic_t op, const qadic_ctx_t ctx)
{
    const slong N = padic_prec(rop);
    const slong d = qadic_ctx_degree(ctx);

    /* N(p^v u) = p^{dv} N(u) */

    if (qadic_is_zero(op) || d * op->val >= N)
    {
        padic_zero(rop);
    }
    else
    {
        _qadic_norm_resultant(padic_unit(rop), op->coeffs, op->length,
                              ctx->a, ctx->j, ctx->len, (&ctx->pctx)->p,
                              N - d * op->val);
        padic_val(rop) = d * op->val;
    }
}
