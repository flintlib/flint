/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2014 Fredrik Johansson
    Copyright (C) 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"
#include "gr.h"
#include "gr_poly.h"

void
_fmpz_mod_poly_div_series(fmpz * Q, const fmpz * A, slong Alen,
    const fmpz * B, slong Blen, slong n, const fmpz_mod_ctx_t ctx)
{
    fmpz_t u, d;

    fmpz_init(d);
    fmpz_init(u);

    if (!fmpz_is_one(B + 0))
    {
       fmpz_gcdinv(d, u, B + 0, fmpz_mod_ctx_modulus(ctx));

       if (!fmpz_is_one(d)) /* check for invertibility */
           flint_throw(FLINT_ERROR, "Impossible inverse in %s\n", __func__);
    }
    else
       fmpz_one(u);

    Alen = FLINT_MIN(Alen, n);
    Blen = FLINT_MIN(Blen, n);

    if (Blen == 1)
    {
        if (fmpz_is_one(B + 0))
            _fmpz_vec_set(Q, A, Alen);
        else
           _fmpz_mod_poly_scalar_mul_fmpz(Q, A, Alen, u, ctx);

        _fmpz_vec_zero(Q + Alen, n - Alen);
    }
    else if (n < 32 || Blen < 20)
    {
        slong i, j;

        if (fmpz_is_one(B + 0))
            fmpz_set(Q + 0, A + 0);
        else
        {
           fmpz_mod_mul(Q + 0, u, A + 0, ctx);
        }

        for (i = 1; i < n; i++)
        {
            fmpz_mul(Q + i, B + 1, Q + i - 1);

            for (j = 2; j < FLINT_MIN(i + 1, Blen); j++)
                fmpz_addmul(Q + i, B + j, Q + i - j);

            if (i < Alen)
                fmpz_sub(Q + i, A + i, Q + i);
            else
                fmpz_neg(Q + i, Q + i);

            if (!fmpz_is_one(B + 0))
                fmpz_mul(Q + i, Q + i, u);

            fmpz_mod_set_fmpz(Q + i, Q + i, ctx);
        }
    }
    else
    {
        gr_ctx_t gr_ctx;
        _gr_ctx_init_fmpz_mod_from_ref(gr_ctx, ctx);
        GR_MUST_SUCCEED(_gr_poly_div_series(Q, A, Alen, B, Blen, n, gr_ctx));
    }

    fmpz_clear(d);
    fmpz_clear(u);
}

void fmpz_mod_poly_div_series(fmpz_mod_poly_t Q, const fmpz_mod_poly_t A,
                    const fmpz_mod_poly_t B, slong n, const fmpz_mod_ctx_t ctx)
{
    slong Alen = FLINT_MIN(A->length, n);
    slong Blen = FLINT_MIN(B->length, n);

    if (Blen == 0)
        flint_throw(FLINT_ERROR, "Division by zero in %s\n", __func__);

    if (Alen == 0)
    {
        fmpz_mod_poly_zero(Q, ctx);
        return;
    }

    if (Q == A || Q == B)
    {
        fmpz_mod_poly_t t;
        fmpz_mod_poly_init2(t, n, ctx);
        _fmpz_mod_poly_div_series(t->coeffs, A->coeffs, Alen,
                                B->coeffs, Blen, n, ctx);
        fmpz_mod_poly_swap(Q, t, ctx);
        fmpz_mod_poly_clear(t, ctx);
    }
    else
    {
        fmpz_mod_poly_fit_length(Q, n, ctx);
        _fmpz_mod_poly_div_series(Q->coeffs, A->coeffs, Alen,
                                B->coeffs, Blen, n, ctx);
    }

    _fmpz_mod_poly_set_length(Q, n);
    _fmpz_mod_poly_normalise(Q);
}
