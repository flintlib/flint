/*
    Copyright (C) 2026 Mael Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpq.h"
#include "fmpz_mpoly_factor.h"
#include "gr_poly.h"

/* The algorithm itself is _fmpz_bpoly_resultant, over the dense bivariate
   type of the fmpz_mpoly module. What is left here is the conversion: a
   polynomial over a polynomial ring over Z has the layout of an fmpz_bpoly
   already, but one over Q has to have its denominators cleared first, so both
   are copied out rather than cast. */

int
_gr_poly_resultant_modular(gr_ptr res, gr_srcptr A, slong lenA,
                           gr_srcptr B, slong lenB, int proved, gr_ctx_t ctx)
{
    const gr_poly_struct * Ax = A;
    const gr_poly_struct * Bx = B;
    gr_poly_struct * resx = res;
    gr_ctx_struct * cctx;
    fmpz_poly_struct * Az, * Bz;
    fmpz_poly_t R;
    fmpz_t da, db, t;
    slong i, k;
    int rational;
    int status;

    if (ctx->which_ring != GR_CTX_GR_POLY)
        return GR_UNABLE;

    cctx = POLYNOMIAL_ELEM_CTX(ctx);

    if (cctx->which_ring == GR_CTX_FMPZ)
        rational = 0;
    else if (cctx->which_ring == GR_CTX_FMPQ)
        rational = 1;
    else
        return GR_UNABLE;

    if (lenB <= 1)
        return _gr_poly_resultant_small(res, A, lenA, B, lenB, ctx);

    if (Ax[lenA - 1].length == 0 || Bx[lenB - 1].length == 0)
        return GR_UNABLE;

    fmpz_init_set_ui(da, 1);
    fmpz_init_set_ui(db, 1);
    fmpz_init(t);

    Az = flint_malloc(lenA * sizeof(fmpz_poly_struct));
    Bz = flint_malloc(lenB * sizeof(fmpz_poly_struct));

    for (i = 0; i < lenA; i++)
        fmpz_poly_init(Az + i);
    for (i = 0; i < lenB; i++)
        fmpz_poly_init(Bz + i);

    if (!rational)
    {
        /* a gr_poly over fmpz has the same layout as an fmpz_poly */
        for (i = 0; i < lenA; i++)
            fmpz_poly_set(Az + i, (const fmpz_poly_struct *) (Ax + i));
        for (i = 0; i < lenB; i++)
            fmpz_poly_set(Bz + i, (const fmpz_poly_struct *) (Bx + i));
    }
    else
    {
        /* clear denominators; res(A/da, B/db) = res(A, B) / (da^n db^m) */
        for (i = 0; i < lenA; i++)
            for (k = 0; k < Ax[i].length; k++)
                fmpz_lcm(da, da, fmpq_denref(((const fmpq *) Ax[i].coeffs) + k));

        for (i = 0; i < lenB; i++)
            for (k = 0; k < Bx[i].length; k++)
                fmpz_lcm(db, db, fmpq_denref(((const fmpq *) Bx[i].coeffs) + k));

        for (i = 0; i < lenA; i++)
        {
            fmpz_poly_fit_length(Az + i, Ax[i].length);
            for (k = 0; k < Ax[i].length; k++)
            {
                const fmpq * c = ((const fmpq *) Ax[i].coeffs) + k;
                fmpz_divexact(t, da, fmpq_denref(c));
                fmpz_mul(Az[i].coeffs + k, fmpq_numref(c), t);
            }
            _fmpz_poly_set_length(Az + i, Ax[i].length);
            _fmpz_poly_normalise(Az + i);
        }

        for (i = 0; i < lenB; i++)
        {
            fmpz_poly_fit_length(Bz + i, Bx[i].length);
            for (k = 0; k < Bx[i].length; k++)
            {
                const fmpq * c = ((const fmpq *) Bx[i].coeffs) + k;
                fmpz_divexact(t, db, fmpq_denref(c));
                fmpz_mul(Bz[i].coeffs + k, fmpq_numref(c), t);
            }
            _fmpz_poly_set_length(Bz + i, Bx[i].length);
            _fmpz_poly_normalise(Bz + i);
        }
    }

    fmpz_poly_init(R);
    status = _fmpz_bpoly_resultant(R, Az, lenA, Bz, lenB, proved)
        ? GR_SUCCESS : GR_UNABLE;

    if (status == GR_SUCCESS)
    {
        gr_poly_fit_length(resx, R->length, cctx);

        if (!rational)
        {
            for (k = 0; k < R->length; k++)
                fmpz_set(((fmpz *) resx->coeffs) + k, R->coeffs + k);
        }
        else
        {
            fmpz_pow_ui(da, da, lenB - 1);
            fmpz_pow_ui(db, db, lenA - 1);
            fmpz_mul(da, da, db);

            for (k = 0; k < R->length; k++)
                fmpq_set_fmpz_frac(((fmpq *) resx->coeffs) + k, R->coeffs + k, da);
        }

        _gr_poly_set_length(resx, R->length, cctx);
        _gr_poly_normalise(resx, cctx);
    }

    fmpz_poly_clear(R);

    for (i = 0; i < lenA; i++)
        fmpz_poly_clear(Az + i);
    for (i = 0; i < lenB; i++)
        fmpz_poly_clear(Bz + i);

    flint_free(Az);
    flint_free(Bz);

    fmpz_clear(da);
    fmpz_clear(db);
    fmpz_clear(t);

    return status;
}

int
gr_poly_resultant_modular(gr_ptr r, const gr_poly_t f,
                          const gr_poly_t g, int proved, gr_ctx_t ctx)
{
    slong len1 = f->length;
    slong len2 = g->length;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    if (len1 == 0 || len2 == 0)
    {
        return gr_zero(r, ctx);
    }

    if (gr_is_zero(GR_ENTRY(f->coeffs, len1 - 1, sz), ctx) != T_FALSE ||
        gr_is_zero(GR_ENTRY(g->coeffs, len2 - 1, sz), ctx) != T_FALSE)
    {
        return GR_UNABLE;
    }

    if (len1 >= len2)
    {
        status |= _gr_poly_resultant_modular(r, f->coeffs, len1, g->coeffs, len2, proved, ctx);
    }
    else
    {
        status |= _gr_poly_resultant_modular(r, g->coeffs, len2, f->coeffs, len1, proved, ctx);

        if (((len1 | len2) & 1) == 0)
            status |= gr_neg(r, r, ctx);
    }

    return status;
}
