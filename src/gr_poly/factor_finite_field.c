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

#include <math.h>
#include "fmpz.h"
#include "fmpz_vec.h"
#include "gr_vec.h"
#include "gr_poly.h"
#include "factor_impl.h"

static int
_gr_poly_factor_no_deflation(gr_ptr c, gr_poly_vec_t fac, fmpz_vec_t exp,
    const gr_poly_t F, int algorithm, gr_ctx_t ctx)
{
    if (algorithm == GR_POLY_FACTOR_ALGORITHM_DEFAULT)
    {
        fmpz_t q;
        slong n = F->length - 1;
        double log2q;
        int status;

        fmpz_init(q);
        status = _gr_poly_factor_ff_info(q, NULL, NULL, ctx);
        log2q = fmpz_dlog(q) * 1.4426950408889634;
        fmpz_clear(q);

        if (status != GR_SUCCESS)
            return status;

        /* Cantor-Zassenhaus is competitive for small q (the Frobenius
           can be computed with a few squarings); Kaltofen-Shoup wins
           as q grows. Empirical crossovers for dense random input
           (nmod): n ~ 125 for q = 2, n ~ 70 for q = 3, n ~ 32 for q = 7,
           n ~ 16 for q = 31, n ~ 8 for q ~ 2^8, n ~ 4 for q ~ 2^16. */
        if (n < 125.0 / pow(log2q, 1.3))
            algorithm = GR_POLY_FACTOR_ALGORITHM_CANTOR_ZASSENHAUS;
        else
            algorithm = GR_POLY_FACTOR_ALGORITHM_KALTOFEN_SHOUP;
    }

    if (algorithm == GR_POLY_FACTOR_ALGORITHM_CANTOR_ZASSENHAUS)
        return gr_poly_factor_cantor_zassenhaus(c, fac, exp, F, ctx);
    else if (algorithm == GR_POLY_FACTOR_ALGORITHM_BERLEKAMP)
        return gr_poly_factor_berlekamp(c, fac, exp, F, ctx);
    else
        return gr_poly_factor_kaltofen_shoup(c, fac, exp, F, ctx);
}

int
_gr_poly_factor_finite_field(gr_ptr c, gr_poly_vec_t fac, fmpz_vec_t exp,
    const gr_poly_t F, int algorithm, gr_ctx_t ctx)
{
    ulong deflation;
    int status = GR_SUCCESS;

    gr_poly_vec_set_length(fac, 0, ctx);
    fmpz_vec_set_length(exp, 0);

    if (F->length == 0)
        return gr_zero(c, ctx);

    if (F->length == 1)
        return gr_set(c, F->coeffs, ctx);

    if (gr_is_zero(gr_poly_coeff_srcptr(F, F->length - 1, ctx), ctx) != T_FALSE)
        return GR_UNABLE;

    deflation = gr_poly_deflation(F, ctx);

    if (deflation == 1)
    {
        status |= _gr_poly_factor_no_deflation(c, fac, exp, F, algorithm, ctx);
    }
    else
    {
        gr_poly_vec_t dfac, ffac;
        fmpz_vec_t dexp, fexp;
        gr_poly_t def, pol;
        gr_ptr cc;
        slong i, j;

        gr_poly_vec_init(dfac, 0, ctx);
        gr_poly_vec_init(ffac, 0, ctx);
        fmpz_vec_init(dexp, 0);
        fmpz_vec_init(fexp, 0);
        gr_poly_init(def, ctx);
        gr_poly_init(pol, ctx);
        cc = gr_heap_init(ctx);

        status |= gr_poly_deflate(def, F, deflation, ctx);
        status |= _gr_poly_factor_no_deflation(c, dfac, dexp, def, algorithm, ctx);

        for (i = 0; i < dfac->length && status == GR_SUCCESS; i++)
        {
            /* the factor of F(x) = def(x^m) corresponding to the factor
               g(x) of def(x) is g(x^m), which may factor further; note
               that g(x^m) has deflation m, so we must not recurse into
               the deflation step */
            status |= gr_poly_inflate(pol, dfac->entries + i, deflation, ctx);
            gr_poly_vec_set_length(ffac, 0, ctx);
            fmpz_vec_set_length(fexp, 0);
            status |= _gr_poly_factor_no_deflation(cc, ffac, fexp, pol, algorithm, ctx);

            for (j = 0; j < ffac->length; j++)
            {
                status |= gr_poly_vec_append(fac, ffac->entries + j, ctx);
                fmpz_vec_append(exp, fexp->entries + j);
                fmpz_mul(exp->entries + exp->length - 1, exp->entries + exp->length - 1, dexp->entries + i);
            }
        }

        gr_poly_vec_clear(dfac, ctx);
        gr_poly_vec_clear(ffac, ctx);
        fmpz_vec_clear(dexp);
        fmpz_vec_clear(fexp);
        gr_poly_clear(def, ctx);
        gr_poly_clear(pol, ctx);
        gr_heap_clear(cc, ctx);
    }

    return status;
}

int
gr_poly_factor_finite_field(gr_ptr c, gr_poly_vec_t fac, fmpz_vec_t exp,
    const gr_poly_t F, gr_ctx_t ctx)
{
    return _gr_poly_factor_finite_field(c, fac, exp, F, GR_POLY_FACTOR_ALGORITHM_DEFAULT, ctx);
}

/* Method with the signature of GR_METHOD_POLY_FACTOR. */
int
_gr_poly_factor_finite_field_method(gr_ptr c, gr_vec_t fac, fmpz_vec_t exp,
    const gr_poly_t F, int FLINT_UNUSED(flags), gr_ctx_t ctx)
{
    gr_ptr cc;
    int status = GR_SUCCESS;

    cc = gr_heap_init(ctx);
    /* gr_vec_t over the polynomial ring is layout-compatible with gr_poly_vec_t */
    status |= gr_poly_factor_finite_field(cc, (gr_poly_vec_struct *) fac, exp, F, ctx);
    /* c is an element of the polynomial ring */
    status |= gr_poly_set_scalar((gr_poly_struct *) c, cc, ctx);
    gr_heap_clear(cc, ctx);

    return status;
}
