/*
    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010, 2022 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "longlong.h"
#include "nmod_poly.h"
#include "nmod_poly_factor.h"

#include "nmod_poly_factor_gr.h"

ulong
nmod_poly_factor(nmod_poly_factor_t res, const nmod_poly_t input)
{
    gr_ctx_t ctx;
    gr_poly_t P;
    gr_poly_vec_t fac;
    fmpz_vec_t exp;
    nmod_t mod = input->mod;
    ulong lc;

    if (input->length <= 1)
    {
        res->num = 0;
        return (input->length == 0) ? 0 : input->coeffs[0];
    }

    _gr_ctx_init_nmod(ctx, &mod);
    GR_MUST_SUCCEED(gr_ctx_set_is_field(ctx, T_TRUE));

    NMOD_POLY_AS_GR(P, input);
    gr_poly_vec_init(fac, 0, ctx);
    fmpz_vec_init(exp, 0);

    GR_MUST_SUCCEED(gr_poly_factor_finite_field(&lc, fac, exp, P, ctx));

    res->num = 0;
    _nmod_poly_factor_set_gr(res, fac, exp, mod, ctx);

    gr_poly_vec_clear(fac, ctx);
    fmpz_vec_clear(exp);
    gr_ctx_clear(ctx);

    return lc;
}

ulong
nmod_poly_factor_with_berlekamp(nmod_poly_factor_t res, const nmod_poly_t input)
{
    ulong lc = nmod_poly_lead(input) == NULL ? 0 : *nmod_poly_lead(input);

    if (input->length <= 1)
    {
        res->num = 0;
        return (input->length == 0) ? 0 : input->coeffs[0];
    }

    res->num = 0;
    _nmod_poly_factor_gr(res, input, GR_POLY_FACTOR_ALGORITHM_BERLEKAMP);
    return lc;
}

ulong
nmod_poly_factor_with_cantor_zassenhaus(nmod_poly_factor_t res, const nmod_poly_t input)
{
    ulong lc = nmod_poly_lead(input) == NULL ? 0 : *nmod_poly_lead(input);

    if (input->length <= 1)
    {
        res->num = 0;
        return (input->length == 0) ? 0 : input->coeffs[0];
    }

    res->num = 0;
    _nmod_poly_factor_gr(res, input, GR_POLY_FACTOR_ALGORITHM_CANTOR_ZASSENHAUS);
    return lc;
}

ulong
nmod_poly_factor_with_kaltofen_shoup(nmod_poly_factor_t res, const nmod_poly_t input)
{
    ulong lc = nmod_poly_lead(input) == NULL ? 0 : *nmod_poly_lead(input);

    if (input->length <= 1)
    {
        res->num = 0;
        return (input->length == 0) ? 0 : input->coeffs[0];
    }

    res->num = 0;
    _nmod_poly_factor_gr(res, input, GR_POLY_FACTOR_ALGORITHM_KALTOFEN_SHOUP);
    return lc;
}
