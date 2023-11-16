/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "gr_poly.h"

int
gr_poly_pow_fmpz(gr_poly_t res,
    const gr_poly_t poly, const fmpz_t exp, gr_ctx_t ctx)
{
    if (fmpz_is_zero(exp))
    {
        return gr_poly_one(res, ctx);
    }
    else if (poly->length == 0)
    {
        if (fmpz_sgn(exp) < 0)
            return GR_DOMAIN;
        else
            return gr_poly_zero(res, ctx);
    }
    else if (poly->length == 1)
    {
        int status;
        gr_poly_fit_length(res, 1, ctx);
        status = gr_pow_fmpz(res->coeffs, poly->coeffs, exp, ctx);
        _gr_poly_set_length(res, 1, ctx);
        _gr_poly_normalise(res, ctx);
        return status;
    }
    else
    {
        if (fmpz_sgn(exp) < 0)
            return GR_DOMAIN;

        if (COEFF_IS_MPZ(*exp))
            return GR_UNABLE;

        return gr_poly_pow_ui(res, poly, *exp, ctx);
    }
}
