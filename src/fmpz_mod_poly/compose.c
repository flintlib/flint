/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"
#include "gr_poly.h"

void _fmpz_mod_poly_compose(fmpz *res, const fmpz *poly1, slong len1,
                                              const fmpz *poly2, slong len2,
                                              const fmpz_mod_ctx_t ctx)
{
    gr_ctx_t gr_ctx;
    _gr_ctx_init_fmpz_mod_from_ref(gr_ctx, ctx);
    GR_MUST_SUCCEED(_gr_poly_compose(res, poly1, len1, poly2, len2, gr_ctx));
}

void fmpz_mod_poly_compose(fmpz_mod_poly_t res,
                    const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2,
                                                      const fmpz_mod_ctx_t ctx)
{
    const slong len1 = poly1->length;
    const slong len2 = poly2->length;

    if (len1 == 0)
    {
        fmpz_mod_poly_zero(res, ctx);
    }
    else if (len1 == 1 || len2 == 0)
    {
        fmpz_mod_poly_set_fmpz(res, poly1->coeffs, ctx);
    }
    else
    {
        const slong lenr = (len1 - 1) * (len2 - 1) + 1;

        if ((res != poly1) && (res != poly2))
        {
            fmpz_mod_poly_fit_length(res, lenr, ctx);
            _fmpz_mod_poly_compose(res->coeffs, poly1->coeffs, len1,
                               poly2->coeffs, len2, ctx);
        }
        else
        {
            fmpz *t = _fmpz_vec_init(lenr);

            _fmpz_mod_poly_compose(t, poly1->coeffs, len1,
                               poly2->coeffs, len2, ctx);
            _fmpz_vec_clear(res->coeffs, res->alloc);
            res->coeffs = t;
            res->alloc  = lenr;
            res->length = lenr;
        }

        _fmpz_mod_poly_set_length(res, lenr);
        _fmpz_mod_poly_normalise(res);
    }
}
