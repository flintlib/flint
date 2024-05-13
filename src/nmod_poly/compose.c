/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2010, 2012 Sebastian Pancratz
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"
#include "gr_poly.h"

void
_nmod_poly_compose(nn_ptr res, nn_srcptr poly1, slong len1,
                               nn_srcptr poly2, slong len2, nmod_t mod)
{
    if (len1 == 1)
        res[0] = poly1[0];
    else if (len2 == 1)
        res[0] = _nmod_poly_evaluate_nmod(poly1, len1, poly2[0], mod);
    else if (len1 <= 7)
        _nmod_poly_compose_horner(res, poly1, len1, poly2, len2, mod);
    else
    {
        /* delegates to Taylor shift, divconquer */
        /* todo: also incorporate the nmod_poly Taylor shift in gr_poly */
        gr_ctx_t ctx;
        _gr_ctx_init_nmod(ctx, &mod);
        GR_MUST_SUCCEED(_gr_poly_compose(res, poly1, len1, poly2, len2, ctx));
    }
}

void nmod_poly_compose(nmod_poly_t res,
                       const nmod_poly_t poly1, const nmod_poly_t poly2)
{
    const slong len1 = poly1->length;
    const slong len2 = poly2->length;

    if (len1 == 0)
    {
        nmod_poly_zero(res);
    }
    else if (len1 == 1 || len2 == 0)
    {
        nmod_poly_fit_length(res, 1);
        res->coeffs[0] = poly1->coeffs[0];
        res->length = (res->coeffs[0] != 0);
    }
    else
    {
        const slong lenr = (len1 - 1) * (len2 - 1) + 1;

        if (res != poly1 && res != poly2)
        {
            nmod_poly_fit_length(res, lenr);
            _nmod_poly_compose(res->coeffs, poly1->coeffs, len1,
                                                   poly2->coeffs, len2, poly1->mod);
        }
        else
        {
            nmod_poly_t t;
            nmod_poly_init2(t, poly1->mod.n, lenr);
            _nmod_poly_compose(t->coeffs, poly1->coeffs, len1,
                                                 poly2->coeffs, len2, poly1->mod);
            nmod_poly_swap(res, t);
            nmod_poly_clear(t);
        }

        res->length = lenr;
        _nmod_poly_normalise(res);
    }
}
