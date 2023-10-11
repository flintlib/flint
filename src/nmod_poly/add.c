/*
    Copyright (C) 2007, David Howden.
    Copyright (C) 2008, 2009, 2010, 2014, 2015 William Hart
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2015 Tommy Hofmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

void
_nmod_poly_add(mp_ptr res, mp_srcptr poly1, slong len1, mp_srcptr poly2,
               slong len2, nmod_t mod)
{
    slong i, min = FLINT_MIN(len1, len2);

    _nmod_vec_add(res, poly1, poly2, min, mod);

    if (poly1 != res)           /* copy any remaining coefficients from poly1 */
        for (i = min; i < len1; i++)
            res[i] = poly1[i];

    if (poly2 != res)           /* copy any remaining coefficients from poly2 */
        for (i = min; i < len2; i++)
            res[i] = poly2[i];
}

void
nmod_poly_add(nmod_poly_t res, const nmod_poly_t poly1,
              const nmod_poly_t poly2)
{
    slong max = FLINT_MAX(poly1->length, poly2->length);

    nmod_poly_fit_length(res, max);

    _nmod_poly_add(res->coeffs, poly1->coeffs, poly1->length, poly2->coeffs,
                   poly2->length, poly1->mod);

    res->length = max;
    _nmod_poly_normalise(res);  /* there may have been cancellation */
}

void nmod_poly_add_series(nmod_poly_t res,
            const nmod_poly_t poly1, const nmod_poly_t poly2, slong n)
{
    slong len1, len2, max = FLINT_MAX(poly1->length, poly2->length);

    if (n < 0)
       n = 0;

    max = FLINT_MIN(max, n);
    len1 = FLINT_MIN(poly1->length, max);
    len2 = FLINT_MIN(poly2->length, max);

    nmod_poly_fit_length(res, max);

    _nmod_poly_add(res->coeffs, poly1->coeffs, len1,
                                    poly2->coeffs, len2, poly1->mod);

    _nmod_poly_set_length(res, max);
    _nmod_poly_normalise(res);
}

void nmod_poly_add_ui(nmod_poly_t res, const nmod_poly_t poly, ulong c)
{

    if (poly->length == 0)
    {
        if (c == 0)
            nmod_poly_zero(res);
        else
        {
            nmod_poly_fit_length(res, 1);
            nmod_poly_set_coeff_ui(res, 0, c);
            _nmod_poly_set_length(res, 1);
        }
    }
    else
    {
        if (c >= poly->mod.n)
            NMOD_RED(c, c, poly->mod);

        nmod_poly_set(res, poly);

        nmod_poly_set_coeff_ui(res, 0, nmod_add(res->coeffs[0], c, poly->mod));

        _nmod_poly_normalise(res);
   }
}

