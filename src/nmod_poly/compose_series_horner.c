/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

void
_nmod_poly_compose_series_horner(mp_ptr res, mp_srcptr poly1, slong len1, 
                            mp_srcptr poly2, slong len2, slong n, nmod_t mod)
{
    if (n == 1)
    {
        res[0] = poly1[0];
    }
    else
    {
        slong i = len1 - 1;
        slong lenr;

        mp_ptr t = _nmod_vec_init(n);

        lenr = len2;
        _nmod_vec_scalar_mul_nmod(res, poly2, len2, poly1[i], mod);
        i--;
        res[0] = nmod_add(res[0], poly1[i], mod);

        while (i > 0)
        {
            i--;
            if (lenr + len2 - 1 < n)
            {
                _nmod_poly_mul(t, res, lenr, poly2, len2, mod);
                lenr = lenr + len2 - 1;
            }
            else
            {
                _nmod_poly_mullow(t, res, lenr, poly2, len2, n, mod);
                lenr = n;
            }
            _nmod_poly_add(res, t, lenr, poly1 + i, 1, mod);
        }

        _nmod_vec_zero(res + lenr, n - lenr);
        _nmod_vec_clear(t);
    }
}

void
nmod_poly_compose_series_horner(nmod_poly_t res, 
                    const nmod_poly_t poly1, const nmod_poly_t poly2, slong n)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;
    slong lenr;

    if (len2 != 0 && poly2->coeffs[0] != 0)
    {
        flint_printf("Exception (nmod_poly_compose_series_horner). Inner polynomial "
               "must have zero constant term.\n");
        flint_abort();
    }

    if (len1 == 0 || n == 0)
    {
        nmod_poly_zero(res);
        return;
    }

    if (len2 == 0 || len1 == 1)
    {
        nmod_poly_fit_length(res, 1);
        res->coeffs[0] = poly1->coeffs[0];
        res->length = 1;
        _nmod_poly_normalise(res);
        return;
    }

    lenr = FLINT_MIN((len1 - 1) * (len2 - 1) + 1, n);
    len1 = FLINT_MIN(len1, lenr);
    len2 = FLINT_MIN(len2, lenr);

    if ((res != poly1) && (res != poly2))
    {
        nmod_poly_fit_length(res, lenr);
        _nmod_poly_compose_series_horner(res->coeffs, poly1->coeffs, len1, 
                                        poly2->coeffs, len2, lenr, res->mod);
        res->length = lenr;
        _nmod_poly_normalise(res);
    }
    else
    {
        nmod_poly_t t;
        nmod_poly_init2_preinv(t, res->mod.n, res->mod.ninv, lenr);
        _nmod_poly_compose_series_horner(t->coeffs, poly1->coeffs, len1,
                                        poly2->coeffs, len2, lenr, res->mod);
        t->length = lenr;
        _nmod_poly_normalise(t);
        nmod_poly_swap(res, t);
        nmod_poly_clear(t);
    }
}
