/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2013 Martin Lee

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_vec.h"
#include "nmod_poly.h"

void _nmod_poly_mulmod_preinv(nn_ptr res, nn_srcptr poly1, slong len1,
                            nn_srcptr poly2, slong len2, nn_srcptr f,
                            slong lenf, nn_srcptr finv, slong lenfinv, nmod_t mod)
{
    nn_ptr T, Q;
    slong lenT, lenQ;

    lenT = len1 + len2 - 1;
    lenQ = lenT - lenf + 1;

    T = _nmod_vec_init(lenT + lenQ);
    Q = T + lenT;

    if (len1 >= len2)
        _nmod_poly_mul(T, poly1, len1, poly2, len2, mod);
    else
        _nmod_poly_mul(T, poly2, len2, poly1, len1, mod);

    _nmod_poly_divrem_newton_n_preinv(Q, res, T, lenT, f, lenf,
                                      finv, lenfinv, mod);
    _nmod_vec_clear(T);
}

void
nmod_poly_mulmod_preinv(nmod_poly_t res, const nmod_poly_t poly1,
                        const nmod_poly_t poly2, const nmod_poly_t f,
                        const nmod_poly_t finv)
{
    slong len1, len2, lenf;
    nn_ptr fcoeffs;

    lenf = f->length;
    len1 = poly1->length;
    len2 = poly2->length;

    if (lenf == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (nmod_poly_mulmod_preinv). Divide by zero.\n");
    }

    if (lenf <= len1 || lenf <= len2)
    {
        flint_throw(FLINT_ERROR, "Exception (nmod_poly_mulmod_preinv). Input larger than modulus.\n");
    }

    if (lenf == 1 || len1 == 0 || len2 == 0)
    {
        nmod_poly_zero(res);
        return;
    }

    if (len1 + len2 - lenf > 0)
    {
        if (f == res)
        {
            fcoeffs = (nn_ptr) flint_malloc(sizeof(ulong) * lenf);
            _nmod_vec_set(fcoeffs, f->coeffs, lenf);
        } else
            fcoeffs = f->coeffs;

        nmod_poly_fit_length(res, lenf - 1);
        _nmod_poly_mulmod_preinv(res->coeffs, poly1->coeffs, len1,
                                     poly2->coeffs, len2, fcoeffs, lenf,
                                         finv->coeffs, finv->length, res->mod);
        if (f == res)
            flint_free(fcoeffs);

        res->length = lenf - 1;
        _nmod_poly_normalise(res);
    } else
    {
        nmod_poly_mul(res, poly1, poly2);
    }
}
