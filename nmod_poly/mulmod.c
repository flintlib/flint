/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

void _nmod_poly_mulmod(mp_ptr res, mp_srcptr poly1, len_t len1, 
                             mp_srcptr poly2, len_t len2, mp_srcptr f,
                            len_t lenf, nmod_t mod)
{
    mp_ptr T, Q;
    len_t lenT, lenQ;

    lenT = len1 + len2 - 1;
    lenQ = lenT - lenf + 1;

    T = _nmod_vec_init(lenT + lenQ);
    Q = T + lenT;

    if (len1 >= len2)
        _nmod_poly_mul(T, poly1, len1, poly2, len2, mod);
    else
        _nmod_poly_mul(T, poly2, len2, poly1, len1, mod);

    _nmod_poly_divrem(Q, res, T, lenT, f, lenf, mod);
    _nmod_vec_clear(T);
}

void
nmod_poly_mulmod(nmod_poly_t res,
    const nmod_poly_t poly1, const nmod_poly_t poly2, const nmod_poly_t f)
{
    len_t len1, len2, lenf;
    mp_ptr fcoeffs;

    lenf = f->length;
    len1 = poly1->length;
    len2 = poly2->length;

    if (lenf == 0)
    {
        printf("Exception (nmod_poly_mulmod). Divide by zero.\n");
        abort();
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
            fcoeffs = flint_malloc(sizeof(mp_limb_t) * lenf);
            _nmod_vec_set(fcoeffs, f->coeffs, lenf);
        }
        else
            fcoeffs = f->coeffs;

        nmod_poly_fit_length(res, lenf - 1);
        _nmod_poly_mulmod(res->coeffs, poly1->coeffs, len1,
                                       poly2->coeffs, len2,
                                       fcoeffs, lenf,
                                       res->mod);
        if (f == res)
            flint_free(fcoeffs);

        res->length = lenf - 1;
        _nmod_poly_normalise(res);
    }
    else
    {
        nmod_poly_mul(res, poly1, poly2);
    }
}
