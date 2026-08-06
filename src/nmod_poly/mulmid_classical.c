/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

/* Assumes poly1 and poly2 are not length 0 and 0 <= nhi < nlo <= len1 + len2 - 1 */
void
_nmod_poly_mulmid_classical(nn_ptr res, nn_srcptr poly1, slong len1,
                            nn_srcptr poly2, slong len2, slong nlo, slong nhi, nmod_t mod)
{
    slong i, n1, n2;
    int squaring;
    ulong c;

    len1 = FLINT_MIN(len1, nhi);
    len2 = FLINT_MIN(len2, nhi);

    if (len1 < len2)
    {
        FLINT_SWAP(nn_srcptr, poly1, poly2);
        FLINT_SWAP(slong, len1, len2);
    }

    if (nhi == 1)
    {
        res[0] = nmod_mul(poly1[0], poly2[0], mod);
        return;
    }

    if (len2 == 1)
    {
        _nmod_vec_scalar_mul_nmod(res, poly1 + nlo, nhi - nlo, poly2[0], mod);
        return;
    }

    squaring = (poly1 == poly2 && len1 == len2);

    dot_params_t params = _nmod_vec_dot_params(len2, mod);

    if (squaring)
    {
        for (i = nlo; i < nhi; i++)
        {
            n1 = FLINT_MAX(0, i - len1 + 1);
            n2 = FLINT_MIN(len1 - 1, (i + 1) / 2 - 1);

            c = _nmod_vec_dot_rev(poly1 + n1, poly1 + i - n2, n2 - n1 + 1, mod, params);
            c = nmod_add(c, c, mod);

            if (i % 2 == 0 && i / 2 < len1)
                NMOD_ADDMUL(c, poly1[i / 2], poly1[i / 2], mod);

            res[i - nlo] = c;
        }
    }
    else
    {
        for (i = nlo; i < nhi; i++)
        {
            n1 = FLINT_MIN(len1 - 1, i);
            n2 = FLINT_MIN(len2 - 1, i);

            res[i - nlo] = _nmod_vec_dot_rev(poly1 + i - n2,
                                       poly2 + i - n1,
                                       n1 + n2 - i + 1, mod, params);
        }
    }
}

void
nmod_poly_mulmid_classical(nmod_poly_t res,
                           const nmod_poly_t poly1, const nmod_poly_t poly2,
                           slong nlo, slong nhi)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;
    slong len;

    if (len1 == 0 || len2 == 0 || nlo >= FLINT_MIN(nhi, len1 + len2 - 1))
    {
        nmod_poly_zero(res);
        return;
    }

    nhi = FLINT_MIN(nhi, len1 + len2 - 1);
    len = nhi - nlo;

    if (res == poly1 || res == poly2)
    {
        nmod_poly_t temp;
        nmod_poly_init2_preinv(temp, poly1->mod.n, poly1->mod.ninv, len);
        _nmod_poly_mulmid_classical(temp->coeffs, poly1->coeffs,
                                    len1, poly2->coeffs, len2, nlo, nhi, poly1->mod);
        nmod_poly_swap(res, temp);
        nmod_poly_clear(temp);
    }
    else
    {
        nmod_poly_fit_length(res, len);
        _nmod_poly_mulmid_classical(res->coeffs, poly1->coeffs,
                                    len1, poly2->coeffs, len2, nlo, nhi, poly1->mod);
    }

    res->length = len;
    _nmod_poly_normalise(res);
}
