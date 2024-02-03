/*
    Copyright (C) 2023 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_vec.h"
#include "nmod_poly.h"
#include "nmod_poly_mat.h"

void nmod_poly_mat_set_trunc(nmod_poly_mat_t res, const nmod_poly_mat_t pmat, long len)
{
    nmod_poly_struct * res_ij;
    nmod_poly_struct * pmat_ij;

    if (pmat == res)
    {
        for (slong i = 0; i < pmat->r; i++)
        {
            for (slong j = 0; j < pmat->c; j++)
            {
                res_ij = nmod_poly_mat_entry(res, i, j);
                if (res_ij->length > len)
                {
                    res_ij->length = len;
                    _nmod_poly_normalise(res_ij);
                }
            }
        }
    }
    else
    {
        slong rlen;
        for (slong i = 0; i < pmat->r; i++)
        {
            for (slong j = 0; j < pmat->c; j++)
            {
                res_ij = nmod_poly_mat_entry(res, i, j);
                pmat_ij = nmod_poly_mat_entry(pmat, i, j);

                rlen = FLINT_MIN(len, pmat_ij->length);
                while (rlen > 0 && pmat_ij->coeffs[rlen - 1] == 0)
                    rlen--;

                nmod_poly_fit_length(res_ij, rlen);
                _nmod_vec_set(res_ij->coeffs, pmat_ij->coeffs, rlen);
                _nmod_poly_set_length(res_ij, rlen);
            }
        }
    }
}
