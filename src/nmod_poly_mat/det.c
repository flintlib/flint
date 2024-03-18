/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "perm.h"
#include "nmod_vec.h"
#include "nmod_mat.h"
#include "nmod_poly.h"
#include "nmod_poly_mat.h"

void
nmod_poly_mat_det_fflu(nmod_poly_t det, const nmod_poly_mat_t A)
{
    slong n = nmod_poly_mat_nrows(A);

    if (n == 0)
        nmod_poly_one(det);
    else
    {
        nmod_poly_mat_t tmp;
        slong * perm;
        nmod_poly_mat_init_set(tmp, A);
        perm = _perm_init(n);

        nmod_poly_mat_fflu(tmp, det, perm, tmp, 1);
        if (_perm_parity(perm, n))
            nmod_poly_neg(det, det);

        _perm_clear(perm);
        nmod_poly_mat_clear(tmp);
    }
}

void
nmod_poly_mat_det_interpolate(nmod_poly_t det, const nmod_poly_mat_t A)
{
    slong i, l, n, len;

    nmod_mat_t X;
    mp_ptr x, d;

    n = A->r;

    if (n == 0)
    {
        nmod_poly_one(det);
        return;
    }

    l = nmod_poly_mat_max_length(A);

    if (l == 0)
    {
        nmod_poly_zero(det);
        return;
    }

    /* Bound degree based on Laplace expansion */
    len = n*(l - 1) + 1;

    /* Not enough points to interpolate */
    if (len > nmod_poly_mat_modulus(A))
    {
        nmod_poly_mat_det_fflu(det, A);
        return;
    }

    x = _nmod_vec_init(len);
    d = _nmod_vec_init(len);
    nmod_mat_init(X, n, n, nmod_poly_mat_modulus(A));

    for (i = 0; i < len; i++)
    {
        x[i] = i;
        nmod_poly_mat_evaluate_nmod(X, A, x[i]);
        d[i] = nmod_mat_det(X);
    }

    nmod_poly_interpolate_nmod_vec(det, x, d, len);

    _nmod_vec_clear(x);
    _nmod_vec_clear(d);
    nmod_mat_clear(X);
}

void
nmod_poly_mat_det(nmod_poly_t det, const nmod_poly_mat_t A)
{
    slong n = A->r;

    if (n == 0)
    {
        nmod_poly_one(det);
    }
    else if (n == 1)
    {
        nmod_poly_set(det, nmod_poly_mat_entry(A, 0, 0));
    }
    else if (n == 2)
    {
        nmod_poly_t tmp;
        nmod_poly_init(tmp, nmod_poly_mat_modulus(A));
        nmod_poly_mul(det, nmod_poly_mat_entry(A, 0, 0),
                           nmod_poly_mat_entry(A, 1, 1));
        nmod_poly_mul(tmp, nmod_poly_mat_entry(A, 0, 1),
                           nmod_poly_mat_entry(A, 1, 0));
        nmod_poly_sub(det, det, tmp);
        nmod_poly_clear(tmp);
    }
    else if (n < 15)  /* should be entry sensitive too */
    {
        nmod_poly_mat_det_fflu(det, A);
    }
    else
    {
        nmod_poly_mat_det_interpolate(det, A);
    }
}
