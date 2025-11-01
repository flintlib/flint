/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* #include "flint.h" */
#include "nmod_poly.h"
#include "nmod_vec.h"
#include "stdio.h"

void _nmod_poly_divrem_xnmc(nn_ptr RQ, nn_srcptr A, slong len, ulong n, ulong c, nmod_t mod)
{
    /* assumes len >= n */
    slong i;
    ulong j, r, val;

    if (RQ != A)
        for (j = 0; j < n; j++)
            RQ[len-n+j] = A[len-n+j];

    r = len % n;
    i = len - r - n;  /* multiple of n, >= 0 by assumption */

    for (j = 0; j < r; j++)
    {
        val = nmod_mul(RQ[i+n+j], c, mod);
        RQ[i+j] = n_addmod(val, A[i+j], mod.n);
    }

    i -= n;
    while (i >= 0)
    {
        for (j = 0; j < n; j++)
        {
            val = nmod_mul(RQ[i+n+j], c, mod);
            RQ[i+j] = n_addmod(val, A[i+j], mod.n);
        }
        i -= n;
    }
}

void nmod_poly_divrem_xnmc(nmod_poly_t Q, nmod_poly_t R, nmod_poly_t A, ulong n, ulong c)
{
    const ulong len = A->length;

    if (len <= n)
    {
        nmod_poly_zero(Q);
        nmod_poly_set(R, A);
        return;
    }

    nn_ptr RQ = _nmod_vec_init(len);

    /* perform division */
    _nmod_poly_divrem_xnmc(RQ, A->coeffs, len, n, c, A->mod);

    /* copy remainder R */
    nmod_poly_fit_length(R, n);
    _nmod_vec_set(R->coeffs, RQ, n);
    _nmod_poly_set_length(R, n);
    _nmod_poly_normalise(R);

    /* copy quotient Q */
    nmod_poly_fit_length(Q, len - n);
    _nmod_vec_set(Q->coeffs, RQ + n, len - n);
    _nmod_poly_set_length(Q, len - n);

    _nmod_vec_clear(RQ);
}
