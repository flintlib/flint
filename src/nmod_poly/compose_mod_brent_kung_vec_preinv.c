/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2014 Martin Lee
    Copyright (C) 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "nmod_mat.h"

void
_nmod_poly_compose_mod_brent_kung_vec_preinv(nmod_poly_struct * res,
                const nmod_poly_struct * polys, slong lenpolys, slong l,
                         mp_srcptr g, slong glen, mp_srcptr poly, slong len,
                                   mp_srcptr polyinv, slong leninv, nmod_t mod)
{
    nmod_mat_t A, B, C;
    mp_ptr t, h;
    slong i, j, k, n, m, len2 = l, len1;

    n = len - 1;

    m = n_sqrt(n*len2) + 1;

    h = _nmod_vec_init(n);
    t = _nmod_vec_init(n);

    k = len/m + 1;

    nmod_mat_init(A, m, n, mod.n);
    nmod_mat_init(B, k*len2, m, mod.n);
    nmod_mat_init(C, k*len2, n, mod.n);

    /* Set rows of B to the segments of polys */
    for (j = 0; j < len2; j++)
    {
        len1 = (polys + j)->length;

        for (i = 0; i < len1/m; i++)
            _nmod_vec_set(B->rows[i + j*k], (polys + j)->coeffs + i*m, m);

        _nmod_vec_set(B->rows[i + j*k], (polys + j)->coeffs + i*m, len1%m);
    }

    /* Set rows of A to powers of last element of polys */
    _nmod_poly_powers_mod_preinv_naive(A->rows, g, glen,
                                           m, poly, len, polyinv, leninv, mod);

    nmod_mat_mul(C, B, A);

    /* Evaluate block composition using the Horner scheme */
    if (n == 1)
    {
        h[0] = n_mulmod2_preinv(A->rows[m - 1][0],
                                               A->rows[1][0], mod.n, mod.ninv);
    } else
    {
        _nmod_poly_mulmod_preinv(h, A->rows[m - 1], n, A->rows[1], n, poly,
                                                    len, polyinv, leninv, mod);
    }

    for (j = 0; j < len2; j++)
    {
        _nmod_vec_set((res + j)->coeffs, C->rows[(j + 1)*k - 1], n);

        if (n == 1)
        {
            for (i = 2; i <= k; i++)
            {
                t[0] = n_mulmod2_preinv(res[j].coeffs[0],
                                                        h[0], mod.n, mod.ninv);
                res[j].coeffs[0] = n_addmod(t[0],
                                             C->rows[(j + 1)*k - i][0], mod.n);
            }
        } else
        {
            for (i = 2; i <= k; i++)
            {
                _nmod_poly_mulmod_preinv(t, res[j].coeffs,
                                     n, h, n, poly, len, polyinv, leninv, mod);
                _nmod_poly_add(res[j].coeffs, t, n,
                                               C->rows[(j + 1)*k - i], n, mod);
            }
        }
    }

    _nmod_vec_clear(h);
    _nmod_vec_clear(t);

    nmod_mat_clear(A);
    nmod_mat_clear(B);
    nmod_mat_clear(C);
}

void
nmod_poly_compose_mod_brent_kung_vec_preinv(nmod_poly_struct * res,
                     const nmod_poly_struct * polys, slong len1, slong n,
        const nmod_poly_t g, const nmod_poly_t poly, const nmod_poly_t polyinv)
{
    slong len2 = poly->length;
    slong len3, i;

    for (i = 0; i < len1; i++)
    {
        len3 = (polys + i)->length;

        if (len3 >= len2)
        {
            flint_throw(FLINT_ERROR, "(nmod_poly_compose_mod_brent_kung_vec_preinv): "
                 "The degree of the first polynomial must be smaller than that of the modulus\n");
        }
    }

    if (n > len1)
    {
        flint_throw(FLINT_ERROR, "(nmod_poly_compose_mod_brent_kung_vec_preinv): "
                "n is larger than the length of polys\n");
    }

    if (n == 0)
        return;

    if (len2 == 1)
    {
        for (i = 0; i < n; i++)
            nmod_poly_zero(res + i);

        return;
    }

    if (len2 == 2)
    {
        for (i = 0; i < n; i++)
            nmod_poly_set(res + i, polys + i);

        return;
    }

    for (i = 0; i < n; i++)
    {
        nmod_poly_fit_length(res + i, len2 - 1);
        _nmod_poly_set_length(res + i, len2 - 1);
    }

    _nmod_poly_compose_mod_brent_kung_vec_preinv(res, polys, len1, n,
                g->coeffs, g->length, poly->coeffs, len2, polyinv->coeffs,
                                                   polyinv->length, poly->mod);

    for (i = 0; i < n; i++)
        _nmod_poly_normalise(res + i);
}
