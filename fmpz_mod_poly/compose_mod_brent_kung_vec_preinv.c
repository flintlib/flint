/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2014 Martin Lee
    Copyright (C) 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz_vec.h"
#include "fmpz_mod_poly.h"
#include "fmpz_mat.h"
#include "ulong_extras.h"

void
_fmpz_mod_poly_compose_mod_brent_kung_vec_preinv(fmpz_mod_poly_struct * res,
                                                const fmpz_mod_poly_struct *
                                                polys, slong lenpolys, slong l,
                                                const fmpz * g, slong glen,
						const fmpz * poly, slong len,
                                                const fmpz * polyinv,
                                                slong leninv, const fmpz_t p)
{
    fmpz_mat_t A, B, C;
    fmpz *t, *h;
    slong i, j, k, n, m, len2 = l, len1;

    n = len - 1;

    m = n_sqrt(n * len2) + 1;

    h = _fmpz_vec_init(n);
    t = _fmpz_vec_init(n);

    k = len / m + 1;

    fmpz_mat_init(A, m, n);
    fmpz_mat_init(B, k * len2, m);
    fmpz_mat_init(C, k * len2, n);

    /* Set rows of B to the segments of polys */
    for (j = 0; j < len2; j++)
    {
        len1 = (polys + j)->length;
        for (i = 0; i < len1 / m; i++)
            _fmpz_vec_set(B->rows[i + j * k], (polys + j)->coeffs + i * m, m);
        _fmpz_vec_set(B->rows[i + j * k], (polys + j)->coeffs + i * m,
                      len1 % m);
    }

    /* Set rows of A to powers of last element of polys */
    _fmpz_mod_poly_powers_mod_preinv_naive(A->rows, g, glen,
                                            m, poly, len, polyinv, leninv, p);

    fmpz_mat_mul(C, B, A);
    for (i = 0; i < k * len2; i++)
        for (j = 0; j < n; j++)
            fmpz_mod(C->rows[i] + j, C->rows[i] + j, p);

    /* Evaluate block composition using the Horner scheme */
    if (n == 1)
    {
        fmpz_mul(h + 0, A->rows[m - 1] + 0, A->rows[1] + 0);
	fmpz_mod(h + 0, h + 0, p);
    } else
    {                                                                                            _fmpz_mod_poly_mulmod_preinv(h, A->rows[m - 1], n, A->rows[1], n, poly,
                                 len, polyinv, leninv, p);
    }

    for (j = 0; j < len2; j++)
    {
        _fmpz_vec_set((res + j)->coeffs, C->rows[(j + 1) * k - 1], n);

	if (n == 1)
        {
            for (i = 2; i <= k; i++)
            {
               fmpz_mul(t + 0, res[j].coeffs + 0, h + 0);
	       fmpz_add(res[j].coeffs + 0, t + 0, C->rows[(j + 1)*k - i] + 0);
	       fmpz_mod(res[j].coeffs + 0, res[j].coeffs + 0, p);
            }
        } else
        {
            for (i = 2; i <= k; i++)
            {
                _fmpz_mod_poly_mulmod_preinv(t, res[j].coeffs, n, h, n, poly,
                                         len, polyinv, leninv, p);
                _fmpz_mod_poly_add(res[j].coeffs, t, n,
                               C->rows[(j + 1)*k - i], n, p);
            }
	}
    }

    _fmpz_vec_clear(h, n);
    _fmpz_vec_clear(t, n);

    fmpz_mat_clear(A);
    fmpz_mat_clear(B);
    fmpz_mat_clear(C);
}

void fmpz_mod_poly_compose_mod_brent_kung_vec_preinv(
    fmpz_mod_poly_struct * res,
    const fmpz_mod_poly_struct * polys,
    slong len1,
    slong n,
    const fmpz_mod_poly_t g,
    const fmpz_mod_poly_t poly,
    const fmpz_mod_poly_t polyinv,
    const fmpz_mod_ctx_t ctx)
{
    slong len2 = poly->length;
    slong len3, i;

    for (i = 0; i < len1; i++)
    {
        len3 = (polys + i)->length;
        if (len3 >= len2)
        {
            flint_printf
                ("Exception (fmpz_mod_poly_compose_mod_brent_kung_vec_preinv)."
                 "The degree of the first polynomial must be smaller than that of the "
                 " modulus\n");
            flint_abort();
        }
    }

    if (n > len1)
    {
        flint_printf
            ("Exception (fmpz_mod_poly_compose_mod_brent_kung_vec_preinv)."
             "n is larger than the length of polys\n");
        flint_abort();
    }

    if (n == 0)
        return;

    if (len2 == 1)
    {
        for (i = 0; i < n; i++)
            fmpz_mod_poly_zero(res + i, ctx);

        return;
    }

    if (len2 == 2)
    {
        for (i = 0; i < n; i++)
            fmpz_mod_poly_set(res + i, polys + i, ctx);

        return;
    }

    for (i = 0; i < n; i++)
    {
        fmpz_mod_poly_fit_length(res + i, len2 - 1, ctx);
	    _fmpz_mod_poly_set_length(res + i, len2 - 1);
    }

    _fmpz_mod_poly_compose_mod_brent_kung_vec_preinv(res, polys, len1, n,
                                                     g->coeffs, g->length,
                                                     poly->coeffs, len2,
                                                     polyinv->coeffs,
                                                     polyinv->length,
                                                     fmpz_mod_ctx_modulus(ctx));

    for (i = 0; i < n; i++)
        _fmpz_mod_poly_normalise(res + i);
}
