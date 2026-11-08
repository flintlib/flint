/*
    Copyright (C) 2016 Vincent Delecroix

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

void
_nmod_poly_power_sums(nn_ptr res, nn_srcptr poly, slong len, slong n,
                      nmod_t mod)
{
    if (10 * n >= len + 75)
        _nmod_poly_power_sums_schoenhage(res, poly, len, n, mod);
    else
        _nmod_poly_power_sums_naive(res, poly, len, n, mod);
}

void
nmod_poly_power_sums(nmod_poly_t res, const nmod_poly_t poly, slong n)
{
    slong len = poly->length;
    slong i = 0;

    if (len == 0)
    {
        flint_throw(FLINT_ERROR, "(nmod_poly_power_sums_naive): Zero polynomial.\n");
    }
    while (poly->coeffs[i] == 0)
        i++;

    if (n <= 0 || len == 1)
        nmod_poly_zero(res);
    else if (len == i + 1)
    {
        nmod_poly_fit_length(res, 1);
        _nmod_poly_set_length(res, 1);
        NMOD_RED(res->coeffs[0], (ulong) (len - 1), poly->mod);
    }
    else
    {
        if (*nmod_poly_lead(poly) != 1)
        {
            nmod_poly_t t;
            nmod_poly_init_preinv(t, poly->mod.n, poly->mod.ninv);
            nmod_poly_make_monic(t, poly);
            nmod_poly_fit_length(res, n);
            _nmod_poly_power_sums(res->coeffs, t->coeffs + i,
                                  len - i, n, t->mod);
            nmod_poly_clear(t);
        }
        else if (poly == res)
        {
            nmod_poly_t t;
            nmod_poly_init_preinv(t, poly->mod.n, poly->mod.ninv);
            nmod_poly_fit_length(t, n);
            _nmod_poly_power_sums(t->coeffs, poly->coeffs + i,
                                  len - i, n, t->mod);
            nmod_poly_swap(t, res);
            nmod_poly_clear(t);
        }
        else
        {
            nmod_poly_fit_length(res, n);
            _nmod_poly_power_sums(res->coeffs, poly->coeffs + i,
                                  len - i, n, poly->mod);
        }
        if (i)
            NMOD_RED(res->coeffs[0], (ulong) (len - 1), poly->mod);
        _nmod_poly_set_length(res, n);
        _nmod_poly_normalise(res);
    }
}

/* Newton's identities: for the monic poly of degree d = len - 1 with
   coefficients a_j, the power sums are p_0 = d and, for k >= 1,
   p_k = -(k a_(d-k) + sum_(i=1)^(k-1) a_(d-k+i) p_i), with a_j = 0 for j < 0.
   For k <= 13 (at most 12 products a_(d-k+i) p_i) the sum is computed term
   by term: there the parameters and the short-length dispatch of
   _nmod_vec_dot cost more than they save. The longer sums are dot
   products. */
void
_nmod_poly_power_sums_naive(nn_ptr res, nn_srcptr poly, slong len, slong n,
                            nmod_t mod)
{
    slong i, k;
    slong kmax = FLINT_MIN(n, len);
    slong kshort = FLINT_MIN(kmax, 14);
    ulong kmod = 0;     /* k mod n, as nmod_mul needs reduced operands */

    NMOD_RED(res[0], (ulong) (len - 1), mod);
    for (k = 1; k < kshort; k++)
    {
        kmod = (kmod + 1 == mod.n) ? 0 : kmod + 1;
        res[k] = nmod_mul(poly[len - 1 - k], kmod, mod);
        for (i = 1; i < k; i++)
            res[k] =
                nmod_add(res[k], nmod_mul(poly[len - 1 - k + i], res[i], mod),
                         mod);
        res[k] = nmod_neg(res[k], mod);
    }

    if (kshort < kmax || n > len)
    {
        dot_params_t params = _nmod_vec_dot_params(len - 1, mod);

        for (k = kshort; k < kmax; k++)
        {
            kmod = (kmod + 1 == mod.n) ? 0 : kmod + 1;
            res[k] = _nmod_vec_dot(poly + len - k, res + 1, k - 1, mod, params);
            res[k] = nmod_neg(nmod_addmul(res[k], poly[len - 1 - k], kmod, mod), mod);
        }
        for (k = len; k < n; k++)
        {
            res[k] = _nmod_vec_dot(poly, res + k - len + 1, len - 1, mod, params);
            res[k] = nmod_neg(res[k], mod);
        }
    }
}

void
nmod_poly_power_sums_naive(nmod_poly_t res, const nmod_poly_t poly, slong n)
{
    if (poly->length == 0)
    {
        flint_throw(FLINT_ERROR, "(nmod_poly_power_sums_naive): Zero polynomial.\n");
    }
    else if (n <= 0 || poly->length == 1)
    {
        nmod_poly_zero(res);
    }
    else
    {
        if (*nmod_poly_lead(poly) != 1)
        {
            nmod_poly_t t;
            nmod_poly_init_preinv(t, poly->mod.n, poly->mod.ninv);
            nmod_poly_make_monic(t, poly);
            nmod_poly_fit_length(res, n);
            _nmod_poly_power_sums_naive(res->coeffs, t->coeffs,
                                        t->length, n, t->mod);
            nmod_poly_clear(t);
        }
        else if (poly == res)
        {
            nmod_poly_t t;
            nmod_poly_init_preinv(t, poly->mod.n, poly->mod.ninv);
            nmod_poly_fit_length(t, n);
            _nmod_poly_power_sums_naive(t->coeffs, poly->coeffs,
                                        poly->length, n, t->mod);
            nmod_poly_swap(t, res);
            nmod_poly_clear(t);
        }
        else
        {
            nmod_poly_fit_length(res, n);
            _nmod_poly_power_sums_naive(res->coeffs, poly->coeffs,
                                        poly->length, n, poly->mod);
        }
        _nmod_poly_set_length(res, n);
        _nmod_poly_normalise(res);
    }
}

void
_nmod_poly_power_sums_schoenhage(nn_ptr res, nn_srcptr poly, slong len,
                                 slong n, nmod_t mod)
{
    nn_ptr a, b;

    a = (nn_ptr) flint_malloc((2 * len - 1) * sizeof(ulong));
    b = a + len;

    _nmod_poly_reverse(a, poly, len, len);
    _nmod_poly_derivative(b, poly, len, mod);
    _nmod_poly_reverse(b, b, len - 1, len - 1);

    _nmod_poly_div_series(res, b, len - 1, a, len, n, mod);

    flint_free(a);
}

void
nmod_poly_power_sums_schoenhage(nmod_poly_t res, const nmod_poly_t poly,
                                slong n)
{
    if (poly->length == 0)
    {
        flint_throw(FLINT_ERROR, "(nmod_poly_power_sums_schoenhage): Zero polynomial.\n");
    }
    else if ((n <= 0) || (poly->length == 1))
    {
        nmod_poly_zero(res);
    }
    else
    {
        if (*nmod_poly_lead(poly) != 1)
        {
            nmod_poly_t t;
            nmod_poly_init_preinv(t, poly->mod.n, poly->mod.ninv);
            nmod_poly_make_monic(t, poly);
            nmod_poly_fit_length(res, n);
            _nmod_poly_power_sums_schoenhage(res->coeffs, t->coeffs,
                                             t->length, n, t->mod);
            nmod_poly_clear(t);
        }
        else if (poly == res)
        {
            nmod_poly_t t;
            nmod_poly_init_preinv(t, poly->mod.n, poly->mod.ninv);
            nmod_poly_fit_length(t, n);
            _nmod_poly_power_sums_schoenhage(t->coeffs, poly->coeffs,
                                             poly->length, n, t->mod);
            nmod_poly_swap(t, res);
            nmod_poly_clear(t);
        }
        else
        {
            nmod_poly_fit_length(res, n);
            _nmod_poly_power_sums_schoenhage(res->coeffs, poly->coeffs,
                                             poly->length, n, poly->mod);
        }
        _nmod_poly_set_length(res, n);
        _nmod_poly_normalise(res);
    }
}
