/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "nmod_poly_mini.h"
#include "flint-impl.h"

void
_nmod_poly_exp_series_monomial_ui(ulong_ptr res, ulong coeff, ulong power,
                                        slong n, nmod_t mod)
{
    slong k, r;
    ulong rfac;
    ulong a;

    r = (n - 1) / power;
    rfac = n_factorial_mod2_preinv(r, mod.n, mod.ninv);
    rfac = n_invmod(rfac, mod.n);

    if (power > 1)
        _NMOD_VEC_ZERO(res, n);

    res[0] = UWORD(1);

    if (coeff == UWORD(1))
    {
        a = rfac;
        for (k = r; k >= 1; k--)
        {
            res[k * power] = a;
            a = n_mulmod2_preinv(a, k, mod.n, mod.ninv);
        }
    }
    else
    {
        a = coeff;
        for (k = power; k < n; k += power)
        {
            res[k] = a;
            a = n_mulmod2_preinv(a, coeff, mod.n, mod.ninv);
        }

        a = rfac;
        for (k = r; k >= 1; k--)
        {
            res[k * power] = n_mulmod2_preinv(res[k * power],
                a, mod.n, mod.ninv);
            a = n_mulmod2_preinv(a, k, mod.n, mod.ninv);
        }
    }
}

void
nmod_poly_exp_series_monomial_ui(nmod_poly_t res, ulong coeff,
                                    ulong power, slong n)
{
    if (n == 0)
    {
        nmod_poly_zero(res);
        return;
    }

    if (coeff == UWORD(0))
    {
        nmod_poly_fit_length(res, 1);
        res->coeffs[0] = UWORD(1);
        res->length = 1;
        return;
    }

    if (power == 0)
        flint_throw(FLINT_ERROR, "Constant term != 0 in nmod_poly_exp_series_monomial_ui\n");

    if (coeff != UWORD(1))
        coeff = n_mod2_preinv(coeff, res->mod.n, res->mod.ninv);

    if (n == 1 || power >= n)
    {
        nmod_poly_fit_length(res, 1);
        res->coeffs[0] = UWORD(1);
        res->length = 1;
    }

    nmod_poly_fit_length(res, n);
    _nmod_poly_exp_series_monomial_ui(res->coeffs, coeff, power, n, res->mod);
    res->length = n;
    _nmod_poly_normalise(res);
}
