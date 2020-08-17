/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

void
_nmod_poly_log_series_monomial_ui(mp_ptr res, mp_limb_t coeff, ulong power,
                                        slong n, nmod_t mod)
{
    slong j, k, rlen;
    mp_limb_t a;

    _nmod_vec_zero(res, n);

    if (power >= n)
        return;

    rlen = (n - 1) / power;
    a = coeff;
    coeff = n_negmod(coeff, mod.n);

    /* Construct geometric series */
    if (coeff == UWORD(1))
    {
        for (j = 0; j < rlen; j++)
            res[j] = a;
    }
    else if (a == UWORD(1))
    {
        for (j = 0; j < rlen; j++)
            res[j] = (j % 2) ? coeff : a;
    }
    else
    {
        for (j = 0; j < rlen; j++)
        {
            res[j] = a;
            a = n_mulmod2_preinv(a, coeff, mod.n, mod.ninv);
        }
    }

    /* Integrate */
    _nmod_poly_integral(res, res, rlen + 1, mod);

    /* Expand */
    if (power != 1)
    {
        for (j = rlen * power + 1; j < n; j++)
            res[j] = UWORD(0);
        for (j = rlen; j > 0; j--)
        {
            res[j * power] = res[j];
            for (k = power; k > 0; k--)
                res[j * power - k] = UWORD(0);
        }
    }
}

void
nmod_poly_log_series_monomial_ui(nmod_poly_t res, mp_limb_t coeff,
                                    ulong power, slong n)
{
    if (power == 0)
    {
        flint_printf("Exception (nmod_poly_log_series_monomial_ui). \n"
               "Constant term != 1.\n");
        flint_abort();
    }

    if (coeff != UWORD(1))
        coeff = n_mod2_preinv(coeff, res->mod.n, res->mod.ninv);

    if (n <= 1 || coeff == UWORD(0))
    {
        nmod_poly_zero(res);
        return;
    }

    nmod_poly_fit_length(res, n);
    _nmod_poly_log_series_monomial_ui(res->coeffs, coeff, power, n, res->mod);
    res->length = n;
    _nmod_poly_normalise(res);
}
