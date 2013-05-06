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

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

void
_nmod_poly_log_series(mp_ptr res, mp_srcptr f, long n, nmod_t mod)
{
    mp_ptr f_diff;
    mp_ptr f_inv;

    f_diff = _nmod_vec_init(n);
    f_inv = _nmod_vec_init(n);

    _nmod_poly_derivative(f_diff, f, n, mod); f_diff[n-1] = 0UL;
    _nmod_poly_inv_series(f_inv, f, n, mod);
    _nmod_poly_mullow(res, f_diff, n - 1, f_inv, n - 1, n - 1, mod);
    _nmod_poly_integral(res, res, n, mod);

    _nmod_vec_clear(f_diff);
    _nmod_vec_clear(f_inv);
}

void
nmod_poly_log_series(nmod_poly_t res, const nmod_poly_t f, long n)
{
    mp_ptr f_coeffs;
    long k;
    long flen = f->length;

    if (flen < 1 || f->coeffs[0] != 1UL)
    {
        printf("Exception (nmod_poly_log_series). Constant term != 1.\n");
        abort();
    }

    if (flen == 1 || n < 2)
    {
        nmod_poly_zero(res);
        return;
    }

    nmod_poly_fit_length(res, n);

    /* Efficiently handle monomials */
    for (k = 1; f->coeffs[k] == 0UL && k < n - 1; k++);
    if (k == flen - 1 || k == n - 1)
    {
        flen = FLINT_MIN(flen, n);
        _nmod_poly_log_series_monomial_ui(res->coeffs,
            f->coeffs[flen-1], flen - 1, n, res->mod);
    }
    else
    {
        if (flen < n)
        {
            f_coeffs = _nmod_vec_init(n);
            flint_mpn_copyi(f_coeffs, f->coeffs, flen);
            flint_mpn_zero(f_coeffs + flen, n - flen);
        }
        else
            f_coeffs = f->coeffs;

        _nmod_poly_log_series(res->coeffs, f_coeffs, n, res->mod);

        if (flen < n)
            _nmod_vec_clear(f_coeffs);
    }

    res->length = n;
	_nmod_poly_normalise(res);
}
