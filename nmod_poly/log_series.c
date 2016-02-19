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

    Copyright (C) 2011, 2016 Fredrik Johansson

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

void
_nmod_poly_log_series(mp_ptr res, mp_srcptr f, slong flen, slong n, nmod_t mod)
{
    flen = FLINT_MIN(flen, n);

    if (flen == 1)
    {
        res[0] = 1;
        _nmod_vec_zero(res + 1, n - 1);
    }
    else
    {
        mp_ptr tmp = _nmod_vec_init(2 * n);
        _nmod_poly_derivative(tmp, f, flen, mod);
        _nmod_poly_div_series(tmp + n, tmp, flen - 1, f,
            FLINT_MIN(flen, n - 1), n - 1, mod);
        _nmod_poly_integral(res, tmp + n, n, mod);
        _nmod_vec_clear(tmp);
    }
}

void
nmod_poly_log_series(nmod_poly_t res, const nmod_poly_t f, slong n)
{
    slong flen = f->length;

    if (flen < 1 || f->coeffs[0] != UWORD(1))
    {
        flint_printf("Exception (nmod_poly_log_series). Constant term != 1.\n");
        abort();
    }

    if (flen == 1 || n < 2)
    {
        nmod_poly_zero(res);
        return;
    }

    nmod_poly_fit_length(res, n);
    _nmod_poly_log_series(res->coeffs, f->coeffs, f->length, n, res->mod);
    res->length = n;
	_nmod_poly_normalise(res);
}

