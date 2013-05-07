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

    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2010 William Hart
    Copyright (C) 2012 Lina Kulakova

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz_vec.h"
#include "fmpz_mod_poly.h"

void
_fmpz_mod_poly_pow_trunc(fmpz * res, fmpz * poly,
                         ulong e, long trunc, const fmpz_t p)
{
    _fmpz_mod_poly_pow_trunc_binexp(res, poly, e, trunc, p);
}

void
fmpz_mod_poly_pow_trunc(fmpz_mod_poly_t res,
                  const fmpz_mod_poly_t poly, ulong e, long trunc)
{
    const long len = poly->length;
    fmpz * q;
    int qcopy = 0;

    if (len < 2 || e < 3UL || trunc == 0)
    {
        if (len == 0 || trunc == 0)
            fmpz_mod_poly_zero(res);
        else if (len == 1)
        {
            fmpz_mod_poly_fit_length(res, 1);
            fmpz_powm_ui(res->coeffs, poly->coeffs, e, &res->p);
            _fmpz_mod_poly_set_length(res, 1);
            _fmpz_mod_poly_normalise(res);
        }
        else if (e == 0UL)
        {
            fmpz_mod_poly_set_coeff_ui(res, 0, 1UL);
            _fmpz_mod_poly_set_length(res, 1);
            _fmpz_mod_poly_normalise(res);
        }
        else if (e == 1UL)
        {
            fmpz_mod_poly_set(res, poly);
            fmpz_mod_poly_truncate(res, trunc);
        }
        else  /* e == 2UL */
            fmpz_mod_poly_mullow(res, poly, poly, trunc);

        return;
    }

    if (poly->length < trunc)
    {
        q = _fmpz_vec_init(trunc);
        _fmpz_vec_set(q, poly->coeffs, poly->length);
        _fmpz_vec_zero(q + poly->length, trunc - poly->length);
        qcopy = 1;
    } else
        q = poly->coeffs;

    if (res != poly || qcopy)
    {
        fmpz_mod_poly_fit_length(res, trunc);
        _fmpz_mod_poly_pow_trunc(res->coeffs, q, e, trunc, &poly->p);
    }
    else
    {
        fmpz_mod_poly_t t;
        fmpz_mod_poly_init2(t, &poly->p, trunc);
        _fmpz_mod_poly_pow_trunc(t->coeffs, q, e, trunc, &poly->p);
        fmpz_mod_poly_swap(res, t);
        fmpz_mod_poly_clear(t);
    }

    if (qcopy)
        _fmpz_vec_clear(q, trunc);

    _fmpz_mod_poly_set_length(res, trunc);
    _fmpz_mod_poly_normalise(res);
}

