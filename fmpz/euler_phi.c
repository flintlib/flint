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

    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#include "fmpz.h"

void
fmpz_factor_euler_phi(fmpz_t res, const fmpz_factor_t fac)
{
    fmpz_t t;
    slong i;

    fmpz_init(t);

    fmpz_one(res);

    for (i = 0; i < fac->num; i++)
    {
        fmpz_sub_ui(t, fac->p + i, 1);
        fmpz_mul(res, res, t);

        if (fac->exp[i] != 1)
        {
            fmpz_pow_ui(t, fac->p + i, fac->exp[i] - 1);
            fmpz_mul(res, res, t);
        }
    }

    fmpz_clear(t);
}

void
fmpz_euler_phi(fmpz_t res, const fmpz_t n)
{
    fmpz_factor_t fac;

    if (fmpz_sgn(n) <= 0)
    {
        fmpz_zero(res);
        return;
    }

    if (fmpz_abs_fits_ui(n))
    {
        fmpz_set_ui(res, n_euler_phi(fmpz_get_ui(n)));
        return;
    }

    fmpz_factor_init(fac);
    fmpz_factor(fac, n);
    fmpz_factor_euler_phi(res, fac);
    fmpz_factor_clear(fac);
}

