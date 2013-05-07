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

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_factor.h"
#include "arith.h"
#include "ulong_extras.h"


void arith_euler_phi(fmpz_t res, const fmpz_t n)
{
    fmpz_factor_t factors;
    fmpz_t t;
    ulong exp;
    long i;

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

    fmpz_factor_init(factors);
    fmpz_factor(factors, n);
    fmpz_one(res);

    fmpz_init(t);
    for (i = 0; i < factors->num; i++)
    {
        fmpz_sub_ui(t, factors->p + i, 1UL);
        fmpz_mul(res, res, t);
        exp = factors->exp[i];
        if (exp != 1)
        {
            fmpz_pow_ui(t, factors->p + i, exp - 1UL);
            fmpz_mul(res, res, t);
        }
    }

    fmpz_clear(t);
    fmpz_factor_clear(factors);
}
