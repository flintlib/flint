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
#include "fmpz_vec.h"
#include "fmpz_factor.h"
#include "arith.h"
#include "ulong_extras.h"

/* note: destroys factors! */
void
_arith_divisor_sigma(fmpz_t res, const fmpz_factor_t factors, ulong k)
{
    long i;
    fmpz * p;
    fmpz_t r;

    fmpz_one(res);

    if (factors->num == 0)
        return;

    fmpz_init(r);

    if (k == 0)
    {
        for (i = 0; i < factors->num; i++)
        {
            fmpz_set_ui(r, factors->exp[i] + 1UL);
            fmpz_mul(res, res, r);
        }
        return;
    }
    else
    {
        for (i = 0; i < factors->num; i++)
        {
            p = factors->p + i;
            fmpz_set(p, factors->p + i);
            fmpz_pow_ui(p, p, k);
            fmpz_pow_ui(r, p, factors->exp[i]  + 1UL);
            fmpz_sub_ui(r, r, 1UL);
            fmpz_sub_ui(p, p, 1UL);
            fmpz_divexact(p, r, p);
        }

        _fmpz_vec_prod(res, factors->p, factors->num);
    }

    fmpz_clear(r);
}

void
arith_divisor_sigma(fmpz_t res, const fmpz_t n, ulong k)
{
    fmpz_factor_t factors;

    if (fmpz_is_zero(n))
    {
        fmpz_zero(res);
        return;
    }

    fmpz_factor_init(factors);
    fmpz_factor(factors, n);
    _arith_divisor_sigma(res, factors, k);
    fmpz_factor_clear(factors);
}
