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
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "arith.h"
#include "ulong_extras.h"


void fmpz_divisor_sigma(fmpz_t res, ulong n, ulong k)
{
    int i;
    fmpz_t p, r;

    if (n == 0)
    {
        fmpz_zero(res);
        return;
    }

    n_factor_t factors;
    n_factor_init(&factors);
    n_factor(&factors, n, 0);

    fmpz_init(p);
    fmpz_init(r);
    fmpz_set_ui(res, 1UL);

    // TODO: use a balanced product in large cases,
    // speedup for small n.
    for (i = 0; i<factors.num; i++)
    {
        if (k == 0)
            fmpz_mul_ui(res, res, factors.exp[i]+1);
        else
        {
            fmpz_set_ui(p, factors.p[i]);
            fmpz_pow_ui(p, p, k);
            fmpz_pow_ui(r, p, factors.exp[i]+1);
            fmpz_sub_ui(r, r, 1);
            fmpz_sub_ui(p, p, 1);
            fmpz_divexact(p, r, p);
            fmpz_mul(res, res, p);
        }
    }
    fmpz_clear(p);
    fmpz_clear(r);
}
