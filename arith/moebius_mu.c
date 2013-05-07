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


int arith_moebius_mu(const fmpz_t n)
{
    fmpz_factor_t factors;
    long i;
    int mu;

    if (fmpz_abs_fits_ui(n))
        return n_moebius_mu(fmpz_get_ui(n));

    fmpz_factor_init(factors);
    fmpz_factor(factors, n);

    mu = 1;
    for (i = 0; i < factors->num; i++)
    {
        if (factors->exp[i] != 1UL)
        {
            mu = 0;
            break;
        }
    }

    if (factors->num % 2)
        mu = -mu;

    fmpz_factor_clear(factors);
    return mu;
}
