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

#include <stdio.h>
#include "flint.h"
#include "fmpz.h"
#include "arith.h"
#include "ulong_extras.h"

void fmpz_bernoulli_denom(fmpz_t den, ulong n)
{
    long i;
    mp_limb_t p;

    if (n % 2 == 1 || n == 0)
    {
        fmpz_set_ui(den, (n == 1) ? 2 : 1);
        return;
    }

    n_compute_primes(n);

    fmpz_set_ui(den, 6UL);
    for (i = 2; i < n; i++)
    {
        p = flint_primes[i];
        if (p - 1 > n)
            break;
        if (n % (p - 1) == 0)
            fmpz_mul_ui(den, den, p);
    }
}
