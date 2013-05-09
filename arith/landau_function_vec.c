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

#include <math.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "arith.h"
#include "ulong_extras.h"
#include "longlong.h"


void arith_landau_function_vec(fmpz * res, len_t len)
{
    mp_limb_t p, pmax;
    mp_limb_t pk, pkhi;
    fmpz_t a;
    ulong k, n;

    if (len < 1)
        return;

    for (k = 0; k < len; k++)
        fmpz_one(res + k);

    pmax = 1.328 * sqrt(len*log(len) + 1);

    fmpz_init(a);

    for (p = 2UL; p <= pmax; p = n_nextprime(p, 0))
    {
        for (n = len - 1; n >= p; n--)
        {
            pk = p;
            pkhi = 0UL;

            for (k = 1; k <= len; k++)
            {
                if (pk > n || pkhi)
                    break;

                fmpz_mul_ui(a, res + n - pk, pk);
                if (fmpz_cmp(res + n, a) < 0)
                    fmpz_set(res + n, a);

                umul_ppmm(pkhi, pk, pk, p);
            }
        }
    }

    fmpz_clear(a);
}
