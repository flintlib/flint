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
#include <math.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "arith.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"


void
arith_number_of_partitions_nmod_vec(mp_ptr res, len_t len, nmod_t mod)
{
    mp_ptr tmp;
    mp_limb_t r;
    len_t k, n;

    r = mod.n - 1UL;

    if (len < 1)
        return;

    tmp = _nmod_vec_init(len);
    _nmod_vec_zero(tmp, len);

    tmp[0] = 1UL;

    for (n = k = 1; n + 4*k + 2 < len; k += 2)
    {
        tmp[n] = r;
        tmp[n + k] = r;
        tmp[n + 3*k + 1] = 1UL;
        tmp[n + 4*k + 2] = 1UL;
        n += 6*k + 5;
    }

    if (n < len) tmp[n] = r;
    if (n + k < len) tmp[n + k] = r;
    if (n + 3*k + 1 < len) tmp[n + 3*k + 1] = 1L;

    _nmod_poly_inv_series(res, tmp, len, mod);

    _nmod_vec_clear(tmp);
}
