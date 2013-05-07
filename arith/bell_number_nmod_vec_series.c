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

#include <gmp.h>
#include "flint.h"
#include "arith.h"
#include "ulong_extras.h"

void
arith_bell_number_nmod_vec_series(mp_ptr res, long n, nmod_t mod)
{
    mp_limb_t fac, c;
    mp_ptr tmp;
    long k;

    if (n < 1)
        return;

    tmp = flint_malloc(sizeof(mp_limb_t) * n);

    /* Divide by factorials */
    fac = n_factorial_mod2_preinv(n-1, mod.n, mod.ninv);
    c = n_invmod(fac, mod.n);

    for (k = n - 1; k > 0; k--)
    {
        tmp[k] = c;
        c = n_mulmod2_preinv(c, k, mod.n, mod.ninv);
    }
    tmp[0] = 0UL;

    _nmod_poly_exp_series(res, tmp, n, mod);

    /* Multiply by factorials */
    c = 1UL;
    for (k = 1; k < n; k++)
    {
        c = n_mulmod2_preinv(c, k, mod.n, mod.ninv);
        res[k] = n_mulmod2_preinv(res[k], c, mod.n, mod.ninv);
    }

    flint_free(tmp);
}
