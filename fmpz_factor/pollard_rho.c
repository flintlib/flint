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

    Copyright (C) 2015 Kushagra Singh

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"

int
fmpz_factor_pollard_rho(fmpz_t p_factor, flint_rand_t state, const fmpz_t n, 
                        mp_limb_t max_tries, mp_limb_t max_iters)
{

    fmpz_t a, x, maxval;
    int ret;

    fmpz_init(a);
    fmpz_init(x);
    fmpz_init(maxval);
    ret = 0;

    fmpz_sub_ui(maxval, n, 3);
    fmpz_randm(a, state, maxval);
    fmpz_add_ui(a, a, 1);
    fmpz_sub_ui(maxval, n, 1);
    fmpz_randm(x, state, maxval);
    fmpz_add_ui(x, x, 1);

    while (max_tries--)
    {
        ret = fmpz_factor_pollard_rho_ac(p_factor, n, a, x, max_iters);

        if (ret == 1)
            break;
    }

    fmpz_clear(a);
    fmpz_clear(x);
    fmpz_clear(maxval);

    return ret;
}