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

/* Sets f to (x^2 + a) % n */

void
sq_and_add_a(fmpz_t f, fmpz_t x, fmpz_t a, fmpz_t n)
{
    fmpz_powm_ui(f, x, 2, n);
    fmpz_add(f, f, a);
    fmpz_mod(f, f, n);
}

int
fmpz_factor_pollard_rho_brent(fmpz_t n, fmpz_t p_factor)
{

    fmpz_t x, y, a, m, q, gcdval, ys, subval, minval;
    mp_limb_t iter, i, k;
    int ret;

    FLINT_TEST_INIT(state);

    fmpz_init(x);
    fmpz_init(y);
    fmpz_init(ys);
    fmpz_init(a);
    fmpz_init(m);
    fmpz_init(subval);
    fmpz_init(minval);
    fmpz_init_set_ui(q, 1);
    fmpz_init_set_ui(gcdval, 1);

    fmpz_randm(a, state, n - 3);
    fmpz_randm(y, state, n - 1);
    fmpz_randm(m, state, n - 1);

    iter = 1;
    i = 0;

    do {
        fmpz_set(x, y);

        for (i = 0; i < iter; i++)
            sq_and_add_a(y, y, a, n);

        k = 0;
        do {
            fmpz_set_ui(minval, iter - k);   /* minval = min(m, iter - k) */
            
            if (fmpz_cmp(m, minval) < 0)
                fmpz_set(minval, m);
            
            fmpz_set(ys, y);
            
            for (i = 0; i < fmpz_get_ui(minval); i++)
            {
                sq_and_add_a(y, y, a, n);
                fmpz_sub(subval, x, y);
                fmpz_abs(subval, subval);
                fmpz_mul(q, q, subval);
                fmpz_mod(q, q, n);
            }
            fmpz_gcd(gcdval, q, n);
            k += fmpz_get_ui(m);
        } while ((k < iter) && !fmpz_cmp_ui(gcdval, 1));
        iter *= 2;
    } while(!fmpz_cmp_ui(gcdval, 1));   /* while gcdval == 1 */

    if (!fmpz_cmp(gcdval, n))
    {
        do {
            sq_and_add_a(ys, ys, a, n);
            fmpz_sub(subval, x, ys);
            fmpz_abs(subval, subval);
            fmpz_gcd(gcdval, subval, n);
        } while (!fmpz_cmp_ui(gcdval, 1));
    }

    if (!fmpz_cmp(gcdval, n))
        ret = 0;
    else
    {
        ret = 1;
        fmpz_set(p_factor, gcdval);
    }

    fmpz_clear(x);
    fmpz_clear(y);
    fmpz_clear(ys);
    fmpz_clear(a);
    fmpz_clear(m);
    fmpz_clear(subval);
    fmpz_clear(minval);
    fmpz_clear(q);
    fmpz_clear(gcdval);

    return ret;
}