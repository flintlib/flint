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

/* Select Montgomery Elliptic Curve given a sigma
   (Suyama's parameterization) 
   Returns 1 in case factor is found while selecting
   the curev. */

void
fmpz_ec_select_curve(fmpz_t f, fmpz_t x, fmpz_t a, fmpz_t sig, fmpz_t n)
{
    fmpz_t u, v, w, y, z;
    int ret;
    ret = 0;

    fmpz_init_set(u, sig);
    fmpz_init(v);
    fmpz_init(w);
    fmpz_init(y);
    fmpz_init(z);

    fmpz_mul_2exp(v, u, 2);     /* v = sig * 4 */
    fmpz_mod(v, v, n);
    fmpz_mul(w, u, u);
    fmpz_mod(w, w, n);
    fmpz_sub_ui(u, w, 5);       /* u = sig^2 - 5 */

    fmpz_mul(w, u, u);          /* w = u * u */
    fmpz_mod(w, w, n);
    fmpz_mul(x, w, u);          /* x = u * u * u */
    fmpz_mod(x, x, n);

    fmpz_mul(w, v, v);          /* w = v * v */
    fmpz_mod(w, w, n);
    fmpz_mul(z, w, v);          /* z = v * v * v */
    fmpz_mod(z, z, n);

    fmpz_mul(w, x, v);
    fmpz_mod(w, w, n);
    fmpz_mul_2exp(y, w, 2);
    fmpz_mod(y, y, n);
    fmpz_mul_ui(w, u, 3);
    fmpz_mod(w, w, n);
    fmpz_sub(u, v, u);

    fmpz_add(v, v, w);
    if (fmpz_cmp(v, n) >= 0)
        fmpz_sub(v, v, n);
    fmpz_mul(w, u, u);
    fmpz_mod(w, w ,n);
    fmpz_mul(u, u, w);
    fmpz_mod(u, u, n);
    fmpz_mul(a, u, v);
    fmpz_mod(a, a, n);

    fmpz_mul(v, y, z);
    fmpz_mod(v, v, n);
    fmpz_gcdinv(f, u, v, n);

    if (fmpz_is_one(f) == 0)
    {
        ret = 1;
        goto cleanup;
    }

    fmpz_mul(v, u, y);
    fmpz_mod(v, v, n);
    fmpz_mul(x, x, v);
    fmpz_mod(x, x, n);

    fmpz_mul(v, u, z);
    fmpz_mod(v, v, n);
    fmpz_mul(w, a, v);
    fmpz_mod(w, w, n);
    fmpz_sub_ui(a, w, 2);


    cleanup:

    fmpz_clear(u);
    fmpz_clear(v);
    fmpz_clear(w);   
    fmpz_clear(y);   
    fmpz_clear(z); 

    return ret; 

}