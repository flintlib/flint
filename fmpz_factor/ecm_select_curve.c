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

/* Also selects initial point Q0 [x0 :: z0]  (z0 = 1) */

int
fmpz_factor_ecm_select_curve(fmpz_t f, fmpz_t sig, fmpz_t n, ecm_t ecm_inf)
{
    fmpz_set(ecm_inf->u, sig);

    fmpz_mul_2exp(ecm_inf->v, ecm_inf->u, 2);     /* v = sig * 4 */
    fmpz_mod(ecm_inf->v, ecm_inf->v, n);
    fmpz_mul(ecm_inf->w, ecm_inf->u, ecm_inf->u);
    fmpz_mod(ecm_inf->w, ecm_inf->w, n);
    fmpz_sub_ui(ecm_inf->u, ecm_inf->w, 5);       /* u = sig^2 - 5 */

    fmpz_mul(ecm_inf->w, ecm_inf->u, ecm_inf->u);          /* w = u * u */
    fmpz_mod(ecm_inf->w, ecm_inf->w, n);
    fmpz_mul(ecm_inf->x, ecm_inf->w, ecm_inf->u);          /* x = u * u * u */
    fmpz_mod(ecm_inf->x, ecm_inf->x, n);

    fmpz_mul(ecm_inf->w, ecm_inf->v, ecm_inf->v);          /* w = v * v */
    fmpz_mod(ecm_inf->w, ecm_inf->w, n);
    fmpz_mul(ecm_inf->z, ecm_inf->w, ecm_inf->v);          /* z = v * v * v */
    fmpz_mod(ecm_inf->z, ecm_inf->z, n);

    fmpz_mul(ecm_inf->w, ecm_inf->x, ecm_inf->v);
    fmpz_mod(ecm_inf->w, ecm_inf->w, n);
    fmpz_mul_2exp(ecm_inf->t, ecm_inf->w, 2);
    fmpz_mod(ecm_inf->t, ecm_inf->t, n);
    fmpz_mul_ui(ecm_inf->w, ecm_inf->u, 3);
    fmpz_mod(ecm_inf->w, ecm_inf->w, n);
    fmpz_sub(ecm_inf->u, ecm_inf->v, ecm_inf->u);

    fmpz_add(ecm_inf->v, ecm_inf->v, ecm_inf->w);
    if (fmpz_cmp(ecm_inf->v, n) >= 0)
        fmpz_sub(ecm_inf->v, ecm_inf->v, n);
    fmpz_mul(ecm_inf->w, ecm_inf->u, ecm_inf->u);
    fmpz_mod(ecm_inf->w, ecm_inf->w ,n);
    fmpz_mul(ecm_inf->u, ecm_inf->u, ecm_inf->w);
    fmpz_mod(ecm_inf->u, ecm_inf->u, n);
    fmpz_mul(ecm_inf->a24, ecm_inf->u, ecm_inf->v);
    fmpz_mod(ecm_inf->a24, ecm_inf->a24, n);

    fmpz_mul(ecm_inf->v, ecm_inf->t, ecm_inf->z);
    fmpz_mod(ecm_inf->v, ecm_inf->v, n);
    fmpz_gcdinv(f, ecm_inf->u, ecm_inf->v, n);

    if (fmpz_is_one(f) == 0)
        return 1;

    fmpz_mul(ecm_inf->v, ecm_inf->u, ecm_inf->t);
    fmpz_mod(ecm_inf->v, ecm_inf->v, n);
    fmpz_mul(ecm_inf->x, ecm_inf->x, ecm_inf->v);
    fmpz_mod(ecm_inf->x, ecm_inf->x, n);

    fmpz_mul(ecm_inf->v, ecm_inf->u, ecm_inf->z);
    fmpz_mod(ecm_inf->v, ecm_inf->v, n);
    fmpz_mul(ecm_inf->w, ecm_inf->a24, ecm_inf->v);
    fmpz_mod(ecm_inf->w, ecm_inf->w, n);
    fmpz_sub_ui(ecm_inf->a24, ecm_inf->w, 2);

    fmpz_add_ui(ecm_inf->a24, ecm_inf->a24, 2);
    fmpz_fdiv_q_2exp(ecm_inf->a24, ecm_inf->a24, 2);
    fmpz_set_ui(ecm_inf->z, 1);

    return 0;
}
