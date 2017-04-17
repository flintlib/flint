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

    Authored 2015 by Daniel S. Roche; US Government work in the public domain. 

******************************************************************************/

#include "fmpz_sparse.h"

void fmpz_sparse_rem_cyc_mod_diverse(nmod_poly_t res, 
    const fmpz_sparse_t poly, ulong a, ulong e, ulong q)
{
    slong i;
    ulong q_inv = n_preinvert_limb(q);
    FLINT_ASSERT (n_powmod2_preinv(a, q-1, q, q_inv) == 1UL);
    if (nmod_poly_modulus(res) == q) nmod_poly_zero(res);
    else {
        /* polynomial has the wrong modulus, must be re-initialized */
        nmod_poly_clear(res);
        nmod_poly_init2(res, q, e+1); /* initializes allocated size */
    }
    nmod_poly_set_coeff_ui(res, e, 1); /* set to x^e, to reserve zeros */
    for (i=0; i < poly->length; ++i) {
        ulong rese = fmpz_fdiv_ui(poly->expons + i, e);
        ulong resc = n_powmod2_preinv(
                a, fmpz_fdiv_ui(poly->expons+i, q-1), q, q_inv);
        resc = n_mulmod2_preinv(
                resc, fmpz_fdiv_ui(poly->coeffs+i, q), q, q_inv);
        resc = n_addmod(
                nmod_poly_get_coeff_ui(res, rese), resc, q);
        nmod_poly_set_coeff_ui(res, rese, resc);
    }
    nmod_poly_set_coeff_ui(res, e, 0); /* set x^e coeff back to 0 */
}
