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

    Copyright (C) 2013 William Hart
   
******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

void
_fmpz_poly_preinvert(fmpz * B_inv, const fmpz * B, slong n)
{
   fmpz * X2n = _fmpz_vec_init(2*n - 1);

   fmpz_set_ui(X2n + 2*n - 2, 1);
   _fmpz_poly_div(B_inv, X2n, 2*n - 1, B, n);

   _fmpz_vec_clear(X2n, 2*n - 1);

   _fmpz_poly_reverse(B_inv, B_inv, n, n);
}

void
fmpz_poly_preinvert(fmpz_poly_t B_inv, const fmpz_poly_t B)
{
    slong n = B->length;
    fmpz_poly_t temp;
    fmpz * Binv_coeffs;

    if (n == 0)
    {
        flint_printf("Exception (fmpz_poly_preinvert). Division by zero.\n");
        abort();
    }

    if (B == B_inv)
    {
       fmpz_poly_init2(temp, n);
       Binv_coeffs = temp->coeffs;
    } else
    {
       fmpz_poly_fit_length(B_inv, n);
       Binv_coeffs = B_inv->coeffs;
    }

    _fmpz_poly_preinvert(Binv_coeffs, B->coeffs, n);


    if (B == B_inv)
    {
       _fmpz_poly_set_length(temp, n);
       fmpz_poly_swap(B_inv, temp);
       fmpz_poly_clear(temp);
    } else
       _fmpz_poly_set_length(B_inv, n);

    /* no need to normalise */
}
