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

void fmpz_sparse_evaluate_mod(fmpz_t res, const fmpz_sparse_t poly, 
    const fmpz_t a, const fmpz_t m)
{
  fmpz_t powdiff;
  fmpz_t apow;
  fmpz_t temp;
  slong i = poly->length - 1;

  fmpz_zero(res);

  if (poly->length == 0) return;

  fmpz_init_set(powdiff, poly->expons+i);
  fmpz_init_set_ui(apow, UWORD(1));
  fmpz_init2(temp, fmpz_size(m));

  while (1) {
    fmpz_powm(temp, a, powdiff, m);
    fmpz_mul(apow, apow, temp);
    fmpz_mod(apow, apow, m);
    fmpz_addmul(res, poly->coeffs+i, apow);

    if (i == 0) break;

    --i;
    fmpz_sub (powdiff, poly->expons+i, poly->expons+(i+1));
  }

  fmpz_mod(res, res, m);

  fmpz_clear(powdiff);
  fmpz_clear(apow);
  fmpz_clear(temp);
}
