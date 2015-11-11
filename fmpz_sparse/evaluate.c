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

#include "fmpz.h"
#include "fmpz_sparse.h"

void fmpz_sparse_evaluate(fmpz_t res, const fmpz_sparse_t f, const fmpz_t a)
{
  ulong powdiff;
  fmpz_t apow;
  fmpz_t temp;
  slong i;

  fmpz_zero(res);

  if (f->length == 0) return;

  if (fmpz_is_zero(a)) return;

  if (fmpz_is_one(a))
  {
    /* special case for f(1) */
    for (i=0; i<f->length; ++i) fmpz_add(res, res, f->coeffs+i);
    return;
  }

  else if (fmpz_is_pm1(a))
  {
    /* special case for f(-1) */
    for (i=0; i<f->length; ++i)
    {
      if (fmpz_is_even(f->expons+i)) fmpz_add(res, res, f->coeffs+i);
      else fmpz_sub(res, res, f->coeffs+i);
    }
    return;
  }

  i = f->length - 1;
  FLINT_ASSERT(fmpz_abs_fits_ui(f->expons+0));
  FLINT_ASSERT(fmpz_sgn(f->expons+i) == 1);

  powdiff = fmpz_get_ui(f->expons+i);
  fmpz_init_set_ui(apow, UWORD(1));
  fmpz_init(temp);

  while (1) {
    fmpz_pow_ui(temp, a, powdiff);
    fmpz_mul(apow, apow, temp);
    fmpz_addmul(res, f->coeffs+i, apow);

    if (i == 0) break;

    --i;
    powdiff = fmpz_get_ui(f->expons+i) - fmpz_get_ui(f->expons+(i+1));
  }

  fmpz_clear(apow);
  fmpz_clear(temp);
}
