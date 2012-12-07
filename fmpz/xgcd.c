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

    Copyright (C) 2012 William Hart

******************************************************************************/

#include <mpir.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

void fmpz_xgcd(fmpz_t d, fmpz_t a, fmpz_t b, const fmpz_t f, const fmpz_t g)
{
   fmpz_t t;
   fmpz_init(t);

   if (fmpz_cmp(f, g) < 0)
   {
      fmpz_gcdinv(d, a, f, g);
      fmpz_mul(t, a, f);
      fmpz_sub(t, d, t);
      fmpz_divexact(b, t, g);
   } else if (fmpz_cmp(g, f) < 0)
   {
      fmpz_gcdinv(d, b, g, f);
      fmpz_mul(t, b, g);
      fmpz_sub(t, d, t);
      fmpz_divexact(a, t, f);
   } else /* f == g */
   {
      fmpz_set(d, f);
      fmpz_set_ui(a, 1);
      fmpz_set_si(b, 0);
   }

   fmpz_clear(t);
}
