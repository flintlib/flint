/*============================================================================

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

===============================================================================*/
/****************************************************************************

   Copyright (C) 2005-2008 Damien Stehle.
   Copyright (C) 2007 David Cade.
   Copyright (C) 2010 William Hart
   
*****************************************************************************/

#include <stdlib.h>
#include <mpir.h>
#include <mpfr.h>
#include "flint.h"
#include "mpfr_vec.h"
#include "mpfr_mat.h"
#include "fmpz_mat.h"

void GSO_mpfr(mpfr_mat_t r, mpfr_t max, mpfr_mat_t mu, fmpz_mat_t G, ulong a, ulong kappa, ulong zeroes)
{
  ulong i, j;
  mpfr_t tmp, rtmp;

  mpfr_init2(tmp, mu->prec);
  mpfr_init2(rtmp, mu->prec);
  
  mpfr_set_ui(max, 0, GMP_RNDN);
  
  for (i = a; i < kappa; i++)
  {
      if (i >= zeroes + 2)
	  {
	     mpfr_set_z(rtmp, COEFF_TO_PTR(G->rows[kappa][i]), GMP_RNDN);
		 _mpfr_vec_scalar_product(tmp, mu->rows[i] + zeroes + 1, r->rows[kappa] + zeroes + 1, i - zeroes - 1);
		 mpfr_sub(r->rows[kappa] + i, rtmp, tmp, GMP_RNDN);
	  } else
	     mpfr_set_z(r->rows[kappa] + i, COEFF_TO_PTR(G->rows[kappa][i]), GMP_RNDN);
      
      mpfr_div(mu->rows[kappa] + i, r->rows[kappa] + i, r->rows[i] + i, GMP_RNDN);
      if (mpfr_cmpabs(mu->rows[kappa] + i, max) > 0)
		 mpfr_set(max, mu->rows[kappa] + i, GMP_RNDN);
   }

   mpfr_abs(max, max, GMP_RNDN);

   mpfr_clear(tmp);
   mpfr_clear(rtmp);
}


