/*
    Copyright (C) 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <math.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

double _fmpz_poly_evaluate_horner_d_2exp(slong * exp, const fmpz * poly, slong n, double d)
{
   slong i, size_p = FLINT_ABS(_fmpz_vec_max_bits(poly, n));
   ulong vbits = ceil(fabs(log(fabs(d))/log(2.0)));
   ulong prec = vbits*(n - 1) + size_p + FLINT_BIT_COUNT(n);
   mpf_t mpf_d, mpf_coeff, output;
   double res;
   
   if (d == 0)
      return fmpz_get_d(poly + 0);

   mpf_set_default_prec(prec);
   
   mpf_init(mpf_coeff);   
   mpf_init(output);   
   mpf_init(mpf_d);

   fmpz_get_mpf(output, poly + n - 1);
   
   mpf_set_d(mpf_d, d); /* set fval to mpf from the double val */
   
   for (i = n - 2; i >= 0; i--)
   {
      mpf_mul(output, output, mpf_d);
      fmpz_get_mpf(mpf_coeff, poly + i);    
      mpf_add(output, output, mpf_coeff);
   }

   res = mpf_get_d_2exp(exp, output);
   
   if (mpf_sgn(output) < 0 && res >= 0.0)
      res = -res; /* work around bug in earlier versions of GMP/MPIR */

   mpf_clear(mpf_coeff);
   mpf_clear(output);
   mpf_clear(mpf_d);

   return res;
}

double fmpz_poly_evaluate_horner_d_2exp(slong * exp, const fmpz_poly_t poly, double d)
{
   if (poly->length == 0)
   {
      *exp = 0;
      return 0.0;
   }

   return _fmpz_poly_evaluate_horner_d_2exp(exp, poly->coeffs,
                                                              poly->length, d);
}
