/*
    Copyright (C) 2010, 2016, 2020 William Hart
    Copyright (C) 2010 Andy Novocin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

#define CLD_EPS 0.00000001

static double _log2(double n)  
{  
    return log(n)/log(2);  
}

static int _d_cmp_2exp(double a, slong a_exp, double b, slong b_exp)
{
   long t;

   if (a_exp == 0)
   {
      if (b_exp == 0)
      {
         if (a > 1.5*b)
            return 2;
         else if (b > 1.5*a)
            return -2;
         else if (a >= b)
            return 1;
         else
            return -1;
      }

      t = 1 + (long) _log2(a);

      if (t >= b_exp + 2) /* a, a_exp >= 2*(b, b_exp) */ 
         return 2;
      else if (b_exp >= t + 2) /* b, b_exp >= 2*(a, a_exp) */
         return -2;
      else /* divide both values by 4 and try again */
         return _d_cmp_2exp(a/4, 0, b*pow(2.0, (double) b_exp - 2.0), 0);
   } else if (b_exp == 0)
      return -_d_cmp_2exp(b, b_exp, a, a_exp);
   else /* neither a_exp not b_exp is zero */
   {
      if (a_exp >= b_exp + 2) /* a, a_exp >= 2*(b, b_exp) */ 
         return 2;
      else if (b_exp >= a_exp + 2) /* b, b_exp >= 2*(a, a_exp) */
         return -2;
      else if (a_exp >= b_exp) /* shift exponents so one is 0 */
         return _d_cmp_2exp(a, a_exp - b_exp, b, 0);
      else
         return -_d_cmp_2exp(b, b_exp - a_exp, a, 0);
   }
}

void fmpz_poly_CLD_bound(fmpz_t res, const fmpz_poly_t f, slong n)
{
   /*
      Given: f = a_0 + ... + a_N x^N and n in {0, 1, ..., N - 1}
      Let B_1(r) = 1/(r^{n+1})(|a_0| + ... + |a_n|r^n)
      Let B_2(r) = 1/(r^{n+1})(|a_{n + 1}|r^{n+1} + ... + |a_N|r^N)
      Find r > 0 such that Max{B_1(r), B_2(r)} is minimized
      Return N*Max{B_1(r), B_2(r)}

   */

   /*
      We let hi(x)    = 1/x^(n+1) * (|a_{n+1}|*x^(n+1) + ... + |a_N|*x^N)
      so that hi(r)   = |a_{n+1}| + ... + |a_N|*r^(N-n-1) = B_2(r)

      and let lo(x)   = |a_n|*x + ... + |a_0|*x^(n+1)
      so that lo(1/r) = |a_n|/r + ... + |a_0|/r^{n+1}     = B_1(r)
   */

   /*
      We notionally begin by letting lo_eval = B_1(0) which has the
      value +INF. And we let hi_eval = B_2(0) which has the value
      |a_{n+1}|.

      We see that lo_eval is too large, so we set step = 1 and add
      step to the power of the evaluation point for B_1. We keep adding
      step until lo_eval is not too large. We then halve step and switch
      its direction. We keep doing this until we have refined the value
      sufficiently.
   */

   fmpz_poly_t lo, hi;
   double rpow = 0, hi_eval, lo_eval, max_eval;
   double step = 1.0;
   slong size_f = FLINT_ABS(fmpz_poly_max_bits(f));
   slong hi_exp = 0, lo_exp = 0;
   slong rbits, max_exp, rexp = 0;
   int too_much = 0;

   fmpz_poly_init(lo);
   fmpz_poly_init(hi);

   /* compute lo(x) and hi(x) */
   fmpz_poly_set_trunc(lo, f, n + 1);
   fmpz_poly_scalar_abs(lo, lo);
   fmpz_poly_reverse(lo, lo, n + 2);

   fmpz_poly_shift_right(hi, f, n + 1);
   fmpz_poly_scalar_abs(hi, hi);

   /* refine guess */
   while (1)
   {
      /* compute approx max_exp that evaluation will use */
      rbits = fabs(rpow);
      max_exp = rbits*hi->length + size_f + 1;

      if (rbits > 200 || rexp != 0) /* r itself has a large exponent */
      {
         /* r is really 2^rpow * 2^rexp */
         double r = pow(2.0, rpow);
         hi_eval = fmpz_poly_evaluate_horner_d_2exp2(&hi_exp, hi, r, rexp);
         lo_eval = fmpz_poly_evaluate_horner_d_2exp2(&lo_exp, lo, 1/r, -rexp);
      }  /* if max exponent may overwhelm a double (with safety margin) */
      else if (max_exp > 950 || too_much) /* result of eval has large exponent */
      {
         double r = pow(2.0, rpow);
         hi_eval = fmpz_poly_evaluate_horner_d_2exp(&hi_exp, hi, r);
         lo_eval = fmpz_poly_evaluate_horner_d_2exp(&lo_exp, lo, 1/r);
      }
      else /* everything can be handled with doubles */
      {
         double r = pow(2.0, rpow);
         hi_eval = fmpz_poly_evaluate_horner_d(hi, r);
         lo_eval = fmpz_poly_evaluate_horner_d(lo, 1/r);
         hi_exp = lo_exp = 0;
      }

      /* if doubles suffice */
      if (hi_exp == 0 && lo_exp == 0)
      {
         if (1.5*lo_eval < hi_eval)
         {
            /* lo_eval is too small */
            if (step >= 0.0)
               step = -step/2.0;

            rpow += step;
         } else if (lo_eval > 1.5*hi_eval)
         {
            /* lo_eval is too big */
            if (step < 0.0)
               step = -step/2.0;

            rpow += step;
         } else if (hi_eval != hi_eval || lo_eval != lo_eval)
         {
            /* doubles were insufficient after all */
            too_much = 1; 
         } else /* we are done */
         {
            if (hi_eval > lo_eval)
               max_eval = hi_eval;
            else
               max_eval = lo_eval;

            /* we adjust due to rounding errors in horner's evaluation */
            fmpz_set_d(res, max_eval*(f->length - 1)*(1.0 + f->length*ldexp(2.0, -51)));

            goto cleanup;
         }
      } else
      {
         /* 
            too big for doubles alone, so we represent in d*2^exp format

            _d_cmp_2exp will return 2 when 2*lo < hi, -2 when 2*hi < lo
            1 when lo <= hi and -1 when lo > hi
         */
         int cmp;

         /* just compare exponents, because it's probably not going to converge otherwise */
         if (rbits > 200 || rexp != 0)
         {
             if ((hi_exp + _log2(hi_eval)) > 1.01 * (lo_exp + _log2(lo_eval)))
                cmp = 2;
             else if ((lo_exp + _log2(lo_eval)) > 1.01 * (hi_exp + _log2(hi_eval)))
                cmp = -2;
             else if ((hi_exp + _log2(hi_eval)) >= (lo_exp + _log2(lo_eval)))
                cmp = 1;
             else
                cmp = -1;
         }
         else
         {
             cmp = _d_cmp_2exp(hi_eval, hi_exp, lo_eval, lo_exp);
         }

         if (fabs(step) < CLD_EPS) /* give up if step is now very small */
         {
            cmp = cmp == 2 ? 1 : -1;
         }

         switch (cmp)
         {
         case 2:
            if (step >= 0.0)
               step = -step/2.0;
            rpow += step;
            break;
         case -2:
            if (step < 0.0)
               step = -step/2.0;
            rpow += step;
            break;
         case 1:
            /* we adjust due to rounding errors in horner's evaluation */
            fmpz_set_d_2exp(res, hi_eval*(f->length - 1)*(1.0 + f->length*ldexp(2.0, -51)), hi_exp);
            goto cleanup;
         case -1:
            /* we adjust due to rounding errors in horner's evaluation */
            fmpz_set_d_2exp(res, lo_eval*(f->length - 1)*(1.0 + f->length*ldexp(2.0, -51)), lo_exp);
            goto cleanup;
         }
      }

      if (rpow > 1000.0)
      {
         rpow -= 1000.0;
         rexp += 1000;
      } else if (rpow < -1000.0)
      {
         rpow += 1000.0;
         rexp -= 1000;
      }
   }

cleanup:

   fmpz_poly_clear(lo);
   fmpz_poly_clear(hi);
}
