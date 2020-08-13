/*
    Copyright (C) 2017 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <stdlib.h>
#include "flint.h"
#include "mpoly.h"

int mpoly_monomial_exists1(slong * index, const ulong * poly_exps,
                                      const ulong exp, slong len, ulong maskhi)
{
   slong n = len;
   slong i = 0;

   if ((exp^maskhi) > (poly_exps[0]^maskhi)) /* greater than first term */
   {
      (*index) = 0;
      return 0;
   }

   while (n > 1) /* do binary search */
   {
      slong half = n/2;

      /* if in first half */
      if ((exp^maskhi) > (poly_exps[i + half]^maskhi))
         n = half;
      else /* in second half */
      {
         n -= half;
         i += half;
      }
   }

   /* if equal to term at index i */
   if (exp == poly_exps[i])
   {
      (*index) = i;
      return 1;
   } else /* less than term at index i, but doesn't exist */
   {
      (*index) = i + 1;
      return 0;
   } 
}

int mpoly_monomial_exists(slong * index, const ulong * poly_exps,
                  const ulong * exp, slong len, slong N, const ulong * cmpmask)
{
   slong n = len;
   slong i = 0;

   if (len == 0) /* no terms to search */
   {
      (*index) = 0;
      return 0;
   }

   /* specialised version if exponent vectors are one word */
   if (N == 1)
      return mpoly_monomial_exists1(index, poly_exps, *exp, len, cmpmask[0]);

   if (mpoly_monomial_gt(exp, poly_exps, N, cmpmask)) /* greater than first term */
   {
      (*index) = 0;
      return 0;
   }

   while (n > 1) /* do binary search */
   {
      slong half = n/2;

      /* if in first half */
      if (mpoly_monomial_gt(exp, poly_exps + (i + half)*N, N, cmpmask))
         n = half;
      else /* in second half */
      {
         n -= half;
         i += half;
      }
   }

   /* if equal to term at index i */
   if (mpoly_monomial_equal(exp, poly_exps + i*N, N))
   {
      (*index) = i;
      return 1;
   } else /* less than term at index i, but doesn't exist */
   {
      (*index) = i + 1;
      return 0;
   } 
}
