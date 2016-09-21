/*
    Copyright (C) 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"

#if FLINT64

void _fmpz_mpoly_set_monomial_8_64(ulong * exp1, const ulong * exp2, 
                                       ulong degree, slong n, int deg, int rev)
{
   slong q, r, j, k1 = 0, k2 = 0;
   ulong v;

   q = n/8;
   r = n%8;

   if (rev)
   {
      for (j = 0; j < q; j++)
      {
         if (deg && j == 0)
            v = (degree << 56);
         else
            v = (exp2[k2++] << 56);
         v += (exp2[k2++] << 48);
         v += (exp2[k2++] << 40);
         v += (exp2[k2++] << 32);
         v += (exp2[k2++] << 24);
         v += (exp2[k2++] << 16);
         v += (exp2[k2++] << 8);
         v += (exp2[k2++]);
         exp1[k1++] = v;
      }    

      if (r > 0)
      {
         if (deg && q == 0)
            v = (degree << 56);
         else
            v = (exp2[k2++] << 56);

         for (j = 1; j < r; j++)
            v += (exp2[k2++] << (56 - 8*j));

         exp1[k1++] = v;
      }
   } else
   {
      k2 = q;
      if (r > 0)
      {
         if (deg && q == 0)
            v = (degree << 56);
         else
            v = (exp2[k2--] << 56);

         for (j = 1; j < r; j++)
            v += (exp2[k2--] << (56 - 8*j));

         exp1[k1++] = v;
      }

      for (j = q - 1; j >= 0; j--)
      {
         if (deg && j == 0)
            v = (degree << 56);
         else
            v = (exp2[k2--] << 56);
         v += (exp2[k2--] << 48);
         v += (exp2[k2--] << 40);
         v += (exp2[k2--] << 32);
         v += (exp2[k2--] << 24);
         v += (exp2[k2--] << 16);
         v += (exp2[k2--] << 8);
         v += (exp2[k2--]);

         exp1[k1++] = v;
      }
   }
}

void _fmpz_mpoly_set_monomial_16_64(ulong * exp1, const ulong * exp2, 
                                       ulong degree, slong n, int deg, int rev)
{
   slong q, r, j, k1 = 0, k2 = 0;
   ulong v;

   q = n/4;
   r = n%4;

   if (rev)
   {
      for (j = 0; j < q; j++)
      {
         if (deg && j == 0)
            v = (degree << 48);
         else
            v = (exp2[k2++] << 48);
         v += (exp2[k2++] << 32);
         v += (exp2[k2++] << 16);
         v += (exp2[k2++]);
         exp1[k1++] = v;
      }    

      if (r > 0)
      {
         if (deg && q == 0)
            v = (degree << 48);
         else
            v = (exp2[k2++] << 48);

         for (j = 1; j < r; j++)
            v += (exp2[k2++] << (48 - 16*j));

         exp1[k1++] = v;
      }
   } else
   {
      k2 = q;
      if (r > 0)
      {
         if (deg && q == 0)
            v = (degree << 48);
         else
            v = (exp2[k2--] << 48);

         for (j = 1; j < r; j++)
            v += (exp2[k2--] << (48 - 16*j));

         exp1[k1++] = v;
      }

      for (j = q - 1; j >= 0; j--)
      {
         if (deg && j == 0)
            v = (degree << 48);
         else
            v = (exp2[k2--] << 48);
         v += (exp2[k2--] << 32);
         v += (exp2[k2--] << 16);
         v += (exp2[k2--]);

         exp1[k1++] = v;
      }
   }
}

void _fmpz_mpoly_set_monomial_32_64(ulong * exp1, const ulong * exp2,
                                       ulong degree, slong n, int deg, int rev)
{
   slong q, r, j, k1 = 0, k2 = 0;
   ulong v;

   q = n/2;
   r = n%2;

   if (rev)
   {
      for (j = 0; j < q; j++)
      {
         if (deg && j == 0)
            v = (degree << 32);
         else
            v = (exp2[k2++] << 32);
         v += (exp2[k2++]);
         exp1[k1++] = v;
      }    

      if (r != 0)
      {
         if (deg && q == 0)
            exp1[k1++] = (degree << 32);
         else
            v = (exp2[k2++] << 32);
      }
   } else
   {
      k2 = q;
      if (r != 0)
      {
         if (deg && q == 0)
            exp1[k1++] = (degree << 32);
         else
            exp1[k1++] = (exp2[k2--] << 32);
      }

      for (j = q - 1; j >= 0; j--)
      {
         if (deg && j == 0)
            v = (degree << 32);
         else
            v = (exp2[k2--] << 32);
         v += (exp2[k2--]);

         exp1[k1++] = v;
      }
   }
}

void _fmpz_mpoly_set_monomial_64_64(ulong * exp1, const ulong * exp2,
                                       ulong degree, slong n, int deg, int rev)
{
   slong k1 = 0;
   
   if (rev)
   {
      if (deg)
         exp1[k1++] = degree;

      for ( ; k1 < n; k1++)
         exp1[k1] = exp2[k1 - deg];
   } else
   {
      if (deg)
         exp1[k1++] = degree;

      for ( ; k1 < n; k1++)
         exp1[k1] = exp2[n - k1 - deg - 1];      
   }
}

#else

void _fmpz_mpoly_set_monomial_8_32(ulong * exp1, const ulong * exp2,
                                       ulong degree, slong n, int deg, int rev)
{
   slong q, r, j, k1 = 0, k2 = 0;
   ulong v;

   q = n/4;
   r = n%4;

   if (rev)
   {
      for (j = 0; j < q; j++)
      {
         if (deg && j == 0)
            v = (degree << 24);
         else
            v = (exp2[k2++] << 24);
         v += (exp2[k2++] << 16);
         v += (exp2[k2++] << 8);
         v += (exp2[k2++]);
         exp1[k1++] = v;
      }    

      if (r > 0)
      {
         if (deg && q == 0)
            v = (degree << 24);
         else
            v = (exp2[k2++] << 24);

         for (j = 1; j < r; j++)
            v += (exp2[k2++] << (24 - 8*j));

         exp1[k1++] = v;
      }
   } else
   {
      k2 = q;
      if (r > 0)
      {
         if (deg && q == 0)
            v = (degree << 24);
         else
            v = (exp2[k2--] << 24);

         for (j = 1; j < r; j++)
            v += (exp2[k2--] << (24 - 8*j));

         exp1[k1++] = v;
      }

      for (j = q - 1; j >= 0; j--)
      {
         if (deg && j == 0)
            v = (degree << 24);
         else
            v = (exp2[k2--] << 24);
         v += (exp2[k2--] << 16);
         v += (exp2[k2--] << 8);
         v += (exp2[k2--]);

         exp1[k1++] = v;
      }
   }
}

void _fmpz_mpoly_set_monomial_16_32(ulong * exp1, const ulong * exp2,
                                       ulong degree, slong n, int deg, int rev)
{
   slong q, r, j, k1 = 0, k2 = 0;
   ulong v;

   q = n/2;
   r = n%2;

   if (rev)
   {
      for (j = 0; j < q; j++)
      {
         if (deg && j == 0)
            v = (degree << 16);
         else
            v = (exp2[k2++] << 16);
         v += (exp2[k2++]);
         exp1[k1++] = v;
      }    

      if (r != 0)
      {
         if (deg && q == 0)
            exp1[k1++] = (degree << 16);
         else
            v = (exp2[k2++] << 16);
      }
   } else
   {
      k2 = q;
      if (r != 0)
      {
         if (deg && q == 0)
            exp1[k1++] = (degree << 16);
         else
            exp1[k1++] = (exp2[k2--] << 16);
      }

      for (j = q - 1; j >= 0; j--)
      {
         if (deg && j == 0)
            v = (degree << 16);
         else
            v = (exp2[k2--] << 16);
         v += (exp2[k2--]);

         exp1[k1++] = v;
      }
   }
}

void _fmpz_mpoly_set_monomial_32_32(ulong * exp1, const ulong * exp2,
                                       ulong degree, slong n, int deg, int rev)
{
   slong k1 = 0;
   
   if (rev)
   {
      if (deg)
         exp1[k1++] = degree;

      for ( ; k1 < n; k1++)
         exp1[k1] = exp2[k1 - deg];
   } else
   {
      if (deg)
         exp1[k1++] = degree;

      for ( ; k1 < n; k1++)
         exp1[k1] = exp2[n - k1 - deg - 1];      
   }
}

#endif

void _fmpz_mpoly_set_monomial(ulong * exp1, const ulong * exp2,
                                         slong bits, slong n, int deg, int rev)
{
   slong i;
   ulong degree = 0;

   if (deg)
   {
      for (i = 0; i < n - deg; i++)
         degree += exp2[i];
   }
   
#if FLINT64
   switch (bits)
   {
   case 8:
      _fmpz_mpoly_set_monomial_8_64(exp1, exp2, degree, n, deg, rev);
   break;
   case 16:
      _fmpz_mpoly_set_monomial_16_64(exp1, exp2, degree, n, deg, rev);
   break;
   case 32:
      _fmpz_mpoly_set_monomial_32_64(exp1, exp2, degree, n, deg, rev);
   break;
   case 64:
      _fmpz_mpoly_set_monomial_64_64(exp1, exp2, degree, n, deg, rev);
   break;
   }
#else
   switch (poly->bits)
   {
   case 8:
      _fmpz_mpoly_set_monomial_8_32(exp1, exp2, degree, n, deg, rev);
   break;
   case 16:
      _fmpz_mpoly_set_monomial_16_32(exp1, exp2, degree, n, deg, rev);
   break;
   case 32:
      _fmpz_mpoly_set_monomial_32_32(exp1, exp2, degree, n, deg, rev);
   break;
   }
#endif
}

void fmpz_mpoly_set_monomial(fmpz_mpoly_t poly, 
                       slong n, const ulong * exps, const fmpz_mpoly_ctx_t ctx)
{
   slong i, bits, max_bits, N;
   ulong maxdeg = 0;
   ulong * ptr;
   int deg, rev;

   degrev_from_ord(deg, rev, ctx->ord);

   if (deg)
   {
      for (i = 0; i < ctx->n - 1; i++)
         maxdeg += exps[i];
   } else
   {
      for (i = 0; i < ctx->n; i++)
      {
         if (exps[i] > maxdeg)
            maxdeg = exps[i];
      }
   }

   max_bits = poly->bits;
   bits = FLINT_BIT_COUNT(maxdeg);
   
   while (bits >= max_bits)
      max_bits *= 2;

   ptr = _fmpz_mpoly_unpack_monomials(max_bits, poly->exps, 
                                             poly->bits, ctx->n, poly->length);

   if (ptr != poly->exps)
   {
      flint_free(poly->exps);   
      poly->exps = ptr;
      poly->bits = max_bits;
   }

   fmpz_mpoly_fit_length(poly, n + 1, ctx);

   N = (poly->bits*ctx->n - 1)/FLINT_BITS + 1;
   
   _fmpz_mpoly_set_monomial(poly->exps + n*N, exps, poly->bits,
                                                             ctx->n, deg, rev);

   if (n + 1 > poly->length)
      _fmpz_mpoly_set_length(poly, n + 1, ctx);
}
