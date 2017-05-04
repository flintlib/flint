/*
    Copyright (C) 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "mpoly.h"

#if FLINT64

void mpoly_get_monomial_8_64(ulong * exp1, const ulong * exp2,
                                                     slong n, int deg, int rev)
{
   slong q, r, j, k1 = 0, k2 = 0;
   ulong v;

   q = n/8;
   r = n%8;

   if (rev)
   {
      for (j = 0; j < q; j++)
      {
         v = exp2[k2++];
         if (!deg || j > 0)
            exp1[k1++] = (v >> 56);
         exp1[k1++] = ((v >> 48) & UWORD(0xFF));
         exp1[k1++] = ((v >> 40) & UWORD(0xFF));
         exp1[k1++] = ((v >> 32) & UWORD(0xFF));
         exp1[k1++] = ((v >> 24) & UWORD(0xFF));
         exp1[k1++] = ((v >> 16) & UWORD(0xFF));
         exp1[k1++] = ((v >> 8) & UWORD(0xFF));
         exp1[k1++] = (v & UWORD(0xFF));
      }    

      if (r > 0)
      {
         v = exp2[k2++];

         for (j = (deg && q == 0); j < r; j++)
            exp1[k1++] = ((v >> (56 - 8*j)) & UWORD(0xFF));
      }
   } else
   {
      k2 = q;
      if (r > 0)
      {
         v = exp2[k2--];

         for (j = r - 1; j >= (deg && q == 0); j--)
            exp1[k1++] = ((v >> (56 - 8*j)) & UWORD(0xFF));
      }

      for (j = q - 1; j >= 0; j--)
      {
         v = exp2[k2--];
         exp1[k1++] = (v & UWORD(0xFF));
         exp1[k1++] = ((v >> 8) & UWORD(0xFF));
         exp1[k1++] = ((v >> 16) & UWORD(0xFF));
         exp1[k1++] = ((v >> 24) & UWORD(0xFF));
         exp1[k1++] = ((v >> 32) & UWORD(0xFF));
         exp1[k1++] = ((v >> 40) & UWORD(0xFF));
         exp1[k1++] = ((v >> 48) & UWORD(0xFF));
         if (!deg || j > 0)
            exp1[k1++] = (v >> 56);
      }
   }
}

void mpoly_get_monomial_16_64(ulong * exp1, const ulong * exp2,
                                                     slong n, int deg, int rev)
{
   slong q, r, j, k1 = 0, k2 = 0;
   ulong v;

   q = n/4;
   r = n%4;

   if (rev)
   {
      for (j = 0; j < q; j++)
      {
         v = exp2[k2++];
         if (!deg || j > 0)
            exp1[k1++] = (v >> 48);
         exp1[k1++] = ((v >> 32) & UWORD(0xFFFF));
         exp1[k1++] = ((v >> 16) & UWORD(0xFFFF));
         exp1[k1++] = (v & UWORD(0xFFFF));
      }    

      if (r > 0)
      {
         v = exp2[k2++];

         for (j = (deg && q == 0); j < r; j++)
            exp1[k1++] = ((v >> (48 - 16*j)) & UWORD(0xFFFF));
      }
   } else
   {
      k2 = q;
      if (r > 0)
      {
         v = exp2[k2--];

         for (j = r - 1; j >= (deg && q == 0); j--)
            exp1[k1++] = ((v >> (48 - 16*j)) & UWORD(0xFFFF));
      }

      for (j = q - 1; j >= 0; j--)
      {
         v = exp2[k2--];
         exp1[k1++] = (v & UWORD(0xFFFF));
         exp1[k1++] = ((v >> 16) & UWORD(0xFFFF));
         exp1[k1++] = ((v >> 32) & UWORD(0xFFFF));
         if (!deg || j > 0)
            exp1[k1++] = (v >> 48);
      }
   }
}

void mpoly_get_monomial_32_64(ulong * exp1, const ulong * exp2,
                                                     slong n, int deg, int rev)
{
   slong q, r, j, k1 = 0, k2 = 0;
   ulong v;

   q = n/2;
   r = n%2;

   if (rev)
   {
      for (j = 0; j < q; j++)
      {
         v = exp2[k2++];
         if (!deg || j > 0)
            exp1[k1++] = (v >> 32);
         exp1[k1++] = (v & UWORD(0xFFFFFFFF));
      }    

      if (r != 0 && (!deg || q != 0))
         exp1[k1++] = (exp2[k2++] >> 32);
   } else
   {
      k2 = q;
      if (r > 0 && (!deg || q != 0))
         exp1[k1++] = (exp2[k2--] >> 32);

      for (j = q - 1; j >= 0; j--)
      {
         v = exp2[k2--];
         exp1[k1++] = (v & UWORD(0xFFFFFFFF));
         if (!deg || j > 0)
            exp1[k1++] = (v >> 32);
      }
   }
}

void mpoly_get_monomial_64_64(ulong * exp1, const ulong * exp2,
                                                     slong n, int deg, int rev)
{
   slong i;

   if (rev)
   {
      for (i = deg; i < n; i++)
         exp1[i - deg] = exp2[i];
   } else
   {
      for (i = n - 1; i >= deg; i--)
         exp1[n - i - 1] = exp2[i];
   }
}

#else

void mpoly_get_monomial_8_32(ulong * exp1, const ulong * exp2,
                                                     slong n, int deg, int rev)
{
   slong q, r, j, k1 = 0, k2 = 0;
   ulong v;

   q = n/4;
   r = n%4;

   if (rev)
   {
      for (j = 0; j < q; j++)
      {
         v = exp2[k2++];
         if (!deg || j > 0)
            exp1[k1++] = (v >> 24);
         exp1[k1++] = ((v >> 16) & UWORD(0xFF));
         exp1[k1++] = ((v >> 8) & UWORD(0xFF));
         exp1[k1++] = (v & UWORD(0xFF));
      }    

      if (r > 0)
      {
         v = exp2[k2++];

         for (j = (deg && q == 0); j < r; j++)
            exp1[k1++] = ((v >> (24 - 8*j)) & UWORD(0xFF));
      }
   } else
   {
      k2 = q;
      if (r > 0)
      {
         v = exp2[k2--];

         for (j = r - 1; j >= (deg && q == 0); j--)
            exp1[k1++] = ((v >> (24 - 8*j)) & UWORD(0xFF));
      }

      for (j = q - 1; j >= 0; j--)
      {
         v = exp2[k2--];
         exp1[k1++] = (v & UWORD(0xFF));
         exp1[k1++] = ((v >> 8) & UWORD(0xFF));
         exp1[k1++] = ((v >> 16) & UWORD(0xFF));
         if (!deg || j > 0)
            exp1[k1++] = (v >> 24);
      }
   }
}

void mpoly_get_monomial_16_32(ulong * exp1, const ulong * exp2,
                                                     slong n, int deg, int rev)
{
   slong q, r, j, k1 = 0, k2 = 0;
   ulong v;

   q = n/2;
   r = n%2;

   if (rev)
   {
      for (j = 0; j < q; j++)
      {
         v = exp2[k2++];
         if (!deg || j > 0)
            exp1[k1++] = (v >> 16);
         exp1[k1++] = (v & UWORD(0xFFFF));
      }    

      if (r != 0 && (!deg || q != 0))
         exp1[k1++] = (exp2[k2++] >> 16);
   } else
   {
      k2 = q;
      if (r > 0 && (!deg || q != 0))
         exp1[k1++] = (exp2[k2--] >> 16);

      for (j = q - 1; j >= 0; j--)
      {
         v = exp2[k2--];
         exp1[k1++] = (v & UWORD(0xFFFF));
         if (!deg || j > 0)
            exp1[k1++] = (v >> 16);
      }
   }
}

void mpoly_get_monomial_32_32(ulong * exp1, const ulong * exp2,
                                                     slong n, int deg, int rev)
{
   slong i;

   if (rev)
   {
      for (i = deg; i < n; i++)
         exp1[i - deg] = exp2[i];
   } else
   {
      for (i = n - 1; i >= deg; i--)
         exp1[n - i - 1] = exp2[i];
   }
}

#endif

void mpoly_get_monomial(ulong * exps, const ulong * poly_exps,
                                         slong bits, slong n, int deg, int rev)
{
#if FLINT64
   switch (bits)
   {
   case 8:
      mpoly_get_monomial_8_64(exps, poly_exps, n, deg, rev);
   break;
   case 16:
      mpoly_get_monomial_16_64(exps, poly_exps, n, deg, rev);
   break;
   case 32:
      mpoly_get_monomial_32_64(exps, poly_exps, n, deg, rev);
   break;
   case 64:
      mpoly_get_monomial_64_64(exps, poly_exps, n, deg, rev);
   break;
   }
#else
   switch (bits)
   {
   case 8:
      mpoly_get_monomial_8_32(exps, poly_exps, n, deg, rev);
   break;
   case 16:
      mpoly_get_monomial_16_32(exps, poly_exps, n, deg, rev);
   break;
   case 32:
      mpoly_get_monomial_32_32(exps, poly_exps, n, deg, rev);
   break;
   }
#endif
}

