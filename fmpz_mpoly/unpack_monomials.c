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

void _fmpz_mpoly_unpack_monomials_8to64(ulong * exps1, const ulong * exps2, 
                                                            slong n, slong len)
{
   slong i, j, q, r, k1 = 0, k2 = 0;
   ulong v;

   q = n/8;
   r = n%8;

   for (i = 0; i < len; i++)
   {
      for (j = 0; j < q; j++)
      {
         v = exps2[k2++];
         exps1[k1++] = (v >> 56);
         exps1[k1++] = ((v >> 48) & UWORD(0xFF));
         exps1[k1++] = ((v >> 40) & UWORD(0xFF));
         exps1[k1++] = ((v >> 32) & UWORD(0xFF));
         exps1[k1++] = ((v >> 24) & UWORD(0xFF));
         exps1[k1++] = ((v >> 16) & UWORD(0xFF));
         exps1[k1++] = ((v >> 8) & UWORD(0xFF));
         exps1[k1++] = (v & UWORD(0xFF));
      }
      
      if (r > 0)
      {
         v = exps2[k2++];

         for (j = 0; j < r; j++)
            exps1[k1++] = ((v >> (56 - 8*j)) & UWORD(0xFF));
      }
   }
}

void _fmpz_mpoly_unpack_monomials_16to64(ulong * exps1, const ulong * exps2, 
                                                            slong n, slong len)
{
   slong i, j, q, r, k1 = 0, k2 = 0;
   ulong v;

   q = n/4;
   r = n%4;

   for (i = 0; i < len; i++)
   {
      for (j = 0; j < q; j++)
      {
         v = exps2[k2++];
         exps1[k1++] = (v >> 48);
         exps1[k1++] = ((v >> 32) & UWORD(0xFFFF));
         exps1[k1++] = ((v >> 16) & UWORD(0xFFFF));
         exps1[k1++] = (v & UWORD(0xFFFF));
      }
      
      if (r > 0)
      {
         v = exps2[k2++];

         for (j = 0; j < r; j++)
            exps1[k1++] = ((v >> (48 - 16*j)) & UWORD(0xFFFF));
      }
   }
}

void _fmpz_mpoly_unpack_monomials_32to64(ulong * exps1, const ulong * exps2, 
                                                            slong n, slong len)
{
   slong i, j, q, r, k1 = 0, k2 = 0;
   ulong v;

   q = n/2;
   r = n%2;

   for (i = 0; i < len; i++)
   {
      for (j = 0; j < q; j++)
      {
         v = exps2[k2++];
         exps1[k1++] = (v >> 32);
         exps1[k1++] = (v & UWORD(0xFFFFFFFF));
      }
      
      if (r != 0)
         exps1[k1++] = (exps2[k2++] >> 32);
   }
}

void _fmpz_mpoly_unpack_monomials_8to32(ulong * exps1, const ulong * exps2, 
                                                            slong n, slong len)
{
   slong i, j, q, r, r2, k1 = 0, k2 = 0;
   ulong v;

   q = n/8;
   r = n%8;

   for (i = 0; i < len; i++)
   {
      for (j = 0; j < q; j++)
      {
         v = exps2[k2++];
         exps1[k1++] = ((v >> 56) << 32) + ((v >> 48) & UWORD(0xFF));
         exps1[k1++] = ((v >> 8) & UWORD(0xFF00000000)) + ((v >> 32) & UWORD(0xFF));
         exps1[k1++] = ((v << 8) & UWORD(0xFF00000000)) + ((v >> 16) & UWORD(0xFF));
         exps1[k1++] = ((v << 24) & UWORD(0xFF00000000)) + (v & UWORD(0xFF));
      }
      
      r2 = r;

      if (r2 > 0)
      {
         v = exps2[k2++];
      
         if (r2 >= 4)
         {
            exps1[k1++] = ((v >> 56) << 32) + ((v >> 48) & UWORD(0xFF));
            exps1[k1++] = ((v >> 8) & UWORD(0xFF00000000)) + ((v >> 32) & UWORD(0xFF));
            r2 -= 4;
            v <<= 32;
         }

         for (j = 0; j + 1 < r2; j += 2)
            exps1[k1++] = ((v >> (24 - 8*j)) & UWORD(0xFF00000000)) + ((v >> (56 - 8*j)) & UWORD(0xFF));

         if (r2 & 1)
            exps1[k1++] = ((v >> (24 - 8*j)) & UWORD(0xFF00000000));
      }
   }
}

void _fmpz_mpoly_unpack_monomials_16to32(ulong * exps1, const ulong * exps2, 
                                                            slong n, slong len)
{
   slong i, j, q, r, k1 = 0, k2 = 0;
   ulong v;

   q = n/4;
   r = n%4;

   for (i = 0; i < len; i++)
   {
      for (j = 0; j < q; j++)
      {
         v = exps2[k2++];
         exps1[k1++] = ((v >> 48) << 32) + ((v >> 32) & UWORD(0xFFFF));
         exps1[k1++] = ((v << 16) & UWORD(0xFFFF00000000)) + (v & UWORD(0xFFFF));
      }
      
      if (r > 0)
      {
         v = exps2[k2++];
      
         if (r >= 2)
         {
            exps1[k1++] = ((v >> 48) << 32) + ((v >> 32) & UWORD(0xFFFF));
            v <<= 32;
         }

         if (r & 1)
            exps1[k1++] = (v >> 16);
      }
   }
}

void _fmpz_mpoly_unpack_monomials_8to16(ulong * exps1, const ulong * exps2, 
                                                            slong n, slong len)
{
   slong i, j, q, r, r2, k1 = 0, k2 = 0;
   ulong v;

   q = n/8;
   r = n%8;

   for (i = 0; i < len; i++)
   {
      for (j = 0; j < q; j++)
      {
         v = exps2[k2++];
         exps1[k1++] = ((v >> 56) << 48) + ((v >> 16) & UWORD(0xFF00000000))
                    + ((v >> 24) & UWORD(0xFF0000)) + ((v >> 32) & UWORD(0xFF));
         exps1[k1++] = ((v << 24) & UWORD(0xFF000000000000)) + ((v << 16) & UWORD(0xFF00000000))
                    + ((v << 8) & UWORD(0xFF0000)) + (v & UWORD(0xFF));
      }
      
      r2 = r;

      if (r2 > 0)
      {
         v = exps2[k2++];
      
         if (r2 >= 4)
         {
            exps1[k1++] = ((v >> 56) << 48) + ((v >> 16) & UWORD(0xFF00000000))
                    + ((v >> 24) & UWORD(0xFF0000)) + ((v >> 32) & UWORD(0xFF));
            v <<= 32;
            r2 -= 4;
         }

         if (r2 > 0)
         {
            exps1[k1] = ((v >> 56) << 48);

            if (r2 > 1)
            {
               exps1[k1] += ((v >> 16) & UWORD(0xFF00000000));

               if (r2 > 2)
                  exps1[k1] += ((v >> 24) & UWORD(0xFF0000));
            }
            
            k1++;
         }
      }
   }
}

#else

void _fmpz_mpoly_unpack_monomials_8to32(ulong * exps1, const ulong * exps2, 
                                                            slong n, slong len)
{
   slong i, j, q, r, k1 = 0, k2 = 0;
   ulong v;

   q = n/4;
   r = n%4;

   for (i = 0; i < len; i++)
   {
      for (j = 0; j < q; j++)
      {
         v = exps2[k2++];
         exps1[k1++] = (v >> 24);
         exps1[k1++] = ((v >> 16) & UWORD(0xFF));
         exps1[k1++] = ((v >> 8) & UWORD(0xFF));
         exps1[k1++] = (v & UWORD(0xFF));
      }
      
      if (r > 0)
      {
         v = exps2[k2++];

         for (j = 0; j < r; j++)
            exps1[k1++] = ((v >> (24 - 8*j)) & UWORD(0xFF));
      }
   }
}

void _fmpz_mpoly_unpack_monomials_16to32(ulong * exps1, const ulong * exps2, 
                                                            slong n, slong len)
{
   slong i, j, q, r, k1 = 0, k2 = 0;
   ulong v;

   q = n/2;
   r = n%2;

   for (i = 0; i < len; i++)
   {
      for (j = 0; j < q; j++)
      {
         v = exps2[k2++];
         exps1[k1++] = ((v >> 16);
         exps1[k1++] = (v & UWORD(0xFFFF));
      }
      
      if (r != 0)
         exps1[k1++] = (exps2[k2++] >> 16);
   }
}

void _fmpz_mpoly_unpack_monomials_8to16(ulong * exps1, const ulong * exps2, 
                                                            slong n, slong len)
{
   slong i, j, q, r, r2, k1 = 0, k2 = 0;
   ulong v;

   q = n/4;
   r = n%4;

   for (i = 0; i < len; i++)
   {
      for (j = 0; j < q; j++)
      {
         v = exps2[k2++];
         exps1[k1++] = ((v >> 24) << 16) + ((v >> 16) & UWORD(0xFF));
         exps1[k1++] = ((v << 8) & UWORD(0xFF0000)) + (v & UWORD(0xFF));
      }
      
      r2 = r;

      if (r2 > 0)
      {
         v = exps2[k2++];
      
         if (r2 >= 2)
         {
            exps1[k1++] = ((v >> 24) << 16) + ((v >> 16) & UWORD(0xFF0000));

            r2 -= 2;
            v <<= 16;
         }

         if (r2 > 0)
            exps1[k1++] = (v >> 8);
      }
   }
}

#endif

ulong * _fmpz_mpoly_unpack_monomials(slong bits1, const ulong * exps2, 
                                               slong bits2, slong n, slong len)
{
   ulong * exps1;
   slong N;

   if (bits1 == bits2)
      return (ulong *) exps2;

   N = (n*bits1 - 1)/FLINT_BITS + 1; /* number of words per exponent vector */

   exps1 = (ulong *) flint_malloc(N*len*sizeof(ulong));

#if FLINT64
   if (bits1 == 64)
   {
      if (bits2 == 8)
         _fmpz_mpoly_unpack_monomials_8to64(exps1, exps2, n, len);
      else if (bits2 == 16)
         _fmpz_mpoly_unpack_monomials_16to64(exps1, exps2, n, len);
      else /* bits2 == 32 */
         _fmpz_mpoly_unpack_monomials_32to64(exps1, exps2, n, len);
   } else 
#endif
   if (bits1 == 32)
   {
      if (bits2 == 8)
         _fmpz_mpoly_unpack_monomials_8to32(exps1, exps2, n, len);
      else /* bits2 == 16 */
         _fmpz_mpoly_unpack_monomials_16to32(exps1, exps2, n, len);
   } else  /* bits1 == 16, bits2 = 8 */
         _fmpz_mpoly_unpack_monomials_8to16(exps1, exps2, n, len);

   return exps1;
}
