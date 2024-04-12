/* mpn_get_d -- limbs to double conversion.

   THE FUNCTIONS IN THIS FILE ARE FOR INTERNAL USE ONLY.  THEY'RE ALMOST
   CERTAIN TO BE SUBJECT TO INCOMPATIBLE CHANGES OR DISAPPEAR COMPLETELY IN
   FUTURE GNU MP RELEASES.

Copyright 2003, 2004 Free Software Foundation, Inc.
Copyright 2021 William Hart (adaption to Flint)

This file is part of the GNU MP Library.

The GNU MP Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 2.1 of the License, or (at your
option) any later version.

The GNU MP Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the GNU MP Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
MA 02110-1301, USA. */

#define CONST_1024	      (1024)
#define CONST_NEG_1023	      (-1023)
#define CONST_NEG_1022_SUB_53 (-1022 - 53)

/* Return the value {ptr,size}*2^exp, and negative if sign<0.
   Must have size>=1, and a non-zero high limb ptr[size-1].

   {ptr,size} is truncated towards zero.  This is consistent with other gmp
   conversions, like mpz_set_f or mpz_set_q, and is easy to implement and
   test.

   In the past conversions had attempted (imperfectly) to let the hardware
   float rounding mode take effect, but that gets tricky since multiple
   roundings need to be avoided, or taken into account, and denorms mean the
   effective precision of the mantissa is not constant.	 (For reference,
   mpz_get_d on IEEE systems was ok, except it operated on the absolute
   value.  mpf_get_d and mpq_get_d suffered from multiple roundings and from
   not always using enough bits to get the rounding right.)

   It's felt that GMP is not primarily concerned with hardware floats, and
   really isn't enhanced by getting involved with hardware rounding modes
   (which could even be some weird unknown style), so something unambiguous
   and straightforward is best.


   The IEEE code below is the usual case, it knows either a 32-bit or 64-bit
   limb and is done with shifts and masks.  The 64-bit case in particular
   should come out nice and compact.

   The generic code works one bit at a time, which will be quite slow, but
   should support any binary-based "double" and be safe against any rounding
   mode.  Note in particular it works on IEEE systems too.


   Traps:

   Hardware traps for overflow to infinity, underflow to zero, or
   unsupported denorms may or may not be taken.	 The IEEE code works bitwise
   and so probably won't trigger them, the generic code works by float
   operations and so probably will.  This difference might be thought less
   than ideal, but again its felt straightforward code is better than trying
   to get intimate with hardware exceptions (of perhaps unknown nature).


   Not done:

   mpz_get_d in the past handled size==1 with a cast limb->double.  This
   might still be worthwhile there (for up to the mantissa many bits), but
   for mpn_get_d here, the cost of applying "exp" to the resulting exponent
   would probably use up any benefit a cast may have over bit twiddling.
   Also, if the exponent is pushed into denorm range then bit twiddling is
   the only option, to ensure the desired truncation is obtained.


   Other:

   For reference, note that HPPA 8000, 8200, 8500 and 8600 trap FCNV,UDW,DBL
   to the kernel for values >= 2^63.  This makes it slow, and worse the
   Linux kernel (what versions?) apparently uses untested code in its trap
   handling routines, and gets the sign wrong.  We don't use such a limb to
   double cast, neither in the IEEE or generic code.  */

#include "mpn_extras.h"

/* assumes unsigned int is at least 32 bits */
#if defined (FLINT_BIG_ENDIAN) && FLINT_BIG_ENDIAN == 1
union ieee_double_extract
{
  struct
    {
      unsigned int sig:1;
      unsigned int exp:11;
      unsigned int manh:20;
      unsigned int manl:32;
    } s;
  double d;
};
#else
union ieee_double_extract
{
  struct
    {
      unsigned int manl:32;
      unsigned int manh:20;
      unsigned int exp:11;
      unsigned int sig:1;
    } s;
  double d;
};
#endif

double flint_mpn_get_d(mp_srcptr ptr, mp_size_t size, mp_size_t sign, long exp)
{
  FLINT_ASSERT (size >= 0);
  FLINT_ASSERT (size == 0 || ptr[size-1] != 0);

  if (size == 0)
    return 0.0;

  /* Adjust exp to a radix point just above {ptr,size}, guarding against
     overflow.	After this exp can of course be reduced to anywhere within
     the {ptr,size} region without underflow.  */
  if ((unsigned long) (FLINT_BITS * size)
		> (unsigned long) (LONG_MAX - exp))
    {
	goto ieee_infinity;

      /* generic */
      exp = LONG_MAX;
    }
  else
    {
      exp += FLINT_BITS * size;
    }

#define ONE_LIMB    (FLINT_BITS == 64)
#define TWO_LIMBS   (FLINT_BITS == 32)

    {
      union ieee_double_extract	 u;
      mp_limb_t	 m0, m1, m2, rmask;
      int	 lshift, rshift;

      m0 = ptr[size-1];			    /* high limb */
      m1 = (size >= 2 ? ptr[size-2] : 0);   /* second highest limb */
      lshift = flint_clz(m0);

      /* relative to just under high non-zero bit */
      exp -= lshift + 1;

      if (ONE_LIMB)
	{
	  /* lshift to have high of m0 non-zero, and collapse nails */
	  rshift = FLINT_BITS - lshift;
	  rmask = lshift == 0 ? 0 : UWORD_MAX;
	  m0 = (m0 << lshift) | ((m1 >> rshift) & rmask);

	  /* rshift back to have bit 53 of m0 the high non-zero */
	  m0 >>= 11;
	}
      else /* TWO_LIMBS */
	{
	  m2 = (size >= 3 ? ptr[size-3] : 0);  /* third highest limb */

	  /* lshift to have high of m0:m1 non-zero, collapse nails from m0 */
	  rshift = FLINT_BITS - lshift;
	  rmask = (lshift == 0 ? 0 : UWORD_MAX);
	  m0 = (m0 << lshift) | ((m1 >> rshift) & rmask);
	  m1 = (m1 << lshift) | ((m2 >> rshift) & rmask);

	  /* rshift back to have bit 53 of m0:m1 the high non-zero */
	  m1 = (m1 >> 11) | (m0 << (FLINT_BITS-11));
	  m0 >>= 11;
	}

      if (exp >= CONST_1024)
	{
	  /* overflow, return infinity */
	ieee_infinity:
	  m0 = 0;
	  m1 = 0;
	  exp = 1024;
	}
      else if (exp <= CONST_NEG_1023)
	{
	  if (exp <= CONST_NEG_1022_SUB_53)
	    return 0.0;	 /* denorm underflows to zero */

	  rshift = -1022 - exp;
	  FLINT_ASSERT (rshift > 0 && rshift < 53);
	  if (ONE_LIMB)
	    {
	      m0 >>= rshift;
	    }
	  else /* TWO_LIMBS */
	    {
	      if (rshift >= 32)
		{
		  m1 = m0;
		  m0 = 0;
		  rshift -= 32;
		}
	      lshift = FLINT_BITS - rshift;
	      m1 = (m1 >> rshift) | (rshift == 0 ? 0 : m0 << lshift);
	      m0 >>= rshift;
	    }
	  exp = -1023;
	}

      if (ONE_LIMB)
	{
#if FLINT_BITS > 32	/* avoid compiler warning about big shift */
	  u.s.manh = m0 >> 32;
#endif
	  u.s.manl = m0;
	}
      else /* TWO_LIMBS */
	{
	  u.s.manh = m0;
	  u.s.manl = m1;
	}

      u.s.exp = exp + 1023;
      u.s.sig = (sign < 0);
      return u.d;
    }
}
