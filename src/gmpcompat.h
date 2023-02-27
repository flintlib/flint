/*
    Copyright (C) 2013 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <limits.h>

#ifndef GMP_COMPAT_H
#define GMP_COMPAT_H

#define FLINT_MPZ_REALLOC(z, len)       \
    ((len) > ((z)->_mp_alloc)           \
        ? (mp_ptr) _mpz_realloc(z, len) \
        : ((z)->_mp_d))

#define FLINT_MPZ_PTR_SWAP(a, b)    \
  do {                              \
    mpz_ptr __tmp = (a);            \
    (a) = (b);                      \
    (b) = __tmp;                    \
  } while (0)

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

static __inline__
double flint_mpn_get_d (mp_srcptr ptr, mp_size_t size, mp_size_t sign, long exp)
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
      count_leading_zeros (lshift, m0);

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

#if WORD_MAX != LONG_MAX && !defined(__MPIR_VERSION)

#define FLINT_MOCK_MPZ_UI(xxx, yyy) \
   __mpz_struct (xxx)[1] = {{ 1, 0, NULL }}; \
   do { \
      (xxx)->_mp_d = (mp_ptr) &(yyy); \
      if ((yyy) > 0) \
         (xxx)->_mp_size = 1; \
   }while (0)

#define FLINT_MOCK_MPZ_SI(xxx, yyy) \
   __mpz_struct (xxx)[1] = {{ 1, 0, NULL }}; \
   do { \
      (xxx)->_mp_d = (mp_ptr) &(yyy); \
      if ((yyy) < 0) (xxx)->_mp_size = -1, (yyy) = -(yyy); \
      else (xxx)->_mp_size = (yyy) != 0; \
   } while (0)

#define flint_mpz_get_si(xxx) \
   ((xxx)->_mp_size == 0 ? (slong) WORD(0) : \
   ((xxx)->_mp_size > 0 ? (slong) (xxx)->_mp_d[0] : (slong) -(xxx)->_mp_d[0]))

#define flint_mpz_get_ui(xxx) \
   ((xxx)->_mp_size == 0 ? UWORD(0) : (xxx)->_mp_d[0])

static __inline__
void flint_mpz_set_si(mpz_ptr r, slong s)
{
   /* GMP 6.2 lazily performs allocation, deal with that if necessary
      (in older GMP versions, this code is simply never triggered) */
   if (r->_mp_alloc == 0)
   {
      r->_mp_d = (mp_ptr) flint_malloc(sizeof(mp_limb_t));
      r->_mp_alloc = 1;
   }

   if (s < 0) {
      r->_mp_size = -1;
      r->_mp_d[0] = -s;
   } else {
      r->_mp_size = s != 0;
      r->_mp_d[0] = s;
   }
}

static __inline__
void flint_mpz_set_ui(mpz_ptr r, ulong u)
{
   /* GMP 6.2 lazily performs allocation, deal with that if necessary
      (in older GMP versions, this code is simply never triggered) */
   if (r->_mp_alloc == 0)
   {
      r->_mp_d = (mp_ptr) flint_malloc(sizeof(mp_limb_t));
      r->_mp_alloc = 1;
   }

   r->_mp_d[0] = u; 
   r->_mp_size = u != 0;
}

static __inline__
void flint_mpz_init_set_si(mpz_ptr r, slong s)
{
   r->_mp_d = (mp_ptr) flint_malloc(sizeof(mp_limb_t));
   r->_mp_alloc = 1;

   if (s < 0) {
      r->_mp_size = -1;
      r->_mp_d[0] = -s;
   } else {
      r->_mp_size = s != 0;
      r->_mp_d[0] = s;
   }
}

static __inline__
void flint_mpz_init_set_ui(mpz_ptr r, ulong u)
{
   r->_mp_d = (mp_ptr) flint_malloc(sizeof(mp_limb_t));
   r->_mp_alloc = 1;

   r->_mp_d[0] = u; 
   r->_mp_size = u != 0;
}

static __inline__
void flint_mpz_add_ui(mpz_ptr a, mpz_srcptr b, ulong c)
{
   FLINT_MOCK_MPZ_UI(tc, c);
   
   mpz_add(a, b, tc);
}

static __inline__
void flint_mpz_sub_ui(mpz_ptr a, mpz_srcptr b, ulong c)
{
   FLINT_MOCK_MPZ_UI(tc, c);
   
   mpz_sub(a, b, tc);
}

static __inline__
void flint_mpz_mul_ui(mpz_ptr a, mpz_srcptr b, ulong c)
{
   FLINT_MOCK_MPZ_UI(tc, c);
   
   mpz_mul(a, b, tc);
}

static __inline__
void flint_mpz_mul_si(mpz_ptr a, mpz_srcptr b, slong c)
{
   FLINT_MOCK_MPZ_SI(tc, c);
   
   mpz_mul(a, b, tc);
}

static __inline__
void flint_mpz_addmul_ui(mpz_ptr a, mpz_srcptr b, ulong c)
{
   FLINT_MOCK_MPZ_UI(tc, c);
   
   mpz_addmul(a, b, tc);
}

static __inline__
void flint_mpz_submul_ui(mpz_ptr a, mpz_srcptr b, ulong c)
{
   FLINT_MOCK_MPZ_UI(tc, c);
   
   mpz_submul(a, b, tc);
}

static __inline__
void flint_mpz_ui_sub(mpz_ptr a, ulong b, mpz_srcptr c)
{
   FLINT_MOCK_MPZ_UI(tb, b);
   
   mpz_sub(a, tb, c);
}

static __inline__
void flint_mpz_ui_pow_ui(mpz_ptr a, ulong b, ulong c)
{
   FLINT_MOCK_MPZ_UI(tb, b);
   
   mpz_pow_ui(a, tb, c);
}

static __inline__
void flint_mpz_divexact_ui(mpz_ptr a, mpz_srcptr b, ulong c)
{
   FLINT_MOCK_MPZ_UI(tc, c);
   
   mpz_divexact(a, b, tc);
}

static __inline__
int flint_mpz_divisible_ui_p(mpz_srcptr a, ulong c)
{
   FLINT_MOCK_MPZ_UI(tc, c);
   
   return mpz_divisible_p(a, tc);
}

static __inline__
int flint_mpz_congruent_ui_p(mpz_srcptr a, ulong b, ulong c)
{
   FLINT_MOCK_MPZ_UI(tc, c);
   {
      FLINT_MOCK_MPZ_UI(tb, b);
   
      return mpz_congruent_p(a, tb, tc);
   }
}

static __inline__
int flint_mpz_cmp_ui(mpz_srcptr a, ulong c)
{
   FLINT_MOCK_MPZ_UI(tc, c);
   
   return mpz_cmp(a, tc);
}

static __inline__
int flint_mpz_cmp_si(mpz_srcptr a, slong c)
{
   FLINT_MOCK_MPZ_SI(tc, c);
   
   return mpz_cmp(a, tc);
}

static __inline__
int flint_mpq_cmp_si(mpq_srcptr a, slong b, ulong c)
{
   mpq_t tq;
   int res;
   FLINT_MOCK_MPZ_SI(tb, b);
   {
      FLINT_MOCK_MPZ_UI(tc, c);
      mpq_init(tq);
      mpq_set_num(tq, tb);
      mpq_set_den(tq, tc);

      res = mpq_cmp(a, tq);

      mpq_clear(tq);

      return res;
   }
}

static __inline__
int flint_mpq_cmp_ui(mpq_srcptr a, ulong b, ulong c)
{
   mpq_t tq;
   int res;
   FLINT_MOCK_MPZ_UI(tb, b);
   {
      FLINT_MOCK_MPZ_UI(tc, c);
      mpq_init(tq);
      mpq_set_num(tq, tb);
      mpq_set_den(tq, tc);

      res = mpq_cmp(a, tq);

      mpq_clear(tq);

      return res;
   }
}

static __inline__
void flint_mpq_set_si(mpq_ptr a, slong b, ulong c)
{
   FLINT_MOCK_MPZ_SI(tb, b);
   {
      FLINT_MOCK_MPZ_UI(tc, c);

      mpq_set_num(a, tb);
      mpq_set_den(a, tc);
   }
}

static __inline__
void flint_mpq_set_ui(mpq_ptr a, ulong b, ulong c)
{
   FLINT_MOCK_MPZ_UI(tb, b);
   {
      FLINT_MOCK_MPZ_UI(tc, c);

      mpq_set_num(a, tb);
      mpq_set_den(a, tc);
   }
}

static __inline__
ulong flint_mpz_cdiv_q_ui(mpz_ptr q, mpz_srcptr n, ulong c)
{
   if (FLINT_BIT_COUNT(c) <= 32)
      return (ulong) mpz_cdiv_q_ui(q, n, (unsigned long) c);
   else
   {
      ulong r;
      mpz_t rz;

      mpz_init(rz);
      
      FLINT_MOCK_MPZ_UI(tc, c);

      mpz_cdiv_qr(q, rz, n, tc);

      r = flint_mpz_get_ui(rz);

      mpz_clear(rz);

      return r;
   }
}    

static __inline__
ulong flint_mpz_fdiv_q_ui(mpz_ptr q, mpz_srcptr n, ulong c)
{
   if (FLINT_BIT_COUNT(c) <= 32)
      return (ulong) mpz_fdiv_q_ui(q, n, (unsigned long) c);
   else
   {
      ulong r;
      mpz_t rz;

      mpz_init(rz);
      
      FLINT_MOCK_MPZ_UI(tc, c);

      mpz_fdiv_qr(q, rz, n, tc);

      r = flint_mpz_get_ui(rz);

      mpz_clear(rz);

      return r;
   }
}

static __inline__
ulong flint_mpz_tdiv_q_ui(mpz_ptr q, mpz_srcptr n, ulong c)
{
    if (FLINT_BIT_COUNT(c) <= 32)
      return (ulong) mpz_tdiv_q_ui(q, n, (unsigned long) c);
   else
   {
      ulong r;
      mpz_t rz;

      mpz_init(rz);
      
      FLINT_MOCK_MPZ_UI(tc, c);

      mpz_tdiv_qr(q, rz, n, tc);

      r = flint_mpz_get_ui(rz);

      mpz_clear(rz);

      return r;
   }
}

static __inline__
ulong flint_mpz_cdiv_r_ui(mpz_ptr r, mpz_srcptr n, ulong c)
{
   FLINT_MOCK_MPZ_UI(tc, c);

   mpz_cdiv_r(r, n, tc);

   return flint_mpz_get_ui(r);
}

static __inline__
ulong flint_mpz_fdiv_r_ui(mpz_ptr r, mpz_srcptr n, ulong c)
{
   FLINT_MOCK_MPZ_UI(tc, c);

   mpz_fdiv_r(r, n, tc);

   return flint_mpz_get_ui(r);
}

static __inline__
ulong flint_mpz_tdiv_r_ui(mpz_ptr r, mpz_srcptr n, ulong c)
{
   FLINT_MOCK_MPZ_UI(tc, c);

   mpz_tdiv_r(r, n, tc);

   return flint_mpz_get_ui(r);
}

static __inline__
ulong flint_mpz_cdiv_qr_ui(mpz_ptr q, mpz_ptr r, mpz_srcptr n, ulong c)
{
   FLINT_MOCK_MPZ_UI(tc, c);

   mpz_cdiv_qr(q, r, n, tc);

   return flint_mpz_get_ui(r);
}

static __inline__
ulong flint_mpz_fdiv_qr_ui(mpz_ptr q, mpz_ptr r, mpz_srcptr n, ulong c)
{
   FLINT_MOCK_MPZ_UI(tc, c);

   mpz_fdiv_qr(q, r, n, tc);

   return flint_mpz_get_ui(r);
}

static __inline__
ulong flint_mpz_tdiv_qr_ui(mpz_ptr q, mpz_ptr r, mpz_srcptr n, ulong c)
{
   FLINT_MOCK_MPZ_UI(tc, c);

   mpz_tdiv_qr(q, r, n, tc);

   return flint_mpz_get_ui(r);
}

static __inline__
ulong flint_mpz_cdiv_ui(mpz_srcptr n, ulong c)
{
   mpz_t r;
   ulong res;

   mpz_init(r);
   FLINT_MOCK_MPZ_UI(tc, c);

   mpz_cdiv_r(r, n, tc);

   res = flint_mpz_get_ui(r);

   mpz_clear(r);

   return res;
}

static __inline__
ulong flint_mpz_fdiv_ui(mpz_srcptr n, ulong c)
{
   mpz_t r;
   ulong res;

   mpz_init(r);
   FLINT_MOCK_MPZ_UI(tc, c);

   mpz_fdiv_r(r, n, tc);

   res = flint_mpz_get_ui(r);

   mpz_clear(r);

   return res;
}

static __inline__
ulong flint_mpz_tdiv_ui(mpz_srcptr n, ulong c)
{
   mpz_t r;
   ulong res;

   mpz_init(r);
   FLINT_MOCK_MPZ_UI(tc, c);

   mpz_tdiv_r(r, n, tc);

   res = flint_mpz_get_ui(r);

   mpz_clear(r);

   return res;
}

static __inline__
ulong flint_mpz_mod_ui(mpz_ptr r, mpz_srcptr n, ulong c)
{
   FLINT_MOCK_MPZ_UI(tc, c);

   mpz_fdiv_r(r, n, tc);

   return flint_mpz_get_ui(r);
}

static __inline__
void flint_mpz_powm_ui(mpz_ptr r, mpz_srcptr b, ulong exp, mpz_srcptr mod)
{
   FLINT_MOCK_MPZ_UI(texp, exp);

   mpz_powm(r, b, texp, mod);
}

static __inline__
void flint_mpz_pow_ui(mpz_ptr r, mpz_srcptr b, ulong exp)
{
   if (exp >= (UWORD(1) << 32)) {
      printf("Exception (flint_mpz_pow_ui). Power too large.\n");
      flint_abort();
   }
   
   mpz_pow_ui(r, b, (unsigned long) exp);
}

static __inline__
void flint_mpz_fac_ui(mpz_ptr r, ulong n)
{
   if (n >= (UWORD(1) << 32)) {
      printf("Exception (flint_mpz_fac_ui). Value n too large.\n");
      flint_abort();
   }
   
   mpz_fac_ui(r, (unsigned long) n);
}

static __inline__
void flint_mpz_bin_uiui(mpz_ptr r, ulong n, ulong k)
{
   if (n >= (UWORD(1) << 32)) {
      printf("Exception (flint_mpz_bin_uiui). Value n too large.\n");
      flint_abort();
   }
   
   if (k >= (UWORD(1) << 32)) {
      printf("Exception (flint_mpz_bin_uiui). Value k too large.\n");
      flint_abort();
   }
   
   mpz_bin_uiui(r, (unsigned long) n, (unsigned long) k);
}

static __inline__
void flint_mpz_fib_ui(mpz_ptr r, ulong n)
{
   if (n >= (UWORD(1) << 32)) {
      printf("Exception (flint_mpz_fib_ui). Value n too large.\n");
      flint_abort();
   }
   
   mpz_fib_ui(r, (unsigned long) n);
}

/* mpf_set_si() -- Assign a float from a signed int.

Copyright 1993, 1994, 1995, 2000, 2001, 2002, 2004 Free Software Foundation,
Inc.

Copyright 2015 William Hart.

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

static __inline__
void flint_mpf_set_si (mpf_ptr dest, slong val)
{
  mp_size_t size;
  mp_limb_t vl;

  vl = (mp_limb_t) (val >= 0 ? val : -val);

  dest->_mp_d[0] = vl;
  size = vl != 0;

  dest->_mp_exp = size;
  dest->_mp_size = val >= 0 ? size : -size;
}

/* mpf_set_ui() -- Assign a float from an unsigned int.

Copyright 1993, 1994, 1995, 2001, 2002, 2004 Free Software Foundation, Inc.

Copyright 2015 William Hart

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

static __inline__
void flint_mpf_set_ui(mpf_ptr f, ulong val)
{
  mp_size_t size;

  f->_mp_d[0] = val;
  size = val != 0;

  f->_mp_exp = f->_mp_size = size;
}

/* mpf_get_si -- mpf to long conversion

Copyright 2001, 2002, 2004 Free Software Foundation, Inc.

Copyright 2015 William Hart

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
MA 02110-1301, USA.
*/

static __inline__
slong flint_mpf_get_si (mpf_srcptr f)
{
  mp_exp_t exp;
  mp_size_t size, abs_size;
  mp_srcptr fp;
  mp_limb_t fl;

  exp = f->_mp_exp;
  size = f->_mp_size;
  fp = f->_mp_d;

  /* fraction alone truncates to zero
     this also covers zero, since we have exp==0 for zero */
  if (exp <= 0)
    return WORD(0);

  /* there are some limbs above the radix point */

  fl = 0;
  abs_size = FLINT_ABS(size);
  if (abs_size >= exp)
    fl = fp[abs_size-exp];

  if (size > 0)
    return fl;
  else
    /* this form necessary to correctly handle -0x80..00 */
    return ~ (fl - 1);
}

/* mpf_cmp_ui -- Compare a float with an unsigned integer.

Copyright 1993, 1994, 1995, 1999, 2001, 2002 Free Software Foundation, Inc.

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

static __inline__
int flint_mpf_cmp_ui(mpf_srcptr u, ulong vval)
{
  mp_srcptr up;
  mp_size_t usize;
  mp_exp_t uexp;
  mp_limb_t ulimb;

  uexp = u->_mp_exp;
  usize = u->_mp_size;

  if (usize < 0)
    return -1;
  
  if (vval == 0)
    return usize != 0;

  if (uexp > 1)
    return 1;
  if (uexp < 1)
    return -1;

  up = u->_mp_d;

  ulimb = up[usize - 1];
  usize--;

  if (ulimb > vval)
    return 1;
  else if (ulimb < vval)
    return -1;

  while (*up == 0)
    {
      up++;
      usize--;
    }

  if (usize > 0)
    return 1;

  return 0;
}

/* mpf_fits_s*_p -- test whether an mpf fits a C signed type.

Copyright 2001, 2002 Free Software Foundation, Inc.

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

static __inline__
int flint_mpf_fits_slong_p(mpf_srcptr f)
{
  mp_size_t  fs, fn;
  mp_srcptr  fp;
  mp_exp_t   exp;
  mp_limb_t  fl;

  fs = f->_mp_size;
  if (fs == 0)
    return 1;  /* zero fits */

  exp = f->_mp_exp;
  if (exp < 1)
    return 1;  /* -1 < f < 1 truncates to zero, so fits */

  fp = f->_mp_d;
  fn = FLINT_ABS(fs);

  if (exp == 1)
    {
      fl = fp[fn-1];
    }
  else
    return 0;

  return fl <= (fs >= 0 ? (mp_limb_t) WORD_MAX : - (mp_limb_t) WORD_MIN);
}

#else

#define flint_mpz_get_si mpz_get_si
#define flint_mpz_get_ui mpz_get_ui
#define flint_mpz_set_si mpz_set_si
#define flint_mpz_set_ui mpz_set_ui
#define flint_mpz_init_set_si mpz_init_set_si
#define flint_mpz_init_set_ui mpz_init_set_ui
#define flint_mpz_add_ui mpz_add_ui
#define flint_mpz_sub_ui mpz_sub_ui
#define flint_mpz_mul_si mpz_mul_si
#define flint_mpz_mul_ui mpz_mul_ui
#define flint_mpz_addmul_ui mpz_addmul_ui
#define flint_mpz_submul_ui mpz_submul_ui
#define flint_mpz_ui_sub mpz_ui_sub
#define flint_mpz_ui_pow_ui mpz_ui_pow_ui
#define flint_mpz_cdiv_q_ui mpz_cdiv_q_ui
#define flint_mpz_cdiv_r_ui mpz_cdiv_r_ui
#define flint_mpz_cdiv_qr_ui mpz_cdiv_qr_ui
#define flint_mpz_cdiv_ui mpz_cdiv_ui
#define flint_mpz_fdiv_q_ui mpz_fdiv_q_ui
#define flint_mpz_fdiv_r_ui mpz_fdiv_r_ui
#define flint_mpz_fdiv_qr_ui mpz_fdiv_qr_ui
#define flint_mpz_fdiv_ui mpz_fdiv_ui
#define flint_mpz_tdiv_q_ui mpz_tdiv_q_ui
#define flint_mpz_tdiv_r_ui mpz_tdiv_r_ui
#define flint_mpz_tdiv_qr_ui mpz_tdiv_qr_ui
#define flint_mpz_tdiv_ui mpz_tdiv_ui
#define flint_mpz_mod_ui mpz_mod_ui
#define flint_mpz_divexact_ui mpz_divexact_ui
#define flint_mpz_divisible_ui_p mpz_divisible_ui_p
#define flint_mpz_congruent_ui_p mpz_congruent_ui_p
#define flint_mpz_powm_ui mpz_powm_ui
#define flint_mpz_pow_ui mpz_pow_ui
#define flint_mpz_fac_ui mpz_fac_ui
#define flint_mpz_bin_uiui mpz_bin_uiui
#define flint_mpz_fib_ui mpz_fib_ui
#define flint_mpz_cmp_si mpz_cmp_si
#define flint_mpz_cmp_ui mpz_cmp_ui
#define flint_mpq_cmp_si mpq_cmp_si
#define flint_mpq_cmp_ui mpq_cmp_ui
#define flint_mpq_set_si mpq_set_si
#define flint_mpq_set_ui mpq_set_ui

#define flint_mpf_set_si mpf_set_si
#define flint_mpf_set_ui mpf_set_ui
#define flint_mpf_get_si mpf_get_si
#define flint_mpf_cmp_ui mpf_cmp_ui
#define flint_mpf_fits_slong_p mpf_fits_slong_p

#endif

#if WORD_MAX == LONG_MAX

#define flint_mpf_get_d_2exp mpf_get_d_2exp

#else

/* double mpf_get_d_2exp (signed long int *exp, mpf_t src).

Copyright 2001, 2002, 2003, 2004 Free Software Foundation, Inc.

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

/*
  This function is only required for MPIR >= 3.0.0 and for GMP, but we include
  unconditionally when sizeof(ulong) != sizeof(unsigned long)
*/

static __inline__
double flint_mpf_get_d_2exp(slong * exp2, mpf_srcptr src)
{
  mp_size_t size, abs_size;
  mp_limb_t * ptr;
  int cnt;
  slong exp;

  size = src->_mp_size;
  if (size == 0)
    {
      *exp2 = 0;
      return 0.0;
    }

  ptr = src->_mp_d;
  abs_size = FLINT_ABS(size);
  count_leading_zeros (cnt, ptr[abs_size - 1]);

  exp = src->_mp_exp * FLINT_BITS - cnt;
  *exp2 = exp;

  return flint_mpn_get_d (ptr, abs_size, size,
                    (long) - (abs_size * FLINT_BITS - cnt));
}

#endif

#endif
