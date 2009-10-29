/* 
   Copyright 1991, 1992, 1993, 1994, 1996, 1997, 1999, 2000, 2001, 2002, 2003,
   2004, 2005 Free Software Foundation, Inc.

   Copyright 2009 William Hart

   This file is free software; you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation; either version 2.1 of the License, or (at your
   option) any later version.

   This file is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
   or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
   License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with this file; see the file COPYING.LIB.  If not, write to
   the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
   MA 02110-1301, USA. 
*/

#ifndef FLINT_LONGLONG_H
#define FLINT_LONGLONG_H

// this will eventually be replaced by a flag set by configure
#if __GMP_BITS_PER_MP_LIMB == 64 /* x86 : 64 bit */

#define add_ssaaaa(sh, sl, ah, al, bh, bl)                  \
  __asm__ ("addq %5,%q1\n\tadcq %3,%q0"					      \
	   : "=r" (sh), "=&r" (sl)					                  \
	   : "0"  ((mp_limb_t)(ah)), "rme" ((mp_limb_t)(bh)),		\
	     "%1" ((mp_limb_t)(al)), "rme" ((mp_limb_t)(bl)))

#define sub_ddmmss(sh, sl, ah, al, bh, bl)                  \
  __asm__ ("subq %5,%q1\n\tsbbq %3,%q0"					      \
	   : "=r" (sh), "=&r" (sl)					                  \
	   : "0" ((mp_limb_t)(ah)), "rme" ((mp_limb_t)(bh)),		\
	     "1" ((mp_limb_t)(al)), "rme" ((mp_limb_t)(bl)))

#define umul_ppmm(w1, w0, u, v)                         \
  __asm__ ("mulq %3"							                 \
	   : "=a" (w0), "=d" (w1)					              \
	   : "%0" ((mp_limb_t)(u)), "rm" ((mp_limb_t)(v)))

#define udiv_qrnnd(q, r, n1, n0, dx)                                          \
  __asm__ ("divq %4"		      	                                             \
	   : "=a" (q), "=d" (r)						                                    \
	   : "0" ((mp_limb_t)(n0)), "1" ((mp_limb_t)(n1)), "rm" ((mp_limb_t)(dx)))

#define sdiv_qrnnd(q, r, n1, n0, dx)                                          \
  __asm__ ("idivq %4"		      	                                             \
	   : "=a" (q), "=d" (r)						                                    \
	   : "0" ((mp_limb_t)(n0)), "1" ((mp_limb_t)(n1)), "rm" ((mp_limb_t)(dx)))

/* bsrq destination must be a 64-bit register, hence mp_limb_t for __cbtmp. */
#define count_leading_zeros(count, x)					                  \
  do {									                                    \
    mp_limb_t __cbtmp;							                           \
    FLINT_ASSERT ((x) != 0);							                           \
    __asm__ ("bsrq %1,%0" : "=r" (__cbtmp) : "rm" ((mp_limb_t)(x)));	\
    (count) = __cbtmp ^ 63;						                        \
  } while (0)

/* bsfq destination must be a 64-bit register, "%q0" forces this in case
   count is only an int. */
#define count_trailing_zeros(count, x)					                  \
  do {									                                    \
    FLINT_ASSERT ((x) != 0);							                           \
    __asm__ ("bsfq %1,%q0" : "=r" (count) : "rm" ((mp_limb_t)(x)));	\
  } while (0)

#else /* x86 : 32 bit */

#define add_ssaaaa(sh, sl, ah, al, bh, bl) \
  __asm__ ("addl %5,%k1\n\tadcl %3,%k0"					\
	   : "=r" (sh), "=&r" (sl)					\
	   : "0"  ((mp_limb_t)(ah)), "g" ((mp_limb_t)(bh)),			\
	     "%1" ((mp_limb_t)(al)), "g" ((mp_limb_t)(bl)))

#define sub_ddmmss(sh, sl, ah, al, bh, bl) \
  __asm__ ("subl %5,%k1\n\tsbbl %3,%k0"					\
	   : "=r" (sh), "=&r" (sl)					\
	   : "0" ((mp_limb_t)(ah)), "g" ((mp_limb_t)(bh)),			\
	     "1" ((mp_limb_t)(al)), "g" ((mp_limb_t)(bl)))

#define umul_ppmm(w1, w0, u, v) \
  __asm__ ("mull %3"							\
	   : "=a" (w0), "=d" (w1)					\
	   : "%0" ((mp_limb_t)(u)), "rm" ((mp_limb_t)(v)))

#define udiv_qrnnd(q, r, n1, n0, dx) \
  __asm__ ("divl %4"	\
	   : "=a" (q), "=d" (r)						\
	   : "0" ((mp_limb_t)(n0)), "1" ((mp_limb_t)(n1)), "rm" ((mp_limb_t)(dx)))

#define sdiv_qrnnd(q, r, n1, n0, dx) \
  __asm__ ("idivl %4"	\
	   : "=a" (q), "=d" (r)						\
	   : "0" ((mp_limb_t)(n0)), "1" ((mp_limb_t)(n1)), "rm" ((mp_limb_t)(dx)))

#define count_leading_zeros(count, x)					\
  do {									\
    mp_limb_t __cbtmp;							\
    FLINT_ASSERT ((x) != 0);							\
    __asm__ ("bsrl %1,%0" : "=r" (__cbtmp) : "rm" ((mp_limb_t)(x)));	\
    (count) = __cbtmp ^ 31;						\
  } while (0)

#define count_trailing_zeros(count, x)					\
  do {									\
    FLINT_ASSERT ((x) != 0);							\
    __asm__ ("bsfl %1,%0" : "=r" (count) : "rm" ((mp_limb_t)(x)));	\
  } while (0)

#endif /* x86 */

#define smul_ppmm(w1, w0, u, v)              \
  do {									            \
    mp_limb_t __w1;							      \
    mp_limb_t __xm0 = (u), __xm1 = (v);		\
    umul_ppmm (__w1, w0, __xm0, __xm1);		\
    (w1) = __w1 - (-(__xm0 >> (FLINT_BITS-1)) & __xm1)	\
		- (-(__xm1 >> (FLINT_BITS-1)) & __xm0);		      \
  } while (0)

#define invert_limb(invxl, xl)                        \
   do {                                               \
      mp_limb_t dummy;                                \
      udiv_qrnnd (invxl, dummy, ~(xl), ~(0L), xl);    \
   } while (0)

#if 1 // branch free variant
#define udiv_qrnnd_preinv(q, r, nh, nl, d, di)		         \
  do {									                           \
    mp_limb_t _n2, _n10, _nmask, _nadj, _q1;				      \
    mp_limb_t _xh, _xl;							                  \
    _n2 = (nh);								                     \
    _n10 = (nl);							                        \
    _nmask = (mp_limb_signed_t) (_n10) >> (FLINT_BITS - 1); \
    _nadj = _n10 + (_nmask & (d));					            \
    umul_ppmm (_xh, _xl, di, _n2 - _nmask);				      \
    add_ssaaaa (_xh, _xl, _xh, _xl, _n2, _nadj);		      \
    _q1 = ~_xh;								                     \
    umul_ppmm (_xh, _xl, _q1, d);					            \
    add_ssaaaa (_xh, _xl, _xh, _xl, nh, nl);				      \
    _xh -= (d);					/* xh = 0 or -1 */            \
    (r) = _xl + ((d) & _xh);						               \
    (q) = _xh - _q1;							                     \
  } while (0)

#else // branched variant, slower on K10 for general inputs
#define udiv_qrnnd_preinv(q, r, nh, nl, d, di)				\
  do {									                        \
    mp_limb_t _q, _ql, _r;						               \
    mp_limb_t _xh, _xl;							               \
    FLINT_ASSERT ((d) != 0);							         \
    umul_ppmm (_q, _ql, (nh), (di));					      \
    _q += (nh);	/* Compensate, di is 2^64 too small */ \
    umul_ppmm (_xh, _xl, _q, (d));					         \
    sub_ddmmss (_xh, _r, (nh), (nl), _xh, _xl);				\
    if (_xh != 0)							                     \
      {									                        \
	sub_ddmmss (_xh, _r, _xh, _r, 0, (d));				      \
	_q += 1;							                           \
	if (_xh != 0)							                     \
	  {								                           \
	    _r -= (d);							                     \
	    _q += 1;							                     \
	  }								                           \
      }									                        \
    if (_r >= (d))							                  \
      {									                        \
	_r -= (d);							                        \
	_q += 1;							                           \
      }									                        \
    (r) = _r;								                     \
    (q) = _q;								                     \
  } while (0)
#endif

#endif
