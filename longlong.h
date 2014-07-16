/* 
   Copyright 1991, 1992, 1993, 1994, 1996, 1997, 1999, 2000, 2001, 2002, 2003,
   2004, 2005 Free Software Foundation, Inc.

   Copyright 2009 William Hart
   Copyright 2011 Fredrik Johansson

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

/*
   N.B: This file has been adapted from code found in GMP 4.2.1.
*/

#ifndef FLINT_LONGLONG_H
#define FLINT_LONGLONG_H

#ifdef __cplusplus
 extern "C" {
#endif

/* x86 : 64 bit */
#if (GMP_LIMB_BITS == 64 && defined (__amd64__)) 

#define add_sssaaaaaa(sh, sm, sl, ah, am, al, bh, bm, bl)  \
  __asm__ ("addq %8,%q2\n\tadcq %6,%q1\n\tadcq %4,%q0"     \
       : "=r" (sh), "=r" (sm), "=&r" (sl)                  \
       : "0"  ((mp_limb_t)(ah)), "rme" ((mp_limb_t)(bh)),  \
         "1"  ((mp_limb_t)(am)), "rme" ((mp_limb_t)(bm)),  \
         "2"  ((mp_limb_t)(al)), "rme" ((mp_limb_t)(bl)))  \

#define add_ssaaaa(sh, sl, ah, al, bh, bl)                 \
  __asm__ ("addq %5,%q1\n\tadcq %3,%q0"                    \
       : "=r" (sh), "=&r" (sl)                             \
       : "0"  ((mp_limb_t)(ah)), "rme" ((mp_limb_t)(bh)),  \
         "%1" ((mp_limb_t)(al)), "rme" ((mp_limb_t)(bl)))

#define sub_ddmmss(sh, sl, ah, al, bh, bl)                \
  __asm__ ("subq %5,%q1\n\tsbbq %3,%q0"                   \
       : "=r" (sh), "=&r" (sl)                            \
       : "0" ((mp_limb_t)(ah)), "rme" ((mp_limb_t)(bh)),  \
         "1" ((mp_limb_t)(al)), "rme" ((mp_limb_t)(bl)))

#define umul_ppmm(w1, w0, u, v)                         \
  __asm__ ("mulq %3"                                    \
       : "=a" (w0), "=d" (w1)                           \
       : "%0" ((mp_limb_t)(u)), "rm" ((mp_limb_t)(v)))

#define udiv_qrnnd(q, r, n1, n0, dx)                                            \
  __asm__ ("divq %4"                                                            \
       : "=a" (q), "=d" (r)                                                     \
       : "0" ((mp_limb_t)(n0)), "1" ((mp_limb_t)(n1)), "rm" ((mp_limb_t)(dx)))

#define sdiv_qrnnd(q, r, n1, n0, dx)                                            \
  __asm__ ("idivq %4"                                                           \
       : "=a" (q), "=d" (r)                                                     \
       : "0" ((mp_limb_t)(n0)), "1" ((mp_limb_t)(n1)), "rm" ((mp_limb_t)(dx)))

/* bsrq destination must be a 64-bit register, hence mp_limb_t for __cbtmp. */
#define count_leading_zeros(count, x)                                 \
  do {                                                                \
    mp_limb_t __cbtmp;                                                \
    FLINT_ASSERT ((x) != 0);                                          \
    __asm__ ("bsrq %1,%0" : "=r" (__cbtmp) : "rm" ((mp_limb_t)(x)));  \
    (count) = __cbtmp ^ (mp_limb_t) 63;                               \
  } while (0)

/* bsfq destination must be a 64-bit register, "%q0" forces this in case
   count is only an int. */
#define count_trailing_zeros(count, x)                               \
  do {                                                               \
    FLINT_ASSERT ((x) != 0);                                         \
    __asm__ ("bsfq %1,%q0" : "=r" (count) : "rm" ((mp_limb_t)(x)));  \
  } while (0)

#endif /* x86_64 */

/* x86 : 32 bit */
#if (GMP_LIMB_BITS == 32 && (defined (__i386__) \
   || defined (__i486__) || defined(__amd64__)))

#define add_sssaaaaaa(sh, sm, sl, ah, am, al, bh, bm, bl)  \
  __asm__ ("addl %8,%k2\n\tadcl %6,%k1\n\tadcl %4,%k0"     \
       : "=r" (sh), "=r" (sm), "=&r" (sl)                  \
       : "0"  ((mp_limb_t)(ah)), "g" ((mp_limb_t)(bh)),    \
         "1"  ((mp_limb_t)(am)), "g" ((mp_limb_t)(bm)),    \
         "2"  ((mp_limb_t)(al)), "g" ((mp_limb_t)(bl)))    \

#define add_ssaaaa(sh, sl, ah, al, bh, bl)               \
  __asm__ ("addl %5,%k1\n\tadcl %3,%k0"                  \
       : "=r" (sh), "=&r" (sl)                           \
       : "0"  ((mp_limb_t)(ah)), "g" ((mp_limb_t)(bh)),  \
         "%1" ((mp_limb_t)(al)), "g" ((mp_limb_t)(bl)))

#define sub_ddmmss(sh, sl, ah, al, bh, bl)              \
  __asm__ ("subl %5,%k1\n\tsbbl %3,%k0"                 \
       : "=r" (sh), "=&r" (sl)                          \
       : "0" ((mp_limb_t)(ah)), "g" ((mp_limb_t)(bh)),  \
         "1" ((mp_limb_t)(al)), "g" ((mp_limb_t)(bl)))

#define umul_ppmm(w1, w0, u, v)                         \
  __asm__ ("mull %3"                                    \
       : "=a" (w0), "=d" (w1)                           \
       : "%0" ((mp_limb_t)(u)), "rm" ((mp_limb_t)(v)))

#define udiv_qrnnd(q, r, n1, n0, dx)                                            \
  __asm__ ("divl %4"                                                            \
       : "=a" (q), "=d" (r)                                                     \
       : "0" ((mp_limb_t)(n0)), "1" ((mp_limb_t)(n1)), "rm" ((mp_limb_t)(dx)))

#define sdiv_qrnnd(q, r, n1, n0, dx)                                            \
  __asm__ ("idivl %4"                                                           \
       : "=a" (q), "=d" (r)                                                     \
       : "0" ((mp_limb_t)(n0)), "1" ((mp_limb_t)(n1)), "rm" ((mp_limb_t)(dx)))

#define count_leading_zeros(count, x)                                 \
  do {                                                                \
    mp_limb_t __cbtmp;                                                \
    FLINT_ASSERT ((x) != 0);                                          \
    __asm__ ("bsrl %1,%0" : "=r" (__cbtmp) : "rm" ((mp_limb_t)(x)));  \
    (count) = __cbtmp ^ (mp_limb_t) 31;                               \
  } while (0)

#define count_trailing_zeros(count, x)                              \
  do {                                                              \
    FLINT_ASSERT ((x) != 0);                                        \
    __asm__ ("bsfl %1,%0" : "=r" (count) : "rm" ((mp_limb_t)(x)));  \
  } while (0)

#endif /* x86 */

/* Itanium */
#if (GMP_LIMB_BITS == 64 && defined (__ia64))

/* This form encourages gcc (pre-release 3.4 at least) to emit predicated
   "sub r=r,r" and "sub r=r,r,1", giving a 2 cycle latency.  The generic
   code using "al<bl" arithmetically comes out making an actual 0 or 1 in a
   register, which takes an extra cycle.  */
#define sub_ddmmss(sh, sl, ah, al, bh, bl)      \
  do {                                          \
    mp_limb_t __x;                              \
    __x = (al) - (bl);                          \
    if ((al) < (bl))                            \
      (sh) = (ah) - (bh) - 1;                   \
    else                                        \
      (sh) = (ah) - (bh);                       \
    (sl) = __x;                                 \
  } while (0)

#define count_trailing_zeros(count, x)			   \
  do {									               \
    mp_limb_t __ctz_x = (x);		   		      \
    __asm__ ("popcnt %0 = %1"						   \
	     : "=r" (count)						         \
	     : "r" ((__ctz_x-1) & ~__ctz_x));		   \
  } while (0)

/* Do both product parts in assembly, since that gives better code with
   all gcc versions.  Some callers will just use the upper part, and in
   that situation we waste an instruction, but not any cycles.  */
#define umul_ppmm(ph, pl, m0, m1)                                 \
    __asm__ ("xma.hu %0 = %2, %3, f0\n\txma.l %1 = %2, %3, f0"		\
	     : "=&f" (ph), "=f" (pl)				                       	\
	     : "f" (m0), "f" (m1))

#endif /* Itanium */

/* ARM */
#if (GMP_LIMB_BITS == 32 && defined (__arm__))

#define add_ssaaaa(sh, sl, ah, al, bh, bl)                \
  __asm__ ("adds\t%1, %4, %5\n\tadc\t%0, %2, %3"			 \
	   : "=r" (sh), "=&r" (sl)					                \
	   : "r" (ah), "rI" (bh), "%r" (al), "rI" (bl) : "cc")

/* rsbs is only available in ARM and Thumb modes, not in Thumb2 mode */
#if !defined(__thumb2__)

#define sub_ddmmss(sh, sl, ah, al, bh, bl)                    \
  do {									                             \
    if (__builtin_constant_p (al))					              \
      {									                             \
	if (__builtin_constant_p (ah))					              \
	  __asm__ ("rsbs\t%1, %5, %4\n\trsc\t%0, %3, %2"		     \
		   : "=r" (sh), "=&r" (sl)				                    \
		   : "rI" (ah), "r" (bh), "rI" (al), "r" (bl) : "cc");  \
	else								                                \
	  __asm__ ("rsbs\t%1, %5, %4\n\tsbc\t%0, %2, %3"		     \
		   : "=r" (sh), "=&r" (sl)				                    \
		   : "r" (ah), "rI" (bh), "rI" (al), "r" (bl) : "cc");  \
      }									                             \
    else if (__builtin_constant_p (ah))					        \
      {									                             \
	if (__builtin_constant_p (bl))					              \
	  __asm__ ("subs\t%1, %4, %5\n\trsc\t%0, %3, %2"		     \
		   : "=r" (sh), "=&r" (sl)				                    \
		   : "rI" (ah), "r" (bh), "r" (al), "rI" (bl) : "cc");  \
	else								                                \
	  __asm__ ("rsbs\t%1, %5, %4\n\trsc\t%0, %3, %2"		     \
		   : "=r" (sh), "=&r" (sl)				                    \
		   : "rI" (ah), "r" (bh), "rI" (al), "r" (bl) : "cc");  \
      }									                             \
    else if (__builtin_constant_p (bl))					        \
      {									                             \
	if (__builtin_constant_p (bh))					              \
	  __asm__ ("subs\t%1, %4, %5\n\tsbc\t%0, %2, %3"		     \
		   : "=r" (sh), "=&r" (sl)				                    \
		   : "r" (ah), "rI" (bh), "r" (al), "rI" (bl) : "cc");  \
	else								                                \
	  __asm__ ("subs\t%1, %4, %5\n\trsc\t%0, %3, %2"		     \
		   : "=r" (sh), "=&r" (sl)				                    \
		   : "rI" (ah), "r" (bh), "r" (al), "rI" (bl) : "cc");  \
      }									                             \
    else /* only bh might be a constant */				        \
      __asm__ ("subs\t%1, %4, %5\n\tsbc\t%0, %2, %3"			  \
	       : "=r" (sh), "=&r" (sl)					              \
	       : "r" (ah), "rI" (bh), "r" (al), "rI" (bl) : "cc"); \
    } while (0)

#else

#define sub_ddmmss(sh, sl, ah, al, bh, bl)                    \
  do {                                                         \
     __asm__ ("subs\t%1, %4, %5\n\tsbc\t%0, %2, %3"                     \
              : "=r" (sh), "=&r" (sl)                                   \
              : "r" (ah), "r" (bh), "r" (al), "r" (bl) : "cc");       \
 } while (0)

#endif /* not Thumb2 */

#define umul_ppmm(xh, xl, a, b) \
  __asm__ ("umull %0,%1,%2,%3" : "=&r" (xl), "=&r" (xh) : "r" (a), "r" (b))

#define smul_ppmm(xh, xl, a, b) \
  __asm__ ("smull %0,%1,%2,%3" : "=&r" (xl), "=&r" (xh) : "r" (a), "r" (b))

#endif /* ARM */

/* fallback code */
#if !(defined (__i386__) || defined (__i486__) || defined(__amd64__))

#define __BITS4 (GMP_LIMB_BITS/4)
#define __ll_B ((mp_limb_t) 1 << (GMP_LIMB_BITS / 2))
#define __ll_lowpart(t) ((mp_limb_t) (t) & (__ll_B - 1))
#define __ll_highpart(t) ((mp_limb_t) (t) >> (GMP_LIMB_BITS / 2))
#define __highbit (~(mp_limb_t)0 ^ ((~(mp_limb_t)0) >> 1))

#define NEED_CLZ_TAB

#if !(GMP_LIMB_BITS == 32 && defined (__arm__))
#if !(GMP_LIMB_BITS == 64 && defined (__ia64))

#define umul_ppmm(w1, w0, u, v)				 \
  do {									          \
    mp_limb_t __x0, __x1, __x2, __x3;		 \
    mp_limb_t __ul, __vl, __uh, __vh;		 \
    mp_limb_t __u = (u), __v = (v);			 \
									                \
    __ul = __ll_lowpart (__u);				 \
    __uh = __ll_highpart (__u);				 \
    __vl = __ll_lowpart (__v);				 \
    __vh = __ll_highpart (__v);				 \
									                \
    __x0 = (mp_limb_t) __ul * __vl;			 \
    __x1 = (mp_limb_t) __ul * __vh;			 \
    __x2 = (mp_limb_t) __uh * __vl;			 \
    __x3 = (mp_limb_t) __uh * __vh;			 \
									                \
    __x1 += __ll_highpart (__x0);/* this can't give carry */            \
    __x1 += __x2;		/* but this indeed can */		                     \
    if (__x1 < __x2)		/* did we get it? */			                     \
      __x3 += __ll_B;		/* yes, add it in the proper pos. */         \
									                                             \
    (w1) = __x3 + __ll_highpart (__x1);					                  \
    (w0) = (__x1 << GMP_LIMB_BITS/2) + __ll_lowpart (__x0);		\
  } while (0)

#endif
#endif

#if !(GMP_LIMB_BITS == 32 && defined (__arm__))

#define add_ssaaaa(sh, sl, ah, al, bh, bl) \
  do {									          \
    mp_limb_t __x;								 \
    __x = (al) + (bl);							 \
    (sh) = (ah) + (bh) + (__x < (al));		 \
    (sl) = __x;								    \
  } while (0)

#endif

#define add_sssaaaaaa(sh, sm, sl, ah, am, al, bh, bm, bl)           \
  do {                                                              \
    mp_limb_t __t, __u;                                             \
    add_ssaaaa(__t, sl, (mp_limb_t) 0, al, (mp_limb_t) 0, bl);      \
    add_ssaaaa(__u, sm, (mp_limb_t) 0, am, (mp_limb_t) 0, bm);      \
    add_ssaaaa(sh, sm, ah + bh, sm, __u, __t);                      \
  } while (0)

#if !((GMP_LIMB_BITS == 64 && defined (__ia64)) ||      \
      (GMP_LIMB_BITS == 32 && defined (__arm__)))

#define sub_ddmmss(sh, sl, ah, al, bh, bl) \
  do {									          \
    mp_limb_t __x;								 \
    __x = (al) - (bl);							 \
    (sh) = (ah) - (bh) - ((al) < (bl));    \
    (sl) = __x;								    \
  } while (0)

#endif

#define udiv_qrnnd_int(q, r, n1, n0, d)                                \
  do {									                                      \
    mp_limb_t __d1, __d0, __q1, __q0, __r1, __r0, __m;			        \
									                                            \
    FLINT_ASSERT ((d) != 0);							                       \
    FLINT_ASSERT ((n1) < (d));						                       \
    						                                                  \
    __d1 = __ll_highpart (d);						                          \
    __d0 = __ll_lowpart (d);						                          \
					                                                        \
    __q1 = (n1) / __d1;							                             \
    __r1 = (n1) - __q1 * __d1;						                       \
    __m = __q1 * __d0;							                             \
    __r1 = __r1 * __ll_B | __ll_highpart (n0);				              \
    if (__r1 < __m)							                                \
      {									                                      \
	__q1--, __r1 += (d);						                                \
	if (__r1 >= (d)) /* i.e. we didn't get carry when adding to __r1 */ \
	  if (__r1 < __m)						                                   \
	    __q1--, __r1 += (d);					                             \
      }									                                      \
    __r1 -= __m;							                                   \
									                                            \
    __q0 = __r1 / __d1;							                             \
    __r0 = __r1  - __q0 * __d1;						                       \
    __m = __q0 * __d0;							                             \
    __r0 = __r0 * __ll_B | __ll_lowpart (n0);				              \
    if (__r0 < __m)							                                \
      {									                                      \
	__q0--, __r0 += (d);						                                \
	if (__r0 >= (d))						                                   \
	  if (__r0 < __m)						                                   \
	    __q0--, __r0 += (d);					                             \
      }									                                      \
    __r0 -= __m;							                                   \
									                                            \
    (q) = __q1 * __ll_B | __q0;						                       \
    (r) = __r0;								                                \
  } while (0)

#define count_leading_zeros(count, x)                        \
  do {									                            \
    mp_limb_t __xr = (x);							                \
    mp_limb_t __a;								                   \
									                                  \
    if (GMP_LIMB_BITS == 32)						       \
      {									                            \
	__a = __xr < ((mp_limb_t) 1 << 2*__BITS4)				       \
	  ? (__xr < ((mp_limb_t) 1 << __BITS4) ? 1 : __BITS4 + 1) \
	  : (__xr < ((mp_limb_t) 1 << 3*__BITS4) ? 2*__BITS4 + 1	 \
	  : 3*__BITS4 + 1);						                      \
      }									                            \
    else								                               \
      {									                            \
	for (__a = GMP_LIMB_BITS - 8; __a > 0; __a -= 8) \
	  if (((__xr >> __a) & 0xff) != 0)				             \
	    break;							                            \
	++__a;								                            \
      }									                            \
									                                  \
    (count) = GMP_LIMB_BITS + 1 - __a - __flint_clz_tab[__xr >> __a]; \
  } while (0)

#if !(GMP_LIMB_BITS == 64 && defined (__ia64))

#define count_trailing_zeros(count, x)                 \
  do {									                      \
    mp_limb_t __ctz_x = (x);						          \
    mp_limb_t __ctz_c;							             \
    FLINT_ASSERT (__ctz_x != 0);						       \
    count_leading_zeros (__ctz_c, __ctz_x & -__ctz_x); \
    (count) = GMP_LIMB_BITS - 1 - __ctz_c;	 \
  } while (0)

#endif

#define udiv_qrnnd(q, r, n1, n0, d)                  \
    do {                                             \
       mp_limb_t __norm;                             \
       count_leading_zeros(__norm, (d));             \
       if (__norm)                                   \
       {                                             \
           udiv_qrnnd_int((q), (r), ((n1) << __norm) + ((n0) >> (GMP_LIMB_BITS - __norm)), (n0) << __norm, (d) << __norm); \
          (r) = ((mp_limb_t) (r) >> __norm);         \
       } else                                        \
          udiv_qrnnd_int((q), (r), (n1), (n0), (d)); \
    } while (0)

#define sdiv_qrnnd(q, r, n1, n0, d) \
  do {                              \
    mp_limb_t __n1, __n0, __d;      \
    int __sgn1 = 0, __sgn2 = 0;     \
    if ((n1) & __highbit)           \
    {                               \
       __n0 = -(n0);                \
       __n1 = ~(n1) + (__n0 == 0);  \
       __sgn1 = ~__sgn1;            \
    } else                          \
    {                               \
      __n0 = (n0);                  \
      __n1 = (n1);                  \
    }                               \
    if ((d) & __highbit)            \
    {                               \
        __d = -(d);                 \
        __sgn2 = ~__sgn2;           \
    } else                          \
    {                               \
        __d = (d);                  \
    }                               \
    udiv_qrnnd((q), (r), (mp_limb_t) __n1, (mp_limb_t) __n0, (mp_limb_t) __d); \
    if (__sgn1 ^ __sgn2)            \
    {                               \
        (q) = -(q);                 \
        if (!__sgn2)                \
        {                           \
            (q)--;                  \
            (r) = (__d) - (r);      \
        }                           \
    } else if (__sgn1 && __sgn2)    \
    {                               \
       (q)++;                       \
       (r) = (__d) - (r);           \
    }                               \
  } while (0)

#endif /* non x86 fallback code */

#if !(GMP_LIMB_BITS == 32 && defined (__arm__))

#define smul_ppmm(w1, w0, u, v)                         \
  do {                                                  \
    mp_limb_t __w1;                                     \
    mp_limb_t __xm0 = (u), __xm1 = (v);                 \
    umul_ppmm (__w1, w0, __xm0, __xm1);                 \
    (w1) = __w1 - (-(__xm0 >> (FLINT_BITS-1)) & __xm1)  \
        - (-(__xm1 >> (FLINT_BITS-1)) & __xm0);         \
  } while (0)

#endif

#define invert_limb(invxl, xl)                      \
   do {                                             \
      mp_limb_t dummy;                              \
      udiv_qrnnd (invxl, dummy, ~(xl), ~(WORD(0)), xl);  \
   } while (0)

/*
   Branch free variant
 */
#if 1
#define udiv_qrnnd_preinv(q, r, nh, nl, d, di)               \
  do {                                                       \
    mp_limb_t _n2, _n10, _nmask, _nadj, _q1;                 \
    mp_limb_t _xh, _xl;                                      \
    _n2 = (nh);                                              \
    _n10 = (nl);                                             \
    _nmask = (mp_limb_signed_t) (_n10) >> (FLINT_BITS - 1);  \
    _nadj = _n10 + (_nmask & (d));                           \
    umul_ppmm (_xh, _xl, di, _n2 - _nmask);                  \
    add_ssaaaa (_xh, _xl, _xh, _xl, _n2, _nadj);             \
    _q1 = ~_xh;                                              \
    umul_ppmm (_xh, _xl, _q1, d);                            \
    add_ssaaaa (_xh, _xl, _xh, _xl, nh, nl);                 \
    _xh -= (d);                 /* xh = 0 or -1 */           \
    (r) = _xl + ((d) & _xh);                                 \
    (q) = _xh - _q1;                                         \
  } while (0)

/*
   Branched variant, slower on K10 for general inputs
 */
#else
#define udiv_qrnnd_preinv(q, r, nh, nl, d, di)          \
  do {                                                  \
    mp_limb_t _q, _ql, _r;                              \
    mp_limb_t _xh, _xl;                                 \
    FLINT_ASSERT ((d) != 0);                            \
    umul_ppmm (_q, _ql, (nh), (di));                    \
    _q += (nh); /* Compensate, di is 2^64 too small */  \
    umul_ppmm (_xh, _xl, _q, (d));                      \
    sub_ddmmss (_xh, _r, (nh), (nl), _xh, _xl);         \
    if (_xh != 0)                                       \
    {                                                   \
      sub_ddmmss (_xh, _r, _xh, _r, 0, (d));            \
      _q += 1;                                          \
      if (_xh != 0)                                     \
      {                                                 \
        _r -= (d);                                      \
        _q += 1;                                        \
      }                                                 \
    }                                                   \
    if (_r >= (d))                                      \
    {                                                   \
      _r -= (d);                                        \
      _q += 1;                                          \
    }                                                   \
    (r) = _r;                                           \
    (q) = _q;                                           \
  } while (0)
#endif

#ifdef __cplusplus
}
#endif

#endif
