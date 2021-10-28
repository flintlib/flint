/* 
   Copyright 1991, 1992, 1993, 1994, 1996, 1997, 1999, 2000, 2001, 2002, 2003,
   2004, 2005 Free Software Foundation, Inc.

   Copyright 2009, 2015, 2016 William Hart
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

/* Undefine to make the ifndef logic below for the fallback 
   work even if the symbols are already defined (e.g. by givaro).  */
#undef count_leading_zeros
#undef count_trailing_zeros

/* 1 if we know that the hardware is strongly-ordered */
#define FLINT_KNOW_STRONG_ORDER 0

/* x86 : 64 bit */
#if (GMP_LIMB_BITS == 64 && defined (__amd64__)) 

#undef FLINT_KNOW_STRONG_ORDER
#define FLINT_KNOW_STRONG_ORDER 1

#define add_ssssaaaaaaaa(s3, s2, s1, s0, a3, a2, a1, a0, b3, b2, b1, b0)  \
  __asm__ ("addq %11,%q3\n\tadcq %9,%q2\n\tadcq %7,%q1\n\tadcq %5,%q0"    \
       : "=r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0)                    \
       : "0"  ((mp_limb_t)(a3)), "rme" ((mp_limb_t)(b3)),                 \
         "1"  ((mp_limb_t)(a2)), "rme" ((mp_limb_t)(b2)),                 \
         "2"  ((mp_limb_t)(a1)), "rme" ((mp_limb_t)(b1)),                 \
         "3"  ((mp_limb_t)(a0)), "rme" ((mp_limb_t)(b0)))                 \

#define add_sssaaaaaa(sh, sm, sl, ah, am, al, bh, bm, bl)  \
  __asm__ ("addq %8,%q2\n\tadcq %6,%q1\n\tadcq %4,%q0"     \
       : "=r" (sh), "=&r" (sm), "=&r" (sl)                  \
       : "0"  ((mp_limb_t)(ah)), "rme" ((mp_limb_t)(bh)),  \
         "1"  ((mp_limb_t)(am)), "rme" ((mp_limb_t)(bm)),  \
         "2"  ((mp_limb_t)(al)), "rme" ((mp_limb_t)(bl)))  \

#define sub_dddmmmsss(dh, dm, dl, mh, mm, ml, sh, sm, sl)  \
  __asm__ ("subq %8,%q2\n\tsbbq %6,%q1\n\tsbbq %4,%q0"     \
       : "=r" (dh), "=&r" (dm), "=&r" (dl)                  \
       : "0"  ((mp_limb_t)(mh)), "rme" ((mp_limb_t)(sh)),  \
         "1"  ((mp_limb_t)(mm)), "rme" ((mp_limb_t)(sm)),  \
         "2"  ((mp_limb_t)(ml)), "rme" ((mp_limb_t)(sl)))  \

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

#define smul_ppmm(w1, w0, u, v)                         \
  __asm__ ("imulq %3"                                   \
       : "=a" (w0), "=d" (w1)                           \
       : "%0" ((mp_limb_t)(u)), "rm" ((mp_limb_t)(v)))

#define udiv_qrnnd(q, r, n1, n0, dx)                                            \
  __asm__ volatile ("divq %4"                                                            \
       : "=a" (q), "=d" (r)                                                     \
       : "0" ((mp_limb_t)(n0)), "1" ((mp_limb_t)(n1)), "rm" ((mp_limb_t)(dx)))

#define sdiv_qrnnd(q, r, n1, n0, dx)                                            \
  __asm__ volatile ("idivq %4"                                                           \
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

#define byte_swap(x)                                                 \
  do {                                                               \
    __asm__("bswapq %q0" : "=r"(x) : "0"(x));                         \
  } while (0)

#endif /* x86_64 */

/* x86 : 32 bit */
#if (GMP_LIMB_BITS == 32 && (defined (__i386__) \
   || defined (__i486__) || defined(__amd64__)))

#undef FLINT_KNOW_STRONG_ORDER
#define FLINT_KNOW_STRONG_ORDER 1


#define add_ssssaaaaaaaa(s3, s2, s1, s0, a3, a2, a1, a0, b3, b2, b1, b0)  \
  __asm__ ("addl %11,%k3\n\tadcl %9,%k2\n\tadcl %7,%k1\n\tadcl %5,%k0"    \
       : "=r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0)                    \
       : "0"  ((mp_limb_t)(a3)), "g" ((mp_limb_t)(b3)),                   \
         "1"  ((mp_limb_t)(a2)), "g" ((mp_limb_t)(b2)),                   \
         "2"  ((mp_limb_t)(a1)), "g" ((mp_limb_t)(b1)),                   \
         "3"  ((mp_limb_t)(a0)), "g" ((mp_limb_t)(b0)))                   \

#define add_sssaaaaaa(sh, sm, sl, ah, am, al, bh, bm, bl)  \
  __asm__ ("addl %8,%k2\n\tadcl %6,%k1\n\tadcl %4,%k0"     \
       : "=r" (sh), "=r" (sm), "=&r" (sl)                  \
       : "0"  ((mp_limb_t)(ah)), "g" ((mp_limb_t)(bh)),    \
         "1"  ((mp_limb_t)(am)), "g" ((mp_limb_t)(bm)),    \
         "2"  ((mp_limb_t)(al)), "g" ((mp_limb_t)(bl)))    \

#define sub_dddmmmsss(dh, dm, dl, mh, mm, ml, sh, sm, sl)  \
  __asm__ ("subl %8,%k2\n\tsbbl %6,%k1\n\tsbbl %4,%k0"     \
       : "=r" (dh), "=r" (dm), "=&r" (dl)                  \
       : "0"  ((mp_limb_t)(mh)), "g" ((mp_limb_t)(sh)),    \
         "1"  ((mp_limb_t)(mm)), "g" ((mp_limb_t)(sm)),    \
         "2"  ((mp_limb_t)(ml)), "g" ((mp_limb_t)(sl)))    \

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

#define smul_ppmm(w1, w0, u, v)                         \
  __asm__ ("imull %3"                                   \
       : "=a" (w0), "=d" (w1)                           \
       : "%0" ((mp_limb_t)(u)), "rm" ((mp_limb_t)(v)))

#define udiv_qrnnd(q, r, n1, n0, dx)                                            \
  __asm__ volatile ("divl %4"                                                            \
       : "=a" (q), "=d" (r)                                                     \
       : "0" ((mp_limb_t)(n0)), "1" ((mp_limb_t)(n1)), "rm" ((mp_limb_t)(dx)))

#define sdiv_qrnnd(q, r, n1, n0, dx)                                            \
  __asm__ volatile ("idivl %4"                                                           \
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

#define byte_swap(x)                                                 \
  do {                                                               \
    __asm__("bswap %0" : "=r"(x) : "0"(x));                          \
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


#if defined(__arm64__)

#define umul_ppmm(xh, xl, a, b) \
  __asm__ ("mul %0,%2,%3\numulh %1,%2,%3" : "=&r" (xl), "=&r" (xh) : "r" (a), "r" (b))

#define smul_ppmm(xh, xl, a, b) \
  __asm__ ("mul %0,%2,%3\nsmulh %1,%2,%3" : "=&r" (xl), "=&r" (xh) : "r" (a), "r" (b))

#define add_ssaaaa(sh, sl, ah, al, bh, bl)               \
  __asm__ ("adds %1,%3,%5\n\tadc %0,%2,%4"               \
       : "=r" (sh), "=&r" (sl)                           \
       : "r"  ((mp_limb_t)(ah)), "r" ((mp_limb_t)(al)),  \
         "r" ((mp_limb_t)(bh)), "rI" ((mp_limb_t)(bl))   \
       : "cc")

#define add_sssaaaaaa(sh, sm, sl, ah, am, al, bh, bm, bl)                     \
  __asm__ ("adds %2,%5,%8\n\tadcs %1,%4,%7\n\tadc %0,%3,%6"                   \
       : "=r" (sh), "=&r" (sm), "=&r" (sl)                                    \
       : "r"  ((mp_limb_t)(ah)), "r" ((mp_limb_t)(am)), "r" ((mp_limb_t)(al)),\
         "r" ((mp_limb_t)(bh)), "r" ((mp_limb_t)(bm)), "rI" ((mp_limb_t)(bl)) \
       : "cc")

#define add_ssssaaaaaaaa(s3, s2, s1, s0, a3, a2, a1, a0, b3, b2, b1, b0)      \
  __asm__ ("adds %3,%7,%11\n\tadcs %2,%6,%10\n\tadcs %1,%5,%9\n\tadc %0,%4,%8"\
       : "=r" (s3), "=&r" (s2), "=&r" (s1), "=&r" (s0)                        \
       : "r" ((mp_limb_t)(a3)), "r" ((mp_limb_t)(a2)),                        \
         "r" ((mp_limb_t)(a1)), "r" ((mp_limb_t)(a0)),                        \
         "r" ((mp_limb_t)(b3)), "r" ((mp_limb_t)(b2)),                        \
         "r" ((mp_limb_t)(b1)), "rI" ((mp_limb_t)(b0))                        \
       : "cc")

#define sub_ddmmss(sh, sl, ah, al, bh, bl)               \
  __asm__ ("subs %1,%3,%5\n\tsbc %0,%2,%4"               \
       : "=r" (sh), "=&r" (sl)                           \
       : "r"  ((mp_limb_t)(ah)), "r" ((mp_limb_t)(al)),  \
         "r" ((mp_limb_t)(bh)), "rI" ((mp_limb_t)(bl))   \
       : "cc")

#define sub_dddmmmsss(sh, sm, sl, ah, am, al, bh, bm, bl)                     \
  __asm__ ("subs %2,%5,%8\n\tsbcs %1,%4,%7\n\tsbc %0,%3,%6"                   \
       : "=r" (sh), "=&r" (sm), "=&r" (sl)                                    \
       : "r"  ((mp_limb_t)(ah)), "r" ((mp_limb_t)(am)), "r" ((mp_limb_t)(al)),\
         "r" ((mp_limb_t)(bh)), "r" ((mp_limb_t)(bm)), "rI" ((mp_limb_t)(bl)) \
       : "cc")

#endif

/* fallback code */
#if !(defined (__i386__) || defined (__i486__) || defined(__amd64__))

#define __BITS4 (GMP_LIMB_BITS/4)
#define __ll_B ((mp_limb_t) 1 << (GMP_LIMB_BITS / 2))
#define __ll_lowpart(t) ((mp_limb_t) (t) & (__ll_B - 1))
#define __ll_highpart(t) ((mp_limb_t) (t) >> (GMP_LIMB_BITS / 2))
#define __highbit (~(mp_limb_t)0 ^ ((~(mp_limb_t)0) >> 1))

#ifndef NEED_CLZ_TAB
#define NEED_CLZ_TAB
#endif

#if !(GMP_LIMB_BITS == 32 && defined (__arm__))
#if !(GMP_LIMB_BITS == 64 && defined (__ia64))
#if !defined(__arm64__)

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
#endif

#if !(GMP_LIMB_BITS == 32 && defined (__arm__))
#if !defined(__arm64__)

#define add_ssaaaa(sh, sl, ah, al, bh, bl) \
  do {									          \
    mp_limb_t __x;								 \
    __x = (al) + (bl);							 \
    (sh) = (ah) + (bh) + (__x < (al));		 \
    (sl) = __x;								    \
  } while (0)

#endif
#endif

#if !defined(__arm64__)

#define add_sssaaaaaa(sh, sm, sl, ah, am, al, bh, bm, bl)           \
  do {                                                              \
    mp_limb_t __t, __u;                                             \
    add_ssaaaa(__t, sl, (mp_limb_t) 0, al, (mp_limb_t) 0, bl);      \
    add_ssaaaa(__u, sm, (mp_limb_t) 0, am, (mp_limb_t) 0, bm);      \
    add_ssaaaa(sh, sm, ah + bh, sm, __u, __t);                      \
  } while (0)

#define add_ssssaaaaaaaa(s3, s2, s1, s0, a3, a2, a1, a0, b3, b2, b1, b0)       \
  do {                                                                         \
    mp_limb_t __tt;                                                            \
    add_sssaaaaaa(__tt, s1, s0, (mp_limb_t) 0, a1, a0, (mp_limb_t) 0, b1, b0); \
    add_ssaaaa(s3, s2, a3, a2, b3, b2);                                        \
    add_ssaaaa(s3, s2, s3, s2, (mp_limb_t) 0, __tt);                           \
  } while (0)

#endif


#if !((GMP_LIMB_BITS == 64 && defined (__ia64)) ||      \
      (GMP_LIMB_BITS == 32 && defined (__arm__)))
#if !defined(__arm64__)

#define sub_ddmmss(sh, sl, ah, al, bh, bl) \
  do {									          \
    mp_limb_t __x;								 \
    __x = (al) - (bl);							 \
    (sh) = (ah) - (bh) - ((al) < (bl));    \
    (sl) = __x;								    \
  } while (0)

#endif
#endif

#if !defined(__arm64__)
#define sub_dddmmmsss(dh, dm, dl, mh, mm, ml, sh, sm, sl)           \
  do {                                                              \
    mp_limb_t __t, __u;                                             \
    sub_ddmmss(__t, dl, (mp_limb_t) 0, ml, (mp_limb_t) 0, sl);      \
    sub_ddmmss(__u, dm, (mp_limb_t) 0, mm, (mp_limb_t) 0, sm);      \
    sub_ddmmss(dh, dm, (mh) - (sh), dm, -__u, -__t);                \
  } while (0)
#endif

/* MIPS and ARM - Use clz builtins */
/* NOTE: Apple clang version 12.0.5 miscompiles the count_leading_zeros fallback */
#if (defined (__mips__) || defined (__arm__) || defined (__arm64__))

#ifdef _LONG_LONG_LIMB
#define count_leading_zeros(count,x)            \
  do {                                          \
    FLINT_ASSERT ((x) != 0);                    \
    (count) = __builtin_clzll (x);              \
  } while (0)
#else
#define count_leading_zeros(count,x)            \
  do {                                          \
    FLINT_ASSERT ((x) != 0);                    \
    (count) = __builtin_clzl (x);               \
  } while (0)
#endif

#ifdef _LONG_LONG_LIMB
#define count_trailing_zeros(count,x)           \
  do {                                          \
    FLINT_ASSERT ((x) != 0);                    \
    (count) = __builtin_ctzll (x);              \
  } while (0)
#else
#define count_trailing_zeros(count,x)           \
  do {                                          \
    FLINT_ASSERT ((x) != 0);                    \
    (count) = __builtin_ctzl (x);               \
  } while (0)
#endif

#endif /* MIPS, ARM */

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

#ifndef count_leading_zeros
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
#endif

#if !(GMP_LIMB_BITS == 64 && defined (__ia64))

#ifndef count_trailing_zeros
#define count_trailing_zeros(count, x)                 \
  do {									                      \
    mp_limb_t __ctz_x = (x);						          \
    mp_limb_t __ctz_c;							             \
    FLINT_ASSERT (__ctz_x != 0);						       \
    count_leading_zeros (__ctz_c, __ctz_x & -__ctz_x); \
    (count) = GMP_LIMB_BITS - 1 - __ctz_c;	 \
  } while (0)
#endif

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

#define sdiv_qrnnd(q, r, n1, n0, d)         \
  do {                                      \
    mp_limb_t __n1, __n0, __d;              \
    mp_limb_t __q, __r;                     \
    unsigned int __sgn_n = 0, __sgn_d = 0;  \
    if ((n1) & __highbit)                   \
    {                                       \
       __n0 = -(n0);                        \
       __n1 = ~(n1) + (__n0 == 0);          \
       __sgn_n = ~__sgn_n;                  \
    } else                                  \
    {                                       \
      __n0 = (n0);                          \
      __n1 = (n1);                          \
    }                                       \
    if ((d) & __highbit)                    \
    {                                       \
        __d = -(d);                         \
        __sgn_d = ~__sgn_d;                 \
    } else                                  \
    {                                       \
        __d = (d);                          \
    }                                       \
    udiv_qrnnd(__q, __r, __n1, __n0, __d);  \
    q = (__sgn_n == __sgn_d) ? __q : -__q;  \
    r = (__sgn_n == 0) ? __r : -__r;        \
  } while (0)

#if GMP_LIMB_BITS == 32

#define byte_swap(n)                                           \
  do {                                                         \
      /* swap adjacent bytes */                                \
      n = (((n & 0xff00ff00) >> 8) | ((n & 0x00ff00ff) << 8)); \
      /* swap adjacent words */                                \
      n = ((n >> 16) | (n << 16));                             \
  } while (0)

#else

#define byte_swap(n)                                                             \
  do {                                                                           \
      /* swap adjacent bytes */                                                  \
      n = (((n & 0xff00ff00ff00ff00) >> 8) | ((n & 0x00ff00ff00ff00ff) << 8));   \
      /* swap adjacent words */                                                  \
      n = (((n & 0xffff0000ffff0000) >> 16) | ((n & 0x0000ffff0000ffff) << 16)); \
      /* swap adjacent double words */                                           \
      n = ((n >> 32) | (n << 32));                                               \
  } while (0)

#endif /* 64 bits */

#endif /* non x86 fallback code */

/* smul_ppm is defined previously for 32bit arm and for all x86 */
#if !( (GMP_LIMB_BITS == 32 && defined (__arm__))                              \
       || defined (__i386__) || defined (__i486__) || defined(__amd64__))

#if !defined(__arm64__)

#define smul_ppmm(w1, w0, u, v)                         \
  do {                                                  \
    mp_limb_t __w1;                                     \
    mp_limb_t __xm0 = (u), __xm1 = (v);                 \
    umul_ppmm (__w1, w0, __xm0, __xm1);                 \
    (w1) = __w1 - (-(__xm0 >> (FLINT_BITS-1)) & __xm1)  \
        - (-(__xm1 >> (FLINT_BITS-1)) & __xm0);         \
  } while (0)

#endif

#endif

#if GMP_LIMB_BITS == 32

/* rec_word_tab[i] = div(2^24 - 2^14 + 2^9, 2^9 + i) */
static const int rec_word_tab[512] = {
   32737, 32673, 32609, 32546, 32483, 32420, 32357, 32295, 32233, 32171, 32109, 32048, 31987, 31926, 31865, 31805,
   31744, 31684, 31625, 31565, 31506, 31447, 31388, 31329, 31271, 31212, 31154, 31097, 31039, 30982, 30924, 30868,
   30811, 30754, 30698, 30642, 30586, 30530, 30475, 30419, 30364, 30309, 30255, 30200, 30146, 30092, 30038, 29984,
   29930, 29877, 29824, 29771, 29718, 29666, 29613, 29561, 29509, 29457, 29405, 29354, 29303, 29251, 29200, 29150,
   29099, 29049, 28998, 28948, 28898, 28849, 28799, 28750, 28700, 28651, 28602, 28554, 28505, 28457, 28409, 28360,
   28313, 28265, 28217, 28170, 28123, 28075, 28029, 27982, 27935, 27889, 27842, 27796, 27750, 27704, 27658, 27613,
   27568, 27522, 27477, 27432, 27387, 27343, 27298, 27254, 27209, 27165, 27121, 27078, 27034, 26990, 26947, 26904,
   26861, 26818, 26775, 26732, 26690, 26647, 26605, 26563, 26521, 26479, 26437, 26395, 26354, 26312, 26271, 26230,
   26189, 26148, 26108, 26067, 26026, 25986, 25946, 25906, 25866, 25826, 25786, 25747, 25707, 25668, 25628, 25589,
   25550, 25511, 25473, 25434, 25395, 25357, 25319, 25281, 25242, 25205, 25167, 25129, 25091, 25054, 25016, 24979,
   24942, 24905, 24868, 24831, 24794, 24758, 24721, 24685, 24649, 24612, 24576, 24540, 24504, 24469, 24433, 24397,
   24362, 24327, 24291, 24256, 24221, 24186, 24151, 24117, 24082, 24047, 24013, 23979, 23944, 23910, 23876, 23842,
   23808, 23774, 23741, 23707, 23674, 23640, 23607, 23574, 23541, 23508, 23475, 23442, 23409, 23377, 23344, 23312,
   23279, 23247, 23215, 23183, 23151, 23119, 23087, 23055, 23023, 22992, 22960, 22929, 22898, 22866, 22835, 22804,
   22773, 22742, 22711, 22681, 22650, 22619, 22589, 22559, 22528, 22498, 22468, 22438, 22408, 22378, 22348, 22318,
   22289, 22259, 22229, 22200, 22171, 22141, 22112, 22083, 22054, 22025, 21996, 21967, 21938, 21910, 21881, 21853,
   21824, 21796, 21767, 21739, 21711, 21683, 21655, 21627, 21599, 21571, 21544, 21516, 21488, 21461, 21433, 21406,
   21379, 21352, 21324, 21297, 21270, 21243, 21216, 21190, 21163, 21136, 21110, 21083, 21056, 21030, 21004, 20977,
   20951, 20925, 20899, 20873, 20847, 20821, 20795, 20769, 20744, 20718, 20693, 20667, 20642, 20616, 20591, 20566,
   20540, 20515, 20490, 20465, 20440, 20415, 20390, 20366, 20341, 20316, 20292, 20267, 20243, 20218, 20194, 20170,
   20145, 20121, 20097, 20073, 20049, 20025, 20001, 19977, 19953, 19930, 19906, 19882, 19859, 19835, 19812, 19789,
   19765, 19742, 19719, 19696, 19672, 19649, 19626, 19603, 19581, 19558, 19535, 19512, 19489, 19467, 19444, 19422,
   19399, 19377, 19354, 19332, 19310, 19288, 19265, 19243, 19221, 19199, 19177, 19155, 19133, 19112, 19090, 19068,
   19046, 19025, 19003, 18982, 18960, 18939, 18917, 18896, 18875, 18854, 18832, 18811, 18790, 18769, 18748, 18727,
   18706, 18686, 18665, 18644, 18623, 18603, 18582, 18561, 18541, 18520, 18500, 18479, 18459, 18439, 18419, 18398,
   18378, 18358, 18338, 18318, 18298, 18278, 18258, 18238, 18218, 18199, 18179, 18159, 18139, 18120, 18100, 18081,
   18061, 18042, 18022, 18003, 17984, 17964, 17945, 17926, 17907, 17888, 17869, 17850, 17831, 17812, 17793, 17774,
   17755, 17736, 17718, 17699, 17680, 17662, 17643, 17624, 17606, 17587, 17569, 17551, 17532, 17514, 17496, 17477,
   17459, 17441, 17423, 17405, 17387, 17369, 17351, 17333, 17315, 17297, 17279, 17261, 17244, 17226, 17208, 17191,
   17173, 17155, 17138, 17120, 17103, 17085, 17068, 17051, 17033, 17016, 16999, 16982, 16964, 16947, 16930, 16913,
   16896, 16879, 16862, 16845, 16828, 16811, 16794, 16778, 16761, 16744, 16727, 16711, 16694, 16677, 16661, 16644,
   16628, 16611, 16595, 16578, 16562, 16546, 16529, 16513, 16497, 16481, 16464, 16448, 16432, 16416, 16400, 16384
};

#define invert_limb(dinv, d)                                      \
   do {                                                           \
      mp_limb_t _v0, _v1, _v2, _d21, _e, _m0;                     \
      FLINT_ASSERT(((d) & (UWORD(1)<<(GMP_LIMB_BITS - 1))) != 0); \
      _v0 = rec_word_tab[((d) >> 22) & 0x1FF];                    \
      _d21 = ((d) >> 11) + 1;                                     \
      _m0 = _v0*_v0;                                              \
      umul_ppmm(_v1, _e, _m0, _d21);                              \
      _v1 = (_v0 << 4) - _v1 - 1;                                 \
      _e = -_v1*((d) >> 1);                                       \
      _m0 = -((d) & (mp_limb_t) 1);                               \
      _e -= ((_v1 - (_v1 >> 1)) & _m0);                           \
      umul_ppmm(_v2, _m0, _v1, _e);                               \
      _v2 = (_v1 << 15) + (_v2 >> 1);                             \
      umul_ppmm(_v0, _d21, _v2, (d));                             \
      add_ssaaaa(_v0, _d21, _v0, _d21, (mp_limb_t) 0, (d));       \
      (dinv) = _v2 - (_v0 + (d));                                 \
   } while(0)

#else

/* rec_word_tab[i] = div(2^19 - 3*2^8, 2^8 + i) */
static const int rec_word_tab[256] = {
   2045, 2037, 2029, 2021, 2013, 2005, 1998, 1990, 1983, 1975, 1968, 1960, 1953, 1946, 1938, 1931,
   1924, 1917, 1910, 1903, 1896, 1889, 1883, 1876, 1869, 1863, 1856, 1849, 1843, 1836, 1830, 1824,
   1817, 1811, 1805, 1799, 1792, 1786, 1780, 1774, 1768, 1762, 1756, 1750, 1745, 1739, 1733, 1727,
   1722, 1716, 1710, 1705, 1699, 1694, 1688, 1683, 1677, 1672, 1667, 1661, 1656, 1651, 1646, 1641,
   1636, 1630, 1625, 1620, 1615, 1610, 1605, 1600, 1596, 1591, 1586, 1581, 1576, 1572, 1567, 1562,
   1558, 1553, 1548, 1544, 1539, 1535, 1530, 1526, 1521, 1517, 1513, 1508, 1504, 1500, 1495, 1491,
   1487, 1483, 1478, 1474, 1470, 1466, 1462, 1458, 1454, 1450, 1446, 1442, 1438, 1434, 1430, 1426,
   1422, 1418, 1414, 1411, 1407, 1403, 1399, 1396, 1392, 1388, 1384, 1381, 1377, 1374, 1370, 1366,
   1363, 1359, 1356, 1352, 1349, 1345, 1342, 1338, 1335, 1332, 1328, 1325, 1322, 1318, 1315, 1312,
   1308, 1305, 1302, 1299, 1295, 1292, 1289, 1286, 1283, 1280, 1276, 1273, 1270, 1267, 1264, 1261,
   1258, 1255, 1252, 1249, 1246, 1243, 1240, 1237, 1234, 1231, 1228, 1226, 1223, 1220, 1217, 1214,
   1211, 1209, 1206, 1203, 1200, 1197, 1195, 1192, 1189, 1187, 1184, 1181, 1179, 1176, 1173, 1171,
   1168, 1165, 1163, 1160, 1158, 1155, 1153, 1150, 1148, 1145, 1143, 1140, 1138, 1135, 1133, 1130,
   1128, 1125, 1123, 1121, 1118, 1116, 1113, 1111, 1109, 1106, 1104, 1102, 1099, 1097, 1095, 1092,
   1090, 1088, 1086, 1083, 1081, 1079, 1077, 1074, 1072, 1070, 1068, 1066, 1064, 1061, 1059, 1057,
   1055, 1053, 1051, 1049, 1047, 1044, 1042, 1040, 1038, 1036, 1034, 1032, 1030, 1028, 1026, 1024
};

#define invert_limb(dinv, d)                                      \
   do {                                                           \
      mp_limb_t _v0, _v2, _d40, _e, _m0;                          \
      FLINT_ASSERT(((d) & (UWORD(1)<<(GMP_LIMB_BITS - 1))) != 0); \
      _d40 = ((d) >> 24) + 1;                                     \
      _v0 = rec_word_tab[((d) >> 55) & 0xFF];                     \
      _v0 = (_v0 << 11) - ((_v0*_v0*_d40) >> 40) - 1;             \
      _v2 = ((_v0*((((mp_limb_t) 1) << 60) - _v0*_d40)) >> 47);   \
      _v2 += (_v0 << 13);                                         \
      _e = -_v2*((d) >> 1);                                       \
      _m0 = -((d) & (mp_limb_t) 1);                               \
      _e -= ((_v2 - (_v2 >> 1)) & _m0);                           \
      umul_ppmm(_v0, _d40, _v2, _e);                              \
      _v2 = (_v2 << 31) + (_v0 >> 1);                             \
      umul_ppmm(_v0, _d40, _v2, (d));                             \
      add_ssaaaa(_v0, _d40, _v0, _d40, (mp_limb_t) 0, (d));       \
      (dinv) = _v2 - (_v0 + (d));                                 \
   } while (0)

#endif /* 64 bits */

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
