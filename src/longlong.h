/*
   Copyright 1991, 1992, 1993, 1994, 1996, 1997, 1999, 2000, 2001, 2002, 2003,
   2004, 2005 Free Software Foundation, Inc.

   Copyright 2009, 2015, 2016 William Hart
   Copyright 2011 Fredrik Johansson
   Copyright 2023 Albin Ahlbäck

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

#define count_leading_zeros _Pragma("GCC error \"'count_leading_zeros' is deprecated. Use 'flint_clz' instead.\"")
#define count_trailing_zeros _Pragma("GCC error \"'count_trailing_zeros' is deprecated. Use 'flint_ctz' instead.\"")

#ifdef FLINT_HAS_POPCNT
# ifndef _LONG_LONG_LIMB
#  define flint_popcount __builtin_popcountl
# else
#  define flint_popcount __builtin_popcountll
# endif
#endif

#ifdef FLINT_HAS_CLZ
# ifndef _LONG_LONG_LIMB
#  define flint_clz __builtin_clzl
# else
#  define flint_clz __builtin_clzll
# endif
#endif

#ifdef FLINT_HAS_CTZ
# ifndef _LONG_LONG_LIMB
#  define flint_ctz __builtin_ctzl
# else
#  define flint_ctz __builtin_ctzll
# endif
#endif

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

#ifndef FLINT_HAS_CLZ
# define flint_clz flint_clz
static __inline__ flint_bitcnt_t flint_clz(mp_limb_t x)
{
    mp_limb_t count;
    FLINT_ASSERT(x != 0);
    __asm__("bsrq %1,%0" : "=r" (count) : "rm" (x));
    return count ^ (mp_limb_t) 63;
}
#endif

#ifndef FLINT_HAS_CTZ
# define flint_ctz flint_ctz
static __inline__ flint_bitcnt_t flint_ctz(mp_limb_t x)
{
    mp_limb_t count;
    FLINT_ASSERT(x != 0);
    __asm__("bsfq %1,%0" : "=r" (count) : "rm" (x));
    return count;
}
#endif

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

#ifndef FLINT_HAS_CLZ
# define flint_clz flint_clz
static __inline__ flint_bitcnt_t flint_clz(mp_limb_t x)
{
    mp_limb_t count;
    FLINT_ASSERT(x != 0);
    __asm__("bsrl %1,%0" : "=r" (count) : "rm" (x));
    return count ^ (mp_limb_t) 31;
}
#endif

#ifndef FLINT_HAS_CTZ
# define flint_ctz flint_ctz
static __inline__ flint_bitcnt_t flint_ctz(mp_limb_t x)
{
    mp_limb_t count;
    FLINT_ASSERT(x != 0);
    __asm__("bsfl %1,%0" : "=r" (count) : "rm" (x));
    return count;
}
#endif

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
    if ((mp_limb_t) (al) < (mp_limb_t) (bl))    \
      (sh) = (ah) - (bh) - 1;                   \
    else                                        \
      (sh) = (ah) - (bh);                       \
    (sl) = __x;                                 \
  } while (0)

#ifndef FLINT_HAS_CTZ
# define flint_ctz(x)                       \
({                                          \
    mp_limb_t __ctz_x = (x);                \
    mp_limb_t count;                        \
    __asm__ ("popcnt %0 = %1"               \
            : "=r" (count)                  \
            : "r" ((__ctz_x-1) & ~__ctz_x));\
    count;                                  \
})
#endif

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
    (sh) = (ah) - (bh) - ((mp_limb_t) (al) < (mp_limb_t) (bl));    \
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

FLINT_DLL extern const unsigned char __flint_clz_tab[128];

#ifndef flint_clz
# define flint_clz flint_clz
static __inline__ flint_bitcnt_t flint_clz(mp_limb_t x)
{
    mp_limb_t __a, __xr = x;
    if (GMP_LIMB_BITS == 32)
        __a = __xr < ((mp_limb_t) 1 << 2*__BITS4)
        ? (__xr < ((mp_limb_t) 1 << __BITS4) ? 1 : __BITS4 + 1)
        : (__xr < ((mp_limb_t) 1 << 3*__BITS4) ? 2*__BITS4 + 1
        : 3*__BITS4 + 1);
    else
    {
        for (__a = GMP_LIMB_BITS - 8; __a > 0; __a -= 8)
            if (((__xr >> __a) & 0xff) != 0)
                break;
        ++__a;
    }
    return GMP_LIMB_BITS + 1 - __a - __flint_clz_tab[__xr >> __a];
}
#endif

#if !(GMP_LIMB_BITS == 64 && defined (__ia64))
# ifndef flint_ctz
#  define flint_ctz flint_ctz
static __inline__ flint_bitcnt_t flint_ctz(mp_limb_t x)
{
    mp_limb_t __ctz_x = (x);
    FLINT_ASSERT (__ctz_x != 0);
    return GMP_LIMB_BITS - 1 - flint_clz(__ctz_x & -__ctz_x);
}
# endif
#endif

#define udiv_qrnnd(q, r, n1, n0, d)                  \
    do {                                             \
       mp_limb_t __norm = flint_clz(d);              \
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
