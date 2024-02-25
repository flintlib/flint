/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h> /* For memcpy */
#include "mpn_extras.h"

/*

We will now define what we consider as a high multiplication, and how we go from
basecase to bigger cases where Toom-Cook or Schönhage-Strassen multiplication is
required to outperform `flint_mpn_mul'.

Let {a, n} and {b, n} be two n-limbed (positive) integers, of which the product
of a and b is {c, 2 n}. For some applications, only the higher part is required,
and so we only calculate an "approximation" of the most significant n + 1 limbs.
With a sloppy wording, what this means is that we only calculate the part of the
multiplication that contributes to the most significant n + 1 limbs, so carries
are disregarded. This result in that the approximation of c[n - 1, ..., 2 n - 1]
is smaller than the real value, and that the error is at most ~n ULPs in the
least significant limb in the approximation.

With {c, 2 n} denoting the product of {a, n} times {b, n}, let {d, n + 1} denote
the high product. We visualise the high multiplication of two 10-limbed integers
with the following:

                               0 1 2 3 4 5 6 7 8 9
                             0                 h x
                             1               h x x
                             2             h x x x
                             3           h x x x x
                             4         h x x x x x
                             5       h x x x x x x
                             6     h x x x x x x x
                             7   h x x x x x x x x
                             8 h x x x x x x x x x
                             9 x x x x x x x x x x

Here `h' means that only the higher part of this entry was calculated, and `x'
means that the full product of the limbs where calculated.

To utilise multiplication algorithms that exploits symmetries, we divide this
figure into four different parts:

                               0 1 2 3 4 5 6 7 8
                             0          |    h x
                             1          |  h x x
                             2          |h x x x
                             3  _ _ _ _h|x_x_x_x
     n = 9:                  4       h x|x x x x
                             5     h x x|x x x x
                             6   h x x x|x x x x
                             7 h x x x x|x x x x
                             8 x x x x x|x x x x

                              0 1 2 3 4 5 6 7 8 9
                            0          |      h x
                            1          |    h x x
                            2          |  h x x x
                            3          |h x x x x
     n = 10:                4  _ _ _ _h|x_x_x_x_x
                            5       h x|x x x x x
                            6     h x x|x x x x x
                            7   h x x x|x x x x x
                            8 h x x x x|x x x x x
                            9 x x x x x|x x x x x

Observe that we have only one multi-limbed full multiplication, two multi-limbed
high multiplications and one single-limbed high multiplication.

*/

#if FLINT_HAVE_NATIVE_MPN_MULHIGH_BASECASE && FLINT_HAVE_NATIVE_2ADD_N_INPLACE

#if !defined(__amd64__)
# error
#endif

/* NOTE: As we will not reuse factors in mulhigh, we utilize mul instead of mulx
 * to save a few bytes. */
#define mulhigh(p, u, v) \
  do { \
    ulong _scr; \
    __asm__("mulq\t%3" \
      : "=a" (_scr), "=d" (p) \
      : "%0" ((ulong)(u)), "rm" ((ulong)(v))); \
  } while (0)

/* NOTE: Assumes no carry */
#define flint_mpn_add_1(rp, x) \
  do { \
    ulong __rp_save = (rp)[0]; \
    (rp)[0] += (x); \
    if (__rp_save > (rp)[0]) \
    { \
      slong __ix = 0; \
      do \
      { \
        __ix++; \
        (rp)[__ix] += 1; \
      } while ((rp)[__ix] == UWORD(0)); \
    } \
  } while (0)

#define RECURSIVE_THRESHOLD 59
#define _RECURSIVE_THRESHOLD 47
#define FALLBACK_THRESHOLD 330

FLINT_STATIC_NOINLINE
mp_limb_t _flint_mpn_mulhigh_rec(mp_ptr rp, mp_srcptr xp, mp_srcptr yp, mp_size_t n, mp_ptr scr)
{
    if (n < _RECURSIVE_THRESHOLD)
        return _flint_mpn_mulhigh_basecase(rp, xp, yp, n);
    else
    {
        mp_size_t np1o2 = (n + 1) / 2;
        mp_size_t no2 = n / 2;
        mp_limb_t c0, c1 = 0, ret;
        mp_ptr hl, hr;

        /* Top left */
        mulhigh(c0, xp[no2 - 1], yp[np1o2 - 1]);

        /* Bottom right */
        _flint_mpn_mul(rp, xp + no2, np1o2, yp + np1o2, no2);

        /* Bottom left */
        hr = scr;
        ret = _flint_mpn_mulhigh_rec(hr, xp + no2, yp, np1o2, hr + np1o2);
        add_ssaaaa(c1, c0, c1, c0, 0, ret);

        /* Top right */
        hl = scr + np1o2;
        hl[np1o2 - 1] = 0;
        ret = _flint_mpn_mulhigh_rec(hl, xp, yp + np1o2, no2, hl + np1o2);
        add_ssaaaa(c1, c0, c1, c0, 0, ret);

        /* Add c1 to rp */
        flint_mpn_add_1(rp, c1);

        /* Add both high multiplications to rp */
        ret = flint_mpn_2add_n_inplace(rp, hr, hl, np1o2);

        /* Add carry from addition to rp */
        flint_mpn_add_1(rp + np1o2, ret);

        return c0;
    }
}

mp_limb_t _flint_mpn_mulhigh(mp_ptr rp, mp_srcptr xp, mp_srcptr yp, mp_size_t n)
{
    FLINT_ASSERT(n > FLINT_MPN_MULHIGH_FUNC_TAB_WIDTH);

    if (n < RECURSIVE_THRESHOLD)
    {
        FLINT_ASSERT(rp != xp && rp != yp);
        return _flint_mpn_mulhigh_basecase(rp, xp, yp, n);
    }
    else if (n < FALLBACK_THRESHOLD)
    {
        mp_limb_t ret;
        mp_ptr scr;

        FLINT_ASSERT(rp != xp && rp != yp);

        scr = flint_malloc(2 * sizeof(mp_limb_t) * n);
        ret = _flint_mpn_mulhigh_rec(rp, xp, yp, n, scr);
        flint_free(scr);

        return ret;
    }
    else
    {
        /* Aliasing is okay */
        mp_ptr tmp;
        mp_limb_t ret;

        tmp = flint_malloc(2 * sizeof(mp_limb_t) * n);

        _flint_mpn_mul_n(tmp, xp, yp, n);
        memcpy(rp, tmp + n, sizeof(mp_limb_t) * n);
        ret = tmp[n - 1];

        flint_free(tmp);

        return ret;
    }
}
#else
typedef int this_file_is_empty;
#endif
