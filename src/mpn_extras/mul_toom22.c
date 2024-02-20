/*
    Copyright 2006-2010, 2012, 2014, 2018, 2020 Free Software Foundation, Inc.
    Contributed to the GNU project by Torbjorn Granlund.

    Copyright 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
   mpn_toom22_mul -- Multiply {ap,an} and {bp,bn} where an >= bn.  Or more
   accurately, bn <= an < 2bn.

   Contributed to the GNU project by Torbjorn Granlund.

   THE FUNCTION IN THIS FILE IS INTERNAL WITH A MUTABLE INTERFACE.  IT IS ONLY
   SAFE TO REACH IT THROUGH DOCUMENTED INTERFACES.  IN FACT, IT IS ALMOST
   GUARANTEED THAT IT WILL CHANGE OR DISAPPEAR IN A FUTURE GNU MP RELEASE.

Copyright 2006-2010, 2012, 2014, 2018, 2020 Free Software Foundation, Inc.

This file is part of the GNU MP Library.

The GNU MP Library is free software; you can redistribute it and/or modify
it under the terms of either:

  * the GNU Lesser General Public License as published by the Free
    Software Foundation; either version 3 of the License, or (at your
    option) any later version.

or

  * the GNU General Public License as published by the Free Software
    Foundation; either version 2 of the License, or (at your option) any
    later version.

or both in parallel, as here.

The GNU MP Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received copies of the GNU General Public License and the
GNU Lesser General Public License along with the GNU MP Library.  If not,
see https://www.gnu.org/licenses/.  */

#include "flint.h"
#include "mpn_extras.h"

/* Unoptimized */
#define MPN_INCR_U(ptr, size, n)   mpn_add_1((ptr), (ptr), (size), (n))

/* Evaluate in: -1, 0, +inf

  <-s--><--n-->
   ____ ______
  |_a1_|___a0_|
   |b1_|___b0_|
   <-t-><--n-->

  v0  =  a0     * b0       #   A(0)*B(0)
  vm1 = (a0- a1)*(b0- b1)  #  A(-1)*B(-1)
  vinf=      a1 *     b1   # A(inf)*B(inf)
*/

#define TOOM22_MUL_N_REC(p, a, b, n, ws) \
    do { \
        if (FLINT_HAVE_MUL_N_FUNC(n)) \
            flint_mpn_mul_n_func_tab[n](p, a, b); \
        else \
            flint_mpn_mul_toom22(p, a, n, b, n, ws); \
    } while (0)

/* Normally, this calls mul_basecase or toom22_mul.  But when when the fraction
   MUL_TOOM33_THRESHOLD / MUL_TOOM22_THRESHOLD is large, an initially small
   relative unbalance will become a larger and larger relative unbalance with
   each recursion (the difference s-t will be invariant over recursive calls).
   Therefore, we need to call toom32_mul.  FIXME: Suppress depending on
   MUL_TOOM33_THRESHOLD / MUL_TOOM22_THRESHOLD and on MUL_TOOM22_THRESHOLD.  */
#define TOOM22_MUL_REC(p, a, an, b, bn, ws) \
    do { \
        if (FLINT_HAVE_MUL_FUNC(an, bn)) \
            flint_mpn_mul_func_tab[an][bn](p, a, b); \
        else if (bn <= 12 || 4 * an >= 5 * bn) \
            flint_mpn_mul(p, a, an, b, bn); \
        else \
            flint_mpn_mul_toom22(p, a, an, b, bn, ws); \
    } while (0)

/*
    else if (4 * an < 5 * bn)                        \
      mpn_toom22_mul (p, a, an, b, bn, ws);                \
    else                                \
      mpn_toom32_mul (p, a, an, b, bn, ws);                \
*/

/* Needed scratch space is 2*(an + k), where k is the recursion depth. */

void
flint_mpn_mul_toom22(mp_ptr pp,
        mp_srcptr ap, mp_size_t an,
        mp_srcptr bp, mp_size_t bn,
        mp_ptr scratch)
{
    mp_size_t n, s, t;
    int vm1_neg;
    mp_limb_t cy, cy2;
    mp_ptr asm1;
    mp_ptr bsm1;
    TMP_INIT;

    TMP_START;
    if (scratch == NULL)
    {
        scratch = TMP_ALLOC(sizeof(mp_limb_t) * (2 * ((an) + FLINT_BITS)));
    }

#define a0  ap
#define a1  (ap + n)
#define b0  bp
#define b1  (bp + n)

    s = an >> 1;
    n = an - s;
    t = bn - n;

    FLINT_ASSERT(an >= bn);
    FLINT_ASSERT(0 < s && s <= n && (n - s) == (an & 1));
    FLINT_ASSERT(0 < t && t <= s);

    asm1 = pp;
    bsm1 = pp + n;

    vm1_neg = 0;

    /* Compute asm1.  */
    if ((an & 1) == 0) /* s == n */
    {
        if (mpn_cmp(a0, a1, n) < 0)
        {
            mpn_sub_n(asm1, a1, a0, n);
            vm1_neg = 1;
        }
        else
        {
            mpn_sub_n(asm1, a0, a1, n);
        }
    }
    else /* n - s == 1 */
    {
        if (a0[s] == 0 && mpn_cmp(a0, a1, s) < 0)
        {
            mpn_sub_n(asm1, a1, a0, s);
            asm1[s] = 0;
            vm1_neg = 1;
        }
        else
        {
            asm1[s] = a0[s] - mpn_sub_n(asm1, a0, a1, s);
        }
    }

    /* Compute bsm1.  */
    if (t == n)
    {
        if (mpn_cmp(b0, b1, n) < 0)
        {
            mpn_sub_n(bsm1, b1, b0, n);
            vm1_neg ^= 1;
        }
        else
        {
            mpn_sub_n(bsm1, b0, b1, n);
        }
    }
    else
    {
        if (flint_mpn_zero_p(b0 + t, n - t) && mpn_cmp(b0, b1, t) < 0)
        {
            mpn_sub_n(bsm1, b1, b0, t);
            flint_mpn_zero(bsm1 + t, n - t);
            vm1_neg ^= 1;
        }
        else
        {
            mpn_sub(bsm1, b0, n, b1, t);
        }
    }

#define v0              pp                /* 2n */
#define vinf           (pp + 2 * n)       /* s+t */
#define vm1            scratch            /* 2n */
#define scratch_out    scratch + 2 * n

    /* vm1, 2n limbs */
    TOOM22_MUL_N_REC (vm1, asm1, bsm1, n, scratch_out);

    if (s > t)
        TOOM22_MUL_REC(vinf, a1, s, b1, t, scratch_out);
    else
        TOOM22_MUL_N_REC(vinf, a1, b1, s, scratch_out);

    /* v0, 2n limbs */
    TOOM22_MUL_N_REC(v0, ap, bp, n, scratch_out);

    /* H(v0) + L(vinf) */
    cy = mpn_add_n(pp + 2 * n, v0 + n, vinf, n);

    /* L(v0) + (H(v0) + L(vinf)) */
    cy2 = cy + mpn_add_n(pp + n, pp + 2 * n, v0, n);

    /* (H(v0) + L(vinf)) + H(vinf) */
    cy += mpn_add(pp + 2 * n, pp + 2 * n, n, vinf + n, s + t - n);

    if (vm1_neg)
    {
        cy += mpn_add_n(pp + n, pp + n, vm1, 2 * n);
    }
    else
    {
        cy -= mpn_sub_n(pp + n, pp + n, vm1, 2 * n);

        if (FLINT_UNLIKELY(cy + 1 == 0)) /* cy is negative */
        {
            /* The total contribution of v0+vinf-vm1 can not be negative. */
#if FLINT_WANT_ASSERT
            /* The borrow in cy stops the propagation of the carry cy2, */
            FLINT_ASSERT(cy2 == 1);
            cy += mpn_add_1(pp + 2 * n, pp + 2 * n, n, cy2);
            FLINT_ASSERT (cy == 0);
#else
            /* we simply fill the area with zeros. */
            flint_mpn_zero(pp + 2 * n, n);
            FLINT_ASSERT(s + t == n || flint_mpn_zero_p(pp + 3 * n, s + t - n));
#endif
            goto cleanup;
        }
    }

    FLINT_ASSERT(cy  <= 2);
    FLINT_ASSERT(cy2 <= 2);

    MPN_INCR_U(pp + 2 * n, s + t, cy2);
    /* if s+t==n, cy is zero, but we should not access pp[3*n] at all. */
    MPN_INCR_U(pp + 3 * n, s + t - n, cy);

cleanup:
    TMP_END;
}
