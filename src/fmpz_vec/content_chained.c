/*
    Copyright 1994, 1996, 2000, 2001, 2009, 2012, 2019 Free Software Foundation, Inc.
    Copyright (C) 2020 Daniel Schultz
    Copyright (C) 2021 Fredrik Johansson
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    Contains code from mpn_gcd_1 in GMP 6.3.0.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"

/* Calculating the GCD is a very costly operation.
  
   Hence, our attempt is to find a single-limbed integer of which we can start
   the process with, because calculating the GCD of {u, n} and {v, 1} is much
   more effective than calculating if of two n-limbed integers.

   Initial strategy to find a single-limbed integer:
    * `ip` will often be small -- start by checking its total size.
    * `ip` may contain factors of two -- check if reducing it leads to a single
      limb.
    * `vp` may contain single-limbed integers at specific places -- imagine it
      representing a polynomial; many polynomials are small at the ends.
*/

#define SIZ(x) ((x)->_mp_size)
#define PTR(x) ((x)->_mp_d)

typedef struct
{
    mp_srcptr d;
    ulong sz;
}
mpn_struct;

typedef mpn_struct * mpn_ptr;

FLINT_STATIC_NOINLINE void heavy_machinery(fmpz_t, const fmpz *, slong, const fmpz_t);

/* In this function we only take care of small inputs. If none is found, we go
 * to heavy_machinery instead. */
void
_fmpz_vec_content_chained(fmpz_t rp, const fmpz * vp, slong vn, const fmpz_t ip)
{
    slong jx = -1, kx;
    ulong gd; /* single limb gcd */
    int exp;
    mpz_srcptr mp;
    slong msz;
    mp_srcptr md;

#if FLINT_WANT_ASSERT
    gd = 0;
#endif

    /* We will assume that vn >= 2 for everything that will come */
    if (vn < WORD(2))
    {
        if (ip == NULL)
            fmpz_abs(rp, vp + 0);
        else
            fmpz_gcd(rp, vp + 0, ip);

        return;
    }

    /* Check if `ip` is a single limb */
    if (ip != NULL)
    {
        fmpz tip = *ip; /* Read-only */

        if (!COEFF_IS_MPZ(tip))
        {
            if (tip == WORD(0))
            {
                ip = NULL;
                goto check_vp_may_be_zero;
            }

            gd = FLINT_ABS(tip);
            goto ip1;
        }
        else
        {
            mp = COEFF_TO_PTR(tip);
            msz = FLINT_ABS(SIZ(mp));

            if (msz == 1)
            {
                gd = PTR(mp)[0];
                goto ip1;
            }
        }
    }

check_vp_may_be_zero: /* Fast checks only. Result can be zero. */
    for (jx = 0; jx < vn; jx++)
    {
        if (!fmpz_is_zero(vp + jx))
            goto check_vp_not_zero;
        else if (!COEFF_IS_MPZ(vp[jx]))
            goto found_small;
    }

    if (ip == NULL)
        fmpz_zero(rp);
    else
        fmpz_abs(rp, ip);

    return;

check_vp_not_zero: /* Fast checks only. Result cannot be zero from this point onwards. */
    for (jx++; jx < vn; jx++)
        if (!fmpz_is_zero(vp + jx) && !COEFF_IS_MPZ(vp[jx]))
            goto found_small;

    /* Alright, bring it on */
    heavy_machinery(rp, vp, vn, ip);
    return;

found_small: /* found small inside vp */
    gd = FLINT_ABS(vp[jx]);

reduce_gd: /* gd is set, but we need to reduce it */
    exp = flint_ctz(gd);

    /* Take care of ip before next step */
    if (ip != NULL)
    {
        slong tsz;
        mp = COEFF_TO_PTR(*ip); /* ip cannot be small */
        tsz = msz = FLINT_ABS(SIZ(mp));
        md = PTR(mp);

        /* Remove trailing zero limbs */
        while (*md == UWORD(0))
        {
            md++;
            msz--;
        }

        if (tsz != msz)
            /* keep exp */;
        else
            exp = FLINT_MIN(exp, flint_ctz(md[0]));

        gd = mpn_gcd_1(md, msz, gd);

        if (gd == UWORD(1))
            goto gd_eq_1;
    }

    goto compute;

ip1: /* We found ip to be small */
    /* gd is set but unreduced and exp is not set */
    FLINT_ASSERT(gd != UWORD(0));
    exp = flint_ctz(gd);
    gd >>= exp;

    if (gd == UWORD(1))
        goto gd_eq_1;

compute: /* gd is set and ip has been taken care of */
    FLINT_ASSERT(gd != UWORD(0));
    FLINT_ASSERT(gd != UWORD(1));
    for (kx = 0; kx < vn; kx++)
    {
        FLINT_ASSERT(gd & 1 == 1);

        /* If kx == jx, we have already taken care of this entry */
        if (fmpz_is_zero(vp + kx) || kx == jx)
            continue;

        if (gd == UWORD(1))
            goto gd_eq_1;

        if (!COEFF_IS_MPZ(vp[kx]))
        {
            ulong td = FLINT_ABS(vp[kx]);
            int c;

            c = flint_ctz(td);
            td >>= c;

            exp = FLINT_MIN(exp, c);

            /* Make td bigger */
            if (td < gd)
                FLINT_SWAP(ulong, td, gd);

            /* Now td >= gd */
            /* If much bigger, reduce via remainder first */
            if ((td >> 16) > gd)
            {
                td %= gd;
                if (td == UWORD(0)) /* td was a multiple of gd */
                    continue;

                c = flint_ctz(td);
                td >>= c;
            }

            gd = mpn_gcd_11(td, gd);
        }
        else
        {
            mp = COEFF_TO_PTR(vp[kx]);
            msz = FLINT_ABS(SIZ(mp));
            md = PTR(mp);

            gd = mpn_gcd_1(md, msz, gd);
        }
    }

    goto fin;

gd_eq_1: /* gd == 1, but we may need to check for common factors of 2 */
    FLINT_ASSERT(gd == UWORD(1));

    if (exp == 0)
        goto fin;

    for (; kx < vn; kx++)
    {
        if (fmpz_is_zero(vp + kx))
            continue;

        if (!COEFF_IS_MPZ(vp[kx]))
        {
            exp = FLINT_MIN(flint_ctz(FLINT_ABS(vp[kx])), exp);

            if (exp == 0)
                goto fin;
        }
        else
        {
            mp = COEFF_TO_PTR(vp[kx]);
            md = PTR(mp);

            if (*md == UWORD(0))
                continue;

            exp = FLINT_MIN(flint_ctz(*md), exp);

            if (exp == 0)
                goto fin;
        }
    }

fin:
    FLINT_ASSERT(gd != UWORD(0));
    fmpz_set_ui(rp, gd << exp);
    return;
}

/* Everything is an mpz */
FLINT_STATIC_NOINLINE
void heavy_machinery(fmpz_t rp, const fmpz * vp, slong vn, const fmpz_t ip)
{
    mpn_ptr mvp;
    ulong limbs = UWORD_MAX, exp;
    slong idx0, idx1;
    ulong sz0 = UWORD_MAX, sz1 = UWORD_MAX;

    ulong a0;

    slong ix;

    if (vn < WORD(2))
        FLINT_UNREACHABLE;

    mvp = flint_malloc(sizeof(mpn_struct) * (vn + (ip != NULL)));

    /* Start by adding entries to mvp */
    if (ip != NULL)
    {
        mpz_srcptr mp;
        mp_srcptr md;
        ulong sz, tz;

        mp = COEFF_TO_PTR(*ip);
        tz = sz = FLINT_ABS(SIZ(mp));
        md = PTR(mp);

        while (*md == UWORD(0))
        {
            md++;
            sz--;
        }

        limbs = tz - sz;
        mvp[0].d = md;
        mvp[0].sz = sz;
        mvp++;

        sz0 = sz;
        idx0 = -1;
    }

    for (ix = 0, vp += vn - 1; ix < vn; vp--)
    {
        mpz_srcptr mp;
        mp_srcptr md;
        ulong sz, tz;

        if (fmpz_is_zero(vp))
        {
            vn--;
            continue;
        }

        mp = COEFF_TO_PTR(*vp);
        tz = sz = FLINT_ABS(SIZ(mp));
        md = PTR(mp);

        while (*md == UWORD(0))
        {
            md++;
            sz--;
        }

        limbs = FLINT_MIN(limbs, tz - sz);
        mvp[ix].d = md;
        mvp[ix].sz = sz;
        ix++;

        if (sz < sz0)
        {
            sz1 = sz0;
            sz0 = sz;
            idx1 = idx0;
            idx0 = ix;
        }
        else if (sz < sz1)
        {
            sz1 = sz;
            idx1 = ix;
        }
    }

    /* Check if we can fast-track the initial GCD */
    if (sz1 == UWORD(1))
    {
        /* 1 = sz0 = sz1 */
        ulong b0;
        ulong bexp;

        a0 = mvp[idx0].d[0];
        b0 = mvp[idx1].d[0];

        exp = flint_ctz(a0);
        bexp = flint_ctz(b0);

        a0 >>= exp;
        b0 >>= bexp;
        exp = FLINT_MIN(exp, bexp);

        if (b0 < a0)
            FLINT_SWAP(ulong, a0, b0);

        /* Now b0 >= a0 */
        /* If much bigger, reduce via remainder first */
        if ((b0 >> 16) > a0)
        {
            b0 %= a0;
            if (b0 == UWORD(0)) /* b0 was a multiple of a0 */
                goto G1;

            b0 >>= flint_ctz(b0);
        }

        a0 = mpn_gcd_11(b0, a0);

        if (a0 == UWORD(1))
            goto a0_eq_1;

        goto G1;
    }
    else if (sz0 == UWORD(1))
    {
        /* 1 = sz0 < sz1 */
        ulong bexp;
        mp_srcptr bd;
        ulong bn;

        a0 = mvp[idx0].d[0];
        bd = mvp[idx1].d;
        bn = mvp[idx1].sz;

        exp = flint_ctz(a0);
        bexp = flint_ctz(bd[0]);

        a0 >>= exp;
        exp = FLINT_MIN(exp, bexp);

        a0 = mpn_gcd_1(bd, bn, a0);

        if (a0 == UWORD(1))
            goto a0_eq_1;

        goto G1;
    }
    else if (sz1 == UWORD(2))
    {
        /* 2 = sz0 = sz1 */
        ulong a1, b0, b1;
        ulong bexp;
        mp_limb_pair_t a;

        a0 = mvp[idx0].d[0];
        a1 = mvp[idx0].d[1];
        b0 = mvp[idx1].d[0];
        b1 = mvp[idx1].d[1];

        exp = flint_ctz(a0);
        bexp = flint_ctz(b0);

        a0 = (a0 >> exp) | (a1 << (FLINT_BITS - exp));
        a1 = (a1 >> exp);
        b0 = (b0 >> bexp) | (b1 << (FLINT_BITS - bexp));
        b1 = (b1 >> bexp);

        a = mpn_gcd_22(b1, b0, a1, a0);
        a1 = a.m2;
        a0 = a.m1;

        if (a1 == UWORD(0))
        {
            if (a0 == UWORD(1))
                goto a0_eq_1;
            else
                goto G1;
        }

        goto G2;
    }
    else if (sz0 == UWORD(2))
    {
        /* 2 = sz0 < sz1 */
        /* Try to reduce `a` to a single limb */
        ulong a1;
        a0 = mvp[idx0].d[0];
        a1 = mvp[idx0].d[1];

        exp = flint_ctz(a0);

        a0 = (a0 >> exp) | (a1 << (FLINT_BITS - exp));
        a1 = (a1 >> exp);

        if (a1 == UWORD(0))
        {
            exp = FLINT_MIN(flint_ctz(mvp[idx1].d[0]), exp);
            a0 = mpn_gcd_1(mvp[idx1].d, mvp[idx1].sz, a0);

            if (a0 == UWORD(1))
                goto a0_eq_1;
            else
                goto G1;
        }
        else
            goto G2;
    }

Gn: /* TODO:
       0. Allocate temporary space for two integers -- one for the output (say
          OP) and one for the input (say IP),
       1. Take the smallest sample we have (given at idx0), and right shift it
          into OP,
       2. Take the next smallest sample and right shift it into IP,
       3. Compute OP <- GCD(IP, OP),
       4. If we are at the end, exit.
       5. If OP can fit in a limb, go to G1.
       6. If OP can fit in two limbs, go to G2.
       7. Take the next sample (no specific order), right shift it into IP and
          go to step 3. */

G2: /* TODO: Do mpn_divrem_2 and continue with mpn_gcd_22 */

G1: for (ix = -1; ix < vn; ix++)
        ;

a0_eq_1: /* a0 == 1, but we may need to check for common factors of 2 */
    FLINT_ASSERT(a0 == UWORD(1));

end:

    flint_free(mvp - (ip != NULL));
}
