/*
    Copyright 1994, 1996, 2000, 2001, 2009, 2012, 2019 Free Software Foundation, Inc.
    Copyright (C) 2020 Daniel Schultz
    Copyright (C) 2021 Fredrik Johansson
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    Contains code from mpn_gcd_1 in GMP 6.3.0 under label `gcd1_limbs_eq_zero`.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

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

/* TODO
   * Do we want to utilize mpn_gcd_22 coupled with mpn_divrem_2? Both functions
     should be available in all new versions of GMP.
   * Do we want to right shift result before pushing it into mpn_gcd_1?
*/

#define MAX_STACK_ALLOC 1024

#define SIZ(x) ((x)->_mp_size)
#define PTR(x) ((x)->_mp_d)

typedef struct
{
    ulong gd;
    int exp;
    ulong limbs;
    slong jx;
    slong kx;
}
hm_t;

typedef struct
{
    mp_srcptr d;
    ulong sz;
}
mpn_struct;

typedef mpn_struct * mpn_ptr;

FLINT_STATIC_NOINLINE void heavy_machinery(fmpz_t, const fmpz *, slong, const fmpz_t);

void
_fmpz_vec_content_chained(fmpz_t rp, const fmpz * vp, slong vn, const fmpz_t ip)
{
    slong jx = -1, kx = 0;
    ulong gd; /* single limb gcd */
    int exp = 0;
    ulong limbs = 0; /* factors of power of two */
    mpz_srcptr mp;
    slong msz;
    mp_srcptr md;

#if FLINT_WANT_ASSERT
    gd = 0;
#endif

    /* Check if `ip` can be represented as a single limb */
    if (ip != NULL)
    {
        fmpz tip = *ip; /* Read-only */

        if (!COEFF_IS_MPZ(tip))
        {
            if (tip == WORD(0))
            {
                ip = NULL;
                goto check_vp;
            }

            gd = FLINT_ABS(tip);
            goto ip1;
        }
        else
        {
            mp = COEFF_TO_PTR(tip);
            msz = FLINT_ABS(SIZ(mp));
            md = PTR(mp);

            if (msz == 1)
            {
                gd = md[0];
                goto ip1;
            }
            else
            {
                limbs = 0;

                while (md[limbs] == UWORD(0))
                    limbs++;

                if (limbs == msz - 1)
                {
                    /* ip = limb * 2^{FLINT_BITS * limbs} */
                    gd = md[limbs];
                    exp = flint_ctz(md[limbs]);

                    if (limbs == UWORD(0))
                        goto gcd1_limbs_eq_zero;
                    else
                        goto gcd1_limbs_nonzero;
                }
                else if (limbs == msz - 2)
                {
                    exp = flint_ctz(md[limbs]);

                    if ((md[msz] >> exp) == UWORD(0))
                    {
                        /* ip = limb * 2^{FLINT_BITS * limbs + exp} */
                        gd = (md[limbs] >> exp) | (md[msz] << (FLINT_BITS - exp));

                        if (limbs == UWORD(0))
                            goto gcd1_limbs_eq_zero;
                        else
                            goto gcd1_limbs_nonzero;
                    }
                }
            }
        }
    }

check_vp: /* Check if we can find any single-limbs in vp. Fast checks only. */
    for (jx = 0; jx < vn; jx++)
        if (!COEFF_IS_MPZ(vp[jx]) && !fmpz_is_zero(vp + jx))
            goto found_small;

    /* Alright, bring it on */
    heavy_machinery(rp, vp, vn, ip);
    return;

found_small: /* found small inside vp */
    limbs = 0;
    gd = FLINT_ABS(vp[jx]);

reduce_gd: /* gd is set, but we need to reduce it */
    exp = flint_ctz(gd);

    /* Take care of ip before next step */
    if (ip != NULL)
    {
        slong tsz;
        mp = COEFF_TO_PTR(*ip);
        tsz = msz = FLINT_ABS(SIZ(mp));
        md = PTR(mp);

        /* Remove trailing zero limbs */
        while (*md == UWORD(0))
        {
            md++;
            msz--;
        }

        /* Set limbs and exp */
        if (tsz - msz > limbs)
            /* Do nothing */;
        else if (tsz - msz == limbs)
        {
            /* Check for smallest exp while not touching limbs. */
            exp = FLINT_MIN(exp, flint_ctz(*md));
        }
        else /* tsz - msz < limbs */
        {
            /* Set limbs and exp according to ip */
            limbs = tsz - msz;
            exp = flint_ctz(*md);
        }

        gd = mpn_gcd_1(md, msz, gd);
    }

gcd1_limbs_nonzero: /* gd is set and ip has been taken care of, limbs != 0 */
    FLINT_ASSERT(limbs != UWORD(0));
    FLINT_ASSERT(gd != UWORD(0));
    FLINT_ASSERT(gd != UWORD(1));
    for (; kx < vn; kx++)
    {
        FLINT_ASSERT(gd & 1 == 1);

        if (gd == UWORD(1))
            goto gd_eq_1;

        /* If kx == jx, we have already taken care of this entry */
        if (fmpz_is_zero(vp + kx) || kx == jx)
            continue;

        if (!COEFF_IS_MPZ(vp[kx]))
        {
            limbs = 0;
            goto gcd1_limbs_eq_zero;
        }
        else
        {
            slong tsz;
            mp = COEFF_TO_PTR(vp[kx]);
            tsz = msz = FLINT_ABS(SIZ(mp));
            md = PTR(mp);

            if (*md != UWORD(0))
            {
                limbs = 0;
                goto gcd1_limbs_eq_zero;
            }

            do
            {
                md++;
                msz--;
            } while (*md == UWORD(0));

            if (tsz - msz > limbs)
                /* Do nothing */;
            else if (tsz - msz == limbs)
            {
                /* Check for smallest exp while not touching limbs. */
                exp = FLINT_MIN(exp, flint_ctz(*md));
            }
            else /* tsz - msz < limbs */
            {
                /* Set both limbs and exp */
                limbs = tsz - msz;
                exp = flint_ctz(*md);
            }

            gd = mpn_gcd_1(md, msz, gd);
        }
    }

    goto small_end;

ip1: /* We found ip to be small */
    /* limbs = 0, exp unset and gd unreduced */
    FLINT_ASSERT(limbs == UWORD(0));
    FLINT_ASSERT(gd != UWORD(0));
    exp = flint_ctz(gd);
    gd >>= exp;

    if (gd == UWORD(1))
        goto gd_eq_1;

gcd1_limbs_eq_zero: /* gd is set, ip has been taken care of and limbs = 0 */
    FLINT_ASSERT(limbs == UWORD(0));
    FLINT_ASSERT(gd != UWORD(0));
    FLINT_ASSERT(gd != UWORD(1));
    for (; kx < vn; kx++)
    {
        FLINT_ASSERT(gd & 1 == 1);

        if (gd == UWORD(1))
            goto gd_eq_1;

        /* If kx == jx, we have already taken care of this entry */
        if (fmpz_is_zero(vp + kx) || kx == jx)
            continue;

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

    goto small_end;

gd_eq_1: /* gd == 1, but we may need to check for common factors of 2 */
    FLINT_ASSERT(gd == UWORD(1));

    if (limbs == UWORD(0) && exp == 0)
    {
        goto small_end;
    }
    else if (limbs == UWORD(0))
    {
limbs_is_zero:
        /* Check for smallest trailing zero count in first limbs */
        for (; kx < vn; kx++)
        {
            if (fmpz_is_zero(vp + kx))
                continue;
            else if (!COEFF_IS_MPZ(vp[kx]))
            {
                exp = FLINT_MIN(flint_ctz(FLINT_ABS(vp[kx])), exp);

                if (exp == 0)
                    goto small_end;
            }
            else
            {
                mp = COEFF_TO_PTR(vp[kx]);
                md = PTR(mp);

                if (*md == UWORD(0))
                    continue;

                exp = FLINT_MIN(flint_ctz(*md), exp);

                if (exp == 0)
                    goto small_end;
            }
        }
    }
    else
    {
        /* Check for trailing zero limbs and trailing zero count */
        /* Go to previous case if we ever hit a zero limb count */
        for (; kx < vn; kx++)
        {
            if (!COEFF_IS_MPZ(vp[kx]))
            {
                limbs = 0;
                goto limbs_is_zero;
            }
            else
            {
                slong tsz;
                mp = COEFF_TO_PTR(vp[kx]);
                md = PTR(mp);

                if (*md != UWORD(0))
                {
                    limbs = 0;
                    goto limbs_is_zero;
                }

                tsz = msz = FLINT_ABS(SIZ(mp));

                do
                {
                    md++;
                    msz--;
                } while (*md == UWORD(0));

                /* Set limbs and exp */
                if (tsz - msz > limbs)
                    continue; /* Do nothing */
                else if (tsz - msz == limbs)
                {
                    /* Check for smallest exp while not touching limbs. */
                    exp = FLINT_MIN(exp, flint_ctz(*md));
                }
                else /* tsz - msz < limbs */
                {
                    /* Set both limbs and exp */
                    limbs = tsz - msz;
                    exp = flint_ctz(*md);
                }
            }
        }
    }

small_end: /* Finished! Set rp = gd * 2^(FLINT_BITS * limbs + exp) */
    FLINT_ASSERT(gd != UWORD(0));

    if (limbs == UWORD(0) && (exp == 0 || (gd >> (FLINT_BITS - exp)) == UWORD(0)))
    {
        fmpz_set_ui(rp, gd << exp);
    }
    else
    {
        mpz_ptr mrp = _fmpz_promote(rp);
        mp_ptr mrpd;

        if (FLINT_ABS(SIZ(mrp)) < limbs + 2)
            mrpd = _mpz_realloc(mrp, limbs + 2);
        else
            mrpd = PTR(mrp);

        for (jx = 0; jx < limbs; jx++)
            mrpd[jx] = 0;

        mrpd[limbs] = gd << exp;
        mrpd[limbs + 1] = gd >> (FLINT_BITS - exp);
        SIZ(mrp) = limbs + 1 + (mrpd[limbs + 1] != UWORD(0));
    }

    return;
}

/* Everything is big, but vp may contain single-limbed integers. */
FLINT_STATIC_NOINLINE
void heavy_machinery(fmpz_t rp, const fmpz * vp, slong vn, const fmpz_t ip)
{
    mpn_ptr mvp;
    mpz_srcptr mp;
    mp_srcptr md;
    ulong tz, sz;
    slong ix;
    ulong limbs = UWORD_MAX, exp;
    slong idx1, idx2;
    ulong sz1 = UWORD_MAX, sz2 = UWORD_MAX;
    ulong g0, g1;
    mp_ptr gd;
    ulong gn;

    mvp = flint_malloc(sizeof(mpn_struct) * (vn + (ip != NULL)));

    if (ip != NULL)
    {
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
    }

    for (ix = 0, vp += vn - 1; ix < vn; ix++, vp--)
    {
        if (fmpz_is_zero(vp))
            continue;

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

        if (sz < sz1)
        {
            sz2 = sz1;
            sz1 = sz;
            idx2 = idx1;
            idx1 = ix;
        }
        else if (sz < sz2)
        {
            sz2 = sz;
            idx2 = ix;
        }
    }

    mvp -= (ip != NULL);

    if (sz2 == UWORD(1))
    {
        /* 1 = sz1 = sz2 */
        g0 = mvp[idx1 + 1].d[0];
        goto G1;
    }
    else if (sz1 == UWORD(1))
    {
        /* 1 = sz1 < sz2 */
        g0 = mvp[idx1 + 1].d[0];
        goto G2;
    }
    else if (sz2 == UWORD(2))
    {
        /* 2 = sz1 = sz2 */
    }
    else if (sz1 == UWORD(2))
    {
        /* 2 = sz1 < sz2 */
    }
    else
    {
        /* 2 < sz1 <= sz2 */
        
    }

Gn:

G2:

G1:

    flint_free(mvp);
}
