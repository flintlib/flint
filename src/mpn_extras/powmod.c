/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

#define EBIT(e, i) (((e)[(i) / FLINT_BITS] >> ((i) % FLINT_BITS)) & 1)

/*
    GMP's Montgomery arithmetic beats Barrett with a precomputed inverse
    for a small modulus, once the exponent is long enough to pay for its
    setup. Entry n of the table is the exponent bit length from which
    mpz_powm wins at n limbs; past the end of the table it never does.
    Measured on x86-64 over varied bases, which matters: timing one base
    repeatedly understates the window code's cost. The ratios are flat at
    each boundary, so the exact cutoffs are not critical. Past 7 limbs the
    two are within a few percent of each other either way, with no
    consistent winner, so there is nothing to gain from dispatching.
*/
#define POWMOD_MPZ_MAX_LIMBS 7

static const unsigned int powmod_mpz_cutoff_tab[POWMOD_MPZ_MAX_LIMBS + 1] =
{
    0,      /* unused */
    0,      /* 1 limb: always, by a factor 1.5 to 4 */
    64,     /* 2 limbs */
    0,      /* 3 limbs: always, we have no unrolled mulmod for this size */
    16,     /* 4 limbs */
    192,    /* 5 limbs */
    192,    /* 6 limbs */
    192     /* 7 limbs */
};

/*
    a^e mod d through mpz_powm. The shift is undone and redone around it;
    both directions are exact, a and d being shifted values by contract.
*/
static void
_powmod_mpz(mp_ptr res, mp_srcptr a, mp_srcptr e, mp_size_t en, mp_size_t n,
        mp_srcptr d, ulong norm)
{
    mpz_t az, dz, ez, rz;
    mp_srcptr ap, dp;
    mp_ptr t;
    mp_size_t an, dn, rn;
    TMP_INIT;

    TMP_START;

    if (norm)
    {
        t = TMP_ALLOC((2 * n) * sizeof(mp_limb_t));
        mpn_rshift(t, a, n, norm);
        mpn_rshift(t + n, d, n, norm);
        ap = t;
        dp = t + n;
    }
    else
    {
        ap = a;
        dp = d;
    }

    /* GMP wants the top limb of a read-only view to be nonzero */
    for (an = n; an > 0 && ap[an - 1] == 0; an--)
        ;
    for (dn = n; dn > 0 && dp[dn - 1] == 0; dn--)
        ;

    mpz_roinit_n(az, ap, an);
    mpz_roinit_n(dz, dp, dn);
    mpz_roinit_n(ez, e, en);

    mpz_init(rz);
    mpz_powm(rz, az, ez, dz);

    rn = rz->_mp_size;
    flint_mpn_copyi(res, rz->_mp_d, rn);
    flint_mpn_zero(res + rn, n - rn);
    mpz_clear(rz);

    if (norm)
        mpn_lshift(res, res, n, norm);

    TMP_END;
}

/* Window width, cut back so that the table of odd powers stays modest. */
static int
_window_size(flint_bitcnt_t ebits, mp_size_t n)
{
    int w;

    if (ebits < 8)
        w = 1;
    else if (ebits < 24)
        w = 2;
    else if (ebits < 70)
        w = 3;
    else if (ebits < 200)
        w = 4;
    else if (ebits < 500)
        w = 5;
    else if (ebits < 1200)
        w = 6;
    else if (ebits < 3000)
        w = 7;
    else
        w = 8;

    while (w > 1 && (((mp_size_t) 1) << (w - 1)) * n > 65536)
        w--;

    return w;
}

void
flint_mpn_powmod_preinvn(mp_ptr res, mp_srcptr a, mp_srcptr e, mp_size_t en,
        mp_size_t n, mp_srcptr d, mp_srcptr dinv, ulong norm)
{
    flint_bitcnt_t ebits;
    slong i, l, k, tabn;
    int w, started = 0;
    ulong val;
    mp_ptr tab, sqr;
    TMP_INIT;

    while (en > 0 && e[en - 1] == 0)
        en--;

    if (en == 0)
    {
        flint_mpn_zero(res, n);
        res[0] = UWORD(1) << norm;
        return;
    }

    ebits = en * FLINT_BITS - flint_clz(e[en - 1]);

    if (n <= POWMOD_MPZ_MAX_LIMBS && ebits >= powmod_mpz_cutoff_tab[n])
    {
        _powmod_mpz(res, a, e, en, n, d, norm);
        return;
    }

    w = _window_size(ebits, n);
    tabn = ((mp_size_t) 1) << (w - 1);

    TMP_START;

    tab = TMP_ALLOC((tabn + 1) * n * sizeof(mp_limb_t));
    sqr = tab + tabn * n;

    /* the odd powers a, a^3, ..., a^(2^w - 1) */
    flint_mpn_copyi(tab, a, n);

    if (tabn > 1)
    {
        flint_mpn_mulmod_preinvn(sqr, a, a, n, d, dinv, norm);

        for (k = 1; k < tabn; k++)
            flint_mpn_mulmod_preinvn(tab + k * n, tab + (k - 1) * n, sqr,
                    n, d, dinv, norm);
    }

    for (i = ebits - 1; i >= 0; )
    {
        if (!EBIT(e, i))
        {
            if (started)
                flint_mpn_mulmod_preinvn(res, res, res, n, d, dinv, norm);

            i--;
            continue;
        }

        /* the longest window ending in a set bit, hence of odd value */
        l = FLINT_MAX(i - w + 1, 0);

        while (!EBIT(e, l))
            l++;

        val = 0;
        for (k = i; k >= l; k--)
            val = 2 * val + EBIT(e, k);

        if (started)
        {
            for (k = 0; k <= i - l; k++)
                flint_mpn_mulmod_preinvn(res, res, res, n, d, dinv, norm);

            flint_mpn_mulmod_preinvn(res, res, tab + ((val - 1) / 2) * n,
                    n, d, dinv, norm);
        }
        else
        {
            /* the leading window, which needs no squarings before it */
            flint_mpn_copyi(res, tab + ((val - 1) / 2) * n, n);
            started = 1;
        }

        i = l - 1;
    }

    TMP_END;
}
