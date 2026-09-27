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
    setup. Entry n of a row is the exponent bit length from which
    mpz_powm wins at n limbs.

    The crossover depends on the base: a base of more than a limb uses the
    full-size window, and a single limb uses _powmod_ui, all the more
    profitably when its cube fits in a limb, as that allows a wider
    window. It also depends on the shift: for norm != 0, which is nearly
    every modulus, flint_mpn_mulmod_preinvn shifts its 2n-limb product on
    each call, which costs it about ten percent, while mpz_powm only
    unshifts once on the way in and out. The base 2 is dealt with apart.

    Measured on x86-64 averaging over moduli as well as bases, which
    matters: how often the Barrett reduction needs its correction steps
    depends on the modulus, and moves the ratio by up to 15% from one
    modulus to the next at a given size, whereas mpz_powm hardly varies.
    Where mpz_powm never wins by more than a couple of percent, the entry
    says never.
*/
#define POWMOD_MPZ_MAX_LIMBS 22
#define NEVER 0xffff

#define POWMOD_BASE_FULL 0      /* more than a limb */
#define POWMOD_BASE_LIMB 1      /* one limb, whose cube does not fit */
#define POWMOD_BASE_SMALL 2     /* one limb, whose cube fits, not 2 */

static const unsigned short
powmod_mpz_cutoff[3][2][POWMOD_MPZ_MAX_LIMBS + 1] =
{
    {   /* full: no shift, then shift */
        { 0, 0, 64, 0, 16, 128, 128, 256, 1024, 128, 256, 128, 128, 128,
          NEVER, NEVER, NEVER, NEVER, NEVER, NEVER, NEVER, NEVER, NEVER },
        { 0, 0, 64, 0, 0, 16, 32, 64, 64, 64, 64, 64, 64, 64,
          256, 256, 256, 256, 256, 256, 256, 256, 256 },
    },
    {   /* limb */
        { 0, 0, 32, 0, 32, 128, 128, 256, NEVER, 1024, 2048, NEVER, NEVER,
          NEVER, NEVER, NEVER, NEVER, NEVER, NEVER, NEVER, NEVER, NEVER,
          NEVER },
        { 0, 0, 32, 0, 0, 64, 64, 128, 128, 128, 1024, 256, 1024, 512,
          NEVER, NEVER, NEVER, NEVER, NEVER, NEVER, NEVER, NEVER, NEVER },
    },
    {   /* small */
        { 0, 0, 128, 32, 128, 1024, NEVER, NEVER, NEVER, NEVER, NEVER,
          NEVER, NEVER, NEVER, NEVER, NEVER, NEVER, NEVER, NEVER, NEVER,
          NEVER, NEVER, NEVER },
        { 0, 0, 64, 16, 32, 512, 1024, 1024, NEVER, NEVER, NEVER, NEVER,
          NEVER, NEVER, NEVER, NEVER, NEVER, NEVER, NEVER, NEVER, NEVER,
          NEVER, NEVER },
    },
};

/*
    mpz_powm has a shortcut for the base 2 which leaves it nothing to set
    up and no full multiplications. It is then ahead at any size for an
    exponent below this, and for any exponent up to the given sizes.
*/
#define POWMOD_MPZ_BASE2_EBITS 128
#define POWMOD_MPZ_BASE2_NORM0_MAX_LIMBS 7
#define POWMOD_MPZ_BASE2_NORMED_MAX_LIMBS 22

static int
_powmod_want_mpz(int kind, mp_size_t n, flint_bitcnt_t ebits, ulong norm)
{
    return n <= POWMOD_MPZ_MAX_LIMBS
        && ebits >= powmod_mpz_cutoff[kind][norm != 0][n];
}

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

/*
    Window width by exponent length, measured over n = 8 to 32 limbs; the
    optimum hardly depends on n, and the top of the range matches GMP's
    own thresholds. The table has 2^(w-1) entries of n limbs, which is
    capped only to rule out an absurd allocation for a huge modulus: one
    step below the optimal width costs no more than a couple of percent.
*/
#define POWMOD_WINDOW_MAX_TAB_LIMBS (WORD(1) << 23)

static int
_window_size(flint_bitcnt_t ebits, mp_size_t n)
{
    static const unsigned int tab[] =
        { 8, 24, 70, 200, 700, 1800, 4600, 11500, 28000 };
    int w;

    for (w = 1; w < 10 && ebits >= tab[w - 1]; w++)
        ;

    while (w > 1 && (((slong) 1) << (w - 1)) * n > POWMOD_WINDOW_MAX_TAB_LIMBS)
        w--;

    return w;
}

/*
    Left to right sliding window over the ebits bits of e. The first
    window is loaded with set_tab and every later one is folded in with
    mul_tab, after the squarings that make room for it; val is the odd
    value of the window, whose table entry is (val - 1) / 2.
*/
#define POWMOD_SLIDING_WINDOW(set_tab, mul_tab, sqr_res)                    \
    do {                                                                    \
        slong i, l, k;                                                      \
        ulong val;                                                          \
        int started = 0;                                                    \
                                                                            \
        for (i = ebits - 1; i >= 0; )                                       \
        {                                                                   \
            if (!EBIT(e, i))                                                \
            {                                                               \
                if (started)                                                \
                    sqr_res;                                                \
                                                                            \
                i--;                                                        \
                continue;                                                   \
            }                                                               \
                                                                            \
            /* the longest window ending in a set bit, hence odd */         \
            l = FLINT_MAX(i - w + 1, 0);                                    \
                                                                            \
            while (!EBIT(e, l))                                             \
                l++;                                                        \
                                                                            \
            val = 0;                                                        \
            for (k = i; k >= l; k--)                                        \
                val = 2 * val + EBIT(e, k);                                 \
                                                                            \
            if (started)                                                    \
            {                                                               \
                for (k = 0; k <= i - l; k++)                                \
                    sqr_res;                                                \
                                                                            \
                mul_tab;                                                    \
            }                                                               \
            else                                                            \
            {                                                               \
                /* the leading window needs no squarings before it */       \
                set_tab;                                                    \
                started = 1;                                                \
            }                                                               \
                                                                            \
            i = l - 1;                                                      \
        }                                                                   \
    } while (0)

/* The general case: a table of the odd powers a, a^3, ..., a^(2^w - 1). */
static void
_powmod_window(mp_ptr res, mp_srcptr a, mp_srcptr e, flint_bitcnt_t ebits,
        mp_size_t n, mp_srcptr d, mp_srcptr dinv, ulong norm, int w)
{
    slong k, tabn = ((slong) 1) << (w - 1);
    mp_ptr tab, sqr;
    TMP_INIT;

    TMP_START;

    tab = TMP_ALLOC((tabn + 1) * n * sizeof(mp_limb_t));
    sqr = tab + tabn * n;

    flint_mpn_copyi(tab, a, n);

    if (tabn > 1)
    {
        flint_mpn_mulmod_preinvn(sqr, a, a, n, d, dinv, norm);

        for (k = 1; k < tabn; k++)
            flint_mpn_mulmod_preinvn(tab + k * n, tab + (k - 1) * n, sqr,
                    n, d, dinv, norm);
    }

    POWMOD_SLIDING_WINDOW(
        flint_mpn_copyi(res, tab + ((val - 1) / 2) * n, n),
        flint_mpn_mulmod_preinvn(res, res, tab + ((val - 1) / 2) * n,
                n, d, dinv, norm),
        flint_mpn_mulmod_preinvn(res, res, res, n, d, dinv, norm));

    TMP_END;
}

/*
    A base that is a single limb before shifting, as in Miller-Rabin. Its
    odd powers b, b^3, ..., b^(2^w - 1) are kept as plain limbs, for as
    wide a window as they fit in, so the table costs nothing and folding
    in a window is a multiplication by a limb followed by a one-limb
    quotient, O(n) instead of a modular multiplication. Only the
    squarings are left at full cost.

    _powmod_ui_tab fills the table and returns w, which is 1 when even b^3
    does not fit; _powmod_ui_cube_fits says so without building it.
*/
#define POWMOD_UI_MAX_WINDOW 6      /* wide enough for 2^63 */

static int
_powmod_ui_tab(ulong * tab, ulong b)
{
    ulong hi, lo, b2;
    slong k, tabn;
    int w = 1;

    tab[0] = b;
    umul_ppmm(hi, b2, b, b);

    while (hi == 0 && w < POWMOD_UI_MAX_WINDOW)
    {
        tabn = ((slong) 1) << (w - 1);

        for (k = tabn; k < 2 * tabn; k++)
        {
            umul_ppmm(hi, lo, tab[k - 1], b2);

            if (hi != 0)
                break;

            tab[k] = lo;
        }

        if (hi == 0)
            w++;
    }

    return w;
}

static int
_powmod_ui_cube_fits(ulong b)
{
    ulong hi, lo;

    umul_ppmm(hi, lo, b, b);

    if (hi != 0)
        return 0;

    umul_ppmm(hi, lo, lo, b);

    return hi == 0;
}

/* Needs n >= 2, for the one-limb quotient. */
static void
_powmod_ui(mp_ptr res, ulong b, mp_srcptr e, flint_bitcnt_t ebits,
        mp_size_t n, mp_srcptr d, mp_srcptr dinv, ulong norm)
{
    ulong tab[((slong) 1) << (POWMOD_UI_MAX_WINDOW - 1)];
    ulong dinv1, q;
    mp_ptr t;
    int w;
    TMP_INIT;

    w = _powmod_ui_tab(tab, b);
    dinv1 = flint_mpn_preinv1(d[n - 1], d[n - 2]);

    TMP_START;

    /* one limb of headroom for the product by a table entry */
    t = TMP_ALLOC((n + 1) * sizeof(mp_limb_t));

    /* a table entry is below 2^FLINT_BITS <= d, so shifting it is exact */
#define POWMOD_UI_SET(c)                                                    \
    do {                                                                    \
        flint_mpn_zero(t, n);                                               \
        t[0] = (c) << norm;                                                 \
        t[1] = (norm == 0) ? 0 : (c) >> (FLINT_BITS - norm);                \
    } while (0)

#define POWMOD_UI_MUL(c)                                                    \
    do {                                                                    \
        t[n] = mpn_mul_1(t, t, n, (c));                                     \
        flint_mpn_divrem_preinv1(&q, t, n + 1, d, n, dinv1);                \
    } while (0)

    POWMOD_SLIDING_WINDOW(
        POWMOD_UI_SET(tab[(val - 1) / 2]),
        POWMOD_UI_MUL(tab[(val - 1) / 2]),
        flint_mpn_mulmod_preinvn(t, t, t, n, d, dinv, norm));

#undef POWMOD_UI_SET
#undef POWMOD_UI_MUL

    flint_mpn_copyi(res, t, n);

    TMP_END;
}

/*
    The base as a single limb before shifting, if it is one. Shifted, it
    is below 2^(FLINT_BITS + norm), which only the two low limbs can hold.
*/
static int
_powmod_small_base(ulong * b, mp_srcptr a, mp_size_t n, ulong norm)
{
    if (n < 2 || !flint_mpn_zero_p(a + 2, n - 2))
        return 0;

    if (norm == 0)
    {
        if (a[1] != 0)
            return 0;

        *b = a[0];
    }
    else
    {
        if ((a[1] >> norm) != 0)
            return 0;

        *b = (a[0] >> norm) | (a[1] << (FLINT_BITS - norm));
    }

    return 1;
}

void
flint_mpn_powmod_preinvn(mp_ptr res, mp_srcptr a, mp_srcptr e, mp_size_t en,
        mp_size_t n, mp_srcptr d, mp_srcptr dinv, ulong norm)
{
    flint_bitcnt_t ebits;
    ulong b;
    int kind;

    while (en > 0 && e[en - 1] == 0)
        en--;

    if (en == 0)
    {
        flint_mpn_zero(res, n);
        res[0] = UWORD(1) << norm;
        return;
    }

    ebits = en * FLINT_BITS - flint_clz(e[en - 1]);

    if (!_powmod_small_base(&b, a, n, norm))
    {
        kind = POWMOD_BASE_FULL;
    }
    else if (b <= 1)
    {
        /* 0^e = 0 and 1^e = 1, e being positive; shifted, a is the answer */
        flint_mpn_copyi(res, a, n);
        return;
    }
    else if (b == 2)
    {
        if (ebits < POWMOD_MPZ_BASE2_EBITS || n <= (norm == 0 ?
                POWMOD_MPZ_BASE2_NORM0_MAX_LIMBS :
                POWMOD_MPZ_BASE2_NORMED_MAX_LIMBS))
            _powmod_mpz(res, a, e, en, n, d, norm);
        else
            _powmod_ui(res, b, e, ebits, n, d, dinv, norm);
        return;
    }
    else
    {
        kind = _powmod_ui_cube_fits(b) ? POWMOD_BASE_SMALL : POWMOD_BASE_LIMB;
    }

    if (_powmod_want_mpz(kind, n, ebits, norm))
        _powmod_mpz(res, a, e, en, n, d, norm);
    else if (kind == POWMOD_BASE_FULL || (kind == POWMOD_BASE_LIMB && n < 4))
        /* below 4 limbs a one-bit window loses to the full-size one */
        _powmod_window(res, a, e, ebits, n, d, dinv, norm,
                _window_size(ebits, n));
    else
        _powmod_ui(res, b, e, ebits, n, d, dinv, norm);
}
