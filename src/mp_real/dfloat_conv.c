/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include <string.h>
#include <stdint.h>
#include <float.h>
#include "flint.h"
#include "mpn_extras.h"
#include "mp_real.h"
#include "impl.h"

/* Conversions between mp_real balls and the double expansions of the
   dfloat module, given as arrays of doubles (the midpoint components)
   and a double radius, without going through arb.

   mp_real_set_dfloat sums the components exactly in a two's complement
   frame spanning their exponents (at most 2^1024 down to 2^-1074, so
   36 limbs on 64-bit), then attaches the radius.

   mp_real_get_dfloat rounds greedily: each component is the double
   nearest to what remains (54 bits read, rounded to 53), subtracted
   exactly; the remainder, the truncated tail of the mantissa and the
   input radius go to the radius, rounded up.  Components that would be
   subnormal are not produced (their value stays in the remainder), and
   a value or radius beyond the double range gives dfloat's whole line,
   zero components and an infinite radius. */

#define DFC_FRAME_LIMBS ((1024 + 1074 + 53 + 2 * FLINT_BITS) / FLINT_BITS + 2)
/* the frame of mp_real_get_dfloat: K + 1 limbs for up to 8 components
   (17 on 32-bit machines) */
#define DFC_GET_LIMBS ((53 * 8 + 2 * FLINT_BITS) / FLINT_BITS + 2)

/* (m, 1 or 2 limbs) 2^r: m < 2^53, 0 <= r < FLINT_BITS, into t; returns
   the number of limbs */
static slong
_dfc_shifted(nn_ptr t, uint64_t m, int r)
{
    slong len;
#if FLINT_BITS == 64
    t[0] = m;
    len = 1;
#else
    t[0] = (ulong) (m & 0xffffffffu);
    t[1] = (ulong) (m >> 32);
    len = 2;
#endif
    if (r != 0)
    {
        t[len] = mpn_lshift(t, t, len, r);
        len++;
    }
    return len;
}

/* a finite double x != 0 as (-1)^neg m 2^e, m < 2^53 an integer (from
   the bit pattern: no frexp) */
FLINT_FORCE_INLINE uint64_t
_dfc_decode(double x, slong * e, int * neg)
{
    uint64_t u, m;
    int ex;
    memcpy(&u, &x, sizeof(u));
    ex = (int) ((u >> 52) & 0x7ff);
    m = u & (((uint64_t) 1 << 52) - 1);
    if (ex != 0)
        m |= (uint64_t) 1 << 52;
    else
        ex = 1;
    *e = ex - 1075;
    *neg = (int) (u >> 63);
    return m;
}

/* res += [+- rad], rad > 0 finite: as an ulp count of the bottom limb
   directly when that is exact or within 2^-20 (rounded up) and the count
   fits a limb (on 32-bit machines a 53-bit mantissa usually does not),
   else through the bound arithmetic (which pads or truncates the
   mantissa) */
static void
_dfc_add_rad(mp_real_t res, double rad)
{
    slong re, sh;
    int neg;
    uint64_t rm = _dfc_decode(rad, &re, &neg), v;

    if (res->size != 0 && res->err == 0)
    {
        sh = re - FLINT_BITS * (res->exp - res->size);
        if (sh >= 0 && sh < 10)
        {
            v = rm << sh;
            if (v <= (uint64_t) UWORD_MAX)
            {
                res->err = (ulong) v;
                return;
            }
        }
        else if (sh < 0 && sh > -33 && (rm >> (-sh)) >= ((uint64_t) 1 << 20))
        {
            v = (rm >> (-sh)) + ((rm & (((uint64_t) 1 << (-sh)) - 1)) != 0);
            if (v <= (uint64_t) UWORD_MAX)
            {
                res->err = (ulong) v;
                return;
            }
        }
    }
    {
        int e;
        double f = frexp(rad, &e);
        slong anc = (e >= 0) ? e / FLINT_BITS : -((-(slong) e + FLINT_BITS - 1) / FLINT_BITS);
        _mp_real_add_error_ulps_at(res, ldexp(f, (int) (e - FLINT_BITS * anc)), anc);
    }
}

void
mp_real_set_dfloat(mp_real_t res, const double * x, slong n, double rad)
{
    ulong A[DFC_FRAME_LIMBS], t[4];
    uint64_t mi[8];
    slong ei[8];
    int ni[8];
    slong i, lo = WORD_MAX, hi = WORD_MIN, L, b, q, tl, nz = 0;
    int neg = 0;

    if (!(rad >= 0.0 && rad <= DBL_MAX))
        flint_throw(FLINT_ERROR, "mp_real_set_dfloat: the radius must be finite and nonnegative\n");
    FLINT_ASSERT(n <= 8);

    for (i = 0; i < n; i++)
    {
        if (!(fabs(x[i]) <= DBL_MAX))
            flint_throw(FLINT_ERROR, "mp_real_set_dfloat: nonfinite component\n");
        if (x[i] != 0.0)
        {
            mi[nz] = _dfc_decode(x[i], ei + nz, ni + nz);
            lo = FLINT_MIN(lo, ei[nz]);
            hi = FLINT_MAX(hi, ei[nz] + 53);
            nz++;
        }
    }

    if (nz == 0)
    {
        mp_real_zero(res);
    }
    else if (nz == 1)
    {
        ulong m = (ulong) mi[0];
#if FLINT_BITS == 64
        _mp_real_set_mpn_2exp(res, &m, 1, ei[0]);
#else
        ulong mm[2] = { (ulong) (mi[0] & 0xffffffffu), (ulong) (mi[0] >> 32) };
        (void) m;
        _mp_real_set_mpn_2exp(res, mm, 2, ei[0]);
#endif
        if (ni[0])
            mp_real_neg(res, res);
    }
    else
    {
        /* the frame: bit 0 at 2^lo, room for the sum of the components
           of magnitude < 2^hi and a sign bit */
        L = (hi - lo + FLINT_BIT_COUNT((ulong) nz) + 1) / FLINT_BITS + 1;
        FLINT_ASSERT(L <= DFC_FRAME_LIMBS);
        flint_mpn_zero(A, L);

        for (i = 0; i < nz; i++)
        {
            b = ei[i] - lo;
            q = b / FLINT_BITS;
            tl = _dfc_shifted(t, mi[i], (int) (b % FLINT_BITS));
            tl = FLINT_MIN(tl, L - q);
            if (!ni[i])
                mpn_add(A + q, A + q, L - q, t, tl);
            else
                mpn_sub(A + q, A + q, L - q, t, tl);
        }

        if (A[L - 1] >> (FLINT_BITS - 1))
        {
            mpn_neg(A, A, L);
            neg = 1;
        }

        _mp_real_set_mpn_2exp(res, A, L, lo);
        if (neg)
            mp_real_neg(res, res);
    }

    if (rad != 0.0)
        _dfc_add_rad(res, rad);
}

/* v 2^e rounded up to a double, v >= 0 an upper bound already: 0 for
   v = 0, at least the smallest subnormal for v > 0, inf beyond the
   range */
static double
_dfc_up(double v, slong e)
{
    if (v == 0.0)
        return 0.0;
    if (e > 1100)
        return HUGE_VAL;
    if (e < -1200)
        return 0x1p-1074;
    v = ldexp(v, (int) e) * (1 + 0x1p-52);
    /* below the normal range ldexp rounds to the subnormal grid, perhaps
       downward: one more step of it */
    return (v < 0x1p-1022) ? v + 0x1p-1074 : v;
}

/* |V| for the two's complement (V, L) into M, returning the bit length
   (0 for V = 0) and the sign */
static slong
_dfc_abs(nn_ptr M, nn_srcptr V, slong L, int * negative)
{
    slong l;

    *negative = (int) (V[L - 1] >> (FLINT_BITS - 1));
    if (*negative)
        mpn_neg(M, V, L);
    else
        flint_mpn_copyi(M, V, L);
    for (l = L; l > 0 && M[l - 1] == 0; l--)
        ;
    if (l == 0)
        return 0;
    return FLINT_BITS * (l - 1) + FLINT_BIT_COUNT(M[l - 1]);
}

/* the bits [bl - k, bl) of (M, L) as an integer, k <= 53 < FLINT_BITS on
   64-bit (two windows of limbs on 32-bit), bits below 0 read as zero */
static uint64_t
_dfc_top_bits(nn_srcptr M, slong bl, int k)
{
    slong pos = bl - k;
    uint64_t r;

#if FLINT_BITS == 64
    if (pos >= 0)
    {
        slong q = pos / 64;
        int b = (int) (pos % 64);
        r = M[q] >> b;
        if (b != 0 && b + k > 64)
            r |= M[q + 1] << (64 - b);
    }
    else
        r = M[0] << (-pos);
    return r & ((k == 64) ? ~(uint64_t) 0 : (((uint64_t) 1 << k) - 1));
#else
    slong j;
    r = 0;
    for (j = 0; j < k; j++)
    {
        slong p = pos + j;
        if (p >= 0 && ((M[p / FLINT_BITS] >> (p % FLINT_BITS)) & 1))
            r |= (uint64_t) 1 << j;
    }
    return r;
#endif
}

void
mp_real_get_dfloat(double * res, double * rad, slong n, const mp_real_t x)
{
    slong K, L, be, i, bl;
    ulong V[DFC_GET_LIMBS], M[DFC_GET_LIMBS], t[4];
    double r = 0.0;
    int neg;

    FLINT_ASSERT(n <= 8);

    for (i = 0; i < n; i++)
        res[i] = 0.0;

    /* the input radius: err ulps at B^(exp - size) */
    if (x->err != 0)
        r = _dfc_up((double) x->err, FLINT_BITS * (x->exp - x->size));

    if (x->size == 0 || n == 0)
    {
        if (x->size != 0)
        {
            /* |x| < B^exp */
            r += _dfc_up(1.0, FLINT_BITS * x->exp);
        }
        if (rad != NULL)
            *rad = (r == 0.0) ? 0.0 : r * (1 + 0x1p-52);
        return;
    }

    /* beyond the double range: the whole line */
    if (FLINT_BITS * (x->exp - 1) > 1024)
    {
        if (rad != NULL)
            *rad = HUGE_VAL;
        return;
    }

    /* the top K limbs of the mantissa, in a two's complement frame V of
       L = K + 1 limbs with bit 0 at 2^be; the tail below goes to the
       radius as 2^be */
    K = FLINT_MIN(x->size, (53 * n + 2 * FLINT_BITS) / FLINT_BITS + 1);
    FLINT_ASSERT(K + 1 <= DFC_GET_LIMBS);
    L = K + 1;
    flint_mpn_copyi(V, x->d + x->size - K, K);
    V[K] = 0;
    be = FLINT_BITS * (x->exp - K);
    for (i = 0; i < x->size - K; i++)
    {
        if (x->d[i] != 0)
        {
            r += _dfc_up(1.0, be);
            break;
        }
    }

    for (i = 0; i < n; i++)
    {
        uint64_t c;
        slong sh, q, tl;
        int vneg;

        bl = _dfc_abs(M, V, L, &vneg);
        if (bl == 0)
            break;
        /* the component is c 2^(be + sh), c the nearest 53-bit integer
           to |V| 2^-sh (54 bits read, the last one the rounding bit) */
        sh = bl - 53;
        if (be + sh < -1074 + 1)
            break;              /* subnormal: left to the radius */
        c = _dfc_top_bits(M, bl, 54);
        c = (c >> 1) + (c & 1);         /* round half up */
        /* c 2^(be + sh) >= 2^1024: overflow */
        if (be + bl > 1024 || (be + bl == 1024 && (c >> 53) != 0))
        {
            for (i = 0; i < n; i++)
                res[i] = 0.0;
            if (rad != NULL)
                *rad = HUGE_VAL;
            return;
        }
        res[i] = ldexp((double) c, (int) (be + sh));
        if (vneg)
            res[i] = -res[i];

        /* V -= (+-) c 2^sh, sh >= -1 here: c 2^sh = (c >> 1) 2^(sh + 1)
           when sh = -1 (c then even unless bl < 53, where V is exact
           and c = |V|) */
        if (sh < 0)
        {
            /* |V| has fewer than 53 bits: c = |V| exactly */
            flint_mpn_zero(V, L);
            break;
        }
        q = sh / FLINT_BITS;
        tl = _dfc_shifted(t, c, (int) (sh % FLINT_BITS));
        tl = FLINT_MIN(tl, L - q);
        if (vneg)
            mpn_add(V + q, V + q, L - q, t, tl);
        else
            mpn_sub(V + q, V + q, L - q, t, tl);
    }

    /* the remainder |V| 2^be */
    bl = _dfc_abs(M, V, L, &neg);
    if (bl != 0)
    {
        uint64_t top = _dfc_top_bits(M, bl, 53);
        r += _dfc_up((double) (top + 1), be + bl - 53);
    }

    if (x->negative)
        for (i = 0; i < n; i++)
            res[i] = -res[i];

    if (!(r <= DBL_MAX))
    {
        for (i = 0; i < n; i++)
            res[i] = 0.0;
        r = HUGE_VAL;
    }
    else if (r != 0.0)
        r *= (1 + 0x1p-52);

    if (rad != NULL)
        *rad = r;
}
