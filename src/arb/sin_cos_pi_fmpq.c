/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "arb.h"
#include "mp_real.h"

/*
    sin(pi x), cos(pi x) and tan(pi x) for rational x = p/q.

    x is first reduced exactly to v/w in [0, 1/4] (an octant and a
    symmetry; w has the odd part of q).  Then, in this order:

    1. Below ARB_PI_FMPQ_SQRT_PREC bits, square roots by arb_sqrt: the
       closed forms for w <= 6, and the half-angle steps for w = 2^k q0,
       q0 = 1, 3, 5, where _pi_fmpq_sqrt_wins (they beat mp_real's
       versions there, whose square roots are slower).
    2. A word-size q (in lowest terms) goes to
       mp_real_sin_cos_pi_ui_div_ui resp. mp_real_tan_pi_ui_div_ui where
       these evaluate it by an algebraic method (closed forms, half
       angles, the Chebyshev iteration: _mp_real_sin_cos_pi_ui_div_ui_
       algebraic), which beats the general kernel; they reduce p/q
       exactly and keep the relative accuracy of every value.
    3. Otherwise the functions of the ball pi v/w (full relative accuracy,
       as v/w is small or away from the zeros), whose evaluation runs the
       mp_real kernel; the tangent as s/c or c/s.
*/

/* the precision from which mp_real's closed forms and half-angle steps
   beat those through arb_sqrt (measured on 64-bit x86: mp_real's are
   1.5-3 times slower at 64-2048 bits, 10-20% at 8192, about even at
   12288, 0.6 times the time at 2^18) */
#define ARB_PI_FMPQ_SQRT_PREC 12288

/* what is wanted: one of sine and cosine, both, or the tangent (as the
   quotient of both) */
#define PI_FMPQ_ONE 0
#define PI_FMPQ_BOTH 1
#define PI_FMPQ_TAN 2

/* whether the square roots through arb_sqrt beat both mp_real and the
   kernel for the denominator w > 6 at prec bits.  For w = 2^k q0 the
   chain takes r = k - 1, k, k + 1 square roots for q0 = 1, 3, 5 (one
   more for a sine alone, one more for the complement of both); it wins
   from the precision of the table up to ARB_PI_FMPQ_SQRT_PREC, by 5-50%,
   and beyond r = 8 (both) resp. 10 (one) loses or is about even.  For w
   = 10 and 12, both values lose to mp_real's quadratic forms from 3072
   resp. 1536 bits, and the tangent to its closed forms for w = 8 and 12
   (and 10 from 768 bits).  Measured on 64-bit x86 through these
   functions against their mp_real and kernel paths, over numerators
   coprime to w. */
static int
_pi_fmpq_sqrt_wins(ulong w, slong prec, int want)
{
    static const short both_min[9] = {768, 768, 768, 768, 1024, 1536, 1536, 2048, 3072};
    static const short one_min[11] = {256, 256, 256, 512, 768, 1024, 1536, 1536, 2048, 3072, 6144};
    slong k, r;
    ulong q0;

    if (w <= 6 || prec >= ARB_PI_FMPQ_SQRT_PREC)
        return 0;

    k = flint_ctz(w);
    q0 = w >> k;
    if (q0 != 1 && q0 != 3 && q0 != 5)
        return 0;

    if (want == PI_FMPQ_TAN && (w == 8 || w == 12))
        return 0;

    if (w == 10)
    {
        if (want == PI_FMPQ_ONE)
            return prec >= 128;
        return prec >= ((want == PI_FMPQ_TAN) ? 128 : 256)
            && prec < ((want == PI_FMPQ_TAN) ? 768 : 3072);
    }

    if (w == 12)
        return (want == PI_FMPQ_ONE) ? (prec >= 256) : (prec >= 512 && prec < 1536);

    r = k + (slong) (q0 / 2) - 1;

    if (want == PI_FMPQ_ONE)
        return r <= 10 && prec >= one_min[r];
    else
        return r <= 8 && prec >= both_min[r];
}

/* whether the square roots apply to v/w (reduced) */
static int
_pi_fmpq_use_sqrt(const fmpz_t w, slong prec, int want)
{
    if (COEFF_IS_MPZ(*w))
        return 0;
    return *w <= 6 || _pi_fmpq_sqrt_wins(*w, prec, want);
}

/* the working limbs for mp_real at prec bits (mp_real adds the guard
   bits for the denominator) */
static slong
_pi_fmpq_limbs(slong prec)
{
    return (prec + 4 + FLINT_BITS - 1) / FLINT_BITS;
}

/* whether to use mp_real for x = p/q (in lowest terms), given the
   reduced denominator w, setting qq = q; want: 1 = sin, 2 = cos, 3 =
   both, 4 = tan */
static int
_pi_fmpq_use_mp_real(ulong * qq, const fmpq_t x, const fmpz_t w, slong prec, int want)
{
    const fmpz * q = fmpq_denref(x);
    ulong u;

    if (COEFF_IS_MPZ(*q) || *q >= (WORD(1) << (FLINT_BITS - 2)))
        return 0;

    u = *q;
    *qq = u;

    /* the rational values 0, +-1/2, +-1 and tan = 0, +-1 directly */
    if (u <= 2 || (want == 1 && u == 6) || (want == 2 && u == 3) || (want == 4 && u == 4))
        return 0;

    /* (below 12 limbs, mp_real's algebraic methods are the closed forms,
       of which only those of the tangent for q = 8, 12 beat the square
       roots or the ball: a quick exit) */
    if (prec < 768 && !(want == 4 && (u == 8 || u == 12)))
        return 0;

    if (prec < ARB_PI_FMPQ_SQRT_PREC && _pi_fmpq_use_sqrt(w, prec,
            (want <= 2) ? PI_FMPQ_ONE : (want == 3) ? PI_FMPQ_BOTH : PI_FMPQ_TAN))
        return 0;

    /* q = 5 at the reduced denominator 10 goes to mp_real's quadratic
       forms where these win (from 12 limbs; below, mp_real would use its
       kernel, and the functions of the ball are faster) */
    if (u <= 6 && prec < 768)
        return 0;

    return (want == 4) ? _mp_real_tan_pi_ui_div_ui_algebraic(u, _pi_fmpq_limbs(prec))
                       : _mp_real_sin_cos_pi_ui_div_ui_algebraic(u, _pi_fmpq_limbs(prec));
}

/* s, c (either NULL) by mp_real for x = p/q, q = qq a word */
static void
_pi_fmpq_mp_real(arb_t s, arb_t c, const fmpq_t x, ulong qq, slong prec)
{
    mp_real_t S, C;
    ulong p;

    /* p mod 2q: the period of sin and cos */
    p = fmpz_fdiv_ui(fmpq_numref(x), 2 * qq);

    mp_real_init(S);
    mp_real_init(C);
    mp_real_sin_cos_pi_ui_div_ui((s != NULL) ? S : NULL, (c != NULL) ? C : NULL, p, qq, _pi_fmpq_limbs(prec));
    if (s != NULL)
    {
        mp_real_get_arb(s, S);
        arb_set_round(s, s, prec);
    }
    if (c != NULL)
    {
        mp_real_get_arb(c, C);
        arb_set_round(c, c, prec);
    }
    mp_real_clear(S);
    mp_real_clear(C);
}

/* The square roots: cos(pi p/q) for 0 <= p/q <= 1/2 in lowest terms (or
   p = 0), q = 2^k q0 with q0 = 1, 3, 5: closed forms for q <= 6, and
   for even q the half-angle step cos(t/2) = sqrt((1 + cos t)/2), resp.
   sqrt((1 - cos(pi - t))/2) for t > pi/2, whose cancellation (1 - cos
   >= 2^(1 - 2 bits(q))) the working precision covers. */
static void
_pi_fmpq_sqrt_cos(arb_t c, ulong p, ulong q, slong prec)
{
    if (p == 0)
    {
        arb_one(c);
    }
    else if (q <= 6)
    {
        if (q == 2)         /* 1/2 */
        {
            arb_zero(c);
        }
        else if (q == 3)    /* 1/3 */
        {
            arb_one(c);
            arb_mul_2exp_si(c, c, -1);
        }
        else if (q == 4)    /* 1/4 */
        {
            arb_sqrt_ui(c, 2, prec);
            arb_mul_2exp_si(c, c, -1);
        }
        else if (q == 5)    /* 1/5, 2/5: (sqrt(5) +- 1)/4 */
        {
            arb_sqrt_ui(c, 5, prec + 3);
            arb_add_si(c, c, (p == 1) ? 1 : -1, prec);
            arb_mul_2exp_si(c, c, -2);
        }
        else                /* 1/6 */
        {
            FLINT_ASSERT(q == 6);
            arb_sqrt_ui(c, 3, prec);
            arb_mul_2exp_si(c, c, -1);
        }
    }
    else
    {
        slong extra = 2 * FLINT_BIT_COUNT(q) + 2;

        /* (odd q > 6 have no square roots: excluded by the dispatch) */
        FLINT_ASSERT(q % 2 == 0);

        if (4 * p <= q)
        {
            _pi_fmpq_sqrt_cos(c, p, q / 2, prec + extra);
            arb_add_ui(c, c, 1, prec + extra);
        }
        else
        {
            _pi_fmpq_sqrt_cos(c, q / 2 - p, q / 2, prec + extra);
            arb_sub_ui(c, c, 1, prec + extra);
            arb_neg(c, c);
        }

        arb_mul_2exp_si(c, c, -1);
        arb_sqrt(c, c, prec);
    }
}

/* sin(pi p/q) = cos(pi (q/2 - p)/q), as above */
static void
_pi_fmpq_sqrt_sin(arb_t s, ulong p, ulong q, slong prec)
{
    if (q % 2 == 0)
    {
        p = q / 2 - p;

        while (p % 2 == 0 && q % 2 == 0)
        {
            p /= 2;
            q /= 2;
        }

        _pi_fmpq_sqrt_cos(s, p, q, prec);
    }
    else
    {
        _pi_fmpq_sqrt_cos(s, q - 2 * p, 2 * q, prec);
    }
}

/* both: the chain for the larger value, the other as sqrt(1 - x^2) */
static void
_pi_fmpq_sqrt_sin_cos(arb_t s, arb_t c, ulong p, ulong q, slong prec)
{
    slong wp;

    if (p == 0)
    {
        arb_one(c);
        arb_zero(s);
        return;
    }

    if (q == 2)
    {
        arb_zero(c);
        arb_one(s);
        return;
    }

    if (q == 4)
    {
        arb_sqrt_ui(c, 2, prec);
        arb_mul_2exp_si(c, c, -1);
        arb_set(s, c);
        return;
    }

    wp = prec + 3;

    /* the formula with less cancellation */
    if (p <= q / 4)
    {
        _pi_fmpq_sqrt_sin(s, p, q, wp);
        arb_mul(c, s, s, wp);
        arb_sub_ui(c, c, 1, wp);
        arb_neg(c, c);
        arb_sqrt(c, c, prec);
    }
    else
    {
        _pi_fmpq_sqrt_cos(c, p, q, wp);
        arb_mul(s, c, c, wp);
        arb_sub_ui(s, s, 1, wp);
        arb_neg(s, s);
        arb_sqrt(s, s, prec);
    }
}

/* the functions of pi v/w for v/w in [0, 1/4] (s and c distinct; want:
   PI_FMPQ_BOTH or PI_FMPQ_TAN) */
static void
_arb_sin_cos_pi_fmpq_oct(arb_t s, arb_t c,
            const fmpz_t v, const fmpz_t w, slong prec, int want)
{
    if (_pi_fmpq_use_sqrt(w, prec, want))
    {
        _pi_fmpq_sqrt_sin_cos(s, c, *v, *w, prec);
    }
    else
    {
        arb_const_pi(s, prec);
        arb_mul_fmpz(s, s, v, prec);
        arb_div_fmpz(s, s, w, prec);
        arb_sin_cos(s, c, s, prec);
    }
}

static void
_arb_sin_pi_fmpq_oct(arb_t s, const fmpz_t v, const fmpz_t w, slong prec)
{
    if (_pi_fmpq_use_sqrt(w, prec, PI_FMPQ_ONE))
    {
        _pi_fmpq_sqrt_sin(s, *v, *w, prec);
    }
    else
    {
        arb_const_pi(s, prec);
        arb_mul_fmpz(s, s, v, prec);
        arb_div_fmpz(s, s, w, prec);
        arb_sin(s, s, prec);
    }
}

static void
_arb_cos_pi_fmpq_oct(arb_t c, const fmpz_t v, const fmpz_t w, slong prec)
{
    if (_pi_fmpq_use_sqrt(w, prec, PI_FMPQ_ONE))
    {
        _pi_fmpq_sqrt_cos(c, *v, *w, prec);
    }
    else
    {
        arb_const_pi(c, prec);
        arb_mul_fmpz(c, c, v, prec);
        arb_div_fmpz(c, c, w, prec);
        arb_cos(c, c, prec);
    }
}

static unsigned int
reduce_octant(fmpz_t v, fmpz_t w, const fmpq_t x)
{
    const fmpz * p = fmpq_numref(x);
    const fmpz * q = fmpq_denref(x);
    unsigned int octant;
    flint_bitcnt_t vval, wval;

    if (*p > COEFF_MIN / 8 &&
        *p < COEFF_MAX / 8 &&
        *q > 0             &&
        *q < COEFF_MAX / 4)
    {
        slong pp, qq, ww, vv, tt;

        pp = *p;
        qq = *q;

        tt = pp * 4;
        ww = tt / qq;
        vv = tt - qq * ww;
        /* compute correct (floor) quotient and remainder */
        if (vv < 0)
        {
            ww--;
            vv += qq;
        }

        octant = ((ulong) ww) % 8;
        ww = qq * 4;
        if (octant % 2 != 0)
            vv = qq - vv;

        if (vv != 0)
        {
            vval = flint_ctz(vv);
            wval = flint_ctz(ww);
            vval = FLINT_MIN(vval, wval);
            vv >>= vval;
            ww >>= vval;
        }

        fmpz_set_si(v, vv);
        fmpz_set_si(w, ww);
    }
    else
    {
        fmpz_mul_2exp(w, p, 2);
        fmpz_fdiv_qr(w, v, w, q);
        octant = fmpz_fdiv_ui(w, 8);
        fmpz_mul_2exp(w, q, 2);

        if (octant % 2 != 0)
            fmpz_sub(v, q, v);

        vval = fmpz_val2(v);
        wval = fmpz_val2(w);
        vval = FLINT_MIN(vval, wval);

        if (vval != 0)
        {
            fmpz_tdiv_q_2exp(v, v, vval);
            fmpz_tdiv_q_2exp(w, w, vval);
        }
    }

    return octant;
}

void
arb_sin_cos_pi_fmpq(arb_t s, arb_t c, const fmpq_t x, slong prec)
{
    fmpz_t v, w;
    unsigned int octant;
    ulong qq;

    fmpz_init(v);
    fmpz_init(w);

    octant = reduce_octant(v, w, x);

    if (_pi_fmpq_use_mp_real(&qq, x, w, prec, 3))
    {
        _pi_fmpq_mp_real(s, c, x, qq, prec);
    }
    else
    {
        if ((octant + 1) % 4 < 2)
            _arb_sin_cos_pi_fmpq_oct(s, c, v, w, prec, PI_FMPQ_BOTH);
        else
            _arb_sin_cos_pi_fmpq_oct(c, s, v, w, prec, PI_FMPQ_BOTH);

        if ((octant + 6) % 8 < 4)  arb_neg(c, c);
        if (octant >= 4)           arb_neg(s, s);
    }

    fmpz_clear(v);
    fmpz_clear(w);
}

void
arb_sin_pi_fmpq(arb_t s, const fmpq_t x, slong prec)
{
    fmpz_t v, w;
    unsigned int octant;
    ulong qq;

    fmpz_init(v);
    fmpz_init(w);

    octant = reduce_octant(v, w, x);

    if (_pi_fmpq_use_mp_real(&qq, x, w, prec, 1))
    {
        _pi_fmpq_mp_real(s, NULL, x, qq, prec);
    }
    else
    {
        if ((octant + 1) % 4 < 2)
            _arb_sin_pi_fmpq_oct(s, v, w, prec);
        else
            _arb_cos_pi_fmpq_oct(s, v, w, prec);

        if (octant >= 4)
            arb_neg(s, s);
    }

    fmpz_clear(v);
    fmpz_clear(w);
}

void
arb_cos_pi_fmpq(arb_t c, const fmpq_t x, slong prec)
{
    fmpz_t v, w;
    unsigned int octant;
    ulong qq;

    fmpz_init(v);
    fmpz_init(w);

    octant = reduce_octant(v, w, x);

    if (_pi_fmpq_use_mp_real(&qq, x, w, prec, 2))
    {
        _pi_fmpq_mp_real(NULL, c, x, qq, prec);
    }
    else
    {
        if ((octant + 1) % 4 < 2)
            _arb_cos_pi_fmpq_oct(c, v, w, prec);
        else
            _arb_sin_pi_fmpq_oct(c, v, w, prec);

        if ((octant + 6) % 8 < 4)
            arb_neg(c, c);
    }

    fmpz_clear(v);
    fmpz_clear(w);
}

void
arb_tan_pi_fmpq(arb_t t, const fmpq_t x, slong prec)
{
    fmpz_t v, w;
    unsigned int octant;
    ulong qq;

    fmpz_init(v);
    fmpz_init(w);

    /* x = 2k + o/4 + v/w (o even) resp. (o + 1)/4 - v/w (o odd), v/w in
       [0, 1/4]: tan(pi x) = T, 1/T, -1/T, -T for o mod 4 = 0, 1, 2, 3,
       T = tan(pi v/w) */
    octant = reduce_octant(v, w, x) % 4;

    if (_pi_fmpq_use_mp_real(&qq, x, w, prec, 4))
    {
        mp_real_t T;
        mp_real_init(T);
        /* (p mod q: the period of tan) */
        if (mp_real_tan_pi_ui_div_ui(T, fmpz_fdiv_ui(fmpq_numref(x), qq), qq, _pi_fmpq_limbs(prec)))
        {
            mp_real_get_arb(t, T);
            arb_set_round(t, t, prec);
        }
        else
            arb_indeterminate(t);
        mp_real_clear(T);
    }
    else if (fmpz_is_zero(v))
    {
        /* (v = 0 only for even o) 0, or a pole at o = 2 */
        if (octant == 0)
            arb_zero(t);
        else
            arb_indeterminate(t);
    }
    else
    {
        /* T or 1/T as s/c or c/s (one division) */
        arb_t c;
        slong wp = prec + 4;

        arb_init(c);
        _arb_sin_cos_pi_fmpq_oct(t, c, v, w, wp, PI_FMPQ_TAN);
        if (octant == 1 || octant == 2)
            arb_div(t, c, t, prec);
        else
            arb_div(t, t, c, prec);
        if (octant >= 2)
            arb_neg(t, t);
        arb_clear(c);
    }

    fmpz_clear(v);
    fmpz_clear(w);
}
