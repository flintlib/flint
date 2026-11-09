/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acf.h"
#include "elem.h"
#include "impl.h"

/*
    The complex elementary functions of acf: in the strict modes each
    part is correctly rounded like the results of the real functions of
    arf, whatever the size of the part against the other; in the relaxed
    modes each part is within an ulp (ARF_RND_FAST) resp. 0.51 ulp
    (ARF_RND_ACCURATE) of the larger part (_acf_elem_set).

    At one and two limbs by the register compositions of acf/elem.h;
    otherwise the compositions of the real functions of arf at p = prec +
    guard bits (not rounded up to whole limbs) in the relaxed mode
    ARF_RND_FAST, with the products and sums rounded down, so that each
    part comes out within K ulps at p bits of its own value (the
    functions within 1 ulp, every other operation within 1 ulp) -- plus,
    for pow, an absolute error from the argument y log x of the
    exponential.  The parts are then rounded by _arf_set_mpn_err,
    recomputed at a higher precision when the window does not decide the
    rounding (the strict modes) or the error bound is too large for the
    relaxed modes (arf/elem.h).
*/

enum
{
    CELEM_EXP,
    CELEM_EXP_PI_I,
    CELEM_SIN_COS,
    CELEM_SIN_COS_PI,
    CELEM_LOG,
    CELEM_POW,
    CELEM_TAN,
    CELEM_COT,
    CELEM_SEC,
    CELEM_CSC,
    CELEM_TAN_PI,
    CELEM_COT_PI,
    CELEM_SEC_PI,
    CELEM_CSC_PI
};

/* the exponent of x (zero for a special x), positive exponents only (0
   below), clamped to ACF_EXP_CLAMP: 2^40 (2^28 on 32-bit machines; the
   bignum exponents being far beyond: such a part is zero within any
   precision, or sends the exponential beyond its range) */
#if FLINT_BITS == 64
#define ACF_EXP_CLAMP (WORD(1) << 40)
#else
#define ACF_EXP_CLAMP (WORD(1) << 28)
#endif

static slong
_acf_exp_clamped(const arf_t x)
{
    if (arf_is_special(x))
        return 0;
    if (COEFF_IS_MPZ(ARF_EXP(x)))
        return (fmpz_sgn(ARF_EXPREF(x)) > 0) ? ACF_EXP_CLAMP : 0;
    return FLINT_MAX(FLINT_MIN(ARF_EXP(x), ACF_EXP_CLAMP), 0);
}

/* the relative error constants of the compositions: exp (exp, sin_cos
   and a product), exp_pi_i (one more for the product pi b), sin_cos
   (as exp), sin_cos_pi, log (the sum of squares and the logarithm, or
   atan2), pow (the exponential of the parts of t = y log x, with the
   absolute error of t accounted separately), tan, cot, sec, csc and
   their pi versions (see the composition).  A product of two factors
   within 1 ulp each, rounded once, is within 4 ulps of its own precision
   (the ulp of the product relative to those of the factors varies by a
   factor of 2); the measured maxima are 4.0 for exp and sin_cos, 3.98 for
   their pi versions, 1.0 for log, and up to 6.6 for the tangent family
   (whose worst-case bound is about 16 ulps at p bits for the pi
   versions) */
static const unsigned char _celem_err[] = { 4, 5, 4, 5, 3, 4, 13, 13, 13, 13, 15, 15, 15, 15 };

/* the parts of f(x) (f(x, y) for pow) at p bits, within the errors of
   _celem_err in ulps at p bits: R1 + R2 i, and for the sine-cosine
   pairs the cosine R3 + R4 i.  (p need not be a multiple of the limb
   size: the precision of the arf functions, which add their own guard
   bits, decides the limbs they work with.)  For pow, ebu[0] and ebu[1]
   receive the exponent of the absolute error of the parts (from that
   of t = y log x) in ulps at p bits of the larger of R1 and R2, or
   WORD_MIN, and *exact the bits of the parts known to be exact (the
   powers of real and imaginary numbers with dyadic values). */
/* T = pi x rounded down at wp bits (T may alias x), for an argument of
   the exponential functions at plim bits; when |x| >= 2^ARF_EXP_MAGLIM(plim)
   those give nan anyway, and T = 4x (beyond the limit as well) spares
   computing pi to the bits of the exponent of x */
static void
_acf_pi_mul(arf_t T, const arf_t x, slong wp, slong plim)
{
    if (!arf_is_special(x) && arf_cmpabs_2exp_si(x, ARF_EXP_MAGLIM(plim)) >= 0)
    {
        arf_mul_2exp_si(T, x, 2);
    }
    else
    {
        arf_t U;
        arf_init(U);
        arf_const_pi(U, wp, ARF_RND_FAST);
        arf_mul(T, U, x, wp, ARF_RND_DOWN);
        arf_clear(U);
    }
}

static void
_acf_elem_inner(arf_t R1, arf_t R2, arf_t R3, arf_t R4, int f, const arf_t a, const arf_t b,
    const arf_t ya, const arf_t yb, slong p, slong * ebu, int * exact)
{
    arf_t E, S, C, T, U;

    arf_init(E);
    arf_init(S);
    arf_init(C);
    arf_init(T);
    arf_init(U);
    ebu[0] = ebu[1] = WORD_MIN;
    *exact = 0;

    switch (f)
    {
        case CELEM_EXP:
            arf_exp(E, a, p, ARF_RND_FAST);
            arf_sin_cos(S, C, b, p, ARF_RND_FAST);
            arf_mul(R1, E, C, p, ARF_RND_DOWN);
            arf_mul(R2, E, S, p, ARF_RND_DOWN);
            break;

        case CELEM_EXP_PI_I:
            /* exp(-pi b) with the product at enough bits for its
               absolute error to lie below 2^-(p + 2) */
            _acf_pi_mul(T, b, p + _acf_exp_clamped(b) + 4, p);
            arf_neg(T, T);
            arf_exp(E, T, p, ARF_RND_FAST);
            arf_sin_cos_pi(S, C, a, p, ARF_RND_FAST);
            arf_mul(R1, E, C, p, ARF_RND_DOWN);
            arf_mul(R2, E, S, p, ARF_RND_DOWN);
            break;

        case CELEM_SIN_COS:
        case CELEM_SIN_COS_PI:
            /* sin(a + bi) = sin a cosh b + i cos a sinh b,
               cos(a + bi) = cos a cosh b - i sin a sinh b */
            if (f == CELEM_SIN_COS)
            {
                arf_sin_cos(S, C, a, p, ARF_RND_FAST);
                arf_sinh_cosh(T, U, b, p, ARF_RND_FAST);
            }
            else
            {
                arf_sin_cos_pi(S, C, a, p, ARF_RND_FAST);
                _acf_pi_mul(T, b, p + _acf_exp_clamped(b) + 4, p);
                arf_sinh_cosh(T, U, T, p, ARF_RND_FAST);
            }
            arf_mul(R1, S, U, p, ARF_RND_DOWN);
            arf_mul(R2, C, T, p, ARF_RND_DOWN);
            arf_mul(R3, C, U, p, ARF_RND_DOWN);
            arf_mul(R4, S, T, p, ARF_RND_DOWN);
            arf_neg(R4, R4);
            break;

        case CELEM_LOG:
            {
                arf_srcptr M = a, mm = b;

                /* log |z| = log(M^2 + m^2) / 2, M the larger part: from
                   the sum of squares S rounded to p + 10 bits, whose
                   relative error 2^-(p + 9) gives log S within 2^-9 /
                   |log S| ulps, less than 0.01 for |S - 1| >= 1/4; near
                   the unit circle (|S - 1| < 1/4) from D = (M^2 - 1) +
                   m^2 formed exactly (unless the terms are too far apart
                   to interact) and log1p */
                if (arf_cmpabs(b, a) > 0)
                {
                    M = b;
                    mm = a;
                }

                arf_sosq(T, a, b, p + 10, ARF_RND_DOWN);

                if ((ARF_EXP(T) == 0 || ARF_EXP(T) == 1) && fabs(arf_get_d(T, ARF_RND_DOWN) - 1.0) < 0.25)
                {
                    slong gap;

                    arf_mul(T, M, M, ARF_PREC_EXACT, ARF_RND_DOWN);
                    arf_sub_ui(T, T, 1, ARF_PREC_EXACT, ARF_RND_DOWN);
                    arf_mul(U, mm, mm, ARF_PREC_EXACT, ARF_RND_DOWN);
                    gap = (arf_is_zero(T) || arf_is_zero(U)) ? 0 :
                        FLINT_ABS(ARF_EXP(T) - ARF_EXP(U));
                    arf_add(T, T, U, (gap <= 4 * p + 128) ? ARF_PREC_EXACT : p + 64, ARF_RND_DOWN);

                    if (arf_cmpabs_2exp_si(T, -1) < 0)
                        arf_log1p(R1, T, p, ARF_RND_FAST);
                    else
                    {
                        arf_add_ui(T, T, 1, p + 64, ARF_RND_DOWN);
                        arf_log(R1, T, p, ARF_RND_FAST);
                    }
                }
                else
                    arf_log(R1, T, p, ARF_RND_FAST);

                arf_mul_2exp_si(R1, R1, -1);
                arf_atan2(R2, b, a, p, ARF_RND_FAST);
            }
            break;

        case CELEM_POW:
            if ((arf_is_zero(a) || arf_is_zero(b)) &&
                (arf_is_zero(yb) || arf_cmpabs_ui(arf_is_zero(b) ? a : b, 1) == 0))
            {
                /* x real or imaginary, arg x = k pi/2, with y real or
                   |x| = 1 (log |x| = 0): x^y = |x|^yr exp(-yi k pi/2)
                   (cos(k yr pi/2) + i sin(k yr pi/2)) -- the exact
                   zeros of the sine and cosine at the half-integers,
                   which the general composition could never certify */
                int k = arf_is_zero(b) ? (arf_sgn(a) > 0 ? 0 : 2) : (arf_sgn(b) > 0 ? 1 : -1);
                int einex, sinex, cinex;

                /* (the strict rounding of these factors reports the
                   exact values: dyadic powers, and the sine and cosine
                   at the half-integers; a product of exact factors or
                   with a zero factor is exact) */
                arf_mul_si(T, ya, k, ARF_PREC_EXACT, ARF_RND_DOWN);
                arf_mul_2exp_si(T, T, -1);
                sinex = arf_sin_cos_pi(S, NULL, T, p, ARF_RND_DOWN);
                cinex = arf_sin_cos_pi(NULL, C, T, p, ARF_RND_DOWN);
                if (arf_is_zero(yb))
                {
                    arf_abs(U, arf_is_zero(b) ? a : b);
                    einex = arf_pow(E, U, ya, p, ARF_RND_DOWN);
                }
                else
                {
                    arf_mul_si(T, yb, -k, ARF_PREC_EXACT, ARF_RND_DOWN);
                    arf_mul_2exp_si(T, T, -1);
                    _acf_pi_mul(T, T, p + _acf_exp_clamped(T) + 4, p);
                    arf_exp(E, T, p, ARF_RND_FAST);
                    einex = 1;
                }
                cinex |= arf_mul(R1, E, C, p, ARF_RND_DOWN);
                sinex |= arf_mul(R2, E, S, p, ARF_RND_DOWN);
                *exact = ((!einex && !cinex) || arf_is_zero(C)) | (((!einex && !sinex) || arf_is_zero(S)) << 1);
            }
            else
            {
                /* t = y log x at p + extra bits, extra covering the
                   sizes of the products (|log |x|| below |E(M)| + 1 and
                   |arg x| below 4, |y| below 2^E(y)) so that the parts
                   of t are within 2^-(p + 2) absolutely */
                slong ea, eb, pt, et;
                int stage;

                /* (the bit size of |E(M)| + 1, unclamped: |log |x||
                   reaches 2^bits(E) for exponents of any size) */
                ea = FLINT_MAX(arf_is_special(a) ? 0 : fmpz_bits(ARF_EXPREF(a)), arf_is_special(b) ? 0 : fmpz_bits(ARF_EXPREF(b))) + 1;
                eb = FLINT_MAX(_acf_exp_clamped(ya), _acf_exp_clamped(yb));

                /* (for a huge y the size of t is found from a first
                   evaluation at p + 64 bits, the logarithm of an x near
                   1 being small: the final precision then from |t|) */
                for (stage = 0; ; stage++)
                {
                    if (eb < (WORD(1) << 16) && ea < (WORD(1) << 16))
                        pt = p + ea + eb + 8;
                    else
                        pt = p + 64;
                    if (stage == 1)
                        pt = p + et + 8;

                    /* (the parts of log x within 3 ulps at pt bits:
                       the composition of acf_log) */
                    _acf_elem_inner(E, S, C, T, CELEM_LOG, a, b, NULL, NULL, pt, ebu, exact);
                    *exact = 0;

                    /* t = (yr L - yi Theta) + i (yr Theta + yi L) */
                    arf_mul(T, ya, E, pt, ARF_RND_DOWN);
                    arf_mul(U, yb, S, pt, ARF_RND_DOWN);
                    arf_sub(T, T, U, pt, ARF_RND_DOWN);
                    arf_mul(U, ya, S, pt, ARF_RND_DOWN);
                    arf_mul(C, yb, E, pt, ARF_RND_DOWN);
                    arf_add(U, U, C, pt, ARF_RND_DOWN);

                    if ((eb < (WORD(1) << 16) && ea < (WORD(1) << 16)) || stage == 1)
                        break;
                    et = FLINT_MAX(_acf_exp_clamped(T), _acf_exp_clamped(U));
                    if (et > ARF_EXP_MAGLIM(p))
                        break;      /* (the exponential overflows: nan) */
                }

                arf_exp(E, T, p, ARF_RND_FAST);
                arf_sin_cos(S, C, U, p, ARF_RND_FAST);
                arf_mul(R1, E, C, p, ARF_RND_DOWN);
                arf_mul(R2, E, S, p, ARF_RND_DOWN);
                /* the absolute error: 2^-(p + 2) through exp and the
                   sine and cosine, times |E|: 2^(E(E) - p - 1), which is
                   2^(E(E) - E(R) - 1) ulps at p bits of the larger part R
                   (the exponents may be bignums: their difference is
                   not) */
                if (!arf_is_special(E))
                {
                    arf_srcptr R = (arf_cmpabs(R1, R2) >= 0) ? R1 : R2;

                    if (!arf_is_zero(R))
                    {
                        fmpz_t d;
                        fmpz_init(d);
                        fmpz_sub(d, ARF_EXPREF(E), ARF_EXPREF(R));
                        ebu[0] = fmpz_fits_si(d) ? fmpz_get_si(d) - 1 : (fmpz_sgn(d) > 0 ? WORD_MAX / 2 : WORD_MIN);
                        fmpz_clear(d);
                    }
                    else
                        ebu[0] = WORD_MAX / 2;
                    ebu[1] = ebu[0];
                }
            }
            break;

        case CELEM_TAN:
        case CELEM_COT:
        case CELEM_SEC:
        case CELEM_CSC:
        case CELEM_TAN_PI:
        case CELEM_COT_PI:
        case CELEM_SEC_PI:
        case CELEM_CSC_PI:
            {
                /* from s = sin a, c = cos a, sh = sinh b, ch = cosh b (of
                   pi a and pi b for the pi versions):

                       tan(a + bi) = (s c + i sh ch) / (c^2 + sh^2)
                       cot(a + bi) = (s c - i sh ch) / (s^2 + sh^2)
                       sec(a + bi) = (c ch + i s sh) / (c^2 + sh^2)
                       csc(a + bi) = (s ch - i c sh) / (s^2 + sh^2)

                   the denominators |cos z|^2 and |sin z|^2 being sums of
                   squares, so that each part is a quotient of accurately
                   known numbers, also when it is exponentially small
                   against the other (large |b|: the real part of
                   tan(a + bi) is about 4 s c e^(-2|b|)).  Errors in ulps
                   at p bits: s, c, sh, ch 1 each (sh and ch 2 for the
                   rounded pi b, whose absolute error stays below
                   2^-(p + 3)), the squares 3 (5), the denominator 4 (6),
                   the numerators 3 (5), the quotients 8 (12), plus the
                   factors of 2 between the ulps of the factors and of
                   the results: K = 13 (15).  The real arguments, whose
                   results may be dyadic (tan(pi/4) = 1), are left to
                   the real functions. */
                int g = (f - CELEM_TAN) & 3, pi = (f >= CELEM_TAN_PI);
                slong pb, eb;

                if (pi)
                    arf_sin_cos_pi(S, C, a, p, ARF_RND_FAST);
                else
                    arf_sin_cos(S, C, a, p, ARF_RND_FAST);

                /* (a huge |b| needs the precision of the exponent of e^b:
                   ARF_EXP_MAGLIM; beyond the retry limit the results are
                   nan) */
                eb = _acf_exp_clamped(b);
                pb = FLINT_MAX(p, FLINT_MIN(eb / 2 + 4, ACF_ELEM_RETRY_LIMIT_BITS));
                if (pi)
                {
                    _acf_pi_mul(T, b, pb + eb + 6, pb);
                    arf_sinh_cosh(E, U, T, pb, ARF_RND_FAST);
                }
                else
                    arf_sinh_cosh(E, U, b, pb, ARF_RND_FAST);

                /* R3 = the denominator (kept for _acf_tan_delta), R4 scratch */
                if (g == 1 || g == 3)
                    arf_mul(R3, S, S, p, ARF_RND_DOWN);
                else
                    arf_mul(R3, C, C, p, ARF_RND_DOWN);
                arf_mul(R4, E, E, p, ARF_RND_DOWN);
                arf_add(R3, R3, R4, p, ARF_RND_DOWN);

                if (g == 0 || g == 1)
                {
                    arf_mul(R1, S, C, p, ARF_RND_DOWN);
                    arf_mul(R2, E, U, p, ARF_RND_DOWN);
                }
                else if (g == 2)
                {
                    arf_mul(R1, C, U, p, ARF_RND_DOWN);
                    arf_mul(R2, S, E, p, ARF_RND_DOWN);
                }
                else
                {
                    arf_mul(R1, S, U, p, ARF_RND_DOWN);
                    arf_mul(R2, C, E, p, ARF_RND_DOWN);
                }
                if (g == 1 || g == 3)
                    arf_neg(R2, R2);

                arf_div(R1, R1, R3, p, ARF_RND_DOWN);
                arf_div(R2, R2, R3, p, ARF_RND_DOWN);
            }
            break;
    }

    arf_clear(E);
    arf_clear(S);
    arf_clear(C);
    arf_clear(T);
    arf_clear(U);
}

/* For large |b| the imaginary part of tan and cot is exponentially close
   to -+1, which the directed modes could only resolve at a precision of
   about 2.9 |b| bits: there |Im| = 1 + delta with

       delta = -(cos 2a + e^(-2|b|)) / (2 (c^2 + sh^2))  (tan)
       delta = (cos 2a - e^(-2|b|)) / (2 (s^2 + sh^2))   (cot)

   (from sh ch - c^2 - sh^2 = -(cos 2a + e^(-2|b|)) / 2 and sh ch - s^2 -
   sh^2 = (cos 2a - e^(-2|b|)) / 2; of pi a and pi b for the pi
   versions), D the denominator of the composition at p bits.  Returns
   1 with delta in res when its sign is certain (the sum exceeding
   2^(4 - p) times the larger term: the errors of the terms and of the
   sum), else 0. */
static int
_acf_tan_delta(arf_t res, int f, const arf_t a, const arf_t b, const arf_t D, slong p)
{
    int g = (f - CELEM_TAN) & 3, pi = (f >= CELEM_TAN_PI), ok = 0;
    slong eb = _acf_exp_clamped(b), pb;
    arf_t C, E, T, U;

    arf_init(C);
    arf_init(E);
    arf_init(T);
    arf_init(U);

    pb = FLINT_MAX(p, FLINT_MIN(eb / 2 + 4, ACF_ELEM_RETRY_LIMIT_BITS));

    arf_mul_2exp_si(T, a, 1);
    if (pi)
        arf_cos_pi(C, T, p, ARF_RND_FAST);
    else
        arf_cos(C, T, p, ARF_RND_FAST);

    /* e^(-2|b|) (of pi b: rounded at pb + eb + 6 bits, as in the
       composition) */
    if (pi)
    {
        _acf_pi_mul(T, b, pb + eb + 6, pb);
    }
    else
        arf_set(T, b);
    arf_abs(T, T);
    arf_mul_2exp_si(T, T, 1);
    arf_neg(T, T);
    arf_exp(E, T, pb, ARF_RND_FAST);

    if (arf_is_finite(C) && arf_is_finite(E) && arf_is_finite(D))
    {
        if (g == 0)
            arf_add(U, C, E, p, ARF_RND_DOWN);
        else
            arf_sub(U, C, E, p, ARF_RND_DOWN);
        arf_set(T, (arf_cmpabs(C, E) >= 0) ? C : E);
        arf_mul_2exp_si(T, T, 4 - p);
        if (arf_cmpabs(U, T) > 0)
        {
            arf_div(res, U, D, p, ARF_RND_DOWN);
            arf_mul_2exp_si(res, res, -1);
            if (g == 0)
                arf_neg(res, res);
            ok = 1;
        }
    }

    arf_clear(C);
    arf_clear(E);
    arf_clear(T);
    arf_clear(U);
    return ok;
}

/* the exponent gap E(Rref) - E(R) >= 0 between a part R and the larger
   part Rref of the same complex number (clamped; 0 when either is zero),
   taken before either is rounded */
slong
_acf_part_gap(const arf_t R, const arf_t Rref)
{
    slong d;

    if (arf_is_special(R) || arf_is_special(Rref))
        return 0;

    if (!COEFF_IS_MPZ(ARF_EXP(R)) && !COEFF_IS_MPZ(ARF_EXP(Rref)))
        d = ARF_EXP(Rref) - ARF_EXP(R);
    else
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_sub(t, ARF_EXPREF(Rref), ARF_EXPREF(R));
        d = fmpz_fits_si(t) ? fmpz_get_si(t) : (fmpz_sgn(t) > 0 ? ACF_PART_GAP_MAX : 0);
        fmpz_clear(t);
    }

    return FLINT_MAX(0, FLINT_MIN(d, ACF_PART_GAP_MAX));
}

/* the arf R (finite, of at most p bits) rounded into res as the frame
   of m = ceil(p / FLINT_BITS) limbs, within K + 2 ulps at p bits of R
   plus 2^ebu ulps at p bits of the larger part Rref of the complex
   number (WORD_MIN for none), D = E(Rref) - E(R) (_acf_part_gap), or
   exactly; an ulp at p bits is 2^sh units of the frame, sh = FLINT_BITS
   m - p.  In the relaxed modes the error bounds are measured against an
   ulp of the larger part: a smaller part (D >= 2) is accepted with an
   error of up to 2^(D - 1) times that allowed for itself (its rounding
   then adds at most half an ulp of the larger part), and a zero R with
   an absolute error below 2^-4 (ARF_RND_FAST) resp. 2^-8
   (ARF_RND_ACCURATE) ulp of the larger part at prec bits becomes zero.  Returns the inexact flag or ARF_ELEM_UNCERTAIN. */
int
_acf_elem_set(arf_t res, const arf_t R, slong D, slong p, ulong K, slong ebu, int exact, slong prec, arf_rnd_t rnd)
{
    nn_srcptr rd;
    nn_ptr y;
    slong rn, eb, status, m, sh, d;
    int relaxed = (rnd >= ARF_RND_FAST);
    TMP_INIT;

    if (exact)
        return arf_set_round(res, R, prec, relaxed ? arf_rnd_relaxed_to_strict(rnd) : rnd);

    if (arf_is_zero(R))
    {
        /* (an exact zero of the composition: a zero value, barring the
           absolute error of pow, which the retries resolve in the strict
           modes) */
        if (ebu == WORD_MIN)
        {
            arf_zero(res);
            return 0;
        }
        if (relaxed && ebu + ((rnd == ARF_RND_FAST) ? 4 : 8) <= p - prec)
        {
            arf_zero(res);
            return 1;
        }
        return ARF_ELEM_UNCERTAIN;
    }

    /* (the absolute error in ulps of R) */
    if (ebu != WORD_MIN)
        ebu = (ebu > WORD_MAX / 4) ? WORD_MAX / 4 : ebu + D;

    d = (relaxed && D >= 2) ? D - 1 : 0;

    if (ebu != WORD_MIN && ebu - d > p)
        return ARF_ELEM_UNCERTAIN;

    m = (p + FLINT_BITS - 1) / FLINT_BITS;
    sh = FLINT_BITS * m - p;

    ARF_GET_MPN_READONLY(rd, rn, R);
    FLINT_ASSERT(rn <= m);
    TMP_START;
    y = TMP_ALLOC((m + 1) * sizeof(ulong));
    _flint_mpn_zero_small(y, m - rn);
    flint_mpn_copyi(y + m - rn, rd, rn);

    eb = ARF_ELEM_EB(K + 2) + sh;
    if (ebu != WORD_MIN)
        eb = FLINT_MAX(eb, ebu + sh) + 1;
    eb -= d;

    if (!COEFF_IS_MPZ(ARF_EXP(R)))
        status = _arf_set_mpn_err(res, y, m, ARF_EXP(R) - FLINT_BITS * m, eb, eb, ARF_SGNBIT(R), prec, rnd);
    else
    {
        fmpz_t e;
        fmpz_init(e);
        fmpz_sub_ui(e, ARF_EXPREF(R), FLINT_BITS * m);
        status = _arf_set_mpn_err(res, y, m, 0, eb, eb, ARF_SGNBIT(R), prec, rnd);
        if (status == 1)
            fmpz_add(ARF_EXPREF(res), ARF_EXPREF(res), e);
        fmpz_clear(e);
    }
    TMP_END;
    return status;
}

#if FLINT_BITS == 64
/* the frame into res (the zero frame, an underflowed product, declined) */
static int
_acf_frame_set(arf_t res, const _acf_frame_struct * f, int two, slong prec, arf_rnd_t rnd)
{
    if (two)
    {
        if ((f->y[3] | f->y[2]) == 0)
            return ARF_ELEM_UNCERTAIN;
        return _arf_small_set_2(res, f->y[3], f->y[2], f->y[1], f->y[0], f->k, f->err, f->sgnbit, prec, rnd);
    }
    else
    {
        if ((f->y[2] | f->y[1]) == 0)
            return ARF_ELEM_UNCERTAIN;
        return _arf_small_set_1(res, f->y[2], f->y[1], f->y[0], f->k, f->err, f->sgnbit, prec, rnd);
    }
}

/* the register paths for prec <= 2 FLINT_BITS and parts of at most two
   limbs with small exponents: the inexact flags of res1 (and res2), or
   ARF_ELEM_UNCERTAIN */
static int
_acf_elem_small(acf_t res1, acf_t res2, int f, const arf_t a, const arf_t b, const arf_t ya, const arf_t yb, slong prec, arf_rnd_t rnd1, arf_rnd_t rnd)
{
    _acf_frame_struct F1, F2, F3 = {{0, 0, 0, 0}, 0, 0, 0}, F4 = {{0, 0, 0, 0}, 0, 0, 0};
    nn_srcptr ad, bd;
    slong an, bn, ea, eb;
    ulong da1, da0, db1, db0;
    int two, asgn = ARF_SGNBIT(a), bsgn = ARF_SGNBIT(b), status, s1, s2, s3, s4;
    arf_t t1, t2, t3, t4;

    if (COEFF_IS_MPZ(ARF_EXP(a)) || COEFF_IS_MPZ(ARF_EXP(b)) || ARF_SIZE(a) > 2 || ARF_SIZE(b) > 2)
        return ARF_ELEM_UNCERTAIN;

    if (f == CELEM_POW)
    {
        /* (one limb only: x and the nonzero parts of y of one limb,
           with small exponents) */
        ulong dc = 0, dd = 0;
        slong ec = 0, ed = 0;

        if (prec > FLINT_BITS || ARF_SIZE(a) > 1 || ARF_SIZE(b) > 1 || ARF_SIZE(ya) > 1 || ARF_SIZE(yb) > 1
            || (!arf_is_zero(ya) && (COEFF_IS_MPZ(ARF_EXP(ya)) || FLINT_ABS(ARF_EXP(ya)) >= (WORD(1) << 40)))
            || (!arf_is_zero(yb) && (COEFF_IS_MPZ(ARF_EXP(yb)) || FLINT_ABS(ARF_EXP(yb)) >= (WORD(1) << 40)))
            || FLINT_ABS(ARF_EXP(a)) >= (WORD(1) << 40) || FLINT_ABS(ARF_EXP(b)) >= (WORD(1) << 40))
            return ARF_ELEM_UNCERTAIN;
        if (!arf_is_zero(ya))
        {
            dc = ARF_NOPTR_D(ya)[0];
            ec = ARF_EXP(ya);
        }
        if (!arf_is_zero(yb))
        {
            dd = ARF_NOPTR_D(yb)[0];
            ed = ARF_EXP(yb);
        }
        if (_acf_small_pow_1(&F1, &F2, ARF_NOPTR_D(a)[0], ARF_EXP(a), asgn, ARF_NOPTR_D(b)[0], ARF_EXP(b), bsgn,
                dc, ec, ARF_SGNBIT(ya), dd, ed, ARF_SGNBIT(yb)) != 0)
            return ARF_ELEM_UNCERTAIN;
        arf_init(t1);
        arf_init(t2);
        s1 = _acf_frame_set(t1, &F1, 0, prec, rnd1);
        s2 = (s1 == ARF_ELEM_UNCERTAIN) ? s1 : _acf_frame_set(t2, &F2, 0, prec, rnd);
        if (s1 == ARF_ELEM_UNCERTAIN || s2 == ARF_ELEM_UNCERTAIN)
            status = ARF_ELEM_UNCERTAIN;
        else
        {
            arf_swap(acf_realref(res1), t1);
            arf_swap(acf_imagref(res1), t2);
            status = s1 | (s2 << 1);
        }
        arf_clear(t1);
        arf_clear(t2);
        return status;
    }

    ea = ARF_EXP(a);
    eb = ARF_EXP(b);
    if (FLINT_ABS(ea) >= (WORD(1) << 40) || FLINT_ABS(eb) >= (WORD(1) << 40))
        return ARF_ELEM_UNCERTAIN;

    ARF_GET_MPN_READONLY(ad, an, a);
    ARF_GET_MPN_READONLY(bd, bn, b);
    two = !(prec <= FLINT_BITS && an == 1 && bn == 1);
    da1 = ad[an - 1];
    da0 = (an == 2) ? ad[0] : 0;
    db1 = bd[bn - 1];
    db0 = (bn == 2) ? bd[0] : 0;

    switch (f)
    {
        case CELEM_EXP:
        case CELEM_EXP_PI_I:
            status = two ? _acf_small_exp_2(&F1, &F2, da1, da0, ea, asgn, db1, db0, eb, bsgn, f == CELEM_EXP_PI_I)
                         : _acf_small_exp_1(&F1, &F2, da1, ea, asgn, db1, eb, bsgn, f == CELEM_EXP_PI_I);
            break;
        case CELEM_SIN_COS:
        case CELEM_SIN_COS_PI:
            status = two ? _acf_small_sin_cos_2(&F1, &F2, (res2 == NULL) ? NULL : &F3, &F4, da1, da0, ea, asgn, db1, db0, eb, bsgn, f == CELEM_SIN_COS_PI)
                         : _acf_small_sin_cos_1(&F1, &F2, (res2 == NULL) ? NULL : &F3, &F4, da1, ea, asgn, db1, eb, bsgn, f == CELEM_SIN_COS_PI);
            break;
        case CELEM_LOG:
            status = two ? _acf_small_log_2(&F1, &F2, da1, da0, ea, asgn, db1, db0, eb, bsgn)
                         : _acf_small_log_1(&F1, &F2, da1, ea, asgn, db1, eb, bsgn);
            break;
        default:
            status = ARF_ELEM_UNCERTAIN;
    }

    if (status != 0)
        return ARF_ELEM_UNCERTAIN;

    /* (into temporaries: res may alias the input, which a retry rereads) */
    arf_init(t1);
    arf_init(t2);
    arf_init(t3);
    arf_init(t4);
    s1 = _acf_frame_set(t1, &F1, two, prec, rnd1);
    s2 = (s1 == ARF_ELEM_UNCERTAIN) ? s1 : _acf_frame_set(t2, &F2, two, prec, rnd);
    s3 = s4 = 0;
    if (res2 != NULL && s2 != ARF_ELEM_UNCERTAIN)
    {
        s3 = _acf_frame_set(t3, &F3, two, prec, rnd);
        s4 = (s3 == ARF_ELEM_UNCERTAIN) ? s3 : _acf_frame_set(t4, &F4, two, prec, rnd);
    }

    if (s1 == ARF_ELEM_UNCERTAIN || s2 == ARF_ELEM_UNCERTAIN || s3 == ARF_ELEM_UNCERTAIN || s4 == ARF_ELEM_UNCERTAIN)
        status = ARF_ELEM_UNCERTAIN;
    else
    {
        arf_swap(acf_realref(res1), t1);
        arf_swap(acf_imagref(res1), t2);
        if (res2 != NULL)
        {
            arf_swap(acf_realref(res2), t3);
            arf_swap(acf_imagref(res2), t4);
        }
        status = s1 | (s2 << 1) | (s3 << 2) | (s4 << 3);
    }

    arf_clear(t1);
    arf_clear(t2);
    arf_clear(t3);
    arf_clear(t4);
    return status;
}
#endif

/* f(x) (f(x, y) for pow) for finite x (and y), which passed the special
   cases of the public functions; res1 (and res2 for the sine-cosine
   pairs) may alias x; returns the inexact flags, bits 0 and 1 for res1
   and 2 and 3 for res2.  The real part of res1 is rounded in rnd1, the
   others in rnd (the rotation of sinh negates that part afterwards:
   its directed modes flipped). */
static int
_acf_elem(acf_t res1, acf_t res2, int f, const acf_t x, const acf_t y, slong prec, arf_rnd_t rnd1, arf_rnd_t rnd)
{
    arf_srcptr a = acf_realref(x), b = acf_imagref(x);
    arf_srcptr ya = (y == NULL) ? NULL : acf_realref(y), yb = (y == NULL) ? NULL : acf_imagref(y);
    arf_t R1, R2, R3, R4;
    slong extra, p, ebu[2], D1, D2;
    ulong K;
    int status, s1, s2, s3, s4, exact;

    if (prec == ARF_PREC_EXACT)
        flint_throw(FLINT_ERROR, "acf elementary function: exact precision not supported\n");

#if FLINT_BITS == 64
    if (f <= CELEM_POW && prec <= 2 * FLINT_BITS && !arf_is_zero(a) && !arf_is_zero(b))
    {
        status = _acf_elem_small(res1, res2, f, a, b, ya, yb, prec, rnd1, rnd);
        if (status != ARF_ELEM_UNCERTAIN)
            return status;
    }
#endif

    K = _celem_err[f];
    arf_init(R1);
    arf_init(R2);
    arf_init(R3);
    arf_init(R4);

    for (extra = 0; ; extra = (extra == 0) ? 32 : 2 * extra)
    {
        /* (the working precision in bits: rounding it up to whole
           limbs would push the arf functions, which add their own
           guard bits, to one more limb) */
        p = _arf_elem_bits(prec, rnd, extra + FLINT_BIT_COUNT(K + 2));

        _acf_elem_inner(R1, R2, R3, R4, f, a, b, ya, yb, p, ebu, &exact);

        if (!arf_is_finite(R1) || !arf_is_finite(R2) || (res2 != NULL && (!arf_is_finite(R3) || !arf_is_finite(R4))))
        {
            arf_nan(acf_realref(res1));
            arf_nan(acf_imagref(res1));
            if (res2 != NULL)
            {
                arf_nan(acf_realref(res2));
                arf_nan(acf_imagref(res2));
            }
            status = 0;
            break;
        }

        /* beyond the limit of the retries (a dyadic value of pow the
           exactness tests missed, say) the approximations are rounded
           as they stand */
        if (p > 4 * prec + ACF_ELEM_RETRY_LIMIT_BITS)
        {
            exact = 15;
            if (rnd1 >= ARF_RND_FAST) rnd1 = arf_rnd_relaxed_to_strict(rnd1);
            if (rnd >= ARF_RND_FAST) rnd = arf_rnd_relaxed_to_strict(rnd);
        }

        /* (into the temporaries themselves, swapped in at the end: res
           may alias x, which a retry rereads) */
        /* (the gaps to the larger parts before either is rounded) */
        D1 = _acf_part_gap(R1, R2);
        D2 = _acf_part_gap(R2, R1);
        s1 = _acf_elem_set(R1, R1, D1, p, K, ebu[0], exact & 1, prec, rnd1);
        /* (tan and cot: the imaginary part -+(1 + delta) from delta
           when the window fails, R3 holding the denominator) */
        if (rnd < ARF_RND_FAST && (f == CELEM_TAN || f == CELEM_COT || f == CELEM_TAN_PI || f == CELEM_COT_PI)
            && !arf_is_zero(R2) && arf_cmpabs_2exp_si(b, 0) >= 0 && arf_cmpabs_2exp_si(R2, -1) > 0)
        {
            s2 = _acf_elem_set(R2, R2, D2, p, K, ebu[1], exact & 2, prec, rnd);
            /* (R2 is left alone when uncertain) */
            if (s2 == ARF_ELEM_UNCERTAIN && _acf_tan_delta(R3, f, a, b, R3, p) && arf_cmpabs_2exp_si(R3, -prec - 3) < 0)
                s2 = _arf_one_plus_tiny(R2, ARF_SGNBIT(R2), arf_sgn(R3) < 0, prec, rnd);
        }
        else if (f == CELEM_LOG)
        {
            /* (the imaginary part of the logarithm is atan2(b, a), a
               real function of exact numbers: rounded by arf, which
               also decides the values near a boundary, arg x near
               b / a for a tiny b / a) */
            s2 = arf_atan2(R2, b, a, prec, rnd);
        }
        else
            s2 = _acf_elem_set(R2, R2, D2, p, K, ebu[1], exact & 2, prec, rnd);
        s3 = s4 = 0;
        if (res2 != NULL)
        {
            D1 = _acf_part_gap(R3, R4);
            D2 = _acf_part_gap(R4, R3);
            s3 = _acf_elem_set(R3, R3, D1, p, K, WORD_MIN, exact & 4, prec, rnd);
            s4 = _acf_elem_set(R4, R4, D2, p, K, WORD_MIN, exact & 8, prec, rnd);
        }

        /* parts left undecided at a working precision well above prec:
           near a rounding boundary, from the expansion of f at the
           nearest exact point (acf/near.c) when x is close to one */
        if ((s1 == ARF_ELEM_UNCERTAIN || s2 == ARF_ELEM_UNCERTAIN || s3 == ARF_ELEM_UNCERTAIN || s4 == ARF_ELEM_UNCERTAIN)
            && p >= prec + 40 && f != CELEM_POW)
        {
            static const int fn1[] = { ACF_NEAR_EXP, ACF_NEAR_EXP_PI_I, ACF_NEAR_SIN, ACF_NEAR_SIN_PI, ACF_NEAR_LOG,
                ACF_NEAR_NONE, ACF_NEAR_TAN, ACF_NEAR_COT, ACF_NEAR_SEC, ACF_NEAR_CSC,
                ACF_NEAR_TAN_PI, ACF_NEAR_COT_PI, ACF_NEAR_SEC_PI, ACF_NEAR_CSC_PI };
            int * st[4] = { &s1, &s2, &s3, &s4 };
            arf_ptr Rp[4] = { R1, R2, R3, R4 };
            arf_t beta;
            int i, side;

            arf_init(beta);
            for (i = 0; i < 4; i++)
            {
                arf_rnd_t ri = (i == 0) ? rnd1 : rnd;
                int fn = (i < 2) ? fn1[f] : (f == CELEM_SIN_COS) ? ACF_NEAR_COS : ACF_NEAR_COS_PI;

                if (*st[i] != ARF_ELEM_UNCERTAIN || ri >= ARF_RND_FAST || arf_is_special(Rp[i]))
                    continue;
                _acf_round_boundary(beta, Rp[i], prec, ri);
                if (_acf_near_side(&side, fn, i & 1, x, beta, 0))
                    *st[i] = _acf_round_side(Rp[i], beta, side, prec, ri);
            }
            arf_clear(beta);
        }

        if (s1 != ARF_ELEM_UNCERTAIN && s2 != ARF_ELEM_UNCERTAIN && s3 != ARF_ELEM_UNCERTAIN && s4 != ARF_ELEM_UNCERTAIN)
        {
            arf_swap(acf_realref(res1), R1);
            arf_swap(acf_imagref(res1), R2);
            if (res2 != NULL)
            {
                arf_swap(acf_realref(res2), R3);
                arf_swap(acf_imagref(res2), R4);
            }
            status = s1 | (s2 << 1) | (s3 << 2) | (s4 << 3);
            break;
        }
    }

    arf_clear(R1);
    arf_clear(R2);
    arf_clear(R3);
    arf_clear(R4);
    return status;
}

/* The public functions: the special values and the real and imaginary
   arguments (the real functions of arf), then the compositions. */

static int
_acf_nan(acf_t res)
{
    arf_nan(acf_realref(res));
    arf_nan(acf_imagref(res));
    return 0;
}

/* whether a part is a nan or an infinity */
#define ACF_PART_NONFINITE(x) (arf_is_special(x) && !arf_is_zero(x))

int
acf_exp(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
{
    arf_srcptr a = acf_realref(x), b = acf_imagref(x);
    int s;

    if (arf_is_zero(b))
    {
        s = arf_exp(acf_realref(res), a, prec, rnd);
        arf_zero(acf_imagref(res));
        return s;
    }
    if (ACF_PART_NONFINITE(a) || ACF_PART_NONFINITE(b))
        return _acf_nan(res);
    if (arf_is_zero(a))
    {
        /* (cos b + i sin b, computed into the right parts whatever the
           aliasing) */
        arf_t s1, c1;
        arf_init(s1);
        arf_init(c1);
        arf_sin_cos(s1, c1, b, prec, rnd);
        arf_swap(acf_realref(res), c1);
        arf_swap(acf_imagref(res), s1);
        arf_clear(s1);
        arf_clear(c1);
        /* (cos b and sin b are irrational for b != 0) */
        return 3;
    }
    return _acf_elem(res, NULL, CELEM_EXP, x, NULL, prec, rnd, rnd);
}

int
acf_exp_pi_i(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
{
    arf_srcptr a = acf_realref(x), b = acf_imagref(x);
    int s;

    if (ACF_PART_NONFINITE(a) || ACF_PART_NONFINITE(b))
        return _acf_nan(res);
    if (arf_is_zero(b))
    {
        arf_t s1, c1;
        arf_init(s1);
        arf_init(c1);
        s = arf_sin_cos_pi(s1, c1, a, prec, rnd);
        arf_swap(acf_realref(res), c1);
        arf_swap(acf_imagref(res), s1);
        arf_clear(s1);
        arf_clear(c1);
        /* (both parts are exact, at the multiples of 1/2, or neither) */
        return (s != 0) ? 3 : 0;
    }
    /* (exp(-pi b) for a = 0 by the composition, whose sine is zero) */
    return _acf_elem(res, NULL, CELEM_EXP_PI_I, x, NULL, prec, rnd, rnd);
}

int
acf_log(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
{
    arf_srcptr a = acf_realref(x), b = acf_imagref(x);
    int s;

    if (ACF_PART_NONFINITE(a) || ACF_PART_NONFINITE(b))
        return _acf_nan(res);

    if (arf_is_zero(b))
    {
        if (arf_sgn(a) >= 0)
        {
            s = arf_log(acf_realref(res), a, prec, rnd);
            arf_zero(acf_imagref(res));
        }
        else
        {
            arf_t t;
            arf_init(t);
            arf_neg(t, a);
            s = arf_log(acf_realref(res), t, prec, rnd);
            s |= arf_const_pi(acf_imagref(res), prec, rnd) << 1;
            arf_clear(t);
        }
        return s;
    }

    if (arf_is_zero(a))
    {
        arf_t t;
        int neg = arf_sgn(b) < 0;
        arf_init(t);
        arf_abs(t, b);
        s = arf_log(acf_realref(res), t, prec, rnd);
        s |= arf_const_pi_mul_2exp(acf_imagref(res), -1, neg, prec, rnd) << 1;
        arf_clear(t);
        return s;
    }

    return _acf_elem(res, NULL, CELEM_LOG, x, NULL, prec, rnd, rnd);
}

/* sin, cos (either NULL) of x, or of pi x; the real part of the sine
   rounded in rnd1 */
static int
_acf_sin_cos(acf_t rs, acf_t rc, const acf_t x, int pi, slong prec, arf_rnd_t rnd1, arf_rnd_t rnd)
{
    arf_srcptr a = acf_realref(x), b = acf_imagref(x);
    int s = 0;

    if (ACF_PART_NONFINITE(a) || ACF_PART_NONFINITE(b))
    {
        if (rs != NULL) _acf_nan(rs);
        if (rc != NULL) _acf_nan(rc);
        return 0;
    }

    if (arf_is_zero(b))
    {
        arf_t s1, c1;
        arf_init(s1);
        arf_init(c1);
        int ss = 0, sc = 0;
        if (pi)
        {
            if (rs != NULL) ss = (arf_sin_cos_pi(s1, NULL, a, prec, rnd1) != 0);
            if (rc != NULL) sc = (arf_sin_cos_pi(NULL, c1, a, prec, rnd) != 0);
        }
        else
        {
            if (rs != NULL) ss = (arf_sin_cos(s1, NULL, a, prec, rnd1) != 0);
            if (rc != NULL) sc = (arf_sin_cos(NULL, c1, a, prec, rnd) != 0);
        }
        /* (the flags of the real parts: bit 0 of each result) */
        s = (rs != NULL) ? (ss | (sc << 2)) : sc;
        if (rs != NULL)
        {
            arf_swap(acf_realref(rs), s1);
            arf_zero(acf_imagref(rs));
        }
        if (rc != NULL)
        {
            arf_swap(acf_realref(rc), c1);
            arf_zero(acf_imagref(rc));
        }
        arf_clear(s1);
        arf_clear(c1);
        return s;
    }

    if (arf_is_zero(a) && !pi)
    {
        /* sin(bi) = i sinh b, cos(bi) = cosh b */
        arf_t s1, c1;
        arf_init(s1);
        arf_init(c1);
        arf_sinh_cosh((rs != NULL) ? s1 : NULL, (rc != NULL) ? c1 : NULL, b, prec, rnd);
        /* (sinh b and cosh b are irrational for b != 0: the imaginary
           part of the sine, the real part of the cosine) */
        s = (rs != NULL) ? (2 | ((rc != NULL) ? 4 : 0)) : 1;
        if (rs != NULL)
        {
            arf_swap(acf_imagref(rs), s1);
            arf_zero(acf_realref(rs));
        }
        if (rc != NULL)
        {
            arf_swap(acf_realref(rc), c1);
            arf_zero(acf_imagref(rc));
        }
        arf_clear(s1);
        arf_clear(c1);
        return s;
    }

    if (rs != NULL && rc != NULL)
        return _acf_elem(rs, rc, pi ? CELEM_SIN_COS_PI : CELEM_SIN_COS, x, NULL, prec, rnd1, rnd);

    {
        /* one of the pair: the compositions compute both */
        acf_t t1, t2;
        acf_init(t1);
        acf_init(t2);
        s = _acf_elem(t1, t2, pi ? CELEM_SIN_COS_PI : CELEM_SIN_COS, x, NULL, prec, rnd1, rnd);
        if (rs != NULL)
        {
            acf_swap(rs, t1);
            s &= 3;
        }
        else
        {
            acf_swap(rc, t2);
            s >>= 2;
        }
        acf_clear(t1);
        acf_clear(t2);
        return s;
    }
}

int acf_sin(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd) { return _acf_sin_cos(res, NULL, x, 0, prec, rnd, rnd); }
int acf_cos(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd) { return _acf_sin_cos(NULL, res, x, 0, prec, rnd, rnd); }
int acf_sin_cos(acf_t res1, acf_t res2, const acf_t x, slong prec, arf_rnd_t rnd) { return _acf_sin_cos(res1, res2, x, 0, prec, rnd, rnd); }
int acf_sin_pi(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd) { return _acf_sin_cos(res, NULL, x, 1, prec, rnd, rnd); }
int acf_cos_pi(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd) { return _acf_sin_cos(NULL, res, x, 1, prec, rnd, rnd); }
int acf_sin_cos_pi(acf_t res1, acf_t res2, const acf_t x, slong prec, arf_rnd_t rnd) { return _acf_sin_cos(res1, res2, x, 1, prec, rnd, rnd); }

/* sinh z = -i sin(iz), cosh z = cos(iz) with iz = -b + ai */
static int
_acf_sinh_cosh(acf_t rs, acf_t rc, const acf_t x, slong prec, arf_rnd_t rnd)
{
    acf_t iz, t;
    int s;

    arf_init_neg_shallow(acf_realref(iz), acf_imagref(x));
    arf_init_set_shallow(acf_imagref(iz), acf_realref(x));

    /* (through a temporary for the sine, whose parts are swapped and
       one negated: -i (u + vi) = v - ui; and since iz aliases x) */
    acf_init(t);
    s = _acf_sin_cos((rs != NULL) ? t : NULL, rc, iz, 0, prec, _acf_rnd_flip(rnd), rnd);
    if (rs != NULL)
    {
        arf_swap(acf_realref(rs), acf_imagref(t));
        arf_swap(acf_imagref(rs), acf_realref(t));
        arf_neg(acf_imagref(rs), acf_imagref(rs));
        s = ((s & 1) << 1) | ((s >> 1) & 1) | (s & 12);
    }
    acf_clear(t);
    return s;
}

int acf_sinh(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd) { return _acf_sinh_cosh(res, NULL, x, prec, rnd); }
int acf_cosh(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd) { return _acf_sinh_cosh(NULL, res, x, prec, rnd); }
int acf_sinh_cosh(acf_t res1, acf_t res2, const acf_t x, slong prec, arf_rnd_t rnd) { return _acf_sinh_cosh(res1, res2, x, prec, rnd); }

/* tan, cot, sec, csc (g = 0, 1, 2, 3) of x or of pi x, the real part
   rounded in rnd1 (the rotations of the hyperbolic functions negate it);
   the real arguments by the real functions */
static int
_acf_tan_family(acf_t res, const acf_t x, int g, int pi, slong prec, arf_rnd_t rnd1, arf_rnd_t rnd)
{
    arf_srcptr a = acf_realref(x), b = acf_imagref(x);
    int s;

    if (ACF_PART_NONFINITE(a) || ACF_PART_NONFINITE(b))
        return _acf_nan(res);

    if (arf_is_zero(b))
    {
        static int (* const fr[2][4])(arf_ptr, arf_srcptr, slong, arf_rnd_t) = {
            { arf_tan, arf_cot, arf_sec, arf_csc },
            { arf_tan_pi, arf_cot_pi, arf_sec_pi, arf_csc_pi } };

        s = fr[pi][g](acf_realref(res), a, prec, rnd1);
        arf_zero(acf_imagref(res));
        return s;
    }

    return _acf_elem(res, NULL, CELEM_TAN + g + 4 * pi, x, NULL, prec, rnd1, rnd);
}

int acf_tan(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd) { return _acf_tan_family(res, x, 0, 0, prec, rnd, rnd); }
int acf_cot(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd) { return _acf_tan_family(res, x, 1, 0, prec, rnd, rnd); }
int acf_sec(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd) { return _acf_tan_family(res, x, 2, 0, prec, rnd, rnd); }
int acf_csc(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd) { return _acf_tan_family(res, x, 3, 0, prec, rnd, rnd); }
int acf_tan_pi(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd) { return _acf_tan_family(res, x, 0, 1, prec, rnd, rnd); }
int acf_cot_pi(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd) { return _acf_tan_family(res, x, 1, 1, prec, rnd, rnd); }
int acf_sec_pi(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd) { return _acf_tan_family(res, x, 2, 1, prec, rnd, rnd); }
int acf_csc_pi(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd) { return _acf_tan_family(res, x, 3, 1, prec, rnd, rnd); }

/* tanh z = -i tan(iz), coth z = -i cot(-iz), sech z = sec(iz),
   csch z = -i csc(-iz) (cot and csc being odd), with iz = -b + ai and
   -iz = b - ai; the real arguments by the real functions */
static int
_acf_tanh_family(acf_t res, const acf_t x, int g, slong prec, arf_rnd_t rnd)
{
    arf_srcptr a = acf_realref(x), b = acf_imagref(x);
    acf_t w, t;
    int s;

    if (ACF_PART_NONFINITE(a) || ACF_PART_NONFINITE(b))
        return _acf_nan(res);

    if (arf_is_zero(b))
    {
        static int (* const fr[4])(arf_ptr, arf_srcptr, slong, arf_rnd_t) = { arf_tanh, arf_coth, arf_sech, arf_csch };

        s = fr[g](acf_realref(res), a, prec, rnd);
        arf_zero(acf_imagref(res));
        return s;
    }

    if (g == 0)
    {
        arf_init_neg_shallow(acf_realref(w), b);
        arf_init_set_shallow(acf_imagref(w), a);
    }
    else if (g == 2)
    {
        arf_init_neg_shallow(acf_realref(w), b);
        arf_init_set_shallow(acf_imagref(w), a);
        return _acf_tan_family(res, w, 2, 0, prec, rnd, rnd);
    }
    else
    {
        arf_init_set_shallow(acf_realref(w), b);
        arf_init_neg_shallow(acf_imagref(w), a);
    }

    /* -i (u + vi) = v - ui: through a temporary (w aliases x) */
    acf_init(t);
    s = _acf_tan_family(t, w, g, 0, prec, _acf_rnd_flip(rnd), rnd);
    arf_swap(acf_realref(res), acf_imagref(t));
    arf_swap(acf_imagref(res), acf_realref(t));
    arf_neg(acf_imagref(res), acf_imagref(res));
    acf_clear(t);
    return ((s & 1) << 1) | ((s >> 1) & 1);
}

int acf_tanh(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd) { return _acf_tanh_family(res, x, 0, prec, rnd); }
int acf_coth(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd) { return _acf_tanh_family(res, x, 1, prec, rnd); }
int acf_sech(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd) { return _acf_tanh_family(res, x, 2, prec, rnd); }
int acf_csch(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd) { return _acf_tanh_family(res, x, 3, prec, rnd); }

/* whether |E(a) - E(b)| > gap for nonzero finite a, b (bignum exponents
   allowed) */
static int
_acf_exp_gap_exceeds(const arf_t a, const arf_t b, slong gap)
{
    fmpz_t d;
    int r;
    fmpz_init(d);
    fmpz_sub(d, ARF_EXPREF(a), ARF_EXPREF(b));
    fmpz_abs(d, d);
    r = (fmpz_cmp_si(d, gap) > 0);
    fmpz_clear(d);
    return r;
}

/* the exact square root of x (nonzero, with dyadic parts) into res, or 0
   when it is not dyadic: sqrt(a + bi) = (q, b / (2q)) for a >= 0 and
   (|b| / (2q), sgn(b) q) for a < 0, q = sqrt((r + |a|) / 2) with
   r = |a + bi|, every root and quotient tested for exactness (a square
   root at the precision of its exact value, the quotient at that of
   the dividend) */
int
_acf_sqrt_exact(acf_t res, const acf_t x)
{
    arf_srcptr a = acf_realref(x), b = acf_imagref(x);
    arf_t r, q, t;
    int ok;

    arf_init(r);
    arf_init(q);
    arf_init(t);

    if (arf_is_zero(b))
    {
        arf_abs(t, a);
        ok = !arf_sqrt(q, t, arf_bits(t) / 2 + 2, ARF_RND_DOWN);
        if (ok)
        {
            if (arf_sgn(a) > 0)
            {
                arf_set(acf_realref(res), q);
                arf_zero(acf_imagref(res));
            }
            else
            {
                arf_zero(acf_realref(res));
                arf_set(acf_imagref(res), q);
            }
        }
    }
    else if (arf_is_zero(a))
    {
        /* sqrt(bi) = sqrt(|b| / 2) (1 +- i) */
        arf_abs(t, b);
        arf_mul_2exp_si(t, t, -1);
        ok = !arf_sqrt(q, t, arf_bits(t) / 2 + 2, ARF_RND_DOWN);
        if (ok)
        {
            arf_set(acf_realref(res), q);
            if (arf_sgn(b) > 0)
                arf_set(acf_imagref(res), q);
            else
                arf_neg(acf_imagref(res), q);
        }
    }
    else if (_acf_exp_gap_exceeds(a, b, 2 * (arf_bits(a) + arf_bits(b)) + 8))
    {
        /* (c + di)^2 = a + bi with dyadic c, d: a = c^2 - d^2 spans from
           the larger to the smaller of c^2 and d^2, b = 2cd lies between,
           so the exponents of a and b differ by at most about bits(a);
           a larger gap is no square (and a^2 + b^2 would be huge) */
        ok = 0;
    }
    else
    {
        arf_sosq(t, a, b, ARF_PREC_EXACT, ARF_RND_DOWN);
        ok = !arf_sqrt(r, t, arf_bits(t) / 2 + 2, ARF_RND_DOWN);
        if (ok)
        {
            arf_abs(t, a);
            arf_add(t, t, r, ARF_PREC_EXACT, ARF_RND_DOWN);
            arf_mul_2exp_si(t, t, -1);
            ok = !arf_sqrt(q, t, arf_bits(t) / 2 + 2, ARF_RND_DOWN);
        }
        if (ok)
        {
            arf_mul_2exp_si(t, q, 1);
            ok = !arf_div(t, b, t, arf_bits(b) + 2, ARF_RND_DOWN);
        }
        if (ok)
        {
            if (arf_sgn(a) > 0)
            {
                arf_set(acf_realref(res), q);
                arf_swap(acf_imagref(res), t);
            }
            else
            {
                arf_abs(acf_realref(res), t);
                if (arf_sgn(b) > 0)
                    arf_set(acf_imagref(res), q);
                else
                    arf_neg(acf_imagref(res), q);
            }
        }
    }

    arf_clear(r);
    arf_clear(q);
    arf_clear(t);
    return ok;
}

/* the size (in bits) of the exact powers computed */
#define ACF_POW_EXACT_MAX_BITS (WORD(1) << 24)

/* the number of bits from the top bit of the larger part of the nonzero
   finite z to the bottom bit of either part, or -1 if it exceeds lim */
static slong
_acf_span_bits(const acf_t z, slong lim)
{
    arf_srcptr a = acf_realref(z), b = acf_imagref(z);
    fmpz_t top, bot, t;
    slong r;

    if (arf_is_zero(a) || arf_is_zero(b))
        return FLINT_MAX(arf_bits(a), arf_bits(b));

    fmpz_init(top);
    fmpz_init(bot);
    fmpz_init(t);
    fmpz_set(top, (fmpz_cmp(ARF_EXPREF(a), ARF_EXPREF(b)) >= 0) ? ARF_EXPREF(a) : ARF_EXPREF(b));
    fmpz_sub_ui(bot, ARF_EXPREF(a), arf_bits(a));
    fmpz_sub_ui(t, ARF_EXPREF(b), arf_bits(b));
    if (fmpz_cmp(t, bot) < 0)
        fmpz_swap(t, bot);
    fmpz_sub(t, top, bot);
    r = (fmpz_cmp_si(t, lim) > 0) ? -1 : fmpz_get_si(t);
    fmpz_clear(top);
    fmpz_clear(bot);
    fmpz_clear(t);
    return r;
}


/* x^y for a real y = M 2^E (M odd, |M| < 2^20, -8 <= E <= 20) and a
   nonzero x, exactly: the 2^-E-th root by repeated exact square roots
   (declined when a root is not dyadic), the power |M| by binary
   powering (declined beyond ACF_POW_EXACT_MAX_BITS), the reciprocal
   for M < 0 (exact for a norm that is a power of two).  The result is
   rounded into res; returns the inexact flags, or -1 when declined. */
static int
_acf_pow_dyadic(acf_t res, const acf_t x, const arf_t y, slong prec, arf_rnd_t rnd)
{
    acf_t z, w;
    fmpz_t M, e;
    slong E, bits, i;
    ulong n;
    int status = -1, neg;

    acf_init(z);
    acf_init(w);
    fmpz_init(M);
    fmpz_init(e);

    /* y = M 2^E with M odd */
    arf_get_fmpz_2exp(M, e, y);
    E = fmpz_get_si(e);
    /* (n = |M| 2^max(E, 0) below 2^24, the bound on the bits of the
       result, so that it fits a word on every machine) */
    if (fmpz_bits(M) > 20 || E < -8 || fmpz_bits(M) + FLINT_MAX(E, 0) > 24)
        goto cleanup;
    neg = fmpz_sgn(M) < 0;
    n = (ulong) FLINT_ABS(fmpz_get_si(M));
    if (E > 0)
    {
        n <<= E;
        E = 0;
    }

    acf_set(z, x);
    for (i = 0; i < -E; i++)
        if (!_acf_sqrt_exact(z, z))
            goto cleanup;

    /* the span of z, from the top bit of the larger part to the bottom
       bit of either part, bounds the spans of the parts of z^k by k times
       it (plus k): decline beyond ACF_POW_EXACT_MAX_BITS */
    bits = _acf_span_bits(z, ACF_POW_EXACT_MAX_BITS);
    if (bits < 0 || (ulong) bits + 1 > ACF_POW_EXACT_MAX_BITS / n)
        goto cleanup;

    /* w = z^n by binary powering */
    arf_one(acf_realref(w));
    arf_zero(acf_imagref(w));
    for (i = FLINT_BIT_COUNT(n) - 1; i >= 0; i--)
    {
        acf_mul(w, w, w, ARF_PREC_EXACT, ARF_RND_DOWN);
        if ((n >> i) & 1)
            acf_mul(w, w, z, ARF_PREC_EXACT, ARF_RND_DOWN);
    }

    if (neg)
    {
        /* 1 / w = conj(w) / |w|^2, exact for a power of two |w|^2 */
        arf_t t;
        arf_init(t);
        arf_sosq(t, acf_realref(w), acf_imagref(w), ARF_PREC_EXACT, ARF_RND_DOWN);
        if (!(ARF_IS_POW2(t)))
        {
            arf_clear(t);
            goto cleanup;
        }
        arf_neg(acf_imagref(w), acf_imagref(w));
        arf_div(acf_realref(w), acf_realref(w), t, ARF_PREC_EXACT, ARF_RND_DOWN);
        arf_div(acf_imagref(w), acf_imagref(w), t, ARF_PREC_EXACT, ARF_RND_DOWN);
        arf_clear(t);
    }

    if (rnd >= ARF_RND_FAST)
        rnd = arf_rnd_relaxed_to_strict(rnd);
    status = acf_set_round(res, w, prec, rnd);

cleanup:
    acf_clear(z);
    acf_clear(w);
    fmpz_clear(M);
    fmpz_clear(e);
    return status;
}

/* x^n for an integer n = y and x = 2^j u v with v = 1 + i and u = i^k
   (the dyadic numbers whose integer powers all have a few bits, other
   than the real and imaginary ones): x^n = 2^E i^T v^r with, for n >= 0
   and q = n div 2, E = j n + q, T = k n + q, r = n mod 2 (from v^2 =
   2i), and for n < 0 with m = -n (from v^-1 = (1 - i) / 2 and 1 - i =
   i^3 v), E = j n - m + (m div 2), T = k n + 3m + (m div 2), r = m mod
   2.  Returns the inexact flags, or -1 when x is not of that form. */
static int
_acf_pow_unit(acf_t res, const acf_t x, const arf_t y, slong prec, arf_rnd_t rnd)
{
    arf_srcptr a = acf_realref(x), b = acf_imagref(x);
    fmpz_t n, m, q, e;
    acf_t w;
    int k, r, neg, status;
    slong t;

    if (!(ARF_IS_POW2(a)) || !(ARF_IS_POW2(b)) || !fmpz_equal(ARF_EXPREF(a), ARF_EXPREF(b)))
        return -1;

    /* |x|^y = 2^((j + 1/2) y) with |y| beyond the limit of exp (y log|x|
       out of range, as for arf_pow): nan, before y becomes an integer
       of that size */
    if (arf_cmpabs_2exp_si(y, ARF_EXP_MAGLIM(prec)) >= 0)
        return _acf_nan(res);

    fmpz_init(n);
    fmpz_init(m);
    fmpz_init(q);
    fmpz_init(e);
    acf_init(w);

    /* u = i^k: (+, +) -> 0, (-, +) -> 1, (-, -) -> 2, (+, -) -> 3 */
    k = (ARF_SGNBIT(a) && ARF_SGNBIT(b)) ? 2 : ARF_SGNBIT(a) ? 1 : ARF_SGNBIT(b) ? 3 : 0;
    arf_get_fmpz(n, y, ARF_RND_DOWN);
    neg = fmpz_sgn(n) < 0;
    fmpz_abs(m, n);
    r = fmpz_is_odd(m);
    fmpz_fdiv_q_2exp(q, m, 1);

    fmpz_sub_ui(e, ARF_EXPREF(a), 1);       /* j */
    fmpz_mul(e, e, n);
    fmpz_add(e, e, q);
    t = (k * (slong) fmpz_fdiv_ui(n, 4) + (slong) fmpz_fdiv_ui(q, 4)) % 4;
    if (neg)
    {
        fmpz_sub(e, e, m);
        t = (t + 3 * (slong) fmpz_fdiv_ui(m, 4)) % 4;
    }

    arf_one(acf_realref(w));
    if (r)
        arf_one(acf_imagref(w));
    else
        arf_zero(acf_imagref(w));
    while (t-- > 0)
    {
        /* i (c + di) = -d + ci */
        arf_swap(acf_realref(w), acf_imagref(w));
        arf_neg(acf_realref(w), acf_realref(w));
    }
    arf_mul_2exp_fmpz(acf_realref(w), acf_realref(w), e);
    arf_mul_2exp_fmpz(acf_imagref(w), acf_imagref(w), e);

    if (rnd >= ARF_RND_FAST)
        rnd = arf_rnd_relaxed_to_strict(rnd);
    status = acf_set_round(res, w, prec, rnd);

    fmpz_clear(n);
    fmpz_clear(m);
    fmpz_clear(q);
    fmpz_clear(e);
    acf_clear(w);
    return status;
}

int
acf_pow(acf_t res, const acf_t x, const acf_t y, slong prec, arf_rnd_t rnd)
{
    arf_srcptr a = acf_realref(x), b = acf_imagref(x), c = acf_realref(y), d = acf_imagref(y);
    int s;

    if (ACF_PART_NONFINITE(a) || ACF_PART_NONFINITE(b) || ACF_PART_NONFINITE(c) || ACF_PART_NONFINITE(d))
        return _acf_nan(res);

    /* x^0 = 1 (also 0^0), 0^y = 0 for Re(y) > 0 and nan otherwise */
    if (arf_is_zero(c) && arf_is_zero(d))
    {
        arf_one(acf_realref(res));
        arf_zero(acf_imagref(res));
        return 0;
    }
    if (arf_is_zero(a) && arf_is_zero(b))
    {
        if (arf_sgn(c) > 0)
        {
            arf_zero(acf_realref(res));
            arf_zero(acf_imagref(res));
            return 0;
        }
        return _acf_nan(res);
    }

    /* real powers of positive reals (correctly rounded by arf_pow, and
       exact where the result is dyadic) */
    if (arf_is_zero(b) && arf_is_zero(d) && arf_sgn(a) > 0)
    {
        s = arf_pow(acf_realref(res), a, c, prec, rnd);
        arf_zero(acf_imagref(res));
        return s;
    }

    /* dyadic powers (integer exponents, or small dyadic ones of exact
       roots) exactly: their rounding, which the compositions could not
       certify, from the exact values */
    if (arf_is_zero(d) && !ARF_IS_SPECIAL(c) && ARF_EXP(c) <= 21 && ARF_EXP(c) >= -8 + 1 && arf_bits(c) <= 20)
    {
        s = _acf_pow_dyadic(res, x, c, prec, rnd);
        if (s >= 0)
            return s;
    }
    if (arf_is_zero(d) && arf_is_int(c))
    {
        s = _acf_pow_unit(res, x, c, prec, rnd);
        if (s >= 0)
            return s;
    }

    /* near the exact cases (x a power of two times a unit, times 1 + h
       with h tiny) and for a tiny y: the expansions of acf/near.c, which
       decide the values close to a rounding boundary exactly */
    if (_acf_pow_near(res, &s, x, y, prec, rnd))
        return s;

    /* exp(y log x) with the principal logarithm (the composition takes
       zero parts: the logarithm of a real or imaginary x has an exact
       imaginary part, a multiple of pi/2) */
    return _acf_elem(res, NULL, CELEM_POW, x, y, prec, rnd, rnd);
}
