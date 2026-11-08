/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Internal definitions for the complex elementary functions of acf,
   shared with nfloat_complex and acb (which include this header). */

#ifndef ACF_ELEM_H
#define ACF_ELEM_H

#include "acf.h"
#include "../arf/elem.h"

/*
    The complex functions are compositions of the real ones:

        exp(a + bi) = exp(a) (cos b + i sin b)
        exp(pi i (a + bi)) = exp(-pi b) (cos(pi a) + i sin(pi a))
        sin(a + bi) = sin a cosh b + i cos a sinh b
        cos(a + bi) = cos a cosh b - i sin a sinh b
        (sin, cos)(pi (a + bi)): the same with pi a and pi b
        log(a + bi) = log(a^2 + b^2) / 2 + i atan2(b, a)
        x^y = exp(y log x)

    with each part a product of two accurately known real quantities,
    or a real function of an exactly known one (a^2 + b^2 is formed
    exactly, so that log keeps its relative accuracy near the unit
    circle).  Each part therefore comes out with its own relative
    accuracy, whatever its size against the other part.

    At one and two limbs (64-bit) the functions below compose the
    register kernels of the real functions (arf/elem.h) on fixed-point
    frames: a frame (y2, y1, y0) resp. (y3, y2, y1, y0) with an
    integral top limb stands for Y 2^(k - 2 FLINT_BITS) resp.
    Y 2^(k - 3 FLINT_BITS) within err units of its lowest limb, which
    _arf_small_set_1 / _2 (arf), _nfloat_small_set_1 / _2 (nfloat) and
    _arb_small_set_1 / _2 (arb) round.  Each part has its own exponent
    k, since a tiny part is represented by its own mantissa (sin b = b
    (1 - delta) for tiny b, say) rather than as a small number in the
    frame of the other.  The functions return ARF_ELEM_UNCERTAIN (with
    the outputs unspecified) for arguments outside their ranges, which
    the callers evaluate at full precision.

    The real factors (the trigonometric, hyperbolic and exponential
    factors, log |z| and atan2), each a real kernel with its argument
    reduction, are compiled once in acf/elem_small.c; the compositions
    (a few frame products) are inline here.
*/

/* one part of a result: (-1)^sgnbit (y[3], y[2], y[1], y[0]) 2^(k - 3
   FLINT_BITS) at two limbs, or (y[2], y[1], y[0]) 2^(k - 2 FLINT_BITS)
   at one limb (y[3] unused), within err units */
typedef struct
{
    ulong y[4];
    slong k;
    ulong err;
    int sgnbit;
}
_acf_frame_struct;

#if FLINT_BITS == 64

/* Fixed-point products ****************************************************/

/* R = X Y for the frames X = (x2, x1, x0), Y = (y2, y1, y0) (an integral
   limb and two fraction limbs each; the product's integral part must
   fit a limb), truncated: within 3 units below the exact product */
FLINT_FORCE_INLINE void
_acf_small_mul_3(ulong * r2, ulong * r1, ulong * r0,
    ulong x2, ulong x1, ulong x0, ulong y2, ulong y1, ulong y0)
{
    ulong h, l, s2, s1, s0;

    umul_ppmm(s1, s0, x1, y1);
    s2 = x2 * y2;
    umul_ppmm(h, l, x2, y1);
    add_ssaaaa(s2, s1, s2, s1, h, l);
    umul_ppmm(h, l, x1, y2);
    add_ssaaaa(s2, s1, s2, s1, h, l);
    umul_ppmm(h, l, x2, y0);
    add_sssaaaaaa(s2, s1, s0, s2, s1, s0, UWORD(0), h, l);
    umul_ppmm(h, l, x0, y2);
    add_sssaaaaaa(s2, s1, s0, s2, s1, s0, UWORD(0), h, l);
    umul_ppmm(h, l, x1, y0);
    add_sssaaaaaa(s2, s1, s0, s2, s1, s0, UWORD(0), UWORD(0), h);
    umul_ppmm(h, l, x0, y1);
    add_sssaaaaaa(s2, s1, s0, s2, s1, s0, UWORD(0), UWORD(0), h);

    *r2 = s2; *r1 = s1; *r0 = s0;
}

/* the same with three fraction limbs: within 6 units below */
FLINT_FORCE_INLINE void
_acf_small_mul_4(ulong * r3, ulong * r2, ulong * r1, ulong * r0,
    ulong x3, ulong x2, ulong x1, ulong x0, ulong y3, ulong y2, ulong y1, ulong y0)
{
    ulong h, l, s3, s2, s1, s0;

    /* weight B^-3: the low words at r0 */
    umul_ppmm(s1, s0, x2, y1);
    s2 = 0;
    umul_ppmm(h, l, x1, y2);
    add_sssaaaaaa(s2, s1, s0, s2, s1, s0, UWORD(0), h, l);
    umul_ppmm(h, l, x3, y0);
    add_sssaaaaaa(s2, s1, s0, s2, s1, s0, UWORD(0), h, l);
    umul_ppmm(h, l, x0, y3);
    add_sssaaaaaa(s2, s1, s0, s2, s1, s0, UWORD(0), h, l);
    /* weight B^-4: the high words at r0 */
    umul_ppmm(h, l, x2, y0);
    add_sssaaaaaa(s2, s1, s0, s2, s1, s0, UWORD(0), UWORD(0), h);
    umul_ppmm(h, l, x0, y2);
    add_sssaaaaaa(s2, s1, s0, s2, s1, s0, UWORD(0), UWORD(0), h);
    umul_ppmm(h, l, x1, y1);
    add_sssaaaaaa(s2, s1, s0, s2, s1, s0, UWORD(0), UWORD(0), h);
    /* weight B^-2: the high words at r2 */
    s3 = x3 * y3;
    umul_ppmm(h, l, x2, y2);
    add_sssaaaaaa(s3, s2, s1, s3, s2, s1, UWORD(0), h, l);
    umul_ppmm(h, l, x3, y1);
    add_sssaaaaaa(s3, s2, s1, s3, s2, s1, UWORD(0), h, l);
    umul_ppmm(h, l, x1, y3);
    add_sssaaaaaa(s3, s2, s1, s3, s2, s1, UWORD(0), h, l);
    /* weight B^-1 */
    umul_ppmm(h, l, x3, y2);
    add_ssaaaa(s3, s2, s3, s2, h, l);
    umul_ppmm(h, l, x2, y3);
    add_ssaaaa(s3, s2, s3, s2, h, l);

    *r3 = s3; *r2 = s2; *r1 = s1; *r0 = s0;
}

/* the error of a product of frames within ex resp. ey units with values
   below XB resp. YB, including the truncation (the cross term of the
   errors is below a unit) */
#define _ACF_SMALL_MUL_ERR_3(ex, XB, ey, YB) ((ex) * (YB) + (ey) * (XB) + 4)
#define _ACF_SMALL_MUL_ERR_4(ex, XB, ey, YB) ((ex) * (YB) + (ey) * (XB) + 7)

/* pi |x| for |x| = d 2^(e - FLINT_BITS) resp. (d1, d0) 2^(e - 2 FLINT_BITS),
   -FLINT_BITS <= e <= FLINT_BITS - 4 (resp. -2 FLINT_BITS <= e), as the
   fixed-point (yi, y2, y1, y0) with three fraction limbs resp. (yi, y3,
   y2, y1, y0) with four, yi < 2^62: 3|x| exactly plus (pi - 3)|x| from
   the floor of pi - 3 at the fraction limbs (its deficit, times |x| <
   2^61, below 2^-131 resp. 2^-195) with the low partial products
   dropped (below 4 resp. 6 units): within 2^61 + 4 units in all */
FLINT_FORCE_INLINE void
_acf_small_pi_times_1(ulong * yi, ulong * y2, ulong * y1, ulong * y0, ulong d, slong e)
{
    nn_srcptr P = _mp_real_const_ptr(MP_REAL_CONST_ID_PI4, 4);
    ulong xi, x1, x0, c2, c1, c0, h, l, s3, s2, s1, s0;

    _arf_small_fixed_1(&xi, &x1, &x0, d, e);

    /* pi - 3 at three fraction limbs: the fraction of 4P */
    c2 = (P[3] << 2) | (P[2] >> (FLINT_BITS - 2));
    c1 = (P[2] << 2) | (P[1] >> (FLINT_BITS - 2));
    c0 = (P[1] << 2) | (P[0] >> (FLINT_BITS - 2));

    /* 3 X */
    s3 = 3 * xi; s2 = 0; s1 = 0; s0 = 0;
    add_ssssaaaaaaaa(s3, s2, s1, s0, s3, s2, s1, s0, UWORD(0), x1, x0, UWORD(0));
    add_ssssaaaaaaaa(s3, s2, s1, s0, s3, s2, s1, s0, UWORD(0), x1, x0, UWORD(0));
    add_ssssaaaaaaaa(s3, s2, s1, s0, s3, s2, s1, s0, UWORD(0), x1, x0, UWORD(0));

    /* (pi - 3) X: the partial products of weight B^0 .. B^-3 */
    umul_ppmm(h, l, xi, c2);
    add_ssaaaa(s3, s2, s3, s2, h, l);
    umul_ppmm(h, l, xi, c1);
    add_sssaaaaaa(s3, s2, s1, s3, s2, s1, UWORD(0), h, l);
    umul_ppmm(h, l, x1, c2);
    add_sssaaaaaa(s3, s2, s1, s3, s2, s1, UWORD(0), h, l);
    umul_ppmm(h, l, xi, c0);
    add_ssssaaaaaaaa(s3, s2, s1, s0, s3, s2, s1, s0, UWORD(0), UWORD(0), h, l);
    umul_ppmm(h, l, x1, c1);
    add_ssssaaaaaaaa(s3, s2, s1, s0, s3, s2, s1, s0, UWORD(0), UWORD(0), h, l);
    umul_ppmm(h, l, x0, c2);
    add_ssssaaaaaaaa(s3, s2, s1, s0, s3, s2, s1, s0, UWORD(0), UWORD(0), h, l);
    umul_ppmm(h, l, x1, c0);
    add_ssssaaaaaaaa(s3, s2, s1, s0, s3, s2, s1, s0, UWORD(0), UWORD(0), UWORD(0), h);
    umul_ppmm(h, l, x0, c1);
    add_ssssaaaaaaaa(s3, s2, s1, s0, s3, s2, s1, s0, UWORD(0), UWORD(0), UWORD(0), h);

    *yi = s3; *y2 = s2; *y1 = s1; *y0 = s0;
}

FLINT_FORCE_INLINE void
_acf_small_pi_times_2(ulong * yi, ulong * y3, ulong * y2, ulong * y1, ulong * y0, ulong d1, ulong d0, slong e)
{
    nn_srcptr P = _mp_real_const_ptr(MP_REAL_CONST_ID_PI4, 5);
    ulong xi, x2, x1, x0, c3, c2, c1, c0, h, l, s4, s3, s2, s1, s0;

    _arf_small_fixed_2(&xi, &x2, &x1, &x0, d1, d0, e);

    c3 = (P[4] << 2) | (P[3] >> (FLINT_BITS - 2));
    c2 = (P[3] << 2) | (P[2] >> (FLINT_BITS - 2));
    c1 = (P[2] << 2) | (P[1] >> (FLINT_BITS - 2));
    c0 = (P[1] << 2) | (P[0] >> (FLINT_BITS - 2));

    s4 = 3 * xi; s3 = 0; s2 = 0; s1 = 0; s0 = 0;
    add_sssssaaaaaaaaaa(s4, s3, s2, s1, s0, s4, s3, s2, s1, s0, UWORD(0), x2, x1, x0, UWORD(0));
    add_sssssaaaaaaaaaa(s4, s3, s2, s1, s0, s4, s3, s2, s1, s0, UWORD(0), x2, x1, x0, UWORD(0));
    add_sssssaaaaaaaaaa(s4, s3, s2, s1, s0, s4, s3, s2, s1, s0, UWORD(0), x2, x1, x0, UWORD(0));

    umul_ppmm(h, l, xi, c3);
    add_ssaaaa(s4, s3, s4, s3, h, l);
    umul_ppmm(h, l, xi, c2);
    add_sssaaaaaa(s4, s3, s2, s4, s3, s2, UWORD(0), h, l);
    umul_ppmm(h, l, x2, c3);
    add_sssaaaaaa(s4, s3, s2, s4, s3, s2, UWORD(0), h, l);
    umul_ppmm(h, l, xi, c1);
    add_ssssaaaaaaaa(s4, s3, s2, s1, s4, s3, s2, s1, UWORD(0), UWORD(0), h, l);
    umul_ppmm(h, l, x2, c2);
    add_ssssaaaaaaaa(s4, s3, s2, s1, s4, s3, s2, s1, UWORD(0), UWORD(0), h, l);
    umul_ppmm(h, l, x1, c3);
    add_ssssaaaaaaaa(s4, s3, s2, s1, s4, s3, s2, s1, UWORD(0), UWORD(0), h, l);
    umul_ppmm(h, l, xi, c0);
    add_sssssaaaaaaaaaa(s4, s3, s2, s1, s0, s4, s3, s2, s1, s0, UWORD(0), UWORD(0), UWORD(0), h, l);
    umul_ppmm(h, l, x2, c1);
    add_sssssaaaaaaaaaa(s4, s3, s2, s1, s0, s4, s3, s2, s1, s0, UWORD(0), UWORD(0), UWORD(0), h, l);
    umul_ppmm(h, l, x1, c2);
    add_sssssaaaaaaaaaa(s4, s3, s2, s1, s0, s4, s3, s2, s1, s0, UWORD(0), UWORD(0), UWORD(0), h, l);
    umul_ppmm(h, l, x0, c3);
    add_sssssaaaaaaaaaa(s4, s3, s2, s1, s0, s4, s3, s2, s1, s0, UWORD(0), UWORD(0), UWORD(0), h, l);
    umul_ppmm(h, l, x2, c0);
    add_sssssaaaaaaaaaa(s4, s3, s2, s1, s0, s4, s3, s2, s1, s0, UWORD(0), UWORD(0), UWORD(0), UWORD(0), h);
    umul_ppmm(h, l, x1, c1);
    add_sssssaaaaaaaaaa(s4, s3, s2, s1, s0, s4, s3, s2, s1, s0, UWORD(0), UWORD(0), UWORD(0), UWORD(0), h);
    umul_ppmm(h, l, x0, c2);
    add_sssssaaaaaaaaaa(s4, s3, s2, s1, s0, s4, s3, s2, s1, s0, UWORD(0), UWORD(0), UWORD(0), UWORD(0), h);

    *yi = s4; *y3 = s3; *y2 = s2; *y1 = s1; *y0 = s0;
}

/* The real factors ********************************************************/

/* tiny arguments: below these exponents the quadratic terms of the
   expansions of the factors (x^2 / 2 and (pi x)^2 / 2 < 2^(2e + 3)) lie
   within two units of the frames of two resp. three fraction limbs */
#define _ACF_SMALL_TINY_EXP_1 (-FLINT_BITS - 1)
#define _ACF_SMALL_TINY_EXP_2 (-3 * FLINT_BITS / 2 - 1)

/* |sin x| = S 2^(ks - 2 FLINT_BITS) and |cos x| = C 2^-(2 FLINT_BITS) for
   the frames S = (s2, s1, s0) and C = (c2, c1, c0) (values below 1)
   within *err units, with the signs *negs of sin x and *negc of cos x,
   for x = (-1)^xsgn d 2^(e - FLINT_BITS), or of pi x (pi = 1).  Tiny
   arguments (e <= _ACF_SMALL_TINY_EXP_1) give sin x = x and cos x = 1
   within two units (the frame of x, or of pi |x| with its lowest limb
   dropped, has ks = e).  Returns 0, or ARF_ELEM_UNCERTAIN when the
   reduction declines. */
int _acf_small_trig_1(ulong * s2, ulong * s1, ulong * s0, slong * ks, int * negs,
    ulong * c2, ulong * c1, ulong * c0, int * negc, ulong * err, ulong d, slong e, int xsgn, int pi);

/* the same at two limbs for x = (-1)^xsgn (d1, d0) 2^(e - 2 FLINT_BITS):
   frames with three fraction limbs */
int _acf_small_trig_2(ulong * s3, ulong * s2, ulong * s1, ulong * s0, slong * ks, int * negs,
    ulong * c3, ulong * c2, ulong * c1, ulong * c0, int * negc, ulong * err,
    ulong d1, ulong d0, slong e, int xsgn, int pi);

/* bounds on the fixed-point values of the frames of the factors (the
   trigonometric and hyperbolic ones reach pi |x| for tiny arguments) */
#define _ACF_SMALL_TRIG_BOUND 4
#define _ACF_SMALL_HYP_BOUND 4
#define _ACF_SMALL_EXP_BOUND 2

/* |sinh x| = S 2^(ks - 2 FLINT_BITS) and cosh x = C 2^(kc - 2 FLINT_BITS)
   (values below _ACF_SMALL_HYP_BOUND) within *err units for |x| = d
   2^(e - FLINT_BITS), or of pi |x| (pi = 1); tiny arguments as for the
   trigonometric functions; e <= MP_REAL_EXP_FP_MAX_EXP - 2 */
int _acf_small_hyp_1(ulong * s2, ulong * s1, ulong * s0, slong * ks,
    ulong * c2, ulong * c1, ulong * c0, slong * kc, ulong * err, ulong d, slong e, int pi);

int _acf_small_hyp_2(ulong * s3, ulong * s2, ulong * s1, ulong * s0, slong * ks,
    ulong * c3, ulong * c2, ulong * c1, ulong * c0, slong * kc, ulong * err, ulong d1, ulong d0, slong e, int pi);

/* E = exp(x) (frame value below 2) for x = (-1)^sgnbit d 2^(e - FLINT_BITS),
   or exp(-pi x) (pi = 1), as (e2, e1, e0) 2^(k - 2 FLINT_BITS) within
   *err units; e <= MP_REAL_EXP_FP_MAX_EXP - 2.  Tiny arguments by the
   linear term (the quadratic one below a unit). */
int _acf_small_exp_factor_1(ulong * e2, ulong * e1, ulong * e0, slong * k, ulong * err, ulong d, slong e, int sgnbit, int pi);

int _acf_small_exp_factor_2(ulong * e3, ulong * e2, ulong * e1, ulong * e0, slong * k, ulong * err,
    ulong d1, ulong d0, slong e, int sgnbit, int pi);

/* The complex functions ***************************************************/

/* exp(a + bi), or exp(pi i (a + bi)) (pi = 1), for a = (-1)^asgn da
   2^(ea - FLINT_BITS) and b likewise: the frames of the real and
   imaginary parts */
FLINT_FORCE_INLINE int
_acf_small_exp_1(_acf_frame_struct * re, _acf_frame_struct * im,
    ulong da, slong ea, int asgn, ulong db, slong eb, int bsgn, int pi)
{
    ulong e2, e1, e0, s2, s1, s0, c2, c1, c0, errE, errT;
    slong kE, ks;
    int negs, negc;

    /* exp(pi i z) = exp(-pi b) (cos(pi a) + i sin(pi a)) */
    if (pi)
    {
        if (_acf_small_exp_factor_1(&e2, &e1, &e0, &kE, &errE, db, eb, bsgn, 1) != 0 ||
            _acf_small_trig_1(&s2, &s1, &s0, &ks, &negs, &c2, &c1, &c0, &negc, &errT, da, ea, asgn, 1) != 0)
            return ARF_ELEM_UNCERTAIN;
    }
    else
    {
        if (_acf_small_exp_factor_1(&e2, &e1, &e0, &kE, &errE, da, ea, asgn, 0) != 0 ||
            _acf_small_trig_1(&s2, &s1, &s0, &ks, &negs, &c2, &c1, &c0, &negc, &errT, db, eb, bsgn, 0) != 0)
            return ARF_ELEM_UNCERTAIN;
    }

    _acf_small_mul_3(re->y + 2, re->y + 1, re->y + 0, e2, e1, e0, c2, c1, c0);
    re->k = kE;
    re->sgnbit = negc;
    _acf_small_mul_3(im->y + 2, im->y + 1, im->y + 0, e2, e1, e0, s2, s1, s0);
    im->k = kE + ks;
    im->sgnbit = negs;
    re->err = im->err = _ACF_SMALL_MUL_ERR_3(errE, _ACF_SMALL_TRIG_BOUND, errT, _ACF_SMALL_EXP_BOUND);
    return 0;
}

FLINT_FORCE_INLINE int
_acf_small_exp_2(_acf_frame_struct * re, _acf_frame_struct * im,
    ulong da1, ulong da0, slong ea, int asgn, ulong db1, ulong db0, slong eb, int bsgn, int pi)
{
    ulong e3, e2, e1, e0, s3, s2, s1, s0, c3, c2, c1, c0, errE, errT;
    slong kE, ks;
    int negs, negc;

    if (pi)
    {
        if (_acf_small_exp_factor_2(&e3, &e2, &e1, &e0, &kE, &errE, db1, db0, eb, bsgn, 1) != 0 ||
            _acf_small_trig_2(&s3, &s2, &s1, &s0, &ks, &negs, &c3, &c2, &c1, &c0, &negc, &errT, da1, da0, ea, asgn, 1) != 0)
            return ARF_ELEM_UNCERTAIN;
    }
    else
    {
        if (_acf_small_exp_factor_2(&e3, &e2, &e1, &e0, &kE, &errE, da1, da0, ea, asgn, 0) != 0 ||
            _acf_small_trig_2(&s3, &s2, &s1, &s0, &ks, &negs, &c3, &c2, &c1, &c0, &negc, &errT, db1, db0, eb, bsgn, 0) != 0)
            return ARF_ELEM_UNCERTAIN;
    }

    _acf_small_mul_4(re->y + 3, re->y + 2, re->y + 1, re->y + 0, e3, e2, e1, e0, c3, c2, c1, c0);
    re->k = kE;
    re->sgnbit = negc;
    _acf_small_mul_4(im->y + 3, im->y + 2, im->y + 1, im->y + 0, e3, e2, e1, e0, s3, s2, s1, s0);
    im->k = kE + ks;
    im->sgnbit = negs;
    re->err = im->err = _ACF_SMALL_MUL_ERR_4(errE, _ACF_SMALL_TRIG_BOUND, errT, _ACF_SMALL_EXP_BOUND);
    return 0;
}

/* sin(a + bi) and cos(a + bi) (either pair of outputs NULL), or of
   pi (a + bi) (pi = 1) */
FLINT_FORCE_INLINE int
_acf_small_sin_cos_1(_acf_frame_struct * sre, _acf_frame_struct * sim, _acf_frame_struct * cre, _acf_frame_struct * cim,
    ulong da, slong ea, int asgn, ulong db, slong eb, int bsgn, int pi)
{
    ulong s2, s1, s0, c2, c1, c0, h2, h1, h0, g2, g1, g0, errT, errH, err;
    slong ks, kh, kg;
    int negs, negc;

    if (_acf_small_trig_1(&s2, &s1, &s0, &ks, &negs, &c2, &c1, &c0, &negc, &errT, da, ea, asgn, pi) != 0 ||
        _acf_small_hyp_1(&h2, &h1, &h0, &kh, &g2, &g1, &g0, &kg, &errH, db, eb, pi) != 0)
        return ARF_ELEM_UNCERTAIN;

    err = _ACF_SMALL_MUL_ERR_3(errT, _ACF_SMALL_HYP_BOUND, errH, _ACF_SMALL_TRIG_BOUND);

    /* sin(a + bi) = sin a cosh b + i cos a sinh b */
    if (sre != NULL)
    {
        _acf_small_mul_3(sre->y + 2, sre->y + 1, sre->y + 0, s2, s1, s0, g2, g1, g0);
        sre->k = ks + kg;
        sre->sgnbit = negs;
        sre->err = err;
        _acf_small_mul_3(sim->y + 2, sim->y + 1, sim->y + 0, c2, c1, c0, h2, h1, h0);
        sim->k = kh;
        sim->sgnbit = negc ^ bsgn;
        sim->err = err;
    }

    /* cos(a + bi) = cos a cosh b - i sin a sinh b */
    if (cre != NULL)
    {
        _acf_small_mul_3(cre->y + 2, cre->y + 1, cre->y + 0, c2, c1, c0, g2, g1, g0);
        cre->k = kg;
        cre->sgnbit = negc;
        cre->err = err;
        _acf_small_mul_3(cim->y + 2, cim->y + 1, cim->y + 0, s2, s1, s0, h2, h1, h0);
        cim->k = ks + kh;
        cim->sgnbit = negs ^ bsgn ^ 1;
        cim->err = err;
    }

    return 0;
}

FLINT_FORCE_INLINE int
_acf_small_sin_cos_2(_acf_frame_struct * sre, _acf_frame_struct * sim, _acf_frame_struct * cre, _acf_frame_struct * cim,
    ulong da1, ulong da0, slong ea, int asgn, ulong db1, ulong db0, slong eb, int bsgn, int pi)
{
    ulong s3, s2, s1, s0, c3, c2, c1, c0, h3, h2, h1, h0, g3, g2, g1, g0, errT, errH, err;
    slong ks, kh, kg;
    int negs, negc;

    if (_acf_small_trig_2(&s3, &s2, &s1, &s0, &ks, &negs, &c3, &c2, &c1, &c0, &negc, &errT, da1, da0, ea, asgn, pi) != 0 ||
        _acf_small_hyp_2(&h3, &h2, &h1, &h0, &kh, &g3, &g2, &g1, &g0, &kg, &errH, db1, db0, eb, pi) != 0)
        return ARF_ELEM_UNCERTAIN;

    err = _ACF_SMALL_MUL_ERR_4(errT, _ACF_SMALL_HYP_BOUND, errH, _ACF_SMALL_TRIG_BOUND);

    if (sre != NULL)
    {
        _acf_small_mul_4(sre->y + 3, sre->y + 2, sre->y + 1, sre->y + 0, s3, s2, s1, s0, g3, g2, g1, g0);
        sre->k = ks + kg;
        sre->sgnbit = negs;
        sre->err = err;
        _acf_small_mul_4(sim->y + 3, sim->y + 2, sim->y + 1, sim->y + 0, c3, c2, c1, c0, h3, h2, h1, h0);
        sim->k = kh;
        sim->sgnbit = negc ^ bsgn;
        sim->err = err;
    }

    if (cre != NULL)
    {
        _acf_small_mul_4(cre->y + 3, cre->y + 2, cre->y + 1, cre->y + 0, c3, c2, c1, c0, g3, g2, g1, g0);
        cre->k = kg;
        cre->sgnbit = negc;
        cre->err = err;
        _acf_small_mul_4(cim->y + 3, cim->y + 2, cim->y + 1, cim->y + 0, s3, s2, s1, s0, h3, h2, h1, h0);
        cim->k = ks + kh;
        cim->sgnbit = negs ^ bsgn ^ 1;
        cim->err = err;
    }

    return 0;
}

/* The logarithm *************************************************************/

/* |log(1 + D)| = R 2^(ed - 2 FLINT_BITS) for the frame R = (r2, r1, r0)
   and D = (-1)^dneg (d1, d0) 2^(ed - 2 FLINT_BITS) (normalised, exact)
   with ed <= -31, by the series to the fourth power (the tail below
   D^5 / 5); *err the error */
FLINT_FORCE_INLINE void
_acf_small_log1p_series_2(ulong * r2, ulong * r1, ulong * r0, ulong * err, ulong d1, ulong d0, slong ed, int dneg)
{
    ulong w2, w1, w0, t1, t0, s2;
    slong sh = -ed;

    /* Horner: W = 1 -+ u (1/2 -+ u (1/3 -+ u / 4)) for u = |D| */
    _arf_small_rshift_3(&s2, &t1, &t0, UWORD(0), d1, d0, sh + 2);
    if (dneg)
        add_ssaaaa(w1, w0, MP_REAL_SMALL_LOG_C3, MP_REAL_SMALL_LOG_C3, t1, t0);
    else
        sub_ddmmss(w1, w0, MP_REAL_SMALL_LOG_C3, MP_REAL_SMALL_LOG_C3, t1, t0);
    _mp_real_mulhi_2x2_sloppy(&t1, &t0, d1, d0, w1, w0);
    _arf_small_rshift_3(&s2, &t1, &t0, UWORD(0), t1, t0, sh);
    if (dneg)
        add_ssaaaa(w1, w0, MP_REAL_SMALL_LOG_C2, UWORD(0), t1, t0);
    else
        sub_ddmmss(w1, w0, MP_REAL_SMALL_LOG_C2, UWORD(0), t1, t0);
    _mp_real_mulhi_2x2_sloppy(&t1, &t0, d1, d0, w1, w0);
    _arf_small_rshift_3(&s2, &t1, &t0, UWORD(0), t1, t0, sh);
    if (dneg)
        add_sssaaaaaa(w2, w1, w0, UWORD(1), UWORD(0), UWORD(0), UWORD(0), t1, t0);
    else
        sub_dddmmmsss(w2, w1, w0, UWORD(1), UWORD(0), UWORD(0), UWORD(0), t1, t0);

    /* R = u W (which may reach 1 in the frame) */
    _acf_small_mul_3(r2, r1, r0, UWORD(0), d1, d0, w2, w1, w0);

    /* the products (below 3 units each, scaled down by u at the inner
       levels), the shifts, the final product and the tail (u^4 / 5
       relative, below 2^(4 ed + 128) units) */
    *err = 16 + ((4 * ed + 126 >= 0) ? (UWORD(1) << (4 * ed + 126)) : 1);
}

/* the same at three limbs, D = (d2, d1, d0) 2^(ed - 3 FLINT_BITS) with
   ed <= -24, by the series to the sixth power (the tail below D^7 / 7) */
FLINT_FORCE_INLINE void
_acf_small_log1p_series_3(ulong * r3, ulong * r2, ulong * r1, ulong * r0, ulong * err, ulong d2, ulong d1, ulong d0, slong ed, int dneg)
{
    ulong w3, w2, w1, w0, t2, t1, t0;
    slong sh = -ed;
    static const ulong c[5][3] = {
        { MP_REAL_SMALL_LOG_C5, MP_REAL_SMALL_LOG_C5, MP_REAL_SMALL_LOG_C5 },         /* 1/5 */
        { MP_REAL_SMALL_LOG_C4, 0, 0 },                                                 /* 1/4 */
        { MP_REAL_SMALL_LOG_C3, MP_REAL_SMALL_LOG_C3, MP_REAL_SMALL_LOG_C3 },         /* 1/3 */
        { MP_REAL_SMALL_LOG_C2, 0, 0 },                                                 /* 1/2 */
        { 0, 0, 0 } };
    int i;

    /* W = 1/5 -+ u/6, then 1/j -+ u W for j = 4, 3, 2, and 1 -+ u W */
    _arf_small_rshift_4(&w3, &t2, &t1, &t0, UWORD(0), d2, d1, d0, sh);
    _mp_real_small_mulhi_3x1(&t2, &t1, &t0, t2, t1, t0, MP_REAL_SMALL_LOG_C6);
    if (dneg)
        add_sssaaaaaa(w2, w1, w0, c[0][0], c[0][1], c[0][2], t2, t1, t0);
    else
        sub_dddmmmsss(w2, w1, w0, c[0][0], c[0][1], c[0][2], t2, t1, t0);

    for (i = 1; i < 4; i++)
    {
        _mp_real_small_mulhi_3x3(&t2, &t1, &t0, d2, d1, d0, w2, w1, w0);
        _arf_small_rshift_4(&w3, &t2, &t1, &t0, UWORD(0), t2, t1, t0, sh);
        if (dneg)
            add_sssaaaaaa(w2, w1, w0, c[i][0], c[i][1], c[i][2], t2, t1, t0);
        else
            sub_dddmmmsss(w2, w1, w0, c[i][0], c[i][1], c[i][2], t2, t1, t0);
    }

    _mp_real_small_mulhi_3x3(&t2, &t1, &t0, d2, d1, d0, w2, w1, w0);
    _arf_small_rshift_4(&w3, &t2, &t1, &t0, UWORD(0), t2, t1, t0, sh);
    if (dneg)
        add_ssssaaaaaaaa(w3, w2, w1, w0, UWORD(1), UWORD(0), UWORD(0), UWORD(0), UWORD(0), t2, t1, t0);
    else
        sub_ddddmmmmssss(w3, w2, w1, w0, UWORD(1), UWORD(0), UWORD(0), UWORD(0), UWORD(0), t2, t1, t0);

    _acf_small_mul_4(r3, r2, r1, r0, UWORD(0), d2, d1, d0, w3, w2, w1, w0);

    /* (the tail u^6 / 7 relative, below 2^(6 ed + 190) units) */
    *err = 24 + ((6 * ed + 190 >= 0) ? (UWORD(1) << (6 * ed + 190)) : 1);
}

/* log(a^2 + b^2) / 2 for a = da 2^(ea - FLINT_BITS), b = db 2^(eb -
   FLINT_BITS) (nonzero mantissas) as the frame re: the sum of squares
   exactly (the smaller square shifted into place, its dropped bits
   below a unit), then the logarithm kernel of arf at two limbs, or
   near 1 the series of log(1 + D) on the exact D = a^2 + b^2 - 1 */
int _acf_small_log_abs_1(_acf_frame_struct * re, ulong da, slong ea, ulong db, slong eb);

int _acf_small_log_abs_2(_acf_frame_struct * re, ulong da1, ulong da0, slong ea, ulong db1, ulong db0, slong eb);

/* atan2(b, a) as the frame im (k = 0, or k = -d for a tiny angle) for
   a = (-1)^asgn da 2^(ea - FLINT_BITS), b = (-1)^bsgn db 2^(eb -
   FLINT_BITS): c pi/4 + sigma atan(t) with t = (num / den) 2^-d < 1 the
   smaller of |b| / |a| and |a| / |b| (arf/elem.h); a quotient far below
   1 (d large, where the kernel of arf declines) through one division
   instead, atan t = t - t^3 / 3 + t^5 / 5 - ... with the terms that
   matter */
int _acf_small_atan2_1(_acf_frame_struct * im, ulong da, slong ea, int asgn, ulong db, slong eb, int bsgn);

int _acf_small_atan2_2(_acf_frame_struct * im, ulong da1, ulong da0, slong ea, int asgn, ulong db1, ulong db0, slong eb, int bsgn);

/* log(a + bi) = log(a^2 + b^2) / 2 + i atan2(b, a) */
FLINT_FORCE_INLINE int
_acf_small_log_1(_acf_frame_struct * re, _acf_frame_struct * im,
    ulong da, slong ea, int asgn, ulong db, slong eb, int bsgn)
{
    if (_acf_small_log_abs_1(re, da, ea, db, eb) != 0)
        return ARF_ELEM_UNCERTAIN;
    return _acf_small_atan2_1(im, da, ea, asgn, db, eb, bsgn);
}

FLINT_FORCE_INLINE int
_acf_small_log_2(_acf_frame_struct * re, _acf_frame_struct * im,
    ulong da1, ulong da0, slong ea, int asgn, ulong db1, ulong db0, slong eb, int bsgn)
{
    if (_acf_small_log_abs_2(re, da1, da0, ea, db1, db0, eb) != 0)
        return ARF_ELEM_UNCERTAIN;
    return _acf_small_atan2_2(im, da1, da0, ea, asgn, db1, db0, eb, bsgn);
}

/* Powers *******************************************************************/

/* the 5-limb Z = P 2^(k + e) for the 4-limb product P = (r2, r1, r0) m
   of a frame (unit 2^(k - 2 FLINT_BITS)) and a one-limb mantissa m
   (value m 2^(e - FLINT_BITS)), the term on the grid of units
   2^-(3 FLINT_BITS); k + e <= FLINT_BITS - 1 */
FLINT_FORCE_INLINE void
_acf_small_pow_term_1(ulong * Z, ulong r2, ulong r1, ulong r0, slong k, ulong m, slong e)
{
    ulong P[5], h, l;
    slong sh = k + e;

    umul_ppmm(P[1], P[0], r0, m);
    umul_ppmm(h, l, r1, m);
    add_ssaaaa(P[2], P[1], h, l, UWORD(0), P[1]);
    umul_ppmm(h, l, r2, m);
    add_ssaaaa(P[3], P[2], h, l, UWORD(0), P[2]);
    P[4] = 0;

    if (sh >= 0)
    {
        if (sh == 0)
            flint_mpn_copyi(Z, P, 5);
        else
            mpn_lshift(Z, P, 5, sh);
    }
    else if (-sh >= 5 * FLINT_BITS)
        flint_mpn_zero(Z, 5);
    else
        _ARF_POW_RSHIFT(Z, P, 5, -sh);
}

/* the signed sum S = (-1)^s1 A + (-1)^s2 B of 5-limb magnitudes, as a
   magnitude with its sign */
FLINT_FORCE_INLINE void
_acf_small_pow_sum_1(ulong * S, int * ssgn, const ulong * A, int s1, const ulong * B, int s2)
{
    if (s1 == s2)
    {
        mpn_add_n(S, A, B, 5);
        *ssgn = s1;
    }
    else if (mpn_cmp(A, B, 5) >= 0)
    {
        mpn_sub_n(S, A, B, 5);
        *ssgn = s1;
    }
    else
    {
        mpn_sub_n(S, B, A, 5);
        *ssgn = s2;
    }
}

/* x^y = exp(y log x) at one limb for x = a + bi (nonzero parts) and
   y = c + di (either part zero or normal; c = (-1)^csgn dc 2^(ec -
   FLINT_BITS), d likewise), the frames of the parts: the logarithm by
   the register path, t = y log x on the grid of 2^-(3 FLINT_BITS)
   exactly from the frames (whose errors, scaled by |y|, make the
   absolute error of t: 2^et units of the grid), exp(Re t) through the
   fixed-point exponential kernel (|Re t| < 2^62) and the sine and
   cosine of Im t (|Im t| < 2^63) through the reduction of two-limb
   mantissas.  Declines (ARF_ELEM_UNCERTAIN) beyond those ranges, for
   |y| >= 2^52 (the error of t would exceed the precision) and when the
   reduction declines. */
FLINT_FORCE_INLINE int
_acf_small_pow_1(_acf_frame_struct * re, _acf_frame_struct * im,
    ulong da, slong ea, int asgn, ulong db, slong eb, int bsgn,
    ulong dc, slong ec, int csgn, ulong dd, slong ed, int dsgn)
{
    _acf_frame_struct FL, FT;
    ulong Z1[5], Z2[5], Z3[5], Z4[5], R[5], W[5];
    ulong e2, e1, e0, s2, s1, s0, g2, g1, g0, c2, c1, c0, v2, v1, v0, errE, errT, ev;
    slong et, kE, e, bits;
    int rsgn, wsgn, a = 0, sneg = 0, ns, nc, negs, negc;
    unsigned int c;

    if ((dc != 0 && ec > 52) || (dd != 0 && ed > 52))
        return ARF_ELEM_UNCERTAIN;
    if (_acf_small_log_1(&FL, &FT, da, ea, asgn, db, eb, bsgn) != 0)
        return ARF_ELEM_UNCERTAIN;

    /* the terms c Lr, d Lt (the real part), c Lt, d Lr (the imaginary
       part) and the exponent of the absolute error of t on the grid */
    et = 0;
    if (dc != 0)
    {
        _acf_small_pow_term_1(Z1, FL.y[2], FL.y[1], FL.y[0], FL.k, dc, ec);
        _acf_small_pow_term_1(Z3, FT.y[2], FT.y[1], FT.y[0], FT.k, dc, ec);
        et = FLINT_MAX(et, (slong) FLINT_BIT_COUNT(FL.err) + FL.k + ec + FLINT_BITS);
        et = FLINT_MAX(et, (slong) FLINT_BIT_COUNT(FT.err) + FT.k + ec + FLINT_BITS);
    }
    else
    {
        flint_mpn_zero(Z1, 5);
        flint_mpn_zero(Z3, 5);
    }
    if (dd != 0)
    {
        _acf_small_pow_term_1(Z2, FT.y[2], FT.y[1], FT.y[0], FT.k, dd, ed);
        _acf_small_pow_term_1(Z4, FL.y[2], FL.y[1], FL.y[0], FL.k, dd, ed);
        et = FLINT_MAX(et, (slong) FLINT_BIT_COUNT(FT.err) + FT.k + ed + FLINT_BITS);
        et = FLINT_MAX(et, (slong) FLINT_BIT_COUNT(FL.err) + FL.k + ed + FLINT_BITS);
    }
    else
    {
        flint_mpn_zero(Z2, 5);
        flint_mpn_zero(Z4, 5);
    }
    et += 2;    /* (four terms and the truncations) */

    /* Re t = c Lr - d Lt, Im t = c Lt + d Lr */
    _acf_small_pow_sum_1(R, &rsgn, Z1, csgn ^ FL.sgnbit, Z2, dsgn ^ FT.sgnbit ^ 1);
    _acf_small_pow_sum_1(W, &wsgn, Z3, csgn ^ FT.sgnbit, Z4, dsgn ^ FL.sgnbit);

    if (R[4] != 0 || R[3] >= (UWORD(1) << 62) || W[4] != 0 || W[3] >= (UWORD(1) << 63))
        return ARF_ELEM_UNCERTAIN;

    /* E = exp(Re t) */
    kE = _arf_small_exp_fixed_1(&e2, &e1, &e0, R[3], R[2], R[1], R[0], rsgn);
    errE = ARF_SMALL_EXP_ERR_2 + ((et >= FLINT_BITS - 1) ? (UWORD(1) << (et - FLINT_BITS + 1)) : 1);

    /* v from Im t: directly for |Im t| < 1, else the reduction of its
       top two limbs (the dropped bits below 2^(e - 2 FLINT_BITS)) */
    ev = 0;
    if (W[3] == 0)
    {
        v2 = W[2]; v1 = W[1]; v0 = W[0];
    }
    else
    {
        ulong d1, d0;
        c = flint_clz(W[3]);
        bits = 4 * FLINT_BITS - c;
        e = bits - 3 * FLINT_BITS;
        if (c == 0)
        {
            d1 = W[3]; d0 = W[2];
        }
        else
        {
            d1 = (W[3] << c) | (W[2] >> (FLINT_BITS - c));
            d0 = (W[2] << c) | (W[1] >> (FLINT_BITS - c));
        }
        if (e > 56 || _arf_small_trig_reduce_2(&a, &sneg, &v2, &v1, &v0, d1, d0, e) != 0)
            return ARF_ELEM_UNCERTAIN;
        ev = 2 + (UWORD(1) << e);
    }
    _mp_real_small_sin_cos_g_3(&s2, &s1, &s0, &g2, &g1, &g0, v2, v1, v0);
    /* the frames of sin v, cos v at two fraction limbs (the third
       dropped), within the kernel's error (below a unit) + 2 + ev,
       plus the error of t */
    errT = 3 + ev + ((et >= FLINT_BITS) ? (UWORD(1) << (et - FLINT_BITS)) : 1);
    sub_dddmmmsss(c2, c1, c0, UWORD(1), UWORD(0), UWORD(0), UWORD(0), g2, g1);
    ARF_SMALL_TRIG_SIGNS(ns, nc, a, sneg, wsgn);
    negs = ((a & 1) ? nc : ns) ^ wsgn;
    negc = (a & 1) ? ns : nc;
    if (a & 1)
    {
        /* |sin| from cos v, |cos| from sin v */
        ulong t2 = c2, t1 = c1, t0 = c0;
        c2 = 0; c1 = s2; c0 = s1;
        s2 = t2; s1 = t1; s0 = t0;
    }
    else
    {
        s0 = s1; s1 = s2; s2 = 0;
    }

    _acf_small_mul_3(re->y + 2, re->y + 1, re->y + 0, e2, e1, e0, c2, c1, c0);
    _acf_small_mul_3(im->y + 2, im->y + 1, im->y + 0, e2, e1, e0, s2, s1, s0);
    re->k = im->k = kE;
    re->sgnbit = negc;
    im->sgnbit = negs;
    re->err = im->err = _ACF_SMALL_MUL_ERR_3(errE, 1, errT, _ACF_SMALL_EXP_BOUND);
    return 0;
}

#endif /* FLINT_BITS == 64 */

#endif
