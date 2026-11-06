/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Template for the complex types dNc and dNcb, instantiated by d1.c ...
   d4.c after template.inc (whose functions it uses), with DC(name) ->
   dNc_name, DCB(name) -> dNcb_name, DCV(name) -> _dNc_name and
   DCBV(name) -> _dNcb_name.

   Only the performance-critical functions are compiled per format: the
   arithmetic (which calls the real arithmetic, inlined where it is
   small, as for N = 1), the vector operations and the dot products
   (which call the real vector operations: an array of len complex
   numbers is also an array of 2 len real numbers, and the dot products
   work on blocks of rearranged copies), and their gr methods. The rest
   (division, square roots, the elementary functions, conversions) is
   implemented once for all formats in complex.c, through the table of
   the real operations of the format defined here (DFV(ops), DBV(ops)),
   and the functions below are thin wrappers of it.

   Every function allows aliasing of the output with any input. */

#define DC_EXACT_ZERO(x) ((x)->d[0] == 0.0)     /* plain: d[0] = 0 only for 0 */

/* for functions that are rarely used or only call expensive ones:
   compiled for size, without inlining the real arithmetic */
#if defined(__GNUC__)
#define DC_COLD __attribute__((cold))
#else
#define DC_COLD
#endif

/* ------------------------------------------------------------------ */
/* Plain complex numbers                                               */
/* ------------------------------------------------------------------ */

FLINT_FORCE_INLINE void DC(_zero_inline)(DC(t) res) { DF(_zero_inline)(&res->re); DF(_zero_inline)(&res->im); }
void DC(zero)(DC(t) res) { DC(_zero_inline)(res); }
FLINT_FORCE_INLINE void DC(_one_inline)(DC(t) res) { DF(_one_inline)(&res->re); DF(_zero_inline)(&res->im); }
void DC(one)(DC(t) res) { DC(_one_inline)(res); }
FLINT_FORCE_INLINE void DC(_onei_inline)(DC(t) res) { DF(_zero_inline)(&res->re); DF(_one_inline)(&res->im); }
void DC(onei)(DC(t) res) { DC(_onei_inline)(res); }
FLINT_FORCE_INLINE void DC(_set_inline)(DC(t) res, DC(srcptr) x) { *res = *x; }
void DC(set)(DC(t) res, DC(srcptr) x) { DC(_set_inline)(res, x); }
FLINT_FORCE_INLINE void DC(_set_d_inline)(DC(t) res, double x) { DF(_set_d_inline)(&res->re, x); DF(_zero_inline)(&res->im); }
void DC(set_d)(DC(t) res, double x) { DC(_set_d_inline)(res, x); }
void DC(set_d_d)(DC(t) res, double x, double y) { DF(_set_d_inline)(&res->re, x); DF(_set_d_inline)(&res->im, y); }
void DC(set_re_im)(DC(t) res, DF(srcptr) x, DF(srcptr) y)
{
    DF(struct) t = *y;
    res->re = *x;
    res->im = t;
}

FLINT_FORCE_INLINE void DC(_neg_inline)(DC(t) res, DC(srcptr) x) { DF(_neg_inline)(&res->re, &x->re); DF(_neg_inline)(&res->im, &x->im); }
void DC(neg)(DC(t) res, DC(srcptr) x) { DC(_neg_inline)(res, x); }
FLINT_FORCE_INLINE void DC(_conj_inline)(DC(t) res, DC(srcptr) x) { res->re = x->re; DF(_neg_inline)(&res->im, &x->im); }
void DC(conj)(DC(t) res, DC(srcptr) x) { DC(_conj_inline)(res, x); }
FLINT_FORCE_INLINE void DC(_add_inline)(DC(t) res, DC(srcptr) x, DC(srcptr) y) { DF(_add_impl)(&res->re, &x->re, &y->re); DF(_add_impl)(&res->im, &x->im, &y->im); }
void DC(add)(DC(t) res, DC(srcptr) x, DC(srcptr) y) { DC(_add_inline)(res, x, y); }
FLINT_FORCE_INLINE void DC(_sub_inline)(DC(t) res, DC(srcptr) x, DC(srcptr) y) { DF(_sub_impl)(&res->re, &x->re, &y->re); DF(_sub_impl)(&res->im, &x->im, &y->im); }
void DC(sub)(DC(t) res, DC(srcptr) x, DC(srcptr) y) { DC(_sub_inline)(res, x, y); }

void DC(mul_re)(DC(t) res, DC(srcptr) x, DF(srcptr) y)
{
    DF(struct) c = *y;
    DF(_mul_impl)(&res->re, &x->re, &c);
    DF(_mul_impl)(&res->im, &x->im, &c);
}

/* (a + b i) i = -b + a i */
void DC(mul_onei)(DC(t) res, DC(srcptr) x)
{
    DF(struct) t = x->re;
    DF(_neg_inline)(&res->re, &x->im);
    res->im = t;
}

FLINT_FORCE_INLINE void
DC(_mul_2exp_si_inline)(DC(t) res, DC(srcptr) x, slong e)
{
    DF(_mul_2exp_si_impl)(&res->re, &x->re, e);
    DF(_mul_2exp_si_impl)(&res->im, &x->im, e);
}
void DC(mul_2exp_si)(DC(t) res, DC(srcptr) x, slong e) { DC(_mul_2exp_si_inline)(res, x, e); }

/* a real (or imaginary) factor takes two real products; otherwise
   (a c - b d) + (a d + b c) i */
DFLOAT_ARITH_IMPL void
DC(_mul_impl)(DC(t) res, DC(srcptr) x, DC(srcptr) y)
{
    DF(t) p, q, r, s;
    if (DC_EXACT_ZERO(&y->im))
    {
        DF(_mul_impl)(r, &x->re, &y->re);
        DF(_mul_impl)(s, &x->im, &y->re);
    }
    else if (DC_EXACT_ZERO(&x->im))
    {
        DF(_mul_impl)(r, &y->re, &x->re);
        DF(_mul_impl)(s, &y->im, &x->re);
    }
    else
    {
        DF(_mul_impl)(p, &x->re, &y->re);
        DF(_mul_impl)(q, &x->im, &y->im);
        DF(_mul_impl)(r, &x->re, &y->im);
        DF(_mul_impl)(s, &x->im, &y->re);
        DF(_add_impl)(s, s, r);
        DF(_sub_impl)(r, p, q);
    }
    res->re = *r;
    res->im = *s;
}
void DC(mul)(DC(t) res, DC(srcptr) x, DC(srcptr) y) { DC(_mul_impl)(res, x, y); }

/* (a + b)(a - b) + 2 a b i (the real part with a relative error) */
DFLOAT_ARITH_IMPL void
DC(_sqr_impl)(DC(t) res, DC(srcptr) x)
{
    DF(t) p, q, r;
    if (DC_EXACT_ZERO(&x->im))
    {
        DF(_sqr_impl)(&res->re, &x->re);
        DF(_zero_inline)(&res->im);
        return;
    }
    DF(_add_impl)(p, &x->re, &x->im);
    DF(_sub_impl)(q, &x->re, &x->im);
    DF(_mul_impl)(r, &x->re, &x->im);
    DF(_mul_impl)(&res->re, p, q);
    DF(_mul_2exp_si_impl)(&res->im, r, 1);
}
void DC(sqr)(DC(t) res, DC(srcptr) x) { DC(_sqr_impl)(res, x); }

int DC(is_zero)(DC(srcptr) x) { return DF(_is_zero_inline)(&x->re) && DF(_is_zero_inline)(&x->im); }
int DC(is_one)(DC(srcptr) x) { return DF(_is_one_inline)(&x->re) && DF(_is_zero_inline)(&x->im); }
int DC(is_real)(DC(srcptr) x) { return DF(_is_zero_inline)(&x->im); }
int DC(is_finite)(DC(srcptr) x) { return DF(_is_finite_inline)(&x->re) && DF(_is_finite_inline)(&x->im); }
int DC(equal)(DC(srcptr) x, DC(srcptr) y) { return DF(_equal_inline)(&x->re, &y->re) && DF(_equal_inline)(&x->im, &y->im); }


/* ---- division, square roots, absolute value ----

   Fast paths, with the real arithmetic of the format called directly:
   a real or imaginary divisor (two real divisions), and operands whose
   nonzero parts have leading components of magnitude in [2^-250, 2^250]
   (DC_MODERATE), for which no square, product of up to three parts or
   reciprocal of these (nor their trailing components) can overflow or
   underflow, so that no rescaling is needed. Everything else (zero,
   nonfinite or extreme operands) goes through the general functions of
   complex.c, which also give the same results up to rounding.

   * x / y = x conj(y) / |y|^2 and 1 / y = conj(y) / |y|^2: for N = 1 by
     two divisions by |y|^2, for N > 1 by one reciprocal and two
     products (cheaper than two divisions at these precisions).

   * sqrt(x) for x = a + bi from m = |x| and w = m + |a| (no
     cancellation): q = sqrt(w / 2) and b / (2q), except for N = 2, where
     v = 1 / sqrt(2w), q = w v and b v (one reciprocal square root and two
     products, cheaper than a square root and a division; for N >= 3 the
     reciprocal square root of the format is neither faster nor as
     accurate). The parts are (q, b / (2q)) for a >= 0 and
     (|b| / (2q), sgn(b) q) for a < 0.

   * 1 / sqrt(x) = conj(sqrt(x)) / m: for N <= 2, v = 1 / sqrt(m (w^2 +
     b^2)) (= 1 / sqrt(2 m^2 w)) with the parts w v and b v; for N >= 3,
     g = 1 / (sqrt(w / 2) m) with the parts (w / 2) g and (b / 2) g. */

/* the real arithmetic in the fast paths: for N = 2 the kernels inlined
   (the results are finite there, so that the treatment of nonfinite
   values in the real functions is not needed), otherwise the functions */
#if N == 2
#define DC_FMUL(r, x, y) DF(_mul_kernel)((r)->d, NULL, (x)->d, (y)->d, 0)
#define DC_FSQR(r, x) DF(_sqr_kernel)((r)->d, NULL, (x)->d, 0)
#define DC_FADD(r, x, y) DF(_add_kernel)((r)->d, NULL, (x)->d, (y)->d, 0)
#define DC_FSUB(r, x, y) DF(_sub_kernel)((r)->d, NULL, (x)->d, (y)->d, 0)
#else
#define DC_FMUL(r, x, y) DF(_mul_impl)(r, x, y)
#define DC_FSQR(r, x) DF(_sqr_impl)(r, x)
#define DC_FADD(r, x, y) DF(_add_impl)(r, x, y)
#define DC_FSUB(r, x, y) DF(_sub_impl)(r, x, y)
#endif

#define DC_MODERATE(v) (fabs(v) >= 0x1p-250 && fabs(v) <= 0x1p250)
#define DC_MODERATE_OR_ZERO(v) ((v) == 0.0 || DC_MODERATE(v))

/* returns 0 if not handled */
FLINT_FORCE_INLINE int
DC(_div_fast)(DC(t) res, DC(srcptr) x, DC(srcptr) y)
{
    DF(t) s, t, p, q;
    double c0 = y->re.d[0], d0 = y->im.d[0];

    if (d0 == 0.0 && c0 != 0.0)
    {
        /* (a + bi) / c */
        DF(_div_inline)(p, &x->re, &y->re);
        DF(_div_inline)(&res->im, &x->im, &y->re);
        res->re = *p;
        return 1;
    }

    if (c0 == 0.0 && d0 != 0.0)
    {
        /* (a + bi) / (di) = b / d - (a / d) i */
        DF(_div_inline)(p, &x->im, &y->im);
        DF(_div_inline)(q, &x->re, &y->im);
        res->re = *p;
        DF(_neg_inline)(&res->im, q);
        return 1;
    }

    if (!DC_MODERATE(c0) || !DC_MODERATE(d0) ||
        !DC_MODERATE_OR_ZERO(x->re.d[0]) || !DC_MODERATE_OR_ZERO(x->im.d[0]))
        return 0;

    /* (a c + b d) + (b c - a d) i over c^2 + d^2 */
    DC_FSQR(s, &y->re);
    DC_FSQR(t, &y->im);
    DC_FADD(s, s, t);
    DC_FMUL(p, &x->re, &y->re);
    DC_FMUL(t, &x->im, &y->im);
    DC_FADD(p, p, t);
    DC_FMUL(q, &x->im, &y->re);
    DC_FMUL(t, &x->re, &y->im);
    DC_FSUB(q, q, t);
    if (N == 1)
    {
        res->re.d[0] = p->d[0] / s->d[0];
        res->im.d[0] = q->d[0] / s->d[0];
    }
    else
    {
        DF(_inv_inline)(t, s);
        DC_FMUL(&res->re, p, t);
        DC_FMUL(&res->im, q, t);
    }
    return 1;
}

FLINT_FORCE_INLINE int
DC(_inv_fast)(DC(t) res, DC(srcptr) x)
{
    DF(t) s, t, p;
    double c0 = x->re.d[0], d0 = x->im.d[0];

    if (d0 == 0.0 && c0 != 0.0)
    {
        DF(_inv_inline)(&res->re, &x->re);
        DF(_zero_inline)(&res->im);
        return 1;
    }

    if (c0 == 0.0 && d0 != 0.0)
    {
        /* 1 / (di) = -i / d */
        DF(_inv_inline)(p, &x->im);
        DF(_zero_inline)(&res->re);
        DF(_neg_inline)(&res->im, p);
        return 1;
    }

    if (!DC_MODERATE(c0) || !DC_MODERATE(d0))
        return 0;

    DC_FSQR(s, &x->re);
    DC_FSQR(t, &x->im);
    DC_FADD(s, s, t);
    if (N == 1)
    {
        res->re.d[0] = c0 / s->d[0];
        res->im.d[0] = -d0 / s->d[0];
    }
    else
    {
        DF(_inv_inline)(t, s);
        DC_FMUL(p, &x->im, t);
        DC_FMUL(&res->re, &x->re, t);
        DF(_neg_inline)(&res->im, p);
    }
    return 1;
}

/* m = |x| and w = m + |a| for x with moderate parts; u = b^2 */
FLINT_FORCE_INLINE void
DC(_abs_w)(DF(t) m, DF(t) w, DF(t) u, DC(srcptr) x)
{
    DF(t) s;
    DC_FSQR(u, &x->im);
    DC_FSQR(s, &x->re);
    DC_FADD(s, s, u);
    DF(_sqrt_inline)(m, s);
    if (x->re.d[0] >= 0.0)
        DC_FADD(w, m, &x->re);
    else
        DC_FSUB(w, m, &x->re);
}

FLINT_FORCE_INLINE int
DC(_sqrt_fast)(DC(t) res, DC(srcptr) x)
{
    DF(t) m, w, u, v, q, r;
    double a0 = x->re.d[0], b0 = x->im.d[0];

    if (!DC_MODERATE(a0) || !DC_MODERATE(b0))
        return 0;

    DC(_abs_w)(m, w, u, x);

    if (N == 1)
    {
        q->d[0] = sqrt(0.5 * w->d[0]);
        r->d[0] = b0 / (2.0 * q->d[0]);
    }
    else if (N == 2)
    {
        DF(_mul_2exp_si_impl)(v, w, 1);
        DF(_rsqrt_inline)(v, v);
        DC_FMUL(q, w, v);
        DC_FMUL(r, &x->im, v);
    }
    else
    {
        DF(_mul_2exp_si_impl)(v, w, -1);
        DF(_sqrt_inline)(q, v);
        DF(_mul_2exp_si_impl)(v, q, 1);
        DF(_div_inline)(r, &x->im, v);
    }

    if (a0 >= 0.0)
    {
        res->re = *q;
        res->im = *r;
    }
    else
    {
        DF(_abs_inline)(&res->re, r);
        if (b0 < 0.0)
            DF(_neg_inline)(&res->im, q);
        else
            res->im = *q;
    }
    return 1;
}

FLINT_FORCE_INLINE int
DC(_rsqrt_fast)(DC(t) res, DC(srcptr) x)
{
    DF(t) m, w, u, v, q, r;
    double a0 = x->re.d[0], b0 = x->im.d[0];

    if (!DC_MODERATE(a0) || !DC_MODERATE(b0))
        return 0;

    DC(_abs_w)(m, w, u, x);
    if (N <= 2)
    {
        DC_FSQR(v, w);
        DC_FADD(v, v, u);
        DC_FMUL(v, v, m);
        DF(_rsqrt_inline)(v, v);
        DC_FMUL(q, w, v);
        DC_FMUL(r, &x->im, v);
    }
    else
    {
        /* with h = w / 2 and g = 1 / (sqrt(h) m): q = h g, r = (b / 2) g */
        DF(_mul_2exp_si_impl)(w, w, -1);
        DF(_sqrt_inline)(v, w);
        DC_FMUL(v, v, m);
        DF(_inv_inline)(v, v);
        DC_FMUL(q, w, v);
        DF(_mul_2exp_si_impl)(u, &x->im, -1);
        DC_FMUL(r, u, v);
    }

    /* (q, -b v) for a >= 0, (|b| v, -sgn(b) q) for a < 0 */
    if (a0 >= 0.0)
    {
        res->re = *q;
        DF(_neg_inline)(&res->im, r);
    }
    else
    {
        DF(_abs_inline)(&res->re, r);
        if (b0 < 0.0)
            res->im = *q;
        else
            DF(_neg_inline)(&res->im, q);
    }
    return 1;
}

FLINT_FORCE_INLINE int
DC(_abs_fast)(DF(t) res, DC(srcptr) x)
{
    DF(t) m, w, u;
    double a0 = x->re.d[0], b0 = x->im.d[0];

    if (b0 == 0.0)
    {
        DF(_abs_inline)(res, &x->re);
        return 1;
    }
    if (a0 == 0.0)
    {
        DF(_abs_inline)(res, &x->im);
        return 1;
    }
    if (!DC_MODERATE(a0) || !DC_MODERATE(b0))
        return 0;

    DC_FSQR(u, &x->im);
    DC_FSQR(w, &x->re);
    DC_FADD(w, w, u);
    DF(_sqrt_inline)(m, w);
    *res = *m;
    return 1;
}

void DC(div)(DC(t) res, DC(srcptr) x, DC(srcptr) y)
{
    if (!DC(_div_fast)(res, x, y))
        GR_IGNORE(_dfloat_complex_div(&DFV(ops), res->re.d, x->re.d, y->re.d, DFLOAT_MODE_AUTO));
}
void DC(inv)(DC(t) res, DC(srcptr) x)
{
    if (!DC(_inv_fast)(res, x))
        GR_IGNORE(_dfloat_complex_inv(&DFV(ops), res->re.d, x->re.d, DFLOAT_MODE_AUTO));
}
void DC(sqrt)(DC(t) res, DC(srcptr) x)
{
    if (!DC(_sqrt_fast)(res, x))
        _dfloat_complex_sqrt(&DFV(ops), res->re.d, x->re.d, DFLOAT_MODE_AUTO);
}
void DC(rsqrt)(DC(t) res, DC(srcptr) x)
{
    if (!DC(_rsqrt_fast)(res, x))
        GR_IGNORE(_dfloat_complex_rsqrt(&DFV(ops), res->re.d, x->re.d, DFLOAT_MODE_AUTO));
}
void DC(abs)(DF(t) res, DC(srcptr) x)
{
    if (!DC(_abs_fast)(res, x))
        _dfloat_complex_abs(&DFV(ops), res->d, x->re.d);
}

/* ---- through complex.c ---- */

void DC(arg)(DF(t) res, DC(srcptr) x) { _dfloat_complex_arg(&DFV(ops), res->d, x->re.d); }
void DC(exp)(DC(t) res, DC(srcptr) x) { _dfloat_complex_exp(&DFV(ops), res->re.d, x->re.d); }
void DC(log)(DC(t) res, DC(srcptr) x) { GR_IGNORE(_dfloat_complex_log(&DFV(ops), res->re.d, x->re.d)); }
void DC(sin_cos)(DC(t) sn, DC(t) cs, DC(srcptr) x) { _dfloat_complex_sin_cos(&DFV(ops), sn->re.d, cs->re.d, x->re.d); }
void DC(sin)(DC(t) res, DC(srcptr) x) { _dfloat_complex_sin_cos(&DFV(ops), res->re.d, NULL, x->re.d); }
void DC(cos)(DC(t) res, DC(srcptr) x) { _dfloat_complex_sin_cos(&DFV(ops), NULL, res->re.d, x->re.d); }
void DC(sinh)(DC(t) res, DC(srcptr) x) { _dfloat_complex_sinh_cosh(&DFV(ops), res->re.d, NULL, x->re.d); }
void DC(cosh)(DC(t) res, DC(srcptr) x) { _dfloat_complex_sinh_cosh(&DFV(ops), NULL, res->re.d, x->re.d); }
void DC(tan)(DC(t) res, DC(srcptr) x) { GR_IGNORE(_dfloat_complex_tan(&DFV(ops), res->re.d, x->re.d)); }
void DC(tanh)(DC(t) res, DC(srcptr) x) { GR_IGNORE(_dfloat_complex_tanh(&DFV(ops), res->re.d, x->re.d)); }
void DC(pow)(DC(t) res, DC(srcptr) x, DC(srcptr) y) { GR_IGNORE(_dfloat_complex_pow(&DFV(ops), res->re.d, x->re.d, y->re.d)); }
void DC(get_acb)(acb_t res, DC(srcptr) x) { _dfloat_complex_get_acb(&DFV(ops), res, x->re.d); }
void DC(set_acb)(DC(t) res, const acb_t x) { _dfloat_complex_set_acb(&DFV(ops), res->re.d, x); }
char * DC(get_str)(DC(srcptr) x, slong digits) { return _dfloat_complex_get_str_ops(&DFV(ops), x->re.d, digits); }
void DC(print)(DC(srcptr) x) { char * s = DC(get_str)(x, 0); flint_printf("%s", s); flint_free(s); }
void DC(randtest)(DC(t) res, flint_rand_t state) { _dfloat_complex_randtest(&DFV(ops), res->re.d, state, 0); }
void DC(randtest_special)(DC(t) res, flint_rand_t state) { _dfloat_complex_randtest(&DFV(ops), res->re.d, state, 1); }


/* ------------------------------------------------------------------ */
/* Complex balls                                                       */
/* ------------------------------------------------------------------ */

/* the real products with a run-time error mode (see template.inc) */
#define DCB_RMUL(r, x, y) DB(_mul_ool)(r, x, y, mode)
#define DCB_RSQR(r, x) DB(_sqr_ool)(r, x, mode)

FLINT_FORCE_INLINE void DCB(_zero_inline)(DCB(t) res) { DB(_zero_inline)(&res->re); DB(_zero_inline)(&res->im); }
void DCB(zero)(DCB(t) res) { DCB(_zero_inline)(res); }
FLINT_FORCE_INLINE void DCB(_one_inline)(DCB(t) res) { DB(_one_inline)(&res->re); DB(_zero_inline)(&res->im); }
void DCB(one)(DCB(t) res) { DCB(_one_inline)(res); }
FLINT_FORCE_INLINE void DCB(_onei_inline)(DCB(t) res) { DB(_zero_inline)(&res->re); DB(_one_inline)(&res->im); }
void DCB(onei)(DCB(t) res) { DCB(_onei_inline)(res); }
void DCB(indeterminate)(DCB(t) res) { DB(_indeterminate_inline)(&res->re); DB(_indeterminate_inline)(&res->im); }
void DCB(set)(DCB(t) res, DCB(srcptr) x) { *res = *x; }
FLINT_FORCE_INLINE void DCB(_set_d_inline)(DCB(t) res, double x) { DB(_set_d_inline)(&res->re, x); DB(_zero_inline)(&res->im); }
void DCB(set_d)(DCB(t) res, double x) { DCB(_set_d_inline)(res, x); }
void DCB(set_d_d)(DCB(t) res, double x, double y) { DB(_set_d_inline)(&res->re, x); DB(_set_d_inline)(&res->im, y); }
void DCB(set_re_im)(DCB(t) res, DB(srcptr) x, DB(srcptr) y)
{
    DB(struct) t = *y;
    res->re = *x;
    res->im = t;
}

FLINT_FORCE_INLINE void DCB(_neg_inline)(DCB(t) res, DCB(srcptr) x) { DB(_neg_inline)(&res->re, &x->re); DB(_neg_inline)(&res->im, &x->im); }
void DCB(neg)(DCB(t) res, DCB(srcptr) x) { DCB(_neg_inline)(res, x); }
FLINT_FORCE_INLINE void
DCB(_conj_inline)(DCB(t) res, DCB(srcptr) x)
{
    res->re = x->re;
    if (FLINT_UNLIKELY(DB(_is_W)(&res->re)))
        DB(_indeterminate_inline)(&res->re);
    DB(_neg_inline)(&res->im, &x->im);
}
void DCB(conj)(DCB(t) res, DCB(srcptr) x) { DCB(_conj_inline)(res, x); }
FLINT_FORCE_INLINE void DCB(_add_inline)(DCB(t) res, DCB(srcptr) x, DCB(srcptr) y) { DB(_add_impl)(&res->re, &x->re, &y->re); DB(_add_impl)(&res->im, &x->im, &y->im); }
void DCB(add)(DCB(t) res, DCB(srcptr) x, DCB(srcptr) y) { DCB(_add_inline)(res, x, y); }
FLINT_FORCE_INLINE void DCB(_sub_inline)(DCB(t) res, DCB(srcptr) x, DCB(srcptr) y) { DB(_sub_impl)(&res->re, &x->re, &y->re); DB(_sub_impl)(&res->im, &x->im, &y->im); }
void DCB(sub)(DCB(t) res, DCB(srcptr) x, DCB(srcptr) y) { DCB(_sub_inline)(res, x, y); }

void DCB(mul_re)(DCB(t) res, DCB(srcptr) x, DB(srcptr) y)
{
    DB(struct) c = *y;
    DB(_mul_impl)(&res->re, &x->re, &c);
    DB(_mul_impl)(&res->im, &x->im, &c);
}

void DCB(mul_onei)(DCB(t) res, DCB(srcptr) x)
{
    DB(struct) t = x->re;
    DB(_neg_inline)(&res->re, &x->im);
    res->im = t;
    if (FLINT_UNLIKELY(DB(_is_W)(&res->im)))
        DB(_indeterminate_inline)(&res->im);
}

FLINT_FORCE_INLINE void
DCB(_mul_2exp_si_inline)(DCB(t) res, DCB(srcptr) x, slong e)
{
    DB(_mul_2exp_si_impl)(&res->re, &x->re, e);
    DB(_mul_2exp_si_impl)(&res->im, &x->im, e);
}
void DCB(mul_2exp_si)(DCB(t) res, DCB(srcptr) x, slong e) { DCB(_mul_2exp_si_inline)(res, x, e); }

/* as the plain version; an exactly real or imaginary factor takes two
   real products (which is also tighter) */
static DFLOAT_NOINLINE void
DCB(_mul_m)(DCB(t) res, DCB(srcptr) x, DCB(srcptr) y, int mode)
{
    DB(t) p, q, r, s;
    if (DB(_is_zero_inline)(&y->im))
    {
        DCB_RMUL(r, &x->re, &y->re);
        DCB_RMUL(s, &x->im, &y->re);
    }
    else if (DB(_is_zero_inline)(&x->im))
    {
        DCB_RMUL(r, &y->re, &x->re);
        DCB_RMUL(s, &y->im, &x->re);
    }
    else if (DB(_is_zero_inline)(&y->re))
    {
        /* (a + b i) (d i) = -b d + a d i */
        DCB_RMUL(p, &x->im, &y->im);
        DCB_RMUL(s, &x->re, &y->im);
        DB(_neg_inline)(r, p);
    }
    else if (DB(_is_zero_inline)(&x->re))
    {
        DCB_RMUL(p, &y->im, &x->im);
        DCB_RMUL(s, &y->re, &x->im);
        DB(_neg_inline)(r, p);
    }
    else
    {
        DCB_RMUL(p, &x->re, &y->re);
        DCB_RMUL(q, &x->im, &y->im);
        DCB_RMUL(r, &x->re, &y->im);
        DCB_RMUL(s, &x->im, &y->re);
        DB(_add_impl)(s, s, r);
        DB(_sub_impl)(r, p, q);
    }
    res->re = *r;
    res->im = *s;
}

void DCB(mul)(DCB(t) res, DCB(srcptr) x, DCB(srcptr) y) { DCB(_mul_m)(res, x, y, DFLOAT_MODE_AUTO); }

/* (a^2 - b^2) + 2 a b i (tighter than (a + b)(a - b) for balls) */
static DFLOAT_NOINLINE void
DCB(_sqr_m)(DCB(t) res, DCB(srcptr) x, int mode)
{
    DB(t) p, q, r;
    if (DB(_is_zero_inline)(&x->im))
    {
        DCB_RSQR(&res->re, &x->re);
        DB(_zero_inline)(&res->im);
        return;
    }
    if (DB(_is_zero_inline)(&x->re))
    {
        DCB_RSQR(p, &x->im);
        DB(_neg_inline)(&res->re, p);
        DB(_zero_inline)(&res->im);
        return;
    }
    DCB_RSQR(p, &x->re);
    DCB_RSQR(q, &x->im);
    DCB_RMUL(r, &x->re, &x->im);
    DB(_sub_impl)(&res->re, p, q);
    DB(_mul_2exp_si_impl)(&res->im, r, 1);
}

void DCB(sqr)(DCB(t) res, DCB(srcptr) x) { DCB(_sqr_m)(res, x, DFLOAT_MODE_AUTO); }

int DCB(is_zero)(DCB(srcptr) x) { return DB(_is_zero_inline)(&x->re) && DB(_is_zero_inline)(&x->im); }
int DCB(is_one)(DCB(srcptr) x)
{
    return x->re.rad == 0.0 && x->re.d[0] == 1.0 && DF(_tail_abs_sum)(x->re.d) == 0.0 && DB(_is_zero_inline)(&x->im);
}
int DCB(is_real)(DCB(srcptr) x) { return DB(_is_zero_inline)(&x->im); }
int DCB(is_exact)(DCB(srcptr) x) { return DB(is_exact)(&x->re) && DB(is_exact)(&x->im); }
int DCB(is_finite)(DCB(srcptr) x) { return DB(is_finite)(&x->re) && DB(is_finite)(&x->im); }
int DCB(contains_zero)(DCB(srcptr) x) { return DB(_contains_zero_inline)(&x->re) && DB(_contains_zero_inline)(&x->im); }
int DCB(contains)(DCB(srcptr) x, DCB(srcptr) y) { return DB(contains)(&x->re, &y->re) && DB(contains)(&x->im, &y->im); }
int DCB(overlaps)(DCB(srcptr) x, DCB(srcptr) y) { return DB(_overlaps_inline)(&x->re, &y->re) && DB(_overlaps_inline)(&x->im, &y->im); }
int DCB(equal)(DCB(srcptr) x, DCB(srcptr) y) { return DB(_equal_inline)(&x->re, &y->re) && DB(_equal_inline)(&x->im, &y->im); }


/* ---- through complex.c ---- */

int DCB(div)(DCB(t) res, DCB(srcptr) x, DCB(srcptr) y) { return _dfloat_complex_div(&DBV(ops), res->re.d, x->re.d, y->re.d, DFLOAT_MODE_AUTO); }
int DCB(inv)(DCB(t) res, DCB(srcptr) x) { return _dfloat_complex_inv(&DBV(ops), res->re.d, x->re.d, DFLOAT_MODE_AUTO); }
void DCB(sqrt)(DCB(t) res, DCB(srcptr) x) { _dfloat_complex_sqrt(&DBV(ops), res->re.d, x->re.d, DFLOAT_MODE_AUTO); }
int DCB(rsqrt)(DCB(t) res, DCB(srcptr) x) { return _dfloat_complex_rsqrt(&DBV(ops), res->re.d, x->re.d, DFLOAT_MODE_AUTO); }
void DCB(abs)(DB(t) res, DCB(srcptr) x) { _dfloat_complex_abs(&DBV(ops), res->d, x->re.d); }
void DCB(arg)(DB(t) res, DCB(srcptr) x) { _dfloat_complex_arg(&DBV(ops), res->d, x->re.d); }
void DCB(exp)(DCB(t) res, DCB(srcptr) x) { _dfloat_complex_exp(&DBV(ops), res->re.d, x->re.d); }
int DCB(log)(DCB(t) res, DCB(srcptr) x) { return _dfloat_complex_log(&DBV(ops), res->re.d, x->re.d); }
void DCB(sin_cos)(DCB(t) sn, DCB(t) cs, DCB(srcptr) x) { _dfloat_complex_sin_cos(&DBV(ops), sn->re.d, cs->re.d, x->re.d); }
void DCB(sin)(DCB(t) res, DCB(srcptr) x) { _dfloat_complex_sin_cos(&DBV(ops), res->re.d, NULL, x->re.d); }
void DCB(cos)(DCB(t) res, DCB(srcptr) x) { _dfloat_complex_sin_cos(&DBV(ops), NULL, res->re.d, x->re.d); }
void DCB(sinh)(DCB(t) res, DCB(srcptr) x) { _dfloat_complex_sinh_cosh(&DBV(ops), res->re.d, NULL, x->re.d); }
void DCB(cosh)(DCB(t) res, DCB(srcptr) x) { _dfloat_complex_sinh_cosh(&DBV(ops), NULL, res->re.d, x->re.d); }
int DCB(tan)(DCB(t) res, DCB(srcptr) x) { return _dfloat_complex_tan(&DBV(ops), res->re.d, x->re.d); }
int DCB(tanh)(DCB(t) res, DCB(srcptr) x) { return _dfloat_complex_tanh(&DBV(ops), res->re.d, x->re.d); }
int DCB(pow)(DCB(t) res, DCB(srcptr) x, DCB(srcptr) y) { return _dfloat_complex_pow(&DBV(ops), res->re.d, x->re.d, y->re.d); }
void DCB(get_acb)(acb_t res, DCB(srcptr) x) { _dfloat_complex_get_acb(&DBV(ops), res, x->re.d); }
void DCB(set_acb)(DCB(t) res, const acb_t x) { _dfloat_complex_set_acb(&DBV(ops), res->re.d, x); }
char * DCB(get_str)(DCB(srcptr) x, slong digits) { return _dfloat_complex_get_str_ops(&DBV(ops), x->re.d, digits); }
void DCB(print)(DCB(srcptr) x) { char * s = DCB(get_str)(x, 0); flint_printf("%s", s); flint_free(s); }
void DCB(randtest)(DCB(t) res, flint_rand_t state) { _dfloat_complex_randtest(&DBV(ops), res->re.d, state, 0); }
void DCB(randtest_special)(DCB(t) res, flint_rand_t state) { _dfloat_complex_randtest(&DBV(ops), res->re.d, state, 1); }

/* ---- vectors ---- */

/* the complex terms per block of the rearranged copies */
#define DC_BLOCK 32

DFLOAT_FUNC_NOINLINE void DCV(vec_add)(DC(ptr) res, DC(srcptr) x, DC(srcptr) y, slong len)
{ DFV(vec_add)((DF(ptr)) res, (DF(srcptr)) x, (DF(srcptr)) y, 2 * len); }

DFLOAT_FUNC_NOINLINE void DCV(vec_sub)(DC(ptr) res, DC(srcptr) x, DC(srcptr) y, slong len)
{ DFV(vec_sub)((DF(ptr)) res, (DF(srcptr)) x, (DF(srcptr)) y, 2 * len); }

DFLOAT_FUNC_NOINLINE void DCV(vec_mul)(DC(ptr) res, DC(srcptr) x, DC(srcptr) y, slong len)
{
    slong i;
    for (i = 0; i < len; i++)
        DC(_mul_impl)(res + i, x + i, y + i);
}

/* ---- SIMD helpers for the complex vector kernels ---- */

#ifdef DFLOAT_HAVE_SIMD

/* the parts of x[0], x[stride], x[2 stride], x[3 stride] in SoA form
   (lane k gets element k) */
FLINT_FORCE_INLINE void
DC(_vload4)(vec4d * re, vec4d * im, DC(srcptr) x, slong stride)
{
    const double * p0 = x[0].re.d, * p1 = x[stride].re.d;
    const double * p2 = x[2 * stride].re.d, * p3 = x[3 * stride].re.d;
    int k;
    UNROLL
    for (k = 0; k < N; k++)
    {
        re[k] = vec4d_set_d4(p0[k], p1[k], p2[k], p3[k]);
        im[k] = vec4d_set_d4(p0[N + k], p1[N + k], p2[N + k], p3[N + k]);
    }
}

FLINT_FORCE_INLINE void
DC(_vstore4)(DC(ptr) x, slong stride, const vec4d * re, const vec4d * im)
{
    double * p0 = x[0].re.d, * p1 = x[stride].re.d;
    double * p2 = x[2 * stride].re.d, * p3 = x[3 * stride].re.d;
    double t[4];
    int k;
    UNROLL
    for (k = 0; k < N; k++)
    {
        vec4d_store_unaligned(t, re[k]);
        p0[k] = t[0]; p1[k] = t[1]; p2[k] = t[2]; p3[k] = t[3];
        vec4d_store_unaligned(t, im[k]);
        p0[N + k] = t[0]; p1[N + k] = t[1]; p2[N + k] = t[2]; p3[N + k] = t[3];
    }
}

FLINT_FORCE_INLINE void
DC(_vbcast)(vec4d * m, const double * c)
{
    int k;
    UNROLL
    for (k = 0; k < N; k++)
        m[k] = vec4d_set_d(c[k]);
}

/* nonzero if some lane of some of the given results is nonfinite (the
   sum of the absolute values of the components is inf or nan) */
FLINT_FORCE_INLINE int
DC(_vany_nonfinite4)(const vec4d * a, const vec4d * b, const vec4d * c, const vec4d * d)
{
    vec4d big = vec4d_set_d(DBL_MAX);
    vec4d m = vec4d_bit_or(vec4d_cmp_nle(DF(_vabs_sum)(a), big), vec4d_cmp_nle(DF(_vabs_sum)(b), big));
    m = vec4d_bit_or(m, vec4d_cmp_nle(DF(_vabs_sum)(c), big));
    m = vec4d_bit_or(m, vec4d_cmp_nle(DF(_vabs_sum)(d), big));
    return VANY(m);
}

/* nonzero if some lane is an exactly real complex number (zero head of
   the imaginary part): DC(_mul_impl) takes its two-product branch for
   such a factor, so the lane formula would not reproduce it */
FLINT_FORCE_INLINE int
DC(_vany_real)(const vec4d * im)
{
    return VANY(vec4d_cmp_eq(im[0], vec4d_zero()));
}

#endif

/* ---- the complex product by a fixed scalar, vectorized ----

   The scalar loops below are the definitions; the SIMD loops run the
   same kernels (in the same order, on the same operands) four elements
   at a time, with the zero head shift of the scalar wrappers done in
   the lanes, so that the results are identical, and redo a block of
   four in scalar only when some lane has a nonfinite value (which the
   scalar code folds into the head) or, for the products below, an
   exactly real factor (for which DC(_mul_impl) takes another branch).
   The products follow
   DC(_mul_impl) exactly: x c = (x_re c_re - x_im c_im) + (x_im c_re +
   x_re c_im) i. */

#ifdef DFLOAT_HAVE_SIMD
/* the plain kernels in the lanes, each result followed by the zero
   head shift of the scalar wrappers (DF(_vshift_zero_head)), so that
   the lanes reproduce DF(_mul_impl), DF(_add_impl) and DF(_sub_impl)
   up to their treatment of nonfinite values, which the callers leave
   to the scalar code */
FLINT_FORCE_INLINE void
DC(_vmul)(vec4d * r, const vec4d * x, const vec4d * y)
{
    DF(_v_mul_kernel)(r, NULL, x, y, 0, DFLOAT_KERNEL_STRONG);
    DF(_vshift_zero_head)(r);
}

FLINT_FORCE_INLINE void
DC(_vaddsub)(vec4d * r, const vec4d * x, const vec4d * y, int negate)
{
    DF(_v_addsub_kernel)(r, NULL, x, y, negate, 0, DFLOAT_KERNEL_STRONG);
    DF(_vshift_zero_head)(r);
}

/* lanes: (tre, tim) = x c for the broadcast scalar (cre, cim) */
FLINT_FORCE_INLINE void
DC(_vmul_scalar4)(vec4d * tre, vec4d * tim, const vec4d * xre, const vec4d * xim,
    const vec4d * cre, const vec4d * cim)
{
    vec4d p[N], q[N], r[N], s[N];
    DC(_vmul)(p, xre, cre);
    DC(_vmul)(q, xim, cim);
    DC(_vmul)(r, xre, cim);
    DC(_vmul)(s, xim, cre);
    DC(_vaddsub)(tim, s, r, 0);
    DC(_vaddsub)(tre, p, q, 1);
}
#endif

DFLOAT_FUNC_NOINLINE void DCV(vec_mul_scalar)(DC(ptr) res, DC(srcptr) x, slong len, DC(srcptr) c)
{
    DC(struct) cc = *c;
    slong i = 0;
    if (DC_EXACT_ZERO(&cc.im))
    {
        DFV(vec_mul_scalar)((DF(ptr)) res, (DF(srcptr)) x, 2 * len, &cc.re);
        return;
    }
#ifdef DFLOAT_HAVE_SIMD
    if (len >= 4)
    {
        vec4d cre[N], cim[N];
        DC(_vbcast)(cre, cc.re.d);
        DC(_vbcast)(cim, cc.im.d);
        for (; i + 4 <= len; i += 4)
        {
            vec4d xre[N], xim[N], tre[N], tim[N];
            int special;
            DC(_vload4)(xre, xim, x + i, 1);
            special = DC(_vany_real)(xim);
            if (FLINT_LIKELY(!special))
            {
                DC(_vmul_scalar4)(tre, tim, xre, xim, cre, cim);
                special = DC(_vany_nonfinite4)(tre, tim, tre, tim);
            }
            if (FLINT_LIKELY(!special))
            {
                DC(_vstore4)(res + i, 1, tre, tim);
            }
            else
            {
                DC(_mul_impl)(res + i, x + i, &cc);
                DC(_mul_impl)(res + i + 1, x + i + 1, &cc);
                DC(_mul_impl)(res + i + 2, x + i + 2, &cc);
                DC(_mul_impl)(res + i + 3, x + i + 3, &cc);
            }
        }
    }
#endif
    for (; i < len; i++)
        DC(_mul_impl)(res + i, x + i, &cc);
}

static void
DCV(_vec_addmul_scalar)(DC(ptr) res, DC(srcptr) x, slong len, DC(srcptr) c, int sub)
{
    DC(struct) cc = *c;
    DC(t) t;
    slong i = 0;
    if (DC_EXACT_ZERO(&cc.im))
    {
        if (sub)
            DFV(vec_submul_scalar)((DF(ptr)) res, (DF(srcptr)) x, 2 * len, &cc.re);
        else
            DFV(vec_addmul_scalar)((DF(ptr)) res, (DF(srcptr)) x, 2 * len, &cc.re);
        return;
    }
#ifdef DFLOAT_HAVE_SIMD
    if (len >= 4)
    {
        vec4d cre[N], cim[N];
        DC(_vbcast)(cre, cc.re.d);
        DC(_vbcast)(cim, cc.im.d);
        for (; i + 4 <= len; i += 4)
        {
            vec4d xre[N], xim[N], rre[N], rim[N], tre[N], tim[N], ore[N], oim[N];
            int special;
            DC(_vload4)(xre, xim, x + i, 1);
            special = DC(_vany_real)(xim);
            if (FLINT_LIKELY(!special))
            {
                DC(_vload4)(rre, rim, res + i, 1);
                DC(_vmul_scalar4)(tre, tim, xre, xim, cre, cim);
                DC(_vaddsub)(ore, rre, tre, sub);
                DC(_vaddsub)(oim, rim, tim, sub);
                special = DC(_vany_nonfinite4)(tre, tim, ore, oim);
            }
            if (FLINT_LIKELY(!special))
            {
                DC(_vstore4)(res + i, 1, ore, oim);
            }
            else
            {
                int k;
                for (k = 0; k < 4; k++)
                {
                    DC(_mul_impl)(t, x + i + k, &cc);
                    if (sub)
                        DC(_sub_inline)(res + i + k, res + i + k, t);
                    else
                        DC(_add_inline)(res + i + k, res + i + k, t);
                }
            }
        }
    }
#endif
    for (; i < len; i++)
    {
        DC(_mul_impl)(t, x + i, &cc);
        if (sub)
            DC(_sub_inline)(res + i, res + i, t);
        else
            DC(_add_inline)(res + i, res + i, t);
    }
}

/* ---- 2 x 2 transformation of a pair of vectors ----

   (u_i, v_i) <- (M0 u_i + M1 v_i, M2 u_i + M3 v_i).  Each part of each
   output is a sum of four real products, computed as such: the
   products and a single renormalization of their sum (_addn_kernel),
   which is both cheaper and more accurate than complex products and
   sums.  The coefficients are, per output part
   (with nMk the negated imaginary part of Mk):

     u_re = M0r u_re + nM0i u_im + M1r v_re + nM1i v_im
     u_im = M0r u_im +  M0i u_re + M1r v_im +  M1i v_re

   and the same with M2, M3 for v.  The terms with an exactly zero
   coefficient are left out, which changes nothing but the sign of a
   zero result (a zero product contributes exact zeros to every order
   sum): with a real diagonal (the plane rotations with real cosine)
   every part is a sum of three products, with a real matrix of two. */

#define DC_APPLY_GENERIC 0
#define DC_APPLY_REAL_DIAG 1
#define DC_APPLY_REAL 2

/* r = a[0] b[0] + ... + a[nt-1] b[nt-1], scalar (2 <= nt <= 4): the
   products, then their sum with a single renormalization */
FLINT_FORCE_INLINE void
DC(_sdotn)(double * r, const double * const * a, const double * const * b, int nt)
{
    double p[4][N], v;
    const double * xs[4];
    int m;
    for (m = 0; m < nt; m++)
    {
        DF(_mul_kernel)(p[m], NULL, a[m], b[m], 0);
        xs[m] = p[m];
    }
    if (FLINT_UNLIKELY(DF(_s_addn_kernel)(r, NULL, xs, nt, 0, DFLOAT_KERNEL_STRONG)))
        DF(_shift_zero_head)(r);
    v = a[0][0] * b[0][0];
    for (m = 1; m < nt; m++)
        v += a[m][0] * b[m][0];
    DF(_fix_nonfinite_v)(r, v);
}

/* the coefficient table: M0r M0i nM0i M1r M1i nM1i M2r ... nM3i */
FLINT_FORCE_INLINE void
DC(_apply_2x2_coeffs)(DF(struct) * cf, DC(srcptr) M)
{
    int j;
    for (j = 0; j < 4; j++)
    {
        cf[3 * j] = M[j].re;
        cf[3 * j + 1] = M[j].im;
        DF(_neg_inline)(cf + 3 * j + 2, &M[j].im);
    }
}

/* The terms of one output part (part = 0..3 for u_re, u_im, v_re,
   v_im) for the given variant, as indices into the coefficient table
   and the operands (0..3 for u_re, u_im, v_re, v_im); returns the
   number of terms.  With constant arguments this folds to constants.
   DC_APPLY_REAL_DIAG: M0 and M3 exactly real; DC_APPLY_REAL: all four
   exactly real (a matrix with only M1, M2 real is rare and runs as
   generic). */
FLINT_FORCE_INLINE int
DC(_apply_2x2_terms)(int * ci, int * xi, int part, int variant)
{
    /* the part reads row r of M (entries Mj, Mk with j = 2 r, k = j + 1),
       whose diagonal entry is Mj for r = 0 and Mk for r = 1 */
    int r = (part >= 2), j = 2 * r, k = j + 1, im = part & 1, n = 0;
    int skip_j, skip_k;
    skip_j = (variant == DC_APPLY_REAL) || (variant == DC_APPLY_REAL_DIAG && r == 0);
    skip_k = (variant == DC_APPLY_REAL) || (variant == DC_APPLY_REAL_DIAG && r == 1);
    /* Mj times u */
    ci[n] = 3 * j; xi[n] = im; n++;
    if (!skip_j)
    {
        ci[n] = 3 * j + (im ? 1 : 2); xi[n] = !im; n++;
    }
    /* Mk times v */
    ci[n] = 3 * k; xi[n] = 2 + im; n++;
    if (!skip_k)
    {
        ci[n] = 3 * k + (im ? 1 : 2); xi[n] = 2 + !im; n++;
    }
    return n;
}

/* one part (0..3: u_re, u_im, v_re, v_im) of the output, scalar */
FLINT_FORCE_INLINE void
DC(_apply_2x2_part)(double * r, const DF(struct) * cf, const double * const * x, int part, int variant)
{
    int ci[4], xi[4], nt, m;
    const double * as[4], * bs[4];
    nt = DC(_apply_2x2_terms)(ci, xi, part, variant);
    for (m = 0; m < nt; m++)
    {
        as[m] = cf[ci[m]].d;
        bs[m] = x[xi[m]];
    }
    DC(_sdotn)(r, as, bs, nt);
}

FLINT_FORCE_INLINE void
DC(_apply_2x2_1)(DC(ptr) u, DC(ptr) v, const DF(struct) * cf, int variant)
{
    DC(struct) a = *u, b = *v;
    const double * x[4];
    x[0] = a.re.d; x[1] = a.im.d; x[2] = b.re.d; x[3] = b.im.d;
    /* (the parts written out, so that each is specialized) */
    DC(_apply_2x2_part)(u->re.d, cf, x, 0, variant);
    DC(_apply_2x2_part)(u->im.d, cf, x, 1, variant);
    DC(_apply_2x2_part)(v->re.d, cf, x, 2, variant);
    DC(_apply_2x2_part)(v->im.d, cf, x, 3, variant);
}

FLINT_FORCE_INLINE void
DC(_apply_2x2_scalar)(DC(ptr) u, DC(ptr) v, slong stride, slong len, const DF(struct) * cf, int variant)
{
    slong i;
    for (i = 0; i < len; i++)
        DC(_apply_2x2_1)(u + i * stride, v + i * stride, cf, variant);
}

#ifdef DFLOAT_HAVE_SIMD
/* one part of the output in the lanes, as DC(_apply_2x2_part) */
FLINT_FORCE_INLINE void
DC(_vapply_2x2_part)(vec4d * r, const vec4d (* c)[N], vec4d (* x)[N], int part, int variant)
{
    int ci[4], xi[4], nt, m;
    vec4d p[4][N];
    const vec4d * xs[4];
    nt = DC(_apply_2x2_terms)(ci, xi, part, variant);
    for (m = 0; m < nt; m++)
    {
        DC(_vmul)(p[m], c[ci[m]], x[xi[m]]);
        xs[m] = p[m];
    }
    DF(_v_addn_kernel)(r, NULL, xs, nt, 0, DFLOAT_KERNEL_STRONG);
    DF(_vshift_zero_head)(r);
}

/* four (or, at the end, fewer) consecutive pairs in the lanes, as
   DC(_apply_2x2_1); with fewer than four the unused lanes compute a
   copy of the first pair and are not stored */
FLINT_FORCE_INLINE void
DC(_apply_2x2_block)(DC(ptr) u, DC(ptr) v, slong stride, slong cnt,
    const vec4d (* c)[N], const DF(struct) * cf, int variant)
{
    vec4d x[4][N], o[4][N];
    int k;

    if (cnt == 4)
    {
        DC(_vload4)(x[0], x[1], u, stride);
        DC(_vload4)(x[2], x[3], v, stride);
    }
    else
    {
        DC(struct) ub[4], vb[4];
        for (k = 0; k < 4; k++)
        {
            ub[k] = u[(k < cnt ? k : 0) * stride];
            vb[k] = v[(k < cnt ? k : 0) * stride];
        }
        DC(_vload4)(x[0], x[1], ub, 1);
        DC(_vload4)(x[2], x[3], vb, 1);
    }

    DC(_vapply_2x2_part)(o[0], c, x, 0, variant);
    DC(_vapply_2x2_part)(o[1], c, x, 1, variant);
    DC(_vapply_2x2_part)(o[2], c, x, 2, variant);
    DC(_vapply_2x2_part)(o[3], c, x, 3, variant);

    if (FLINT_UNLIKELY(DC(_vany_nonfinite4)(o[0], o[1], o[2], o[3])))
    {
        DC(_apply_2x2_scalar)(u, v, stride, cnt, cf, variant);
    }
    else if (cnt == 4)
    {
        DC(_vstore4)(u, stride, o[0], o[1]);
        DC(_vstore4)(v, stride, o[2], o[3]);
    }
    else
    {
        DC(struct) ub[4], vb[4];
        DC(_vstore4)(ub, 1, o[0], o[1]);
        DC(_vstore4)(vb, 1, o[2], o[3]);
        for (k = 0; k < cnt; k++)
        {
            u[k * stride] = ub[k];
            v[k * stride] = vb[k];
        }
    }
}

FLINT_FORCE_INLINE void
DC(_apply_2x2_simd)(DC(ptr) u, DC(ptr) v, slong stride, slong len, const DF(struct) * cf, int variant)
{
    vec4d c[12][N];
    slong i;
    int j;
    for (j = 0; j < 12; j++)
        DC(_vbcast)(c[j], cf[j].d);
    for (i = 0; i + 4 <= len; i += 4)
        DC(_apply_2x2_block)(u + i * stride, v + i * stride, stride, 4, (const vec4d (*)[N]) c, cf, variant);
    if (i < len)
        DC(_apply_2x2_block)(u + i * stride, v + i * stride, stride, len - i, (const vec4d (*)[N]) c, cf, variant);
}
#endif

/* the SIMD code also handles the short vectors: a block with unused
   lanes costs about as much as a single pair in scalar at N = 2, 3 and
   less at N = 4 */
#if N == 1
#define DC_APPLY_2X2_SIMD_MIN 4
#elif N < 4
#define DC_APPLY_2X2_SIMD_MIN 2
#else
#define DC_APPLY_2X2_SIMD_MIN 1
#endif

DFLOAT_FUNC_NOINLINE void
DCV(vec_apply_2x2)(DC(ptr) u, DC(ptr) v, slong stride, slong len, DC(srcptr) M)
{
    DF(struct) cf[12];
    int variant;

    DC(_apply_2x2_coeffs)(cf, M);

    if (DC_EXACT_ZERO(&M[0].im) && DC_EXACT_ZERO(&M[3].im))
        variant = (DC_EXACT_ZERO(&M[1].im) && DC_EXACT_ZERO(&M[2].im)) ? DC_APPLY_REAL : DC_APPLY_REAL_DIAG;
    else
        variant = DC_APPLY_GENERIC;

#ifdef DFLOAT_HAVE_SIMD
    if (len >= DC_APPLY_2X2_SIMD_MIN)
    {
        if (variant == DC_APPLY_GENERIC)
            DC(_apply_2x2_simd)(u, v, stride, len, cf, DC_APPLY_GENERIC);
        else if (variant == DC_APPLY_REAL_DIAG)
            DC(_apply_2x2_simd)(u, v, stride, len, cf, DC_APPLY_REAL_DIAG);
        else
            DC(_apply_2x2_simd)(u, v, stride, len, cf, DC_APPLY_REAL);
        return;
    }
#endif

    if (variant == DC_APPLY_GENERIC)
        DC(_apply_2x2_scalar)(u, v, stride, len, cf, DC_APPLY_GENERIC);
    else if (variant == DC_APPLY_REAL_DIAG)
        DC(_apply_2x2_scalar)(u, v, stride, len, cf, DC_APPLY_REAL_DIAG);
    else
        DC(_apply_2x2_scalar)(u, v, stride, len, cf, DC_APPLY_REAL);
}

#undef DC_APPLY_2X2_SIMD_MIN
#undef DC_APPLY_GENERIC
#undef DC_APPLY_REAL_DIAG
#undef DC_APPLY_REAL

/* ---- squared absolute value ---- */

FLINT_FORCE_INLINE void
DC(_abs2_impl)(DF(t) res, DC(srcptr) x)
{
    if (DC_EXACT_ZERO(&x->im))
        DF(_sqr_impl)(res, &x->re);
    else if (DC_EXACT_ZERO(&x->re))
        DF(_sqr_impl)(res, &x->im);
    else
    {
        DF(t) a, b;
        DF(_sqr_impl)(a, &x->re);
        DF(_sqr_impl)(b, &x->im);
        DF(_add_impl)(res, a, b);
    }
}

void DC(abs2)(DF(t) res, DC(srcptr) x) { DC(_abs2_impl)(res, x); }

DFLOAT_FUNC_NOINLINE void DCV(vec_addmul_scalar)(DC(ptr) res, DC(srcptr) x, slong len, DC(srcptr) c) { DCV(_vec_addmul_scalar)(res, x, len, c, 0); }
DFLOAT_FUNC_NOINLINE void DCV(vec_submul_scalar)(DC(ptr) res, DC(srcptr) x, slong len, DC(srcptr) c) { DCV(_vec_addmul_scalar)(res, x, len, c, 1); }

static DFLOAT_NOINLINE void
DCV(_rdot)(DF(t) res, DF(srcptr) initial, int subtract, DF(srcptr) x, DF(srcptr) y, slong len)
{
    DFV(_vec_dot_dispatch)(res, initial, subtract, x, y, len, 0);
}

/* The dot products: with X the 2 len real numbers of x, the real part
   is the real dot product of X with the copy (Re y_j, -Im y_j) and the
   imaginary part that of X with the copy (Im y_j, Re y_j), computed
   by blocks (y reversed for rev), the result of one block being the
   initial value of the next. */
static DFLOAT_NOINLINE void
DCV(_vec_dot)(DC(t) res, DC(srcptr) initial, int subtract, DC(srcptr) x, DC(srcptr) y, slong len, int rev)
{
    DF(struct) yc[2 * DC_BLOCK], ys[2 * DC_BLOCK];
    DF(t) sre, sim;
    slong i, j, n;
    int k;

#if N == 1
    /* plain doubles: a direct loop (the copies would cost more than the
       products), with two accumulators per part */
    {
        double r0 = 0.0, r1 = 0.0, i0 = 0.0, i1 = 0.0, a, b, c, d;
        slong st = rev ? -1 : 1;
        const DC(struct) * yy = rev ? y + len - 1 : y;
        for (i = 0; i + 2 <= len; i += 2)
        {
            a = x[i].re.d[0]; b = x[i].im.d[0]; c = yy[i * st].re.d[0]; d = yy[i * st].im.d[0];
            r0 += a * c - b * d;
            i0 += a * d + b * c;
            a = x[i + 1].re.d[0]; b = x[i + 1].im.d[0]; c = yy[(i + 1) * st].re.d[0]; d = yy[(i + 1) * st].im.d[0];
            r1 += a * c - b * d;
            i1 += a * d + b * c;
        }
        if (i < len)
        {
            a = x[i].re.d[0]; b = x[i].im.d[0]; c = yy[i * st].re.d[0]; d = yy[i * st].im.d[0];
            r0 += a * c - b * d;
            i0 += a * d + b * c;
        }
        r0 += r1;
        i0 += i1;
        if (subtract)
        {
            r0 = -r0;
            i0 = -i0;
        }
        if (initial != NULL)
        {
            r0 = initial->re.d[0] + r0;
            i0 = initial->im.d[0] + i0;
        }
        res->re.d[0] = r0;
        res->im.d[0] = i0;
        (void) yc; (void) ys; (void) sre; (void) sim; (void) j; (void) n; (void) k;
        return;
    }
#endif

    if (len <= 0)
    {
        if (initial == NULL)
            DC(_zero_inline)(res);
        else
            DC(_set_inline)(res, initial);
        return;
    }

    for (i = 0; i < len; i += n)
    {
        n = FLINT_MIN(DC_BLOCK, len - i);
        for (j = 0; j < n; j++)
        {
            DC(srcptr) yj = rev ? y + len - 1 - i - j : y + i + j;
            yc[2 * j] = yj->re;
            for (k = 0; k < N; k++)
                yc[2 * j + 1].d[k] = -yj->im.d[k];
            ys[2 * j] = yj->im;
            ys[2 * j + 1] = yj->re;
        }
        if (i == 0)
        {
            DCV(_rdot)(sre, (initial == NULL) ? NULL : &initial->re, subtract, (DF(srcptr)) (x + i), yc, 2 * n);
            DCV(_rdot)(sim, (initial == NULL) ? NULL : &initial->im, subtract, (DF(srcptr)) (x + i), ys, 2 * n);
        }
        else
        {
            DCV(_rdot)(sre, sre, subtract, (DF(srcptr)) (x + i), yc, 2 * n);
            DCV(_rdot)(sim, sim, subtract, (DF(srcptr)) (x + i), ys, 2 * n);
        }
    }
    res->re = *sre;
    res->im = *sim;
}

DFLOAT_FUNC_NOINLINE void DCV(vec_dot)(DC(t) res, DC(srcptr) initial, int subtract, DC(srcptr) x, DC(srcptr) y, slong len)
{ DCV(_vec_dot)(res, initial, subtract, x, y, len, 0); }

DFLOAT_FUNC_NOINLINE void DCV(vec_dot_rev)(DC(t) res, DC(srcptr) initial, int subtract, DC(srcptr) x, DC(srcptr) y, slong len)
{ DCV(_vec_dot)(res, initial, subtract, x, y, len, 1); }

/* exp by blocks: the real exponential of the real parts and the sine
   and cosine of the imaginary parts (the vectorized real functions),
   then two real vector products */
DFLOAT_FUNC_NOINLINE void DCV(vec_exp)(DC(ptr) res, DC(srcptr) x, slong len)
{
    DF(struct) a[DC_BLOCK], b[DC_BLOCK], e[DC_BLOCK], s[DC_BLOCK], c[DC_BLOCK];
    slong i, j, n;
    for (i = 0; i < len; i += n)
    {
        n = FLINT_MIN(DC_BLOCK, len - i);
        for (j = 0; j < n; j++)
        {
            a[j] = x[i + j].re;
            b[j] = x[i + j].im;
        }
        DFV(vec_exp)(e, a, n);
        DFV(vec_sin_cos)(s, c, b, n);
        DFV(vec_mul)(a, e, c, n);
        DFV(vec_mul)(b, e, s, n);
        for (j = 0; j < n; j++)
        {
            res[i + j].re = a[j];
            res[i + j].im = b[j];
        }
    }
}

/* ---- vectors ---- */

DFLOAT_FUNC_NOINLINE void DCBV(vec_add)(DCB(ptr) res, DCB(srcptr) x, DCB(srcptr) y, slong len)
{ DBV(vec_add)((DB(ptr)) res, (DB(srcptr)) x, (DB(srcptr)) y, 2 * len); }

DFLOAT_FUNC_NOINLINE void DCBV(vec_sub)(DCB(ptr) res, DCB(srcptr) x, DCB(srcptr) y, slong len)
{ DBV(vec_sub)((DB(ptr)) res, (DB(srcptr)) x, (DB(srcptr)) y, 2 * len); }

static void
DCBV(_vec_mul_m)(DCB(ptr) res, DCB(srcptr) x, DCB(srcptr) y, slong len, int mode)
{
    slong i;
    for (i = 0; i < len; i++)
        DCB(_mul_m)(res + i, x + i, y + i, mode);
}

DFLOAT_FUNC_NOINLINE void DCBV(vec_mul)(DCB(ptr) res, DCB(srcptr) x, DCB(srcptr) y, slong len) { DCBV(_vec_mul_m)(res, x, y, len, DFLOAT_MODE_AUTO); }

static void
DCBV(_vec_mul_scalar_m)(DCB(ptr) res, DCB(srcptr) x, slong len, DCB(srcptr) c, int mode)
{
    DCB(struct) cc = *c;
    slong i;
    if (DB(_is_zero_inline)(&cc.im))
        DBV(_vec_mul_scalar_m)((DB(ptr)) res, (DB(srcptr)) x, 2 * len, &cc.re, mode);
    else
        for (i = 0; i < len; i++)
            DCB(_mul_m)(res + i, x + i, &cc, mode);
}

static void
DCBV(_vec_addmul_scalar_m)(DCB(ptr) res, DCB(srcptr) x, slong len, DCB(srcptr) c, int sub, int mode)
{
    DCB(struct) cc = *c;
    DCB(t) t;
    slong i;
    if (DB(_is_zero_inline)(&cc.im))
    {
        DBV(_vec_addmul_scalar_m)((DB(ptr)) res, (DB(srcptr)) x, 2 * len, &cc.re, sub, mode);
        return;
    }
    for (i = 0; i < len; i++)
    {
        DCB(_mul_m)(t, x + i, &cc, mode);
        if (sub)
            DCB(_sub_inline)(res + i, res + i, t);
        else
            DCB(_add_inline)(res + i, res + i, t);
    }
}

DFLOAT_FUNC_NOINLINE void DCBV(vec_mul_scalar)(DCB(ptr) res, DCB(srcptr) x, slong len, DCB(srcptr) c) { DCBV(_vec_mul_scalar_m)(res, x, len, c, DFLOAT_MODE_AUTO); }
DFLOAT_FUNC_NOINLINE void DCBV(vec_addmul_scalar)(DCB(ptr) res, DCB(srcptr) x, slong len, DCB(srcptr) c) { DCBV(_vec_addmul_scalar_m)(res, x, len, c, 0, DFLOAT_MODE_AUTO); }
DFLOAT_FUNC_NOINLINE void DCBV(vec_submul_scalar)(DCB(ptr) res, DCB(srcptr) x, slong len, DCB(srcptr) c) { DCBV(_vec_addmul_scalar_m)(res, x, len, c, 1, DFLOAT_MODE_AUTO); }

static DFLOAT_NOINLINE void
DCBV(_rdot)(DB(t) res, DB(srcptr) initial, int subtract, DB(srcptr) x, DB(srcptr) y, slong len, int mode)
{
    (void) DBV_vec_dot_fwd(res, initial, subtract, x, y, len, mode);
}

/* the dot products, as the plain ones, with the real ball dot products */
static DFLOAT_NOINLINE void
DCBV(_vec_dot)(DCB(t) res, DCB(srcptr) initial, int subtract, DCB(srcptr) x, DCB(srcptr) y, slong len, int rev, int mode)
{
    DB(struct) yc[2 * DC_BLOCK], ys[2 * DC_BLOCK];
    DB(t) sre, sim;
    slong i, j, n;
    int k;

    if (len <= 0)
    {
        /* (as the real ones: the initial value, normalized) */
        DCBV(_rdot)(&res->re, (initial == NULL) ? NULL : &initial->re, subtract, NULL, NULL, 0, mode);
        DCBV(_rdot)(&res->im, (initial == NULL) ? NULL : &initial->im, subtract, NULL, NULL, 0, mode);
        return;
    }

    for (i = 0; i < len; i += n)
    {
        n = FLINT_MIN(DC_BLOCK, len - i);
        for (j = 0; j < n; j++)
        {
            DCB(srcptr) yj = rev ? y + len - 1 - i - j : y + i + j;
            yc[2 * j] = yj->re;
            for (k = 0; k < N; k++)
                yc[2 * j + 1].d[k] = -yj->im.d[k];
            yc[2 * j + 1].rad = yj->im.rad;
            ys[2 * j] = yj->im;
            ys[2 * j + 1] = yj->re;
        }
        if (i == 0)
        {
            DCBV(_rdot)(sre, (initial == NULL) ? NULL : &initial->re, subtract, (DB(srcptr)) (x + i), yc, 2 * n, mode);
            DCBV(_rdot)(sim, (initial == NULL) ? NULL : &initial->im, subtract, (DB(srcptr)) (x + i), ys, 2 * n, mode);
        }
        else
        {
            DCBV(_rdot)(sre, sre, subtract, (DB(srcptr)) (x + i), yc, 2 * n, mode);
            DCBV(_rdot)(sim, sim, subtract, (DB(srcptr)) (x + i), ys, 2 * n, mode);
        }
    }
    res->re = *sre;
    res->im = *sim;
}

DFLOAT_FUNC_NOINLINE void DCBV(vec_dot)(DCB(t) res, DCB(srcptr) initial, int subtract, DCB(srcptr) x, DCB(srcptr) y, slong len)
{ DCBV(_vec_dot)(res, initial, subtract, x, y, len, 0, DFLOAT_MODE_AUTO); }

DFLOAT_FUNC_NOINLINE void DCBV(vec_dot_rev)(DCB(t) res, DCB(srcptr) initial, int subtract, DCB(srcptr) x, DCB(srcptr) y, slong len)
{ DCBV(_vec_dot)(res, initial, subtract, x, y, len, 1, DFLOAT_MODE_AUTO); }

DFLOAT_FUNC_NOINLINE void DCBV(vec_exp)(DCB(ptr) res, DCB(srcptr) x, slong len)
{
    DB(struct) a[DC_BLOCK], b[DC_BLOCK], e[DC_BLOCK], s[DC_BLOCK], c[DC_BLOCK];
    slong i, j, n;
    for (i = 0; i < len; i += n)
    {
        n = FLINT_MIN(DC_BLOCK, len - i);
        for (j = 0; j < n; j++)
        {
            a[j] = x[i + j].re;
            b[j] = x[i + j].im;
        }
        DBV(vec_exp)(e, a, n);
        DBV(vec_sin_cos)(s, c, b, n);
        DBV(vec_mul)(a, e, c, n);
        DBV(vec_mul)(b, e, s, n);
        for (j = 0; j < n; j++)
        {
            res[i + j].re = a[j];
            res[i + j].im = b[j];
        }
    }
}


/* ------------------------------------------------------------------ */
/* The real operations for complex.c                                   */
/* ------------------------------------------------------------------ */

/* (the plain functions that differ in signature, through wrappers; the
   others are the functions themselves, whose pointer arguments are
   pointers to the doubles of the structs) */
static DC_COLD void DF(_op_mul)(double * r, const double * x, const double * y, int mode) { DF(_mul_impl)((DF(ptr)) r, (DF(srcptr)) x, (DF(srcptr)) y); }
static DC_COLD void DF(_op_sqr)(double * r, const double * x, int mode) { DF(_sqr_impl)((DF(ptr)) r, (DF(srcptr)) x); }
static DC_COLD int DF(_op_div)(double * r, const double * x, const double * y, int mode) { DF(_div_inline)((DF(ptr)) r, (DF(srcptr)) x, (DF(srcptr)) y); return GR_SUCCESS; }
static DC_COLD int DF(_op_sqrt)(double * r, const double * x, int mode) { DF(_sqrt_inline)((DF(ptr)) r, (DF(srcptr)) x); return GR_SUCCESS; }
static DC_COLD int DF(_op_log)(double * r, const double * x) { DF(log)((DF(ptr)) r, (DF(srcptr)) x); return GR_SUCCESS; }
static DC_COLD int DF(_op_log1p)(double * r, const double * x) { DF(log1p)((DF(ptr)) r, (DF(srcptr)) x); return GR_SUCCESS; }
static DC_COLD int DF(_op_tan)(double * r, const double * x) { DF(tan)((DF(ptr)) r, (DF(srcptr)) x); return GR_SUCCESS; }
static DC_COLD int DF(_op_contains_zero)(const double * x) { return x[0] == 0.0; }
static DC_COLD int DF(_op_cmp_all_d)(const double * x, double c, int sign, int strict)
{
    DF(t) t;
    int r;
    if (DF(_is_nan_inline)((DF(srcptr)) x))
        return 0;
    DF(_set_d_inline)(t, c);
    r = DF(_cmp_inline)((DF(srcptr)) x, t);
    return (sign > 0) ? (strict ? r > 0 : r >= 0) : (strict ? r < 0 : r <= 0);
}
static DC_COLD int DF(_op_get_fmpz)(fmpz_t r, const double * x) { return GRF(get_fmpz)(r, (DF(srcptr)) x, NULL); }
static DC_COLD void DF(_op_cmul)(double * r, const double * x, const double * y, int mode) { DC(_mul_impl)((DC(ptr)) r, (DC(srcptr)) x, (DC(srcptr)) y); }
static DC_COLD void DF(_op_csqr)(double * r, const double * x, int mode) { DC(_sqr_impl)((DC(ptr)) r, (DC(srcptr)) x); }

#define DC_OP1(f) ((void (*)(double *)) (f))
#define DC_OP2(f) ((void (*)(double *, const double *)) (f))
#define DC_OP3(f) ((void (*)(double *, const double *, const double *)) (f))

const _dfloat_ops_struct DFV(ops) =
{
    N, 0,
    DC_OP1(DF(zero)),
    DC_OP1(DF(one)),
    (void (*)(double *, double)) DF(set_d),
    (void (*)(double *, const fmpz_t)) DF(set_fmpz),
    DC_OP2(DF(neg)),
    DC_OP2(DF(abs)),
    DC_OP3(DF(add)),
    DC_OP3(DF(sub)),
    DF(_op_mul),
    DF(_op_sqr),
    DF(_op_div),
    DF(_op_sqrt),
    (void (*)(double *, const double *, double)) DF(mul_d),
    (void (*)(double *, const double *, slong)) DF(mul_2exp_si),
    DC_OP2(DF(exp)),
    DC_OP2(DF(expm1)),
    DF(_op_log),
    DF(_op_log1p),
    (void (*)(double *, double *, const double *)) DF(sin_cos),
    DF(_op_tan),
    DC_OP2(DF(tanh)),
    DC_OP3(DF(atan2)),
    DF(_op_contains_zero),
    DF(_op_cmp_all_d),
    DF(_op_get_fmpz),
    (void (*)(arb_t, const double *)) DF(get_arb),
    (void (*)(double *, const arb_t)) DF(set_arb),
    (void (*)(double *, flint_rand_t)) DF(randtest),
    (void (*)(double *, flint_rand_t)) DF(randtest_special),
    DF(_canonicalise_inplace),
    DF(_op_cmul),
    DF(_op_csqr),
};

const _dfloat_ops_struct DBV(ops) =
{
    N, 1,
    DC_OP1(DB(zero)),
    DC_OP1(DB(one)),
    (void (*)(double *, double)) DB(set_d),
    (void (*)(double *, const fmpz_t)) DB(set_fmpz),
    DC_OP2(DB(neg)),
    DC_OP2(DB(abs)),
    DC_OP3(DB(add)),
    DC_OP3(DB(sub)),
    (void (*)(double *, const double *, const double *, int)) DB(_mul_ool),
    (void (*)(double *, const double *, int)) DB(_sqr_ool),
    (int (*)(double *, const double *, const double *, int)) DB(_div_ool),
    (int (*)(double *, const double *, int)) DB(_sqrt_m),
    (void (*)(double *, const double *, double)) DB(mul_d),
    (void (*)(double *, const double *, slong)) DB(mul_2exp_si),
    DC_OP2(DB(exp)),
    DC_OP2(DB(expm1)),
    (int (*)(double *, const double *)) DB(log),
    (int (*)(double *, const double *)) DB(log1p),
    (void (*)(double *, double *, const double *)) DB(sin_cos),
    (int (*)(double *, const double *)) DB(tan),
    DC_OP2(DB(tanh)),
    DC_OP3(DB(atan2)),
    (int (*)(const double *)) DB(contains_zero),
    (int (*)(const double *, double, int, int)) DB(_cmp_all_d),
    (int (*)(fmpz_t, const double *)) DB(_get_fmpz),
    (void (*)(arb_t, const double *)) DB(get_arb),
    (void (*)(double *, const arb_t)) DB(set_arb),
    (void (*)(double *, flint_rand_t)) DB(randtest),
    (void (*)(double *, flint_rand_t)) DB(randtest_special),
    DF(_canonicalise_inplace),
    (void (*)(double *, const double *, const double *, int)) DCB(_mul_m),
    (void (*)(double *, const double *, int)) DCB(_sqr_m),
};

/* ------------------------------------------------------------------ */
/* Generic ring methods (the arithmetic and the vectors; the rest is  */
/* in complex.c)                                                       */
/* ------------------------------------------------------------------ */

#define GRC(name) DCV(gr_##name)
#define GRCB(name) DCBV(gr_##name)

/* the truth of "both" from the truths of the parts */
FLINT_FORCE_INLINE truth_t
DC(_truth_and)(truth_t a, truth_t b)
{
    if (a == T_FALSE || b == T_FALSE)
        return T_FALSE;
    if (a == T_TRUE && b == T_TRUE)
        return T_TRUE;
    return T_UNKNOWN;
}

/* ---- plain ---- */

static void GRC(init)(DC(ptr) x, gr_ctx_t ctx) { DC(_zero_inline)(x); }
static void GRC(clear)(DC(ptr) x, gr_ctx_t ctx) { }
static void GRC(swap)(DC(ptr) x, DC(ptr) y, gr_ctx_t ctx) { DC(struct) t = *x; *x = *y; *y = t; }
static int GRC(zero)(DC(ptr) res, gr_ctx_t ctx) { DC(_zero_inline)(res); return GR_SUCCESS; }
static int GRC(one)(DC(ptr) res, gr_ctx_t ctx) { DC(_one_inline)(res); return GR_SUCCESS; }
static int GRC(neg_one)(DC(ptr) res, gr_ctx_t ctx) { DC(_set_d_inline)(res, -1.0); return GR_SUCCESS; }
static int GRC(i)(DC(ptr) res, gr_ctx_t ctx) { DC(_onei_inline)(res); return GR_SUCCESS; }
static truth_t GRC(is_zero)(DC(srcptr) x, gr_ctx_t ctx) { return DC(_truth_and)(GRF(is_zero)(&x->re, ctx), GRF(is_zero)(&x->im, ctx)); }
static truth_t GRC(is_one)(DC(srcptr) x, gr_ctx_t ctx) { return DC(_truth_and)(GRF(is_one)(&x->re, ctx), GRF(is_zero)(&x->im, ctx)); }
static truth_t GRC(is_neg_one)(DC(srcptr) x, gr_ctx_t ctx) { return DC(_truth_and)(GRF(is_neg_one)(&x->re, ctx), GRF(is_zero)(&x->im, ctx)); }
static truth_t GRC(equal)(DC(srcptr) x, DC(srcptr) y, gr_ctx_t ctx) { return DC(_truth_and)(GRF(equal)(&x->re, &y->re, ctx), GRF(equal)(&x->im, &y->im, ctx)); }
static int GRC(set)(DC(ptr) res, DC(srcptr) x, gr_ctx_t ctx) { DC(_set_inline)(res, x); return GR_SUCCESS; }
static int GRC(neg)(DC(ptr) res, DC(srcptr) x, gr_ctx_t ctx) { DC(_neg_inline)(res, x); return GR_SUCCESS; }
static int GRC(add)(DC(ptr) res, DC(srcptr) x, DC(srcptr) y, gr_ctx_t ctx) { DC(_add_inline)(res, x, y); return GR_SUCCESS; }
static int GRC(sub)(DC(ptr) res, DC(srcptr) x, DC(srcptr) y, gr_ctx_t ctx) { DC(_sub_inline)(res, x, y); return GR_SUCCESS; }
static int GRC(mul)(DC(ptr) res, DC(srcptr) x, DC(srcptr) y, gr_ctx_t ctx) { DC(_mul_impl)(res, x, y); return GR_SUCCESS; }
static int GRC(sqr)(DC(ptr) res, DC(srcptr) x, gr_ctx_t ctx) { DC(_sqr_impl)(res, x); return GR_SUCCESS; }
static int GRC(mul_two)(DC(ptr) res, DC(srcptr) x, gr_ctx_t ctx) { DC(_mul_2exp_si_inline)(res, x, 1); return GR_SUCCESS; }
static int GRC(mul_2exp_si)(DC(ptr) res, DC(srcptr) x, slong e, gr_ctx_t ctx) { DC(_mul_2exp_si_inline)(res, x, e); return GR_SUCCESS; }
static int GRC(conj)(DC(ptr) res, DC(srcptr) x, gr_ctx_t ctx) { DC(_conj_inline)(res, x); return GR_SUCCESS; }
static int GRC(re)(DC(ptr) res, DC(srcptr) x, gr_ctx_t ctx) { res->re = x->re; DF(_zero_inline)(&res->im); return GR_SUCCESS; }

/* (GR_DOMAIN at the exact singularities, as in complex.c, whose plain
   functions return no other status) */
static int GRC(div)(DC(ptr) res, DC(srcptr) x, DC(srcptr) y, gr_ctx_t ctx)
{
    if (DC_EXACT_ZERO(&y->re) && DC_EXACT_ZERO(&y->im))
        return GR_DOMAIN;
    DC(div)(res, x, y);
    return GR_SUCCESS;
}
static int GRC(inv)(DC(ptr) res, DC(srcptr) x, gr_ctx_t ctx)
{
    if (DC_EXACT_ZERO(&x->re) && DC_EXACT_ZERO(&x->im))
        return GR_DOMAIN;
    DC(inv)(res, x);
    return GR_SUCCESS;
}
static int GRC(sqrt)(DC(ptr) res, DC(srcptr) x, gr_ctx_t ctx)
{
    DC(sqrt)(res, x);
    return GR_SUCCESS;
}
static int GRC(rsqrt)(DC(ptr) res, DC(srcptr) x, gr_ctx_t ctx)
{
    if (DC_EXACT_ZERO(&x->re) && DC_EXACT_ZERO(&x->im))
        return GR_DOMAIN;
    DC(rsqrt)(res, x);
    return GR_SUCCESS;
}
static int GRC(abs)(DC(ptr) res, DC(srcptr) x, gr_ctx_t ctx)
{
    DF(t) t;
    DC(abs)(t, x);
    res->re = *t;
    DF(_zero_inline)(&res->im);
    return GR_SUCCESS;
}
static int GRC(im)(DC(ptr) res, DC(srcptr) x, gr_ctx_t ctx) { res->re = x->im; DF(_zero_inline)(&res->im); return GR_SUCCESS; }
static int GRC(abs2)(DC(ptr) res, DC(srcptr) x, gr_ctx_t ctx) { DC(_abs2_impl)(&res->re, x); DF(_zero_inline)(&res->im); return GR_SUCCESS; }

/* only real numbers are ordered (as for acf and nfloat_complex) */
static int GRC(cmp)(int * res, DC(srcptr) x, DC(srcptr) y, gr_ctx_t ctx)
{
    if (!DF(_is_zero_inline)(&x->im) || !DF(_is_zero_inline)(&y->im))
    {
        *res = 0;
        if (DF(_is_nan_inline)(&x->im) || DF(_is_nan_inline)(&y->im))
            return GR_UNABLE;
        return GR_DOMAIN;
    }
    return GRF(cmp)(res, &x->re, &y->re, ctx);
}

/* |x| against |y|: by |x|^2 against |y|^2 (no square root) when the
   larger head of each argument lies in [2^-400, 2^480] or is zero, so
   that the squares neither overflow nor lose their tails to underflow
   at any N (a smaller part whose square underflows is then below the
   precision of the larger one's); otherwise through the scaled abs */
static int GRC(cmpabs)(int * res, DC(srcptr) x, DC(srcptr) y, gr_ctx_t ctx)
{
    double hx = FLINT_MAX(fabs(x->re.d[0]), fabs(x->im.d[0]));
    double hy = FLINT_MAX(fabs(y->re.d[0]), fabs(y->im.d[0]));
    DF(t) a, b;
    if ((hx == 0.0 || (hx >= 0x1p-400 && hx <= 0x1p480)) &&
        (hy == 0.0 || (hy >= 0x1p-400 && hy <= 0x1p480)))
    {
        DC(_abs2_impl)(a, x);
        DC(_abs2_impl)(b, y);
    }
    else
    {
        DC(abs)(a, x);
        DC(abs)(b, y);
    }
    return GRF(cmp)(res, a, b, ctx);
}

static int GRC(vec_apply_2x2)(DC(ptr) u, DC(ptr) v, slong stride, slong len, DC(srcptr) M, gr_ctx_t ctx)
{
    DCV(vec_apply_2x2)(u, v, stride, len, M);
    return GR_SUCCESS;
}

static void GRC(vec_init)(DC(ptr) vec, slong len, gr_ctx_t ctx)
{
    slong i;
    for (i = 0; i < len; i++)
        DC(_zero_inline)(vec + i);
}
static void GRC(vec_clear)(DC(ptr) vec, slong len, gr_ctx_t ctx) { }
static int GRC(vec_set)(DC(ptr) res, DC(srcptr) x, slong len, gr_ctx_t ctx)
{
    if (res != x)
        memcpy(res, x, len * sizeof(DC(struct)));
    return GR_SUCCESS;
}
static int GRC(vec_zero)(DC(ptr) vec, slong len, gr_ctx_t ctx) { GRC(vec_init)(vec, len, ctx); return GR_SUCCESS; }
static int GRC(vec_add)(DC(ptr) res, DC(srcptr) x, DC(srcptr) y, slong len, gr_ctx_t ctx) { DCV(vec_add)(res, x, y, len); return GR_SUCCESS; }
static int GRC(vec_sub)(DC(ptr) res, DC(srcptr) x, DC(srcptr) y, slong len, gr_ctx_t ctx) { DCV(vec_sub)(res, x, y, len); return GR_SUCCESS; }
static int GRC(vec_mul)(DC(ptr) res, DC(srcptr) x, DC(srcptr) y, slong len, gr_ctx_t ctx) { DCV(vec_mul)(res, x, y, len); return GR_SUCCESS; }
static int GRC(vec_mul_scalar)(DC(ptr) res, DC(srcptr) x, slong len, DC(srcptr) c, gr_ctx_t ctx) { DCV(vec_mul_scalar)(res, x, len, c); return GR_SUCCESS; }
static int GRC(vec_addmul_scalar)(DC(ptr) res, DC(srcptr) x, slong len, DC(srcptr) c, gr_ctx_t ctx) { DCV(vec_addmul_scalar)(res, x, len, c); return GR_SUCCESS; }
static int GRC(vec_submul_scalar)(DC(ptr) res, DC(srcptr) x, slong len, DC(srcptr) c, gr_ctx_t ctx) { DCV(vec_submul_scalar)(res, x, len, c); return GR_SUCCESS; }
static int GRC(vec_dot)(DC(ptr) res, DC(srcptr) initial, int subtract, DC(srcptr) x, DC(srcptr) y, slong len, gr_ctx_t ctx) { DCV(_vec_dot)(res, initial, subtract, x, y, len, 0); return GR_SUCCESS; }
static int GRC(vec_dot_rev)(DC(ptr) res, DC(srcptr) initial, int subtract, DC(srcptr) x, DC(srcptr) y, slong len, gr_ctx_t ctx) { DCV(_vec_dot)(res, initial, subtract, x, y, len, 1); return GR_SUCCESS; }
static int GRC(vec_exp)(DC(ptr) res, DC(srcptr) x, slong len, gr_ctx_t ctx) { DCV(vec_exp)(res, x, len); return GR_SUCCESS; }
static int GRC(vec_gather)(DC(ptr) res, DC(srcptr) vec, const slong * idx, slong len, gr_ctx_t ctx)
{
    slong i;
    for (i = 0; i < len; i++)
        res[i] = vec[idx[i]];
    return GR_SUCCESS;
}
static int GRC(vec_scatter)(DC(ptr) vec, const slong * idx, DC(srcptr) src, slong len, gr_ctx_t ctx)
{
    slong i;
    for (i = 0; i < len; i++)
        vec[idx[i]] = src[i];
    return GR_SUCCESS;
}

/* strong variants: the same, then both parts canonicalised */
static DFLOAT_NOINLINE void DC(_canon)(DC(ptr) x) { DF(_canonicalise_inplace)(x->re.d); DF(_canonicalise_inplace)(x->im.d); }
#define GRC_STRONG_UNARY(name) \
static DC_COLD int GRC(name##_strong)(DC(ptr) res, DC(srcptr) x, gr_ctx_t ctx) \
{ DC(name)(res, x); DC(_canon)(res); return GR_SUCCESS; }
#define GRC_STRONG_BINARY(name) \
static DC_COLD int GRC(name##_strong)(DC(ptr) res, DC(srcptr) x, DC(srcptr) y, gr_ctx_t ctx) \
{ DC(name)(res, x, y); DC(_canon)(res); return GR_SUCCESS; }
#define GRC_STRONG_VEC_BINARY(name) \
static DC_COLD int GRC(name##_strong)(DC(ptr) res, DC(srcptr) x, DC(srcptr) y, slong len, gr_ctx_t ctx) \
{ slong i; DCV(name)(res, x, y, len); for (i = 0; i < len; i++) DC(_canon)(res + i); return GR_SUCCESS; }
#define GRC_STRONG_VEC_SCALAR(name) \
static DC_COLD int GRC(name##_strong)(DC(ptr) res, DC(srcptr) x, slong len, DC(srcptr) c, gr_ctx_t ctx) \
{ slong i; DCV(name)(res, x, len, c); for (i = 0; i < len; i++) DC(_canon)(res + i); return GR_SUCCESS; }
#define GRC_STRONG_DOT(name, rev) \
static DC_COLD int GRC(name##_strong)(DC(ptr) res, DC(srcptr) initial, int subtract, DC(srcptr) x, DC(srcptr) y, slong len, gr_ctx_t ctx) \
{ DCV(_vec_dot)(res, initial, subtract, x, y, len, rev); DC(_canon)(res); return GR_SUCCESS; }

GRC_STRONG_BINARY(add)
GRC_STRONG_BINARY(sub)
GRC_STRONG_BINARY(mul)
GRC_STRONG_UNARY(sqr)
GRC_STRONG_VEC_BINARY(vec_add)
GRC_STRONG_VEC_BINARY(vec_sub)
GRC_STRONG_VEC_BINARY(vec_mul)
GRC_STRONG_VEC_SCALAR(vec_mul_scalar)
GRC_STRONG_VEC_SCALAR(vec_addmul_scalar)
GRC_STRONG_VEC_SCALAR(vec_submul_scalar)
GRC_STRONG_DOT(vec_dot, 0)
GRC_STRONG_DOT(vec_dot_rev, 1)
static DC_COLD int GRC(vec_exp_strong)(DC(ptr) res, DC(srcptr) x, slong len, gr_ctx_t ctx)
{ slong i; DCV(vec_exp)(res, x, len); for (i = 0; i < len; i++) DC(_canon)(res + i); return GR_SUCCESS; }
static DC_COLD int GRC(div_strong)(DC(ptr) res, DC(srcptr) x, DC(srcptr) y, gr_ctx_t ctx)
{ int status = GRC(div)(res, x, y, ctx); DC(_canon)(res); return status; }
static DC_COLD int GRC(inv_strong)(DC(ptr) res, DC(srcptr) x, gr_ctx_t ctx)
{ int status = GRC(inv)(res, x, ctx); DC(_canon)(res); return status; }
static DC_COLD int GRC(sqrt_strong)(DC(ptr) res, DC(srcptr) x, gr_ctx_t ctx)
{ int status = GRC(sqrt)(res, x, ctx); DC(_canon)(res); return status; }
static DC_COLD int GRC(rsqrt_strong)(DC(ptr) res, DC(srcptr) x, gr_ctx_t ctx)
{ int status = GRC(rsqrt)(res, x, ctx); DC(_canon)(res); return status; }
static DC_COLD int GRC(abs_strong)(DC(ptr) res, DC(srcptr) x, gr_ctx_t ctx)
{ int status = GRC(abs)(res, x, ctx); DC(_canon)(res); return status; }
static DC_COLD int GRC(abs2_strong)(DC(ptr) res, DC(srcptr) x, gr_ctx_t ctx)
{ GRC(abs2)(res, x, ctx); DC(_canon)(res); return GR_SUCCESS; }
static DC_COLD int GRC(vec_apply_2x2_strong)(DC(ptr) u, DC(ptr) v, slong stride, slong len, DC(srcptr) M, gr_ctx_t ctx)
{
    slong i;
    DCV(vec_apply_2x2)(u, v, stride, len, M);
    for (i = 0; i < len; i++)
    {
        DC(_canon)(u + i * stride);
        DC(_canon)(v + i * stride);
    }
    return GR_SUCCESS;
}

/* ---- LU decomposition: the basecase, as DFV(gr_mat_lu_classical)
   (the elimination column by column, four rows at a time in the SIMD
   lanes, with the products of DC(_vmul_scalar4); the pivot of largest
   |Re| + |Im| of the heads) ---- */

#define DC_ENTRY(M, i, j) (((DC(ptr)) (M)->entries) + (i) * (M)->stride + (j))

static slong
DC(_lu_pivot)(const gr_mat_t M, slong start_row, slong col)
{
    slong i, best = -1;
    double bv = 0.0, v;
    for (i = start_row; i < M->r; i++)
    {
        DC(srcptr) x = DC_ENTRY(M, i, col);
        v = fabs(x->re.d[0]) + fabs(x->im.d[0]);
        if (v != v)
            return i;
        if (v > bv)
        {
            bv = v;
            best = i;
        }
    }
    return best;
}

static void
DC(_lu_eliminate)(gr_mat_t LU, slong row, slong col, slong rank_col, DC(srcptr) d)
{
    slong m = LU->r, n = LU->c, j = row + 1, c;
    DC(struct) e[4];
    DC(t) t;
    DC(srcptr) prow = DC_ENTRY(LU, row, 0);

#ifdef DFLOAT_HAVE_SIMD
    for ( ; j + 4 <= m; j += 4)
    {
        slong k;
        vec4d ere[N], eim[N];
        for (k = 0; k < 4; k++)
            DC(_mul_impl)(e + k, DC_ENTRY(LU, j + k, col), d);
        DC(_vload4)(ere, eim, e, 1);
        for (c = col + 1; c < n; c++)
        {
            vec4d pre[N], pim[N], xre[N], xim[N], tre[N], tim[N], rre[N], rim[N];
            DC(ptr) x = DC_ENTRY(LU, j, c);
            DC(_vbcast)(pre, prow[c].re.d);
            DC(_vbcast)(pim, prow[c].im.d);
            DC(_vload4)(xre, xim, x, LU->stride);
            DC(_vmul_scalar4)(tre, tim, ere, eim, pre, pim);
            DC(_vaddsub)(rre, xre, tre, 1);
            DC(_vaddsub)(rim, xim, tim, 1);
            if (FLINT_LIKELY(!DC(_vany_nonfinite4)(tre, tim, rre, rim)))
            {
                DC(_vstore4)(x, LU->stride, rre, rim);
            }
            else
            {
                for (k = 0; k < 4; k++)
                {
                    DC(_mul_impl)(t, e + k, prow + c);
                    DC(_sub_inline)(DC_ENTRY(LU, j + k, c), DC_ENTRY(LU, j + k, c), t);
                }
            }
        }
        for (k = 0; k < 4; k++)
        {
            DC(_zero_inline)(DC_ENTRY(LU, j + k, col));
            DC(_set_inline)(DC_ENTRY(LU, j + k, rank_col), e + k);
        }
    }
#endif

    for ( ; j < m; j++)
    {
        DC(_mul_impl)(e, DC_ENTRY(LU, j, col), d);
        for (c = col + 1; c < n; c++)
        {
            DC(_mul_impl)(t, e, prow + c);
            DC(_sub_inline)(DC_ENTRY(LU, j, c), DC_ENTRY(LU, j, c), t);
        }
        DC(_zero_inline)(DC_ENTRY(LU, j, col));
        DC(_set_inline)(DC_ENTRY(LU, j, rank_col), e);
    }
}

int DCV(gr_mat_lu_classical)(slong * res_rank, slong * P, gr_mat_t LU, const gr_mat_t A, int rank_check, gr_ctx_t ctx);

int
DCV(gr_mat_lu_classical)(slong * res_rank, slong * P, gr_mat_t LU, const gr_mat_t A, int rank_check, gr_ctx_t ctx)
{
    slong i, m = A->r, n = A->c, r, rank = 0, row = 0, col = 0;
    DC(t) d;
    int status = GR_SUCCESS;

    for (i = 0; i < m; i++)
        P[i] = i;

    if (m == 0 || n == 0)
    {
        *res_rank = 0;
        return GR_SUCCESS;
    }

    status |= gr_mat_set(LU, A, ctx);

    while (row < m && col < n)
    {
        r = DC(_lu_pivot)(LU, row, col);

        if (r == -1)
        {
            if (rank_check)
            {
                rank = 0;
                break;
            }
            col++;
            continue;
        }

        rank++;

        if (r != row)
            _gr_mat_swap_rows(LU, P, row, r, ctx);

        DC(inv)(d, DC_ENTRY(LU, row, col));
        DC(_lu_eliminate)(LU, row, col, rank - 1, d);

        row++;
        col++;
    }

    *res_rank = rank;
    return status;
}

#undef DC_ENTRY

static int GRC(mat_mul)(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
{ return DFV(gr_mat_mul_cplx)(C, A, B, 1, ctx); }

static int GRC(mat_lu)(slong * rank, slong * P, gr_mat_t LU, const gr_mat_t A, int rank_check, gr_ctx_t ctx)
{
    if (A->r < DFLOAT_MAT_LU_CUTOFF || A->c < DFLOAT_MAT_LU_CUTOFF)
        return DCV(gr_mat_lu_classical)(rank, P, LU, A, rank_check, ctx);
    else
        return gr_mat_lu_recursive(rank, P, LU, A, rank_check, ctx);
}
static DC_COLD int GRC(mat_mul_strong)(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
{ return DFV(gr_mat_mul_cplx_strong)(C, A, B, 1, ctx); }

gr_method_tab_input DCV(gr_strong_methods_input)[] =
{
    {GR_METHOD_MAT_MUL,         (gr_funcptr) GRC(mat_mul_strong)},
    {GR_METHOD_ADD,             (gr_funcptr) GRC(add_strong)},
    {GR_METHOD_SUB,             (gr_funcptr) GRC(sub_strong)},
    {GR_METHOD_MUL,             (gr_funcptr) GRC(mul_strong)},
    {GR_METHOD_SQR,             (gr_funcptr) GRC(sqr_strong)},
    {GR_METHOD_VEC_ADD,         (gr_funcptr) GRC(vec_add_strong)},
    {GR_METHOD_VEC_SUB,         (gr_funcptr) GRC(vec_sub_strong)},
    {GR_METHOD_VEC_MUL,         (gr_funcptr) GRC(vec_mul_strong)},
    {GR_METHOD_VEC_MUL_SCALAR,  (gr_funcptr) GRC(vec_mul_scalar_strong)},
    {GR_METHOD_VEC_ADDMUL_SCALAR, (gr_funcptr) GRC(vec_addmul_scalar_strong)},
    {GR_METHOD_VEC_SUBMUL_SCALAR, (gr_funcptr) GRC(vec_submul_scalar_strong)},
    {GR_METHOD_VEC_DOT,         (gr_funcptr) GRC(vec_dot_strong)},
    {GR_METHOD_VEC_DOT_REV,     (gr_funcptr) GRC(vec_dot_rev_strong)},
    {GR_METHOD_VEC_EXP,         (gr_funcptr) GRC(vec_exp_strong)},
    {GR_METHOD_DIV,             (gr_funcptr) GRC(div_strong)},
    {GR_METHOD_INV,             (gr_funcptr) GRC(inv_strong)},
    {GR_METHOD_SQRT,            (gr_funcptr) GRC(sqrt_strong)},
    {GR_METHOD_RSQRT,           (gr_funcptr) GRC(rsqrt_strong)},
    {GR_METHOD_ABS,             (gr_funcptr) GRC(abs_strong)},
    {GR_METHOD_ABS2,            (gr_funcptr) GRC(abs2_strong)},
    {GR_METHOD_VEC_APPLY_2X2,   (gr_funcptr) GRC(vec_apply_2x2_strong)},
    {0,                         (gr_funcptr) NULL},
};

/* (extends the generic table _dfloat_complex_gr_methods_input) */
gr_method_tab_input DCV(gr_methods_input)[] =
{
    {GR_METHOD_CTX_IS_RING,     (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,  (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FIELD,            (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_UNIQUE_FACTORIZATION_DOMAIN, (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FINITE,           (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC, (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_ALGEBRAICALLY_CLOSED, (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_ORDERED_RING,     (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_APPROX_COMMUTATIVE_RING, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_EXACT,            (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_CANONICAL,        (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_HAS_REAL_PREC,       (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_SET_REAL_PREC,       (gr_funcptr) GRF(ctx_set_real_prec)},
    {GR_METHOD_CTX_GET_REAL_PREC,       (gr_funcptr) GRF(ctx_get_real_prec)},
    {GR_METHOD_INIT,            (gr_funcptr) GRC(init)},
    {GR_METHOD_CLEAR,           (gr_funcptr) GRC(clear)},
    {GR_METHOD_SWAP,            (gr_funcptr) GRC(swap)},
    {GR_METHOD_SET_SHALLOW,     (gr_funcptr) GRC(set)},
    {GR_METHOD_ZERO,            (gr_funcptr) GRC(zero)},
    {GR_METHOD_ONE,             (gr_funcptr) GRC(one)},
    {GR_METHOD_NEG_ONE,         (gr_funcptr) GRC(neg_one)},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) GRC(is_zero)},
    {GR_METHOD_IS_ONE,          (gr_funcptr) GRC(is_one)},
    {GR_METHOD_IS_NEG_ONE,      (gr_funcptr) GRC(is_neg_one)},
    {GR_METHOD_EQUAL,           (gr_funcptr) GRC(equal)},
    {GR_METHOD_SET,             (gr_funcptr) GRC(set)},
    {GR_METHOD_NEG,             (gr_funcptr) GRC(neg)},
    {GR_METHOD_ADD,             (gr_funcptr) GRC(add)},
    {GR_METHOD_SUB,             (gr_funcptr) GRC(sub)},
    {GR_METHOD_MUL,             (gr_funcptr) GRC(mul)},
    {GR_METHOD_MUL_TWO,         (gr_funcptr) GRC(mul_two)},
    {GR_METHOD_SQR,             (gr_funcptr) GRC(sqr)},
    {GR_METHOD_MUL_2EXP_SI,     (gr_funcptr) GRC(mul_2exp_si)},
    {GR_METHOD_CONJ,            (gr_funcptr) GRC(conj)},
    {GR_METHOD_RE,              (gr_funcptr) GRC(re)},
    {GR_METHOD_IM,              (gr_funcptr) GRC(im)},
    {GR_METHOD_I,               (gr_funcptr) GRC(i)},
    {GR_METHOD_DIV,             (gr_funcptr) GRC(div)},
    {GR_METHOD_INV,             (gr_funcptr) GRC(inv)},
    {GR_METHOD_SQRT,            (gr_funcptr) GRC(sqrt)},
    {GR_METHOD_RSQRT,           (gr_funcptr) GRC(rsqrt)},
    {GR_METHOD_ABS,             (gr_funcptr) GRC(abs)},
    {GR_METHOD_VEC_INIT,        (gr_funcptr) GRC(vec_init)},
    {GR_METHOD_VEC_CLEAR,       (gr_funcptr) GRC(vec_clear)},
    {GR_METHOD_VEC_SET,         (gr_funcptr) GRC(vec_set)},
    {GR_METHOD_VEC_ZERO,        (gr_funcptr) GRC(vec_zero)},
    {GR_METHOD_VEC_ADD,         (gr_funcptr) GRC(vec_add)},
    {GR_METHOD_VEC_SUB,         (gr_funcptr) GRC(vec_sub)},
    {GR_METHOD_VEC_MUL,         (gr_funcptr) GRC(vec_mul)},
    {GR_METHOD_VEC_MUL_SCALAR,  (gr_funcptr) GRC(vec_mul_scalar)},
    {GR_METHOD_VEC_ADDMUL_SCALAR, (gr_funcptr) GRC(vec_addmul_scalar)},
    {GR_METHOD_VEC_SUBMUL_SCALAR, (gr_funcptr) GRC(vec_submul_scalar)},
    {GR_METHOD_VEC_DOT,         (gr_funcptr) GRC(vec_dot)},
    {GR_METHOD_VEC_DOT_REV,     (gr_funcptr) GRC(vec_dot_rev)},
    {GR_METHOD_VEC_EXP,         (gr_funcptr) GRC(vec_exp)},
    {GR_METHOD_VEC_GATHER,      (gr_funcptr) GRC(vec_gather)},
    {GR_METHOD_VEC_SCATTER,     (gr_funcptr) GRC(vec_scatter)},
    {GR_METHOD_ABS2,            (gr_funcptr) GRC(abs2)},
    {GR_METHOD_CMP,             (gr_funcptr) GRC(cmp)},
    {GR_METHOD_CMPABS,          (gr_funcptr) GRC(cmpabs)},
    {GR_METHOD_VEC_APPLY_2X2,   (gr_funcptr) GRC(vec_apply_2x2)},
    {GR_METHOD_MAT_MUL,         (gr_funcptr) GRC(mat_mul)},
    {GR_METHOD_MAT_LU,          (gr_funcptr) GRC(mat_lu)},
    {GR_METHOD_MAT_NONSINGULAR_SOLVE_TRIL, (gr_funcptr) DFV(gr_mat_nonsingular_solve_tril)},
    {GR_METHOD_MAT_NONSINGULAR_SOLVE_TRIU, (gr_funcptr) DFV(gr_mat_nonsingular_solve_triu)},
    {GR_METHOD_MAT_DET,         (gr_funcptr) gr_mat_det_generic_field},
    {GR_METHOD_MAT_FIND_NONZERO_PIVOT, (gr_funcptr) gr_mat_find_nonzero_pivot_large_abs},
    {0,                         (gr_funcptr) NULL},
};

/* ---- balls ---- */

static void GRCB(init)(DCB(ptr) x, gr_ctx_t ctx) { DCB(_zero_inline)(x); }
static void GRCB(clear)(DCB(ptr) x, gr_ctx_t ctx) { }
static void GRCB(swap)(DCB(ptr) x, DCB(ptr) y, gr_ctx_t ctx) { DCB(struct) t = *x; *x = *y; *y = t; }
static int GRCB(zero)(DCB(ptr) res, gr_ctx_t ctx) { DCB(_zero_inline)(res); return GR_SUCCESS; }
static int GRCB(one)(DCB(ptr) res, gr_ctx_t ctx) { DCB(_one_inline)(res); return GR_SUCCESS; }
static int GRCB(neg_one)(DCB(ptr) res, gr_ctx_t ctx) { DCB(_set_d_inline)(res, -1.0); return GR_SUCCESS; }
static int GRCB(i)(DCB(ptr) res, gr_ctx_t ctx) { DCB(_onei_inline)(res); return GR_SUCCESS; }
static truth_t GRCB(is_zero)(DCB(srcptr) x, gr_ctx_t ctx) { return DC(_truth_and)(GRB(is_zero)(&x->re, ctx), GRB(is_zero)(&x->im, ctx)); }
static truth_t GRCB(is_one)(DCB(srcptr) x, gr_ctx_t ctx) { return DC(_truth_and)(GRB(is_one)(&x->re, ctx), GRB(is_zero)(&x->im, ctx)); }
static truth_t GRCB(is_neg_one)(DCB(srcptr) x, gr_ctx_t ctx) { return DC(_truth_and)(GRB(is_neg_one)(&x->re, ctx), GRB(is_zero)(&x->im, ctx)); }
static truth_t GRCB(equal)(DCB(srcptr) x, DCB(srcptr) y, gr_ctx_t ctx) { return DC(_truth_and)(GRB(equal)(&x->re, &y->re, ctx), GRB(equal)(&x->im, &y->im, ctx)); }
/* (normalizes the whole line, as the real gr set) */
static int GRCB(set)(DCB(ptr) res, DCB(srcptr) x, gr_ctx_t ctx)
{
    GRB(set)(&res->re, &x->re, ctx);
    GRB(set)(&res->im, &x->im, ctx);
    return GR_SUCCESS;
}
static int GRCB(neg)(DCB(ptr) res, DCB(srcptr) x, gr_ctx_t ctx) { DCB(_neg_inline)(res, x); return GR_SUCCESS; }
static int GRCB(add)(DCB(ptr) res, DCB(srcptr) x, DCB(srcptr) y, gr_ctx_t ctx) { DCB(_add_inline)(res, x, y); return GR_SUCCESS; }
static int GRCB(sub)(DCB(ptr) res, DCB(srcptr) x, DCB(srcptr) y, gr_ctx_t ctx) { DCB(_sub_inline)(res, x, y); return GR_SUCCESS; }
static int GRCB(mul)(DCB(ptr) res, DCB(srcptr) x, DCB(srcptr) y, gr_ctx_t ctx) { DCB(_mul_m)(res, x, y, DFLOAT_MODE_AUTO); return GR_SUCCESS; }
static int GRCB(sqr)(DCB(ptr) res, DCB(srcptr) x, gr_ctx_t ctx) { DCB(_sqr_m)(res, x, DFLOAT_MODE_AUTO); return GR_SUCCESS; }
static int GRCB(mul_two)(DCB(ptr) res, DCB(srcptr) x, gr_ctx_t ctx) { DCB(_mul_2exp_si_inline)(res, x, 1); return GR_SUCCESS; }
static int GRCB(mul_2exp_si)(DCB(ptr) res, DCB(srcptr) x, slong e, gr_ctx_t ctx) { DCB(_mul_2exp_si_inline)(res, x, e); return GR_SUCCESS; }
static int GRCB(conj)(DCB(ptr) res, DCB(srcptr) x, gr_ctx_t ctx) { DCB(_conj_inline)(res, x); return GR_SUCCESS; }
static int GRCB(re)(DCB(ptr) res, DCB(srcptr) x, gr_ctx_t ctx) { GRB(set)(&res->re, &x->re, ctx); DB(_zero_inline)(&res->im); return GR_SUCCESS; }
static int GRCB(im)(DCB(ptr) res, DCB(srcptr) x, gr_ctx_t ctx) { GRB(set)(&res->re, &x->im, ctx); DB(_zero_inline)(&res->im); return GR_SUCCESS; }

static void GRCB(vec_init)(DCB(ptr) vec, slong len, gr_ctx_t ctx)
{
    slong i;
    for (i = 0; i < len; i++)
        DCB(_zero_inline)(vec + i);
}
static void GRCB(vec_clear)(DCB(ptr) vec, slong len, gr_ctx_t ctx) { }
static int GRCB(vec_set)(DCB(ptr) res, DCB(srcptr) x, slong len, gr_ctx_t ctx)
{
    if (res != x)
        memcpy(res, x, len * sizeof(DCB(struct)));
    return GR_SUCCESS;
}
static int GRCB(vec_zero)(DCB(ptr) vec, slong len, gr_ctx_t ctx) { GRCB(vec_init)(vec, len, ctx); return GR_SUCCESS; }
static int GRCB(vec_add)(DCB(ptr) res, DCB(srcptr) x, DCB(srcptr) y, slong len, gr_ctx_t ctx) { DCBV(vec_add)(res, x, y, len); return GR_SUCCESS; }
static int GRCB(vec_sub)(DCB(ptr) res, DCB(srcptr) x, DCB(srcptr) y, slong len, gr_ctx_t ctx) { DCBV(vec_sub)(res, x, y, len); return GR_SUCCESS; }
static int GRCB(vec_mul)(DCB(ptr) res, DCB(srcptr) x, DCB(srcptr) y, slong len, gr_ctx_t ctx) { DCBV(_vec_mul_m)(res, x, y, len, DFLOAT_MODE_AUTO); return GR_SUCCESS; }
static int GRCB(vec_mul_scalar)(DCB(ptr) res, DCB(srcptr) x, slong len, DCB(srcptr) c, gr_ctx_t ctx) { DCBV(_vec_mul_scalar_m)(res, x, len, c, DFLOAT_MODE_AUTO); return GR_SUCCESS; }
static int GRCB(vec_addmul_scalar)(DCB(ptr) res, DCB(srcptr) x, slong len, DCB(srcptr) c, gr_ctx_t ctx) { DCBV(_vec_addmul_scalar_m)(res, x, len, c, 0, DFLOAT_MODE_AUTO); return GR_SUCCESS; }
static int GRCB(vec_submul_scalar)(DCB(ptr) res, DCB(srcptr) x, slong len, DCB(srcptr) c, gr_ctx_t ctx) { DCBV(_vec_addmul_scalar_m)(res, x, len, c, 1, DFLOAT_MODE_AUTO); return GR_SUCCESS; }
static int GRCB(vec_dot)(DCB(ptr) res, DCB(srcptr) initial, int subtract, DCB(srcptr) x, DCB(srcptr) y, slong len, gr_ctx_t ctx) { DCBV(_vec_dot)(res, initial, subtract, x, y, len, 0, DFLOAT_MODE_AUTO); return GR_SUCCESS; }
static int GRCB(vec_dot_rev)(DCB(ptr) res, DCB(srcptr) initial, int subtract, DCB(srcptr) x, DCB(srcptr) y, slong len, gr_ctx_t ctx) { DCBV(_vec_dot)(res, initial, subtract, x, y, len, 1, DFLOAT_MODE_AUTO); return GR_SUCCESS; }
static int GRCB(vec_exp)(DCB(ptr) res, DCB(srcptr) x, slong len, gr_ctx_t ctx) { DCBV(vec_exp)(res, x, len); return GR_SUCCESS; }
static int GRCB(vec_gather)(DCB(ptr) res, DCB(srcptr) vec, const slong * idx, slong len, gr_ctx_t ctx)
{
    slong i;
    for (i = 0; i < len; i++)
        res[i] = vec[idx[i]];
    return GR_SUCCESS;
}
static int GRCB(vec_scatter)(DCB(ptr) vec, const slong * idx, DCB(srcptr) src, slong len, gr_ctx_t ctx)
{
    slong i;
    for (i = 0; i < len; i++)
        vec[idx[i]] = src[i];
    return GR_SUCCESS;
}

/* fast variants (DFLOAT_FAST) */
static int GRCB(mul_fast)(DCB(ptr) res, DCB(srcptr) x, DCB(srcptr) y, gr_ctx_t ctx) { DCB(_mul_m)(res, x, y, DFLOAT_MODE_FAST); return GR_SUCCESS; }
static int GRCB(sqr_fast)(DCB(ptr) res, DCB(srcptr) x, gr_ctx_t ctx) { DCB(_sqr_m)(res, x, DFLOAT_MODE_FAST); return GR_SUCCESS; }
static int GRCB(vec_mul_fast)(DCB(ptr) res, DCB(srcptr) x, DCB(srcptr) y, slong len, gr_ctx_t ctx) { DCBV(_vec_mul_m)(res, x, y, len, DFLOAT_MODE_FAST); return GR_SUCCESS; }
static int GRCB(vec_mul_scalar_fast)(DCB(ptr) res, DCB(srcptr) x, slong len, DCB(srcptr) c, gr_ctx_t ctx) { DCBV(_vec_mul_scalar_m)(res, x, len, c, DFLOAT_MODE_FAST); return GR_SUCCESS; }
static int GRCB(vec_addmul_scalar_fast)(DCB(ptr) res, DCB(srcptr) x, slong len, DCB(srcptr) c, gr_ctx_t ctx) { DCBV(_vec_addmul_scalar_m)(res, x, len, c, 0, DFLOAT_MODE_FAST); return GR_SUCCESS; }
static int GRCB(vec_submul_scalar_fast)(DCB(ptr) res, DCB(srcptr) x, slong len, DCB(srcptr) c, gr_ctx_t ctx) { DCBV(_vec_addmul_scalar_m)(res, x, len, c, 1, DFLOAT_MODE_FAST); return GR_SUCCESS; }
static int GRCB(vec_dot_fast)(DCB(ptr) res, DCB(srcptr) initial, int subtract, DCB(srcptr) x, DCB(srcptr) y, slong len, gr_ctx_t ctx) { DCBV(_vec_dot)(res, initial, subtract, x, y, len, 0, DFLOAT_MODE_FAST); return GR_SUCCESS; }
static int GRCB(vec_dot_rev_fast)(DCB(ptr) res, DCB(srcptr) initial, int subtract, DCB(srcptr) x, DCB(srcptr) y, slong len, gr_ctx_t ctx) { DCBV(_vec_dot)(res, initial, subtract, x, y, len, 1, DFLOAT_MODE_FAST); return GR_SUCCESS; }

gr_method_tab_input DCBV(gr_fast_methods_input)[] =
{
    {GR_METHOD_MUL,             (gr_funcptr) GRCB(mul_fast)},
    {GR_METHOD_SQR,             (gr_funcptr) GRCB(sqr_fast)},
    {GR_METHOD_VEC_MUL,         (gr_funcptr) GRCB(vec_mul_fast)},
    {GR_METHOD_VEC_MUL_SCALAR,  (gr_funcptr) GRCB(vec_mul_scalar_fast)},
    {GR_METHOD_VEC_ADDMUL_SCALAR, (gr_funcptr) GRCB(vec_addmul_scalar_fast)},
    {GR_METHOD_VEC_SUBMUL_SCALAR, (gr_funcptr) GRCB(vec_submul_scalar_fast)},
    {GR_METHOD_VEC_DOT,         (gr_funcptr) GRCB(vec_dot_fast)},
    {GR_METHOD_VEC_DOT_REV,     (gr_funcptr) GRCB(vec_dot_rev_fast)},
    {0,                         (gr_funcptr) NULL},
};

/* strong variants */
static DFLOAT_NOINLINE void DCB(_canon)(DCB(ptr) x) { DF(_canonicalise_inplace)(x->re.d); DF(_canonicalise_inplace)(x->im.d); }
#define DCB_MODE(ctx) (DFLOAT_CTX_FAST(ctx) ? DFLOAT_MODE_FAST : DFLOAT_MODE_AUTO)
static DC_COLD int GRCB(add_strong)(DCB(ptr) res, DCB(srcptr) x, DCB(srcptr) y, gr_ctx_t ctx) { DCB(_add_inline)(res, x, y); DCB(_canon)(res); return GR_SUCCESS; }
static DC_COLD int GRCB(sub_strong)(DCB(ptr) res, DCB(srcptr) x, DCB(srcptr) y, gr_ctx_t ctx) { DCB(_sub_inline)(res, x, y); DCB(_canon)(res); return GR_SUCCESS; }
static DC_COLD int GRCB(mul_strong)(DCB(ptr) res, DCB(srcptr) x, DCB(srcptr) y, gr_ctx_t ctx) { DCB(_mul_m)(res, x, y, DCB_MODE(ctx)); DCB(_canon)(res); return GR_SUCCESS; }
static DC_COLD int GRCB(sqr_strong)(DCB(ptr) res, DCB(srcptr) x, gr_ctx_t ctx) { DCB(_sqr_m)(res, x, DCB_MODE(ctx)); DCB(_canon)(res); return GR_SUCCESS; }
static DC_COLD int GRCB(vec_add_strong)(DCB(ptr) res, DCB(srcptr) x, DCB(srcptr) y, slong len, gr_ctx_t ctx)
{ slong i; DCBV(vec_add)(res, x, y, len); for (i = 0; i < len; i++) DCB(_canon)(res + i); return GR_SUCCESS; }
static DC_COLD int GRCB(vec_sub_strong)(DCB(ptr) res, DCB(srcptr) x, DCB(srcptr) y, slong len, gr_ctx_t ctx)
{ slong i; DCBV(vec_sub)(res, x, y, len); for (i = 0; i < len; i++) DCB(_canon)(res + i); return GR_SUCCESS; }
static DC_COLD int GRCB(vec_mul_strong)(DCB(ptr) res, DCB(srcptr) x, DCB(srcptr) y, slong len, gr_ctx_t ctx)
{ slong i; DCBV(_vec_mul_m)(res, x, y, len, DCB_MODE(ctx)); for (i = 0; i < len; i++) DCB(_canon)(res + i); return GR_SUCCESS; }
static DC_COLD int GRCB(vec_mul_scalar_strong)(DCB(ptr) res, DCB(srcptr) x, slong len, DCB(srcptr) c, gr_ctx_t ctx)
{ slong i; DCBV(_vec_mul_scalar_m)(res, x, len, c, DCB_MODE(ctx)); for (i = 0; i < len; i++) DCB(_canon)(res + i); return GR_SUCCESS; }
static DC_COLD int GRCB(vec_addmul_scalar_strong)(DCB(ptr) res, DCB(srcptr) x, slong len, DCB(srcptr) c, gr_ctx_t ctx)
{ slong i; DCBV(_vec_addmul_scalar_m)(res, x, len, c, 0, DCB_MODE(ctx)); for (i = 0; i < len; i++) DCB(_canon)(res + i); return GR_SUCCESS; }
static DC_COLD int GRCB(vec_submul_scalar_strong)(DCB(ptr) res, DCB(srcptr) x, slong len, DCB(srcptr) c, gr_ctx_t ctx)
{ slong i; DCBV(_vec_addmul_scalar_m)(res, x, len, c, 1, DCB_MODE(ctx)); for (i = 0; i < len; i++) DCB(_canon)(res + i); return GR_SUCCESS; }
static DC_COLD int GRCB(vec_dot_strong)(DCB(ptr) res, DCB(srcptr) initial, int subtract, DCB(srcptr) x, DCB(srcptr) y, slong len, gr_ctx_t ctx)
{ DCBV(_vec_dot)(res, initial, subtract, x, y, len, 0, DCB_MODE(ctx)); DCB(_canon)(res); return GR_SUCCESS; }
static DC_COLD int GRCB(vec_dot_rev_strong)(DCB(ptr) res, DCB(srcptr) initial, int subtract, DCB(srcptr) x, DCB(srcptr) y, slong len, gr_ctx_t ctx)
{ DCBV(_vec_dot)(res, initial, subtract, x, y, len, 1, DCB_MODE(ctx)); DCB(_canon)(res); return GR_SUCCESS; }

static DC_COLD int GRCB(vec_exp_strong)(DCB(ptr) res, DCB(srcptr) x, slong len, gr_ctx_t ctx)
{ slong i; DCBV(vec_exp)(res, x, len); for (i = 0; i < len; i++) DCB(_canon)(res + i); return GR_SUCCESS; }

gr_method_tab_input DCBV(gr_strong_methods_input)[] =
{
    {GR_METHOD_ADD,             (gr_funcptr) GRCB(add_strong)},
    {GR_METHOD_SUB,             (gr_funcptr) GRCB(sub_strong)},
    {GR_METHOD_MUL,             (gr_funcptr) GRCB(mul_strong)},
    {GR_METHOD_SQR,             (gr_funcptr) GRCB(sqr_strong)},
    {GR_METHOD_VEC_ADD,         (gr_funcptr) GRCB(vec_add_strong)},
    {GR_METHOD_VEC_SUB,         (gr_funcptr) GRCB(vec_sub_strong)},
    {GR_METHOD_VEC_MUL,         (gr_funcptr) GRCB(vec_mul_strong)},
    {GR_METHOD_VEC_MUL_SCALAR,  (gr_funcptr) GRCB(vec_mul_scalar_strong)},
    {GR_METHOD_VEC_ADDMUL_SCALAR, (gr_funcptr) GRCB(vec_addmul_scalar_strong)},
    {GR_METHOD_VEC_SUBMUL_SCALAR, (gr_funcptr) GRCB(vec_submul_scalar_strong)},
    {GR_METHOD_VEC_DOT,         (gr_funcptr) GRCB(vec_dot_strong)},
    {GR_METHOD_VEC_DOT_REV,     (gr_funcptr) GRCB(vec_dot_rev_strong)},
    {GR_METHOD_VEC_EXP,         (gr_funcptr) GRCB(vec_exp_strong)},
    {0,                         (gr_funcptr) NULL},
};

/* (extends the generic table _dfloat_complex_gr_methods_input) */
gr_method_tab_input DCBV(gr_methods_input)[] =
{
    {GR_METHOD_CTX_IS_RING,     (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,  (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FIELD,            (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_UNIQUE_FACTORIZATION_DOMAIN, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FINITE,           (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC, (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_ALGEBRAICALLY_CLOSED, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_RATIONAL_VECTOR_SPACE, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_REAL_VECTOR_SPACE, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_COMPLEX_VECTOR_SPACE, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_ORDERED_RING,     (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_EXACT,            (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_CANONICAL,        (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_HAS_REAL_PREC,       (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_SET_REAL_PREC,       (gr_funcptr) GRF(ctx_set_real_prec)},
    {GR_METHOD_CTX_GET_REAL_PREC,       (gr_funcptr) GRF(ctx_get_real_prec)},
    {GR_METHOD_INIT,            (gr_funcptr) GRCB(init)},
    {GR_METHOD_CLEAR,           (gr_funcptr) GRCB(clear)},
    {GR_METHOD_SWAP,            (gr_funcptr) GRCB(swap)},
    {GR_METHOD_SET_SHALLOW,     (gr_funcptr) GRCB(set)},
    {GR_METHOD_ZERO,            (gr_funcptr) GRCB(zero)},
    {GR_METHOD_ONE,             (gr_funcptr) GRCB(one)},
    {GR_METHOD_NEG_ONE,         (gr_funcptr) GRCB(neg_one)},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) GRCB(is_zero)},
    {GR_METHOD_IS_ONE,          (gr_funcptr) GRCB(is_one)},
    {GR_METHOD_IS_NEG_ONE,      (gr_funcptr) GRCB(is_neg_one)},
    {GR_METHOD_EQUAL,           (gr_funcptr) GRCB(equal)},
    {GR_METHOD_SET,             (gr_funcptr) GRCB(set)},
    {GR_METHOD_NEG,             (gr_funcptr) GRCB(neg)},
    {GR_METHOD_ADD,             (gr_funcptr) GRCB(add)},
    {GR_METHOD_SUB,             (gr_funcptr) GRCB(sub)},
    {GR_METHOD_MUL,             (gr_funcptr) GRCB(mul)},
    {GR_METHOD_MUL_TWO,         (gr_funcptr) GRCB(mul_two)},
    {GR_METHOD_SQR,             (gr_funcptr) GRCB(sqr)},
    {GR_METHOD_MUL_2EXP_SI,     (gr_funcptr) GRCB(mul_2exp_si)},
    {GR_METHOD_CONJ,            (gr_funcptr) GRCB(conj)},
    {GR_METHOD_RE,              (gr_funcptr) GRCB(re)},
    {GR_METHOD_IM,              (gr_funcptr) GRCB(im)},
    {GR_METHOD_I,               (gr_funcptr) GRCB(i)},
    {GR_METHOD_VEC_INIT,        (gr_funcptr) GRCB(vec_init)},
    {GR_METHOD_VEC_CLEAR,       (gr_funcptr) GRCB(vec_clear)},
    {GR_METHOD_VEC_SET,         (gr_funcptr) GRCB(vec_set)},
    {GR_METHOD_VEC_ZERO,        (gr_funcptr) GRCB(vec_zero)},
    {GR_METHOD_VEC_ADD,         (gr_funcptr) GRCB(vec_add)},
    {GR_METHOD_VEC_SUB,         (gr_funcptr) GRCB(vec_sub)},
    {GR_METHOD_VEC_MUL,         (gr_funcptr) GRCB(vec_mul)},
    {GR_METHOD_VEC_MUL_SCALAR,  (gr_funcptr) GRCB(vec_mul_scalar)},
    {GR_METHOD_VEC_ADDMUL_SCALAR, (gr_funcptr) GRCB(vec_addmul_scalar)},
    {GR_METHOD_VEC_SUBMUL_SCALAR, (gr_funcptr) GRCB(vec_submul_scalar)},
    {GR_METHOD_VEC_DOT,         (gr_funcptr) GRCB(vec_dot)},
    {GR_METHOD_VEC_DOT_REV,     (gr_funcptr) GRCB(vec_dot_rev)},
    {GR_METHOD_VEC_EXP,         (gr_funcptr) GRCB(vec_exp)},
    {GR_METHOD_VEC_GATHER,      (gr_funcptr) GRCB(vec_gather)},
    {GR_METHOD_VEC_SCATTER,     (gr_funcptr) GRCB(vec_scatter)},
    {0,                         (gr_funcptr) NULL},
};
