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

void DC(zero)(DC(t) res) { DF(zero)(&res->re); DF(zero)(&res->im); }
void DC(one)(DC(t) res) { DF(one)(&res->re); DF(zero)(&res->im); }
void DC(onei)(DC(t) res) { DF(zero)(&res->re); DF(one)(&res->im); }
void DC(set)(DC(t) res, DC(srcptr) x) { *res = *x; }
void DC(set_d)(DC(t) res, double x) { DF(set_d)(&res->re, x); DF(zero)(&res->im); }
void DC(set_d_d)(DC(t) res, double x, double y) { DF(set_d)(&res->re, x); DF(set_d)(&res->im, y); }
void DC(set_re_im)(DC(t) res, DF(srcptr) x, DF(srcptr) y)
{
    DF(struct) t = *y;
    res->re = *x;
    res->im = t;
}

void DC(neg)(DC(t) res, DC(srcptr) x) { DF(neg)(&res->re, &x->re); DF(neg)(&res->im, &x->im); }
void DC(conj)(DC(t) res, DC(srcptr) x) { res->re = x->re; DF(neg)(&res->im, &x->im); }
void DC(add)(DC(t) res, DC(srcptr) x, DC(srcptr) y) { DF(add)(&res->re, &x->re, &y->re); DF(add)(&res->im, &x->im, &y->im); }
void DC(sub)(DC(t) res, DC(srcptr) x, DC(srcptr) y) { DF(sub)(&res->re, &x->re, &y->re); DF(sub)(&res->im, &x->im, &y->im); }

void DC(mul_re)(DC(t) res, DC(srcptr) x, DF(srcptr) y)
{
    DF(struct) c = *y;
    DF(mul)(&res->re, &x->re, &c);
    DF(mul)(&res->im, &x->im, &c);
}

/* (a + b i) i = -b + a i */
void DC(mul_onei)(DC(t) res, DC(srcptr) x)
{
    DF(struct) t = x->re;
    DF(neg)(&res->re, &x->im);
    res->im = t;
}

void DC(mul_2exp_si)(DC(t) res, DC(srcptr) x, slong e)
{
    DF(mul_2exp_si)(&res->re, &x->re, e);
    DF(mul_2exp_si)(&res->im, &x->im, e);
}

/* a real (or imaginary) factor takes two real products; otherwise
   (a c - b d) + (a d + b c) i */
void DC(mul)(DC(t) res, DC(srcptr) x, DC(srcptr) y)
{
    DF(t) p, q, r, s;
    if (DC_EXACT_ZERO(&y->im))
    {
        DF(mul)(r, &x->re, &y->re);
        DF(mul)(s, &x->im, &y->re);
    }
    else if (DC_EXACT_ZERO(&x->im))
    {
        DF(mul)(r, &y->re, &x->re);
        DF(mul)(s, &y->im, &x->re);
    }
    else
    {
        DF(mul)(p, &x->re, &y->re);
        DF(mul)(q, &x->im, &y->im);
        DF(mul)(r, &x->re, &y->im);
        DF(mul)(s, &x->im, &y->re);
        DF(add)(s, s, r);
        DF(sub)(r, p, q);
    }
    res->re = *r;
    res->im = *s;
}

/* (a + b)(a - b) + 2 a b i (the real part with a relative error) */
void DC(sqr)(DC(t) res, DC(srcptr) x)
{
    DF(t) p, q, r;
    if (DC_EXACT_ZERO(&x->im))
    {
        DF(sqr)(&res->re, &x->re);
        DF(zero)(&res->im);
        return;
    }
    DF(add)(p, &x->re, &x->im);
    DF(sub)(q, &x->re, &x->im);
    DF(mul)(r, &x->re, &x->im);
    DF(mul)(&res->re, p, q);
    DF(mul_2exp_si)(&res->im, r, 1);
}

int DC(is_zero)(DC(srcptr) x) { return DF(is_zero)(&x->re) && DF(is_zero)(&x->im); }
int DC(is_one)(DC(srcptr) x) { return DF(is_one)(&x->re) && DF(is_zero)(&x->im); }
int DC(is_real)(DC(srcptr) x) { return DF(is_zero)(&x->im); }
int DC(is_finite)(DC(srcptr) x) { return DF(is_finite)(&x->re) && DF(is_finite)(&x->im); }
int DC(equal)(DC(srcptr) x, DC(srcptr) y) { return DF(equal)(&x->re, &y->re) && DF(equal)(&x->im, &y->im); }


/* ---- through complex.c ---- */

void DC(div)(DC(t) res, DC(srcptr) x, DC(srcptr) y) { GR_IGNORE(_dfloat_complex_div(&DFV(ops), res->re.d, x->re.d, y->re.d, DFLOAT_MODE_AUTO)); }
void DC(inv)(DC(t) res, DC(srcptr) x) { GR_IGNORE(_dfloat_complex_inv(&DFV(ops), res->re.d, x->re.d, DFLOAT_MODE_AUTO)); }
void DC(sqrt)(DC(t) res, DC(srcptr) x) { _dfloat_complex_sqrt(&DFV(ops), res->re.d, x->re.d, DFLOAT_MODE_AUTO); }
void DC(rsqrt)(DC(t) res, DC(srcptr) x) { GR_IGNORE(_dfloat_complex_rsqrt(&DFV(ops), res->re.d, x->re.d, DFLOAT_MODE_AUTO)); }
void DC(abs)(DF(t) res, DC(srcptr) x) { _dfloat_complex_abs(&DFV(ops), res->d, x->re.d); }
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

/* the real products with a run-time error mode (see template.inc); the
   square in the fast mode as a product (whose bound is slightly larger),
   which saves an out-of-line copy of the square */
static DC_COLD void DB(_sqr_ool)(DB(t) res, DB(srcptr) x, int mode)
{
    if (mode == DFLOAT_MODE_FAST)
        DB(_mul_ool)(res, x, x, mode);
    else
        DB(sqr)(res, x);
}
#define DCB_RMUL(r, x, y) DB(_mul_ool)(r, x, y, mode)
#define DCB_RSQR(r, x) DB(_sqr_ool)(r, x, mode)

void DCB(zero)(DCB(t) res) { DB(zero)(&res->re); DB(zero)(&res->im); }
void DCB(one)(DCB(t) res) { DB(one)(&res->re); DB(zero)(&res->im); }
void DCB(onei)(DCB(t) res) { DB(zero)(&res->re); DB(one)(&res->im); }
void DCB(indeterminate)(DCB(t) res) { DB(indeterminate)(&res->re); DB(indeterminate)(&res->im); }
void DCB(set)(DCB(t) res, DCB(srcptr) x) { *res = *x; }
void DCB(set_d)(DCB(t) res, double x) { DB(set_d)(&res->re, x); DB(zero)(&res->im); }
void DCB(set_d_d)(DCB(t) res, double x, double y) { DB(set_d)(&res->re, x); DB(set_d)(&res->im, y); }
void DCB(set_re_im)(DCB(t) res, DB(srcptr) x, DB(srcptr) y)
{
    DB(struct) t = *y;
    res->re = *x;
    res->im = t;
}

void DCB(neg)(DCB(t) res, DCB(srcptr) x) { DB(neg)(&res->re, &x->re); DB(neg)(&res->im, &x->im); }
void DCB(conj)(DCB(t) res, DCB(srcptr) x)
{
    res->re = x->re;
    if (FLINT_UNLIKELY(DB(_is_W)(&res->re)))
        DB(indeterminate)(&res->re);
    DB(neg)(&res->im, &x->im);
}
void DCB(add)(DCB(t) res, DCB(srcptr) x, DCB(srcptr) y) { DB(add)(&res->re, &x->re, &y->re); DB(add)(&res->im, &x->im, &y->im); }
void DCB(sub)(DCB(t) res, DCB(srcptr) x, DCB(srcptr) y) { DB(sub)(&res->re, &x->re, &y->re); DB(sub)(&res->im, &x->im, &y->im); }

void DCB(mul_re)(DCB(t) res, DCB(srcptr) x, DB(srcptr) y)
{
    DB(struct) c = *y;
    DB(mul)(&res->re, &x->re, &c);
    DB(mul)(&res->im, &x->im, &c);
}

void DCB(mul_onei)(DCB(t) res, DCB(srcptr) x)
{
    DB(struct) t = x->re;
    DB(neg)(&res->re, &x->im);
    res->im = t;
    if (FLINT_UNLIKELY(DB(_is_W)(&res->im)))
        DB(indeterminate)(&res->im);
}

void DCB(mul_2exp_si)(DCB(t) res, DCB(srcptr) x, slong e)
{
    DB(mul_2exp_si)(&res->re, &x->re, e);
    DB(mul_2exp_si)(&res->im, &x->im, e);
}

/* as the plain version; an exactly real or imaginary factor takes two
   real products (which is also tighter) */
static DFLOAT_NOINLINE void
DCB(_mul_m)(DCB(t) res, DCB(srcptr) x, DCB(srcptr) y, int mode)
{
    DB(t) p, q, r, s;
    if (DB(is_zero)(&y->im))
    {
        DCB_RMUL(r, &x->re, &y->re);
        DCB_RMUL(s, &x->im, &y->re);
    }
    else if (DB(is_zero)(&x->im))
    {
        DCB_RMUL(r, &y->re, &x->re);
        DCB_RMUL(s, &y->im, &x->re);
    }
    else if (DB(is_zero)(&y->re))
    {
        /* (a + b i) (d i) = -b d + a d i */
        DCB_RMUL(p, &x->im, &y->im);
        DCB_RMUL(s, &x->re, &y->im);
        DB(neg)(r, p);
    }
    else if (DB(is_zero)(&x->re))
    {
        DCB_RMUL(p, &y->im, &x->im);
        DCB_RMUL(s, &y->re, &x->im);
        DB(neg)(r, p);
    }
    else
    {
        DCB_RMUL(p, &x->re, &y->re);
        DCB_RMUL(q, &x->im, &y->im);
        DCB_RMUL(r, &x->re, &y->im);
        DCB_RMUL(s, &x->im, &y->re);
        DB(add)(s, s, r);
        DB(sub)(r, p, q);
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
    if (DB(is_zero)(&x->im))
    {
        DCB_RSQR(&res->re, &x->re);
        DB(zero)(&res->im);
        return;
    }
    if (DB(is_zero)(&x->re))
    {
        DCB_RSQR(p, &x->im);
        DB(neg)(&res->re, p);
        DB(zero)(&res->im);
        return;
    }
    DCB_RSQR(p, &x->re);
    DCB_RSQR(q, &x->im);
    DCB_RMUL(r, &x->re, &x->im);
    DB(sub)(&res->re, p, q);
    DB(mul_2exp_si)(&res->im, r, 1);
}

void DCB(sqr)(DCB(t) res, DCB(srcptr) x) { DCB(_sqr_m)(res, x, DFLOAT_MODE_AUTO); }

int DCB(is_zero)(DCB(srcptr) x) { return DB(is_zero)(&x->re) && DB(is_zero)(&x->im); }
int DCB(is_one)(DCB(srcptr) x)
{
    return x->re.rad == 0.0 && x->re.d[0] == 1.0 && DF(_tail_abs_sum)(x->re.d) == 0.0 && DB(is_zero)(&x->im);
}
int DCB(is_real)(DCB(srcptr) x) { return DB(is_zero)(&x->im); }
int DCB(is_exact)(DCB(srcptr) x) { return DB(is_exact)(&x->re) && DB(is_exact)(&x->im); }
int DCB(is_finite)(DCB(srcptr) x) { return DB(is_finite)(&x->re) && DB(is_finite)(&x->im); }
int DCB(contains_zero)(DCB(srcptr) x) { return DB(contains_zero)(&x->re) && DB(contains_zero)(&x->im); }
int DCB(contains)(DCB(srcptr) x, DCB(srcptr) y) { return DB(contains)(&x->re, &y->re) && DB(contains)(&x->im, &y->im); }
int DCB(overlaps)(DCB(srcptr) x, DCB(srcptr) y) { return DB(overlaps)(&x->re, &y->re) && DB(overlaps)(&x->im, &y->im); }
int DCB(equal)(DCB(srcptr) x, DCB(srcptr) y) { return DB(equal)(&x->re, &y->re) && DB(equal)(&x->im, &y->im); }


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
        DC(mul)(res + i, x + i, y + i);
}

DFLOAT_FUNC_NOINLINE void DCV(vec_mul_scalar)(DC(ptr) res, DC(srcptr) x, slong len, DC(srcptr) c)
{
    DC(struct) cc = *c;
    slong i;
    if (DC_EXACT_ZERO(&cc.im))
        DFV(vec_mul_scalar)((DF(ptr)) res, (DF(srcptr)) x, 2 * len, &cc.re);
    else
        for (i = 0; i < len; i++)
            DC(mul)(res + i, x + i, &cc);
}

static void
DCV(_vec_addmul_scalar)(DC(ptr) res, DC(srcptr) x, slong len, DC(srcptr) c, int sub)
{
    DC(struct) cc = *c;
    DC(t) t;
    slong i;
    if (DC_EXACT_ZERO(&cc.im))
    {
        if (sub)
            DFV(vec_submul_scalar)((DF(ptr)) res, (DF(srcptr)) x, 2 * len, &cc.re);
        else
            DFV(vec_addmul_scalar)((DF(ptr)) res, (DF(srcptr)) x, 2 * len, &cc.re);
        return;
    }
    for (i = 0; i < len; i++)
    {
        DC(mul)(t, x + i, &cc);
        if (sub)
            DC(sub)(res + i, res + i, t);
        else
            DC(add)(res + i, res + i, t);
    }
}

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
            DC(zero)(res);
        else
            DC(set)(res, initial);
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
    if (DB(is_zero)(&cc.im))
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
    if (DB(is_zero)(&cc.im))
    {
        DBV(_vec_addmul_scalar_m)((DB(ptr)) res, (DB(srcptr)) x, 2 * len, &cc.re, sub, mode);
        return;
    }
    for (i = 0; i < len; i++)
    {
        DCB(_mul_m)(t, x + i, &cc, mode);
        if (sub)
            DCB(sub)(res + i, res + i, t);
        else
            DCB(add)(res + i, res + i, t);
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
static DC_COLD void DF(_op_mul)(double * r, const double * x, const double * y, int mode) { DF(mul)((DF(ptr)) r, (DF(srcptr)) x, (DF(srcptr)) y); }
static DC_COLD void DF(_op_sqr)(double * r, const double * x, int mode) { DF(sqr)((DF(ptr)) r, (DF(srcptr)) x); }
static DC_COLD int DF(_op_div)(double * r, const double * x, const double * y, int mode) { DF(div)((DF(ptr)) r, (DF(srcptr)) x, (DF(srcptr)) y); return GR_SUCCESS; }
static DC_COLD int DF(_op_sqrt)(double * r, const double * x, int mode) { DF(sqrt)((DF(ptr)) r, (DF(srcptr)) x); return GR_SUCCESS; }
static DC_COLD int DF(_op_log)(double * r, const double * x) { DF(log)((DF(ptr)) r, (DF(srcptr)) x); return GR_SUCCESS; }
static DC_COLD int DF(_op_log1p)(double * r, const double * x) { DF(log1p)((DF(ptr)) r, (DF(srcptr)) x); return GR_SUCCESS; }
static DC_COLD int DF(_op_tan)(double * r, const double * x) { DF(tan)((DF(ptr)) r, (DF(srcptr)) x); return GR_SUCCESS; }
static DC_COLD int DF(_op_contains_zero)(const double * x) { return x[0] == 0.0; }
static DC_COLD int DF(_op_cmp_all_d)(const double * x, double c, int sign, int strict)
{
    DF(t) t;
    int r;
    if (DF(is_nan)((DF(srcptr)) x))
        return 0;
    DF(set_d)(t, c);
    r = DF(cmp)((DF(srcptr)) x, t);
    return (sign > 0) ? (strict ? r > 0 : r >= 0) : (strict ? r < 0 : r <= 0);
}
static DC_COLD int DF(_op_get_fmpz)(fmpz_t r, const double * x) { return GRF(get_fmpz)(r, (DF(srcptr)) x, NULL); }
static DC_COLD void DF(_op_cmul)(double * r, const double * x, const double * y, int mode) { DC(mul)((DC(ptr)) r, (DC(srcptr)) x, (DC(srcptr)) y); }
static DC_COLD void DF(_op_csqr)(double * r, const double * x, int mode) { DC(sqr)((DC(ptr)) r, (DC(srcptr)) x); }

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

static void GRC(init)(DC(ptr) x, gr_ctx_t ctx) { DC(zero)(x); }
static void GRC(clear)(DC(ptr) x, gr_ctx_t ctx) { }
static void GRC(swap)(DC(ptr) x, DC(ptr) y, gr_ctx_t ctx) { DC(struct) t = *x; *x = *y; *y = t; }
static int GRC(zero)(DC(ptr) res, gr_ctx_t ctx) { DC(zero)(res); return GR_SUCCESS; }
static int GRC(one)(DC(ptr) res, gr_ctx_t ctx) { DC(one)(res); return GR_SUCCESS; }
static int GRC(neg_one)(DC(ptr) res, gr_ctx_t ctx) { DC(set_d)(res, -1.0); return GR_SUCCESS; }
static int GRC(i)(DC(ptr) res, gr_ctx_t ctx) { DC(onei)(res); return GR_SUCCESS; }
static truth_t GRC(is_zero)(DC(srcptr) x, gr_ctx_t ctx) { return DC(_truth_and)(GRF(is_zero)(&x->re, ctx), GRF(is_zero)(&x->im, ctx)); }
static truth_t GRC(is_one)(DC(srcptr) x, gr_ctx_t ctx) { return DC(_truth_and)(GRF(is_one)(&x->re, ctx), GRF(is_zero)(&x->im, ctx)); }
static truth_t GRC(is_neg_one)(DC(srcptr) x, gr_ctx_t ctx) { return DC(_truth_and)(GRF(is_neg_one)(&x->re, ctx), GRF(is_zero)(&x->im, ctx)); }
static truth_t GRC(equal)(DC(srcptr) x, DC(srcptr) y, gr_ctx_t ctx) { return DC(_truth_and)(GRF(equal)(&x->re, &y->re, ctx), GRF(equal)(&x->im, &y->im, ctx)); }
static int GRC(set)(DC(ptr) res, DC(srcptr) x, gr_ctx_t ctx) { DC(set)(res, x); return GR_SUCCESS; }
static int GRC(neg)(DC(ptr) res, DC(srcptr) x, gr_ctx_t ctx) { DC(neg)(res, x); return GR_SUCCESS; }
static int GRC(add)(DC(ptr) res, DC(srcptr) x, DC(srcptr) y, gr_ctx_t ctx) { DC(add)(res, x, y); return GR_SUCCESS; }
static int GRC(sub)(DC(ptr) res, DC(srcptr) x, DC(srcptr) y, gr_ctx_t ctx) { DC(sub)(res, x, y); return GR_SUCCESS; }
static int GRC(mul)(DC(ptr) res, DC(srcptr) x, DC(srcptr) y, gr_ctx_t ctx) { DC(mul)(res, x, y); return GR_SUCCESS; }
static int GRC(sqr)(DC(ptr) res, DC(srcptr) x, gr_ctx_t ctx) { DC(sqr)(res, x); return GR_SUCCESS; }
static int GRC(mul_two)(DC(ptr) res, DC(srcptr) x, gr_ctx_t ctx) { DC(mul_2exp_si)(res, x, 1); return GR_SUCCESS; }
static int GRC(mul_2exp_si)(DC(ptr) res, DC(srcptr) x, slong e, gr_ctx_t ctx) { DC(mul_2exp_si)(res, x, e); return GR_SUCCESS; }
static int GRC(conj)(DC(ptr) res, DC(srcptr) x, gr_ctx_t ctx) { DC(conj)(res, x); return GR_SUCCESS; }
static int GRC(re)(DC(ptr) res, DC(srcptr) x, gr_ctx_t ctx) { res->re = x->re; DF(zero)(&res->im); return GR_SUCCESS; }
static int GRC(im)(DC(ptr) res, DC(srcptr) x, gr_ctx_t ctx) { res->re = x->im; DF(zero)(&res->im); return GR_SUCCESS; }

static void GRC(vec_init)(DC(ptr) vec, slong len, gr_ctx_t ctx)
{
    slong i;
    for (i = 0; i < len; i++)
        DC(zero)(vec + i);
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

gr_method_tab_input DCV(gr_strong_methods_input)[] =
{
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
    {0,                         (gr_funcptr) NULL},
};

/* ---- balls ---- */

static void GRCB(init)(DCB(ptr) x, gr_ctx_t ctx) { DCB(zero)(x); }
static void GRCB(clear)(DCB(ptr) x, gr_ctx_t ctx) { }
static void GRCB(swap)(DCB(ptr) x, DCB(ptr) y, gr_ctx_t ctx) { DCB(struct) t = *x; *x = *y; *y = t; }
static int GRCB(zero)(DCB(ptr) res, gr_ctx_t ctx) { DCB(zero)(res); return GR_SUCCESS; }
static int GRCB(one)(DCB(ptr) res, gr_ctx_t ctx) { DCB(one)(res); return GR_SUCCESS; }
static int GRCB(neg_one)(DCB(ptr) res, gr_ctx_t ctx) { DCB(set_d)(res, -1.0); return GR_SUCCESS; }
static int GRCB(i)(DCB(ptr) res, gr_ctx_t ctx) { DCB(onei)(res); return GR_SUCCESS; }
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
static int GRCB(neg)(DCB(ptr) res, DCB(srcptr) x, gr_ctx_t ctx) { DCB(neg)(res, x); return GR_SUCCESS; }
static int GRCB(add)(DCB(ptr) res, DCB(srcptr) x, DCB(srcptr) y, gr_ctx_t ctx) { DCB(add)(res, x, y); return GR_SUCCESS; }
static int GRCB(sub)(DCB(ptr) res, DCB(srcptr) x, DCB(srcptr) y, gr_ctx_t ctx) { DCB(sub)(res, x, y); return GR_SUCCESS; }
static int GRCB(mul)(DCB(ptr) res, DCB(srcptr) x, DCB(srcptr) y, gr_ctx_t ctx) { DCB(_mul_m)(res, x, y, DFLOAT_MODE_AUTO); return GR_SUCCESS; }
static int GRCB(sqr)(DCB(ptr) res, DCB(srcptr) x, gr_ctx_t ctx) { DCB(_sqr_m)(res, x, DFLOAT_MODE_AUTO); return GR_SUCCESS; }
static int GRCB(mul_two)(DCB(ptr) res, DCB(srcptr) x, gr_ctx_t ctx) { DCB(mul_2exp_si)(res, x, 1); return GR_SUCCESS; }
static int GRCB(mul_2exp_si)(DCB(ptr) res, DCB(srcptr) x, slong e, gr_ctx_t ctx) { DCB(mul_2exp_si)(res, x, e); return GR_SUCCESS; }
static int GRCB(conj)(DCB(ptr) res, DCB(srcptr) x, gr_ctx_t ctx) { DCB(conj)(res, x); return GR_SUCCESS; }
static int GRCB(re)(DCB(ptr) res, DCB(srcptr) x, gr_ctx_t ctx) { GRB(set)(&res->re, &x->re, ctx); DB(zero)(&res->im); return GR_SUCCESS; }
static int GRCB(im)(DCB(ptr) res, DCB(srcptr) x, gr_ctx_t ctx) { GRB(set)(&res->re, &x->im, ctx); DB(zero)(&res->im); return GR_SUCCESS; }

static void GRCB(vec_init)(DCB(ptr) vec, slong len, gr_ctx_t ctx)
{
    slong i;
    for (i = 0; i < len; i++)
        DCB(zero)(vec + i);
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
static DC_COLD int GRCB(add_strong)(DCB(ptr) res, DCB(srcptr) x, DCB(srcptr) y, gr_ctx_t ctx) { DCB(add)(res, x, y); DCB(_canon)(res); return GR_SUCCESS; }
static DC_COLD int GRCB(sub_strong)(DCB(ptr) res, DCB(srcptr) x, DCB(srcptr) y, gr_ctx_t ctx) { DCB(sub)(res, x, y); DCB(_canon)(res); return GR_SUCCESS; }
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
