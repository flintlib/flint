/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "padic_radix.h"
#include "impl.h"
#include "longlong.h"
#include "fmpz.h"
#include "gr.h"
#include "gr_mat.h"

/* ------------------------------------------------------------------------- */
/*    Internal helpers                                                       */
/* ------------------------------------------------------------------------- */

/*
    Whether |x| < p^k (i.e. x fits in k digits), determined without counting
    digits: we only look at the limb sitting on the digit boundary and compare
    it against the precomputed bound p^(k mod e) (radix->bpow[...]). When k is a
    multiple of e the boundary is a clean limb boundary and no compare is needed.
    Equivalent to radix_integer_size_digits(x) <= k but cheaper.
*/
static int
_radix_integer_fits_digits(const radix_integer_t x, slong k, const radix_t radix)
{
    slong un = FLINT_ABS(x->size);
    slong limbs = radix_div_exp(k, radix);
    slong rem = k - limbs * (slong) radix->exp;

    if (rem == 0)
        return un <= limbs;
    if (un <= limbs)
        return 1;
    if (un == limbs + 1)
        return x->d[limbs] < radix->bpow[rem];
    return 0;
}

/* x := x mod p^k (nonnegative residue) in place, k >= 1 */
static inline void
_radix_integer_mod_digits_inplace(radix_integer_t x, slong k, const radix_t radix)
{
    slong n = x->size, limbs, rem;

    if (n < 0)
    {
        radix_integer_mod_digits(x, x, k, radix);
        return;
    }

    limbs = radix_div_exp(k, radix);
    if (n <= limbs)
        return;

    rem = k - limbs * (slong) radix->exp;
    if (rem == 0)
    {
        n = limbs;
    }
    else
    {
        n = limbs + 1;
        x->d[limbs] = n_rem_precomp(x->d[limbs], radix->bpow[rem], radix->bpow_div + rem);
    }

    while (n > 0 && x->d[n - 1] == 0)
        n--;
    x->size = n;
}

/*
    Compare the nonnegative residues of two nonnegative radix_integers modulo
    p^k, without allocating. Returns 1 if equal, 0 otherwise. Assumes x and y
    have nonnegative size (which holds for all inexact units and for unsigned
    exact units).
*/
static int
_radix_integer_residue_eq_mod(const radix_integer_t x, const radix_integer_t y,
    slong k, const radix_t radix)
{
    slong e = radix->exp;
    slong limbs, rem, i, xs, ys, sig, loop_to;

    if (k <= 0)
        return 1;

    limbs = k / e;
    rem = k - limbs * e;
    xs = x->size;
    ys = y->size;

    /* both operands are zero in every limb at or above sig, so the residues can
       only differ within the first sig limbs; this keeps the work O(unit size)
       even when k (and hence limbs) is enormous. */
    sig = FLINT_MAX(xs, ys);
    loop_to = FLINT_MIN(limbs, sig);

    for (i = 0; i < loop_to; i++)
    {
        ulong xv = (i < xs) ? x->d[i] : 0;
        ulong yv = (i < ys) ? y->d[i] : 0;
        if (xv != yv)
            return 0;
    }

    if (limbs >= sig)
        return 1;                       /* all higher digits (incl. rem) are 0 == 0 */

    if (rem != 0)
    {
        ulong m = radix->bpow[rem];
        ulong xv = (limbs < xs) ? x->d[limbs] : 0;
        ulong yv = (limbs < ys) ? y->d[limbs] : 0;
        xv = n_rem_precomp(xv, m, radix->bpow_div + rem);
        yv = n_rem_precomp(yv, m, radix->bpow_div + rem);
        if (xv != yv)
            return 0;
    }

    return 1;
}

/*
    Whether the unit integers ux, uy (of either sign) are congruent modulo p^k.
    Bounded by the operands' significant size: a nonnegative unit pads its high
    digits with 0 and a negative unit pads with p-1, so once k passes the content
    a sign mismatch is automatically incongruent and matching signs reduce to the
    content. This avoids materialising a residue of ~k digits.
*/
static int
_radix_units_congruent_mod(const radix_integer_t ux, const radix_integer_t uy,
    slong k, const radix_t radix)
{
    slong e = radix->exp;
    slong content;
    int res;
    radix_integer_t rx, ry;

    if (k <= 0)
        return 1;

    if (ux->size >= 0 && uy->size >= 0)
        return _radix_integer_residue_eq_mod(ux, uy, k, radix);

    content = (FLINT_MAX(FLINT_ABS(ux->size), FLINT_ABS(uy->size)) + 1) * e;

    if (k > content)
    {
        if ((ux->size < 0) != (uy->size < 0))
            return 0;                   /* differ in the padding region within p^k */
        k = content;                    /* identical padding: only the content matters */
    }

    /* k is now bounded by 'content' (a few limbs); compare nonnegative residues */
    radix_integer_init(rx, radix);
    radix_integer_init(ry, radix);
    radix_integer_mod_digits(rx, ux, k, radix);
    radix_integer_mod_digits(ry, uy, k, radix);
    res = radix_integer_equal(rx, ry, radix);
    radix_integer_clear(rx, radix);
    radix_integer_clear(ry, radix);
    return res;
}

/*
    Test whether the nonnegative residue of a nonnegative radix_integer u
    modulo p^k equals 1, without allocating.
*/
static int
_radix_integer_residue_is_one_mod(const radix_integer_t u, slong k, const radix_t radix)
{
    slong e = radix->exp;
    slong limbs, rem, i, us;

    if (k <= 0)
        return 1;                 /* everything is congruent mod p^0 */

    limbs = k / e;
    rem = k - limbs * e;
    us = u->size;

    for (i = 0; i < limbs; i++)
    {
        ulong uv = (i < us) ? u->d[i] : 0;
        ulong ov = (i == 0) ? 1 : 0;
        if (uv != ov)
            return 0;
    }

    if (rem != 0)
    {
        ulong m = radix->bpow[rem];
        ulong uv = (limbs < us) ? u->d[limbs] : 0;
        ulong ov = (limbs == 0) ? (UWORD(1) % m) : 0;
        uv = n_rem_precomp(uv, m, radix->bpow_div + rem);
        if (uv != ov)
            return 0;
    }

    return 1;
}

/*
    Test whether the nonnegative residue of a nonnegative radix_integer u
    modulo p^k equals p^k - 1 (i.e. u == -1 mod p^k: every one of the low k
    p-adic digits is p - 1), without allocating.
*/
static int
_radix_integer_residue_is_neg_one_mod(const radix_integer_t u, slong k, const radix_t radix)
{
    slong e = radix->exp;
    slong limbs, rem, i, us;
    ulong full = LIMB_RADIX(radix) - 1;     /* a full limb of (p-1) digits */

    if (k <= 0)
        return 1;

    limbs = k / e;
    rem = k - limbs * e;
    us = u->size;

    for (i = 0; i < limbs; i++)
    {
        ulong uv = (i < us) ? u->d[i] : 0;
        if (uv != full)
            return 0;
    }

    if (rem != 0)
    {
        ulong m = radix->bpow[rem];
        ulong uv = (limbs < us) ? u->d[limbs] : 0;
        uv = n_rem_precomp(uv, m, radix->bpow_div + rem);
        if (uv != m - 1)
            return 0;
    }

    return 1;
}

/* Pull the p-adic valuation of the unit into the exponent v. */
static void
_padic_radix_canonicalise(padic_radix_t x, gr_ctx_t ctx)
{
    radix_struct * radix = PADIC_RADIX_CTX_RADIX(ctx);
    slong w;

    if (x->u.size == 0)
    {
        x->v = 0;
        return;
    }

    w = radix_integer_valuation_digits(&x->u, radix);

    if (w != 0)
    {
    {
        slong nl = radix_div_exp(w, radix), nd = w - nl * radix->exp, n;
        if (nl != 0)
            radix_integer_rshift_limbs(&x->u, &x->u, nl, radix);
        n = FLINT_ABS(x->u.size);
        if (nd != 0 && n != 0)
        {
            radix_rshift_digits(x->u.d, x->u.d, n, nd, radix);
            n -= (x->u.d[n - 1] == 0);
            x->u.size = (x->u.size > 0) ? n : -n;
        }
    }
        x->v += w;
    }
}

/*
    Target absolute precision for a value of valuation v under the context
    default precisions: min(prec_abs, v + prec_rel), with PREC_INF acting as
    +infinity.
*/
static slong
_padic_radix_target_prec(slong v, gr_ctx_t ctx)
{
    slong a = PADIC_RADIX_CTX_PREC_ABS(ctx);
    slong r = PADIC_RADIX_CTX_PREC_REL(ctx);
    slong t;

    if (r == PADIC_RADIX_PREC_INF)
    {
        t = PADIC_RADIX_PREC_INF;
    }
    else
    {
        /* r is finite and clamped to [0, ERR_MAX]; v is bounded in practice */
        t = v + r;
        if (t < 0 || t > PADIC_RADIX_ERR_MAX)  /* overflow / very large -> treat as inf-ish */
            t = (v >= 0) ? PADIC_RADIX_PREC_INF : t;
    }

    return FLINT_MIN(a, t);
}

/*
    Impose the context precision on a canonical element, truncating when
    required and recording the resulting absolute error. An exact element is
    preserved as exact when its value fits strictly within the target
    precision. Inexact results always use the nonnegative residue.
*/
static void
_padic_radix_reduce(padic_radix_t x, gr_ctx_t ctx)
{
    radix_struct * radix = PADIC_RADIX_CTX_RADIX(ctx);
    slong target, Neff, r;

    if (x->u.size == 0)
    {
        x->v = 0;
        if (x->N != PADIC_RADIX_EXACT)
            x->N = FLINT_MIN(x->N, PADIC_RADIX_ERR_MAX);
        return;
    }

    target = _padic_radix_target_prec(x->v, ctx);
    Neff = (x->N == PADIC_RADIX_EXACT) ? target : FLINT_MIN(x->N, target);

    if (Neff == PADIC_RADIX_PREC_INF)
        return;   /* stays exact */

    if (target == PADIC_RADIX_PREC_INF)
        Neff = _padic_radix_cap_prec(Neff, x->v, ctx);

    r = Neff - x->v;

    if (r <= 0)
    {
        radix_integer_zero(&x->u, radix);
        x->v = 0;
        x->N = FLINT_MIN(Neff, PADIC_RADIX_ERR_MAX);
        return;
    }

    /* exact value that fits: keep it exact (allows signed units) */
    if (x->N == PADIC_RADIX_EXACT && _radix_integer_fits_digits(&x->u, r, radix))
        return;

    /* otherwise the result is inexact: store the nonnegative residue */
    _radix_integer_mod_digits_inplace(&x->u, r, radix);
    if (x->u.size == 0)
        x->v = 0;
    x->N = FLINT_MIN(Neff, PADIC_RADIX_ERR_MAX);
}

/*
    Finish a conversion: a freshly built (possibly signed) canonical element
    with N == EXACT. Handles the unsigned-mode representability of negative
    values and then imposes the context precision.
*/
static int
_padic_radix_finish(padic_radix_t res, gr_ctx_t ctx)
{
    if (res->u.size < 0 && !PADIC_RADIX_CTX_SIGNED(ctx))
    {
        slong target = _padic_radix_target_prec(res->v, ctx);

        if (target == PADIC_RADIX_PREC_INF)
            return GR_UNABLE;   /* no finite residue: cannot represent exactly */

        res->N = target;        /* force inexact; reduce takes the nonneg residue */
    }

    _padic_radix_reduce(res, ctx);
    return GR_SUCCESS;
}


/* ------------------------------------------------------------------------- */
/*    Low-level helpers                                                      */
/* ------------------------------------------------------------------------- */

static inline ulong
_padic_radix_divrem2(ulong * r, ulong hi, ulong lo, nmod_t B)
{
    ulong q, rr;
    if (B.norm == 0)
    {
        udiv_qrnnd_preinv(q, rr, hi, lo, B.n, B.ninv);
        *r = rr;
    }
    else
    {
        udiv_qrnnd_preinv(q, rr, (hi << B.norm) | (lo >> (FLINT_BITS - B.norm)),
                          lo << B.norm, B.n << B.norm, B.ninv);
        *r = rr >> B.norm;
    }
    return q;
}

/* rd[0, n) = (x0 + x1 B)(y0 + y1 B) mod B^n, n in {1, 2}; normalised length */
static inline slong
_padic_radix_mullow_12(nn_ptr rd, ulong x0, ulong x1, ulong y0, ulong y1, slong n, nmod_t B)
{
    ulong hi, lo, r0, r1, q, t2, t1, t0;

    umul_ppmm(hi, lo, x0, y0);

    if (B.n == UWORD(1) << (FLINT_BITS - 1))
    {
        ulong M = B.n - 1;
        r0 = lo & M;
        rd[0] = r0;
        if (n == 1)
            return (r0 != 0);
        q = (hi << 1) | (lo >> (FLINT_BITS - 1));
        r1 = (q + x0 * y1 + x1 * y0) & M;
        rd[1] = r1;
        return (r1 != 0) ? 2 : (r0 != 0);
    }

    if (n == 1)
    {
        NMOD_RED2(r0, hi, lo, B);
        rd[0] = r0;
        return (r0 != 0);
    }

    q = _padic_radix_divrem2(&r0, hi, lo, B);
    t2 = 0; t1 = 0; t0 = q;
    umul_ppmm(hi, lo, x0, y1);
    add_sssaaaaaa(t2, t1, t0, t2, t1, t0, 0, hi, lo);
    umul_ppmm(hi, lo, x1, y0);
    add_sssaaaaaa(t2, t1, t0, t2, t1, t0, 0, hi, lo);
    NMOD_RED3(r1, t2, t1, t0, B);

    rd[0] = r0;
    rd[1] = r1;
    return (r1 != 0) ? 2 : (r0 != 0);
}

/* res->u = sign * (|x| |y| mod B^n), any aliasing */
static inline void
_padic_radix_unit_mullow(padic_radix_t res, const padic_radix_t x, const padic_radix_t y,
    slong n, const radix_struct * radix)
{
    slong xs = x->u.size, ys = y->u.size;
    slong xn = FLINT_ABS(xs), yn = FLINT_ABS(ys), rn;

    if (n <= 2)
    {
        ulong x0 = x->u.d[0], y0 = y->u.d[0];
        ulong x1 = (xn > 1) ? x->u.d[1] : 0, y1 = (yn > 1) ? y->u.d[1] : 0;
        nn_ptr rd = radix_integer_fit_limbs(&res->u, 2, radix);
        rn = _padic_radix_mullow_12(rd, x0, x1, y0, y1, n, radix->B);
        res->u.size = ((xs ^ ys) >= 0) ? rn : -rn;
    }
    else
    {
        radix_integer_mullow_limbs(&res->u, &x->u, &y->u, n, radix);
    }
}

/* whether the nonzero limb L is divisible by p */
static inline int
_padic_radix_divisible_p(ulong L, const radix_struct * radix)
{
    if (radix->b.n == 2)
        return !(L & 1);
    return n_divisible_odd_gm(L, radix->bpow_oddinv[1].a, radix->bpow_oddinv[1].b);
}

/* digit valuation of a nonzero limb */
static inline slong
_padic_radix_val1(ulong L, const radix_struct * radix)
{
    return (radix->b.n == 2) ? (slong) flint_ctz(L) : (slong) _radix_valuation_digits_1(L, radix);
}

/* ------------------------------------------------------------------------- */
/*    Memory management                                                      */
/* ------------------------------------------------------------------------- */

void
padic_radix_init(padic_radix_t res, gr_ctx_t ctx)
{
    radix_integer_init(&res->u, PADIC_RADIX_CTX_RADIX(ctx));
    res->v = 0;
    res->N = PADIC_RADIX_EXACT;
}

void
padic_radix_clear(padic_radix_t res, gr_ctx_t ctx)
{
    radix_integer_clear(&res->u, PADIC_RADIX_CTX_RADIX(ctx));
}

void
padic_radix_swap(padic_radix_t x, padic_radix_t y, gr_ctx_t ctx)
{
    FLINT_SWAP(padic_radix_struct, *x, *y);
}

void
padic_radix_set_shallow(padic_radix_t res, const padic_radix_t x, gr_ctx_t ctx)
{
    *res = *x;
}

/* ------------------------------------------------------------------------- */
/*    Error term access                                                      */
/* ------------------------------------------------------------------------- */

slong
padic_radix_get_error(const padic_radix_t x, gr_ctx_t ctx)
{
    return x->N;
}

truth_t
padic_radix_is_exact(const padic_radix_t x, gr_ctx_t ctx)
{
    return (x->N == PADIC_RADIX_EXACT) ? T_TRUE : T_FALSE;
}

/* ------------------------------------------------------------------------- */
/*    Assignments                                                            */
/* ------------------------------------------------------------------------- */

int
padic_radix_zero(padic_radix_t res, gr_ctx_t ctx)
{
    radix_integer_zero(&res->u, PADIC_RADIX_CTX_RADIX(ctx));
    res->v = 0;
    res->N = PADIC_RADIX_EXACT;
    return GR_SUCCESS;
}

int
padic_radix_one(padic_radix_t res, gr_ctx_t ctx)
{
    radix_integer_one(&res->u, PADIC_RADIX_CTX_RADIX(ctx));
    res->v = 0;
    res->N = PADIC_RADIX_EXACT;
    _padic_radix_reduce(res, ctx);   /* becomes 0 + O(p^N) if target precision <= 0 */
    return GR_SUCCESS;
}

/* Allows parsing 'O(p^n)' in strings */
static int
padic_radix_big_o_base_fmpz(padic_radix_t res, const padic_radix_t x, const fmpz_t exp, gr_ctx_t ctx)
{
    /* p^exp */
    if (x->v == 1 && radix_integer_is_one(&x->u, PADIC_RADIX_CTX_RADIX(ctx)) && res->N == PADIC_RADIX_EXACT)
    {
        if (fmpz_cmp_si(exp, -PADIC_RADIX_ERR_MAX) >= 0)
        {
            padic_radix_zero(res, ctx);
            res->N = fmpz_cmp_si(exp, PADIC_RADIX_ERR_MAX) <= 0 ? *exp : PADIC_RADIX_ERR_MAX;
            return GR_SUCCESS;
        }
    }

    return GR_UNABLE;
}

int
padic_radix_set(padic_radix_t res, const padic_radix_t x, gr_ctx_t ctx)
{
    if (res != x)
    {
        radix_integer_set(&res->u, &x->u, PADIC_RADIX_CTX_RADIX(ctx));
        res->v = x->v;
        res->N = x->N;
    }

    _padic_radix_reduce(res, ctx);
    return GR_SUCCESS;
}

/*
    exact_set_*: store the exact value, ignoring the context precision (the
    result keeps N == EXACT and is never truncated). The unit may be signed
    regardless of the SIGNED flag, since this is the true value. Fails with
    GR_UNABLE only if the valuation is too large to represent (never in
    practice, as it is bounded by the bit length of the argument).
*/
static int
_padic_radix_exact_set_finish(padic_radix_t res, gr_ctx_t ctx)
{
    _padic_radix_canonicalise(res, ctx);

    if (res->u.size != 0
        && (res->v > PADIC_RADIX_ERR_MAX || res->v < -PADIC_RADIX_ERR_MAX))
        return GR_UNABLE;

    return GR_SUCCESS;
}

int
padic_radix_exact_set_fmpz(padic_radix_t res, const fmpz_t c, gr_ctx_t ctx)
{
    radix_integer_set_fmpz(&res->u, c, PADIC_RADIX_CTX_RADIX(ctx));
    res->v = 0;
    res->N = PADIC_RADIX_EXACT;
    return _padic_radix_exact_set_finish(res, ctx);
}

int
padic_radix_exact_set_ui(padic_radix_t res, ulong c, gr_ctx_t ctx)
{
    radix_integer_set_ui(&res->u, c, PADIC_RADIX_CTX_RADIX(ctx));
    res->v = 0;
    res->N = PADIC_RADIX_EXACT;
    return _padic_radix_exact_set_finish(res, ctx);
}

int
padic_radix_exact_set_si(padic_radix_t res, slong c, gr_ctx_t ctx)
{
    radix_integer_set_si(&res->u, c, PADIC_RADIX_CTX_RADIX(ctx));
    res->v = 0;
    res->N = PADIC_RADIX_EXACT;
    return _padic_radix_exact_set_finish(res, ctx);
}

int
padic_radix_set_fmpz(padic_radix_t res, const fmpz_t c, gr_ctx_t ctx)
{
    int status = padic_radix_exact_set_fmpz(res, c, ctx);
    if (status != GR_SUCCESS)
        return status;
    return _padic_radix_finish(res, ctx);        /* round to context precision */
}

int
padic_radix_get_fmpz(fmpz_t res, const padic_radix_t x, gr_ctx_t ctx)
{
    radix_struct * radix = PADIC_RADIX_CTX_RADIX(ctx);

    if (x->N != PADIC_RADIX_EXACT)
        return GR_UNABLE;       /* not a definite integer */

    if (x->v < 0)
        return GR_DOMAIN;       /* genuine p-adic non-integer */

    if (x->u.size == 0)
    {
        fmpz_zero(res);
        return GR_SUCCESS;
    }

    if ((ulong) x->v > ctx->size_limit)
        return GR_UNABLE;

    radix_integer_get_fmpz(res, &x->u, radix);

    if (x->v != 0)
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_ui_pow_ui(t, GR_PADIC_RADIX_CTX(ctx)->p, x->v);
        fmpz_mul(res, res, t);
        fmpz_clear(t);
    }

    return GR_SUCCESS;
}

/* ------------------------------------------------------------------------- */
/*    Arithmetic                                                             */
/* ------------------------------------------------------------------------- */

int
padic_radix_neg(padic_radix_t res, const padic_radix_t x, gr_ctx_t ctx)
{
    radix_struct * radix = PADIC_RADIX_CTX_RADIX(ctx);

    /* inexact x with a 1-2 limb unit: -x mod p^k as B^L - X truncated to
       k = N - v digits (the unit's low limb is nonzero: no borrow chain) */
    if (x->N != PADIC_RADIX_EXACT && x->u.size > 0 && x->u.size <= 2)
    {
        slong e = radix->exp, k, L, rem, sz = x->u.size;
        slong target = _padic_radix_target_prec(x->v, ctx);

        if (x->N <= target && x->N > x->v)
        {
            ulong x0 = x->u.d[0], x1 = (sz > 1) ? x->u.d[1] : 0, s0, s1;
            nn_ptr d;

            k = x->N - x->v;
            L = radix_div_exp(k + e - 1, radix);

            if (L <= 2 && sz <= L)
            {
                s0 = radix->B.n - x0;
                s1 = (L == 2) ? radix->B.n - 1 - x1 : 0;
                rem = k - (L - 1) * e;
                if (rem < e)
                {
                    if (L == 1)
                        s0 = n_rem_precomp(s0, radix->bpow[rem], radix->bpow_div + rem);
                    else
                        s1 = n_rem_precomp(s1, radix->bpow[rem], radix->bpow_div + rem);
                }
                d = radix_integer_fit_limbs(&res->u, 2, radix);
                d[0] = s0;
                d[1] = s1;
                res->u.size = (s1 != 0) ? 2 : 1;
                res->v = x->v;
                res->N = x->N;
                return GR_SUCCESS;
            }
        }
    }

    if (res != x)

        padic_radix_set(res, x, ctx);

    radix_integer_neg(&res->u, &res->u, radix);

    if (res->N == PADIC_RADIX_EXACT)
        return _padic_radix_finish(res, ctx);

    /* inexact: re-normalise to the nonnegative residue at precision N */
    {
        slong r = res->N - res->v;

        if (r <= 0)
        {
            radix_integer_zero(&res->u, radix);
            res->v = 0;
        }
        else
        {
            _radix_integer_mod_digits_inplace(&res->u, r, radix);
            if (res->u.size == 0)
                res->v = 0;
        }
    }

    return GR_SUCCESS;
}

/*
    Finalize a result whose unit is ALREADY canonical (valuation 0, i.e. its
    lowest digit is a p-adic unit) -- as produced by multiplication, inversion,
    division, and addition/subtraction when the operand valuations differ. Skips
    the canonicalisation scan that _padic_radix_finalize performs.
*/
/*
    Exponents out of range: the result cannot be represented. GR_UNABLE is
    returned with res set to 0 + O(p^-ERR_MAX), so that callers which keep
    computing with the output (e.g. generic loops accumulating the status)
    never see exponents outside the representable range.
*/
static int
_padic_radix_unable(padic_radix_t res, gr_ctx_t ctx)
{
    radix_integer_zero(&res->u, PADIC_RADIX_CTX_RADIX(ctx));
    res->v = 0;
    res->N = -PADIC_RADIX_ERR_MAX;
    return GR_UNABLE;
}

static int
_padic_radix_finalize_assume_canonical(padic_radix_t res, gr_ctx_t ctx)
{
    radix_struct * radix = PADIC_RADIX_CTX_RADIX(ctx);

    /*
        Enforce the representability invariant |v| <= ERR_MAX and (N == EXACT
        or |N| <= ERR_MAX). Keeping both exponents within ERR_MAX = WORD_MAX/4
        guarantees that the sums formed by add/sub/mul (v_x + v_y, v_y + N_x,
        ...) never overflow a signed word.
    */
    if (res->N != PADIC_RADIX_EXACT)
    {
        if (res->N < -PADIC_RADIX_ERR_MAX)
            return _padic_radix_unable(res, ctx);   /* precision too low to represent */
        if (res->N > PADIC_RADIX_ERR_MAX)
            res->N = PADIC_RADIX_ERR_MAX;           /* clamp down: always sound */
    }

    if (res->u.size != 0)
    {
        if (res->v > PADIC_RADIX_ERR_MAX)
        {
            /* u * p^v with v > ERR_MAX is == 0 modulo p^ERR_MAX */
            slong newN = (res->N == PADIC_RADIX_EXACT)
                       ? PADIC_RADIX_ERR_MAX
                       : FLINT_MIN(res->N, PADIC_RADIX_ERR_MAX);
            radix_integer_zero(&res->u, radix);
            res->v = 0;
            res->N = newN;                          /* now finite */
        }
        else if (res->v < -PADIC_RADIX_ERR_MAX)
        {
            return _padic_radix_unable(res, ctx);   /* valuation too small to represent */
        }
    }

    if (res->N == PADIC_RADIX_EXACT)
        return _padic_radix_finish(res, ctx);

    _padic_radix_reduce(res, ctx);
    return GR_SUCCESS;
}

int
_padic_radix_finalize(padic_radix_t res, gr_ctx_t ctx)
{
    _padic_radix_canonicalise(res, ctx);
    return _padic_radix_finalize_assume_canonical(res, ctx);
}

/*
    Reference (unoptimized) add/sub: aligns the higher-valuation operand into a
    temporary and lets the finalizer reduce/truncate. Kept as a correctness oracle
    for the optimized _padic_radix_add_sub below (the test suite cross-checks the
    two). res = x + y, with sub != 0 selecting x - y.
*/
int
_padic_radix_add_sub_reference(padic_radix_t res, const padic_radix_t x,
    const padic_radix_t y, int sub, gr_ctx_t ctx)
{
    radix_struct * radix = PADIC_RADIX_CTX_RADIX(ctx);
    slong vx = x->v, vy = y->v, Nx = x->N, Ny = y->N;
    slong vmin = FLINT_MIN(vx, vy);
    slong vmax = FLINT_MAX(vx, vy);
    slong Nz, target, P, Nresult;
    radix_integer_t T;

    /* sum is known modulo p^min(Nx, Ny); EXACT == WORD_MAX acts as +inf */
    Nz = FLINT_MIN(Nx, Ny);

    /*
        Effective absolute precision we actually need to compute: the error
        bound capped by the context rounding precision. Anything at or above
        p^P is invisible in the result, so an operand whose valuation reaches P
        contributes nothing and must NOT be aligned by shifting (the shift
        v_hi - v_min could be astronomically large, e.g. ERR_MAX digits). P ==
        PREC_INF means a genuinely exact computation with no cap.
    */
    target = _padic_radix_target_prec(vmin, ctx);
    if (Nz == PADIC_RADIX_EXACT)
        P = target;
    else if (target == PADIC_RADIX_PREC_INF)
        P = _padic_radix_cap_prec(Nz, vmin, ctx);
    else
        P = FLINT_MIN(Nz, target);

    /*
        Zero operands: the result is the other operand (no shift needed). This
        must be handled before the horizon optimisation below, which assumes
        the lower-valuation operand is nonzero (otherwise the result valuation
        is not v_min and an operand could be dropped incorrectly).
    */
    if (x->u.size == 0 || y->u.size == 0)
    {
        if (x->u.size == 0 && y->u.size == 0)
        {
            radix_integer_zero(&res->u, radix);
            res->v = 0;
        }
        else if (x->u.size == 0)
        {
            if (sub)
                radix_integer_neg(&res->u, &y->u, radix);
            else
                radix_integer_set(&res->u, &y->u, radix);
            res->v = vy;
        }
        else
        {
            radix_integer_set(&res->u, &x->u, radix);
            res->v = vx;
        }
        res->N = Nz;
        /* result is the other operand's (already canonical) unit, or zero */
        return _padic_radix_finalize_assume_canonical(res, ctx);
    }

    /* whole result is below the precision horizon: it is 0 + O(p^P) */
    if (P != PADIC_RADIX_PREC_INF && vmin >= P)
    {
        radix_integer_zero(&res->u, radix);
        res->v = 0;
        res->N = P;                 /* known only modulo p^P (P <= Nz) */
        return _padic_radix_finalize_assume_canonical(res, ctx);   /* zero is canonical */
    }

    /*
        Equal valuations: the sum/difference can be written straight into the
        result unit (no temporary needed, even when res aliases x or y, since
        radix_integer_add/sub handle that). The leading digits may cancel, so
        this is the one case that still needs a canonicalising scan.
    */
    if (vx == vy)
    {
        if (sub)
            radix_integer_sub(&res->u, &x->u, &y->u, radix);
        else
            radix_integer_add(&res->u, &x->u, &y->u, radix);
        res->v = vmin;
        res->N = Nz;
        return _padic_radix_finalize(res, ctx);
    }

    radix_integer_init(T, radix);

    Nresult = Nz;

    if (P != PADIC_RADIX_PREC_INF && vmax >= P)
    {
        /* both operands nonzero, so the result has valuation v_min; the
           higher-valuation operand is negligible modulo p^P and is dropped.
           The dropped digits mean the result is now known only modulo p^P.
           No shift is needed, so this writes the surviving operand straight
           into res->u (radix_integer_set/neg handle res aliasing it). */
        Nresult = P;
        if (vx < vy)
        {
            radix_integer_set(&res->u, &x->u, radix);    /* result ~ x */
        }
        else
        {
            if (sub)
                radix_integer_neg(&res->u, &y->u, radix);    /* x - y, x negligible -> -y */
            else
                radix_integer_set(&res->u, &y->u, radix);
        }
    }
    else if (target == PADIC_RADIX_PREC_INF && !_padic_radix_exact_digits_ok(vmax - vmin, ctx))
    {
        radix_integer_clear(T, radix);
        return GR_UNABLE;
    }
    else if (vx < vy)
    {
        /* T holds the aligned higher-valuation operand; the combination is
           written directly into res->u (safe even when res aliases x or y:
           T is a distinct buffer and the integer add/sub handle in-place). */
        radix_integer_lshift_digits(T, &y->u, vmax - vmin, radix);
        if (sub)
            radix_integer_sub(&res->u, &x->u, T, radix);
        else
            radix_integer_add(&res->u, &x->u, T, radix);
    }
    else
    {
        radix_integer_lshift_digits(T, &x->u, vmax - vmin, radix);
        if (sub)
            radix_integer_sub(&res->u, T, &y->u, radix);
        else
            radix_integer_add(&res->u, T, &y->u, radix);
    }

    radix_integer_clear(T, radix);

    res->v = vmin;
    res->N = Nresult;

    /*
        With differing valuations the low digit of the result is the low digit
        of the lower-valuation operand (a unit), so the result is already
        canonical and needs no canonicalisation.
    */
    return _padic_radix_finalize_assume_canonical(res, ctx);
}

/*
    res = x + y (sub != 0 selects x - y), optimized.

    Identical in behaviour to _padic_radix_add_sub_reference, but the
    differing-valuation alignment is done in place with the fused
    radix_integer_addlsh / radix_integer_sublsh, eliminating the temporary that
    the reference allocates to hold the shifted higher-valuation operand.

    With v_x < v_y the lower operand is x, the result is x_unit +/- (y_unit << d);
    with v_x > v_y the lower operand is y, the result is (x_unit << d) +/- y_unit.
    Writing lo / hi for the lower / higher valuation operand and d = |v_x - v_y|,
    addition is always lo + (hi << d). Subtraction is lo - (hi << d) when x is the
    lower operand, and -(lo - (hi << d)) = (hi << d) - lo when y is, so in the
    latter case we negate the fused difference.
*/
static int
_padic_radix_add_sub(padic_radix_t res, const padic_radix_t x,
    const padic_radix_t y, int sub, gr_ctx_t ctx)
{
    radix_struct * radix = PADIC_RADIX_CTX_RADIX(ctx);
    slong vx = x->v, vy = y->v, Nx = x->N, Ny = y->N;
    slong vmin = FLINT_MIN(vx, vy);
    slong vmax = FLINT_MAX(vx, vy);
    slong Nz, target, P, Nresult;

    Nz = FLINT_MIN(Nx, Ny);

    target = _padic_radix_target_prec(vmin, ctx);
    if (Nz == PADIC_RADIX_EXACT)
        P = target;
    else if (target == PADIC_RADIX_PREC_INF)
        P = _padic_radix_cap_prec(Nz, vmin, ctx);
    else
        P = FLINT_MIN(Nz, target);

    /* zero operands: result is the other operand */
    if (x->u.size == 0 || y->u.size == 0)
    {
        if (x->u.size == 0 && y->u.size == 0)
        {
            radix_integer_zero(&res->u, radix);
            res->v = 0;
        }
        else if (x->u.size == 0)
        {
            if (sub)
                radix_integer_neg(&res->u, &y->u, radix);
            else
                radix_integer_set(&res->u, &y->u, radix);
            res->v = vy;
        }
        else
        {
            radix_integer_set(&res->u, &x->u, radix);
            res->v = vx;
        }
        res->N = Nz;
        return _padic_radix_finalize_assume_canonical(res, ctx);
    }

    /* whole result below the precision horizon: 0 + O(p^P) */
    if (P != PADIC_RADIX_PREC_INF && vmin >= P)
    {
        radix_integer_zero(&res->u, radix);
        res->v = 0;
        res->N = P;
        return _padic_radix_finalize_assume_canonical(res, ctx);
    }

    /* equal valuations: combine straight into res (cancellation possible). */
    if (vx == vy)
    {
        if (sub)
            radix_integer_sub(&res->u, &x->u, &y->u, radix);
        else
            radix_integer_add(&res->u, &x->u, &y->u, radix);
        res->v = vmin;
        res->N = Nz;

        /* Fast finish: inexact and the lowest digit did not cancel (one
           multiplicative divisibility test), so v = vmin and Neff = P. */
        if (Nz != PADIC_RADIX_EXACT && res->u.size != 0
            && !_padic_radix_divisible_p(res->u.d[0], radix))
        {
            _radix_integer_mod_digits_inplace(&res->u, P - vmin, radix);
            res->N = P;
            return GR_SUCCESS;
        }

        return _padic_radix_finalize(res, ctx);
    }

    Nresult = Nz;

    if (P != PADIC_RADIX_PREC_INF && vmax >= P)
    {
        /* higher-valuation operand is negligible modulo p^P; drop it (no shift). */
        Nresult = P;
        if (vx < vy)
        {
            radix_integer_set(&res->u, &x->u, radix);
        }
        else
        {
            if (sub)
                radix_integer_neg(&res->u, &y->u, radix);
            else
                radix_integer_set(&res->u, &y->u, radix);
        }
    }
    else
    {
        /* align in place via the fused shift, no temporary. lo / hi are the
           lower / higher valuation operands; d is the digit gap. */
        const padic_radix_struct * lo = (vx < vy) ? x : y;
        const padic_radix_struct * hi = (vx < vy) ? y : x;
        slong d = vmax - vmin;
        int negate = (sub && vx > vy);    /* x - y with y the lower operand */

        /* in a context without finite precision, a sum beyond the size
           limit of the context */
        if (target == PADIC_RADIX_PREC_INF && !_padic_radix_exact_digits_ok(d, ctx))
            return GR_UNABLE;

        if (sub)
            radix_integer_sublsh(&res->u, &lo->u, &hi->u, d, radix);
        else
            radix_integer_addlsh(&res->u, &lo->u, &hi->u, d, radix);

        if (negate)
            radix_integer_neg(&res->u, &res->u, radix);
    }

    res->v = vmin;
    res->N = Nresult;

    /* Fast finish: inexact, valuation vmin, so Neff = P. */
    if (Nresult != PADIC_RADIX_EXACT)
    {
        _radix_integer_mod_digits_inplace(&res->u, P - vmin, radix);
        res->N = P;
        return GR_SUCCESS;
    }

    /* low digit is the lower-valuation operand's (a unit): already canonical. */
    return _padic_radix_finalize_assume_canonical(res, ctx);
}

int
padic_radix_add(padic_radix_t res, const padic_radix_t x, const padic_radix_t y, gr_ctx_t ctx)
{
    return _padic_radix_add_sub(res, x, y, 0, ctx);
}

int
padic_radix_sub(padic_radix_t res, const padic_radix_t x, const padic_radix_t y, gr_ctx_t ctx)
{
    return _padic_radix_add_sub(res, x, y, 1, ctx);
}

/*
    Error term of a product: N = min(v_y + N_x, v_x + N_y), an exact factor
    contributing nothing. An inexact zero x = 0 + O(p^N_x) lies in p^N_x Z_p,
    so N_x plays the role of its valuation (the stored v = 0 is meaningless):
    e.g. (0 + O(p^N_x)) (0 + O(p^N_y)) = 0 + O(p^(N_x + N_y)). Using v = 0
    would claim min(N_x, N_y), which is too optimistic when both are negative.
*/
static inline slong
_padic_radix_mul_error(const padic_radix_t x, const padic_radix_t y)
{
    slong vx = (x->u.size == 0) ? x->N : x->v;
    slong vy = (y->u.size == 0) ? y->N : y->v;
    slong Nz = PADIC_RADIX_EXACT;

    if (x->N != PADIC_RADIX_EXACT)
        Nz = vy + x->N;
    if (y->N != PADIC_RADIX_EXACT)
        Nz = FLINT_MIN(Nz, vx + y->N);
    return Nz;
}

/*
    Reference (unoptimized) multiplication: forms the full product of the units
    and lets the finalizer reduce/truncate. Kept as a correctness oracle for the
    optimized padic_radix_mul below (the test suite cross-checks the two). The
    unit product of two canonical units is itself canonical, so this uses the
    canonical-assuming finalizer.
*/
int
_padic_radix_mul_reference(padic_radix_t res, const padic_radix_t x, const padic_radix_t y, gr_ctx_t ctx)
{
    radix_struct * radix = PADIC_RADIX_CTX_RADIX(ctx);
    slong vx = x->v, vy = y->v, Nx = x->N, Ny = y->N;
    slong vz, Nz;
    radix_integer_t U;

    /* an exact zero annihilates exactly, even against an inexact factor */
    if ((x->u.size == 0 && Nx == PADIC_RADIX_EXACT)
     || (y->u.size == 0 && Ny == PADIC_RADIX_EXACT))
        return padic_radix_zero(res, ctx);

    vz = vx + vy;
    Nz = _padic_radix_mul_error(x, y);

    radix_integer_init(U, radix);
    radix_integer_mul(U, &x->u, &y->u, radix);

    FLINT_SWAP(radix_integer_struct, res->u, *U);
    res->v = vz;
    res->N = Nz;

    radix_integer_clear(U, radix);

    return _padic_radix_finalize_assume_canonical(res, ctx);
}

/*
    Multiplication core, as padic_radix_mul with: the inline 1-2 limb
    product, a direct finish for inexact products, and (for addmul) an
    optional cap Kcap: when the product is inexact only its digits below
    min(Neff, Kcap) are formed (not normalised digit-wise when nomod is set;
    a placeholder unit 1 if Kcap <= vz); N is always the true Neff.
*/
static int
_padic_radix_mul_capped(padic_radix_t res, const padic_radix_t x, const padic_radix_t y,
    slong Kcap, int nomod, gr_ctx_t ctx)
{
    const radix_struct * radix = PADIC_RADIX_CTX_RADIX(ctx);
    slong e = radix->exp;
    slong vx = x->v, vy = y->v, Nx = x->N, Ny = y->N;
    slong vz, Nz, xn, yn, target;

    if ((x->u.size == 0 && Nx == PADIC_RADIX_EXACT)
     || (y->u.size == 0 && Ny == PADIC_RADIX_EXACT))
        return padic_radix_zero(res, ctx);

    vz = vx + vy;
    Nz = _padic_radix_mul_error(x, y);

    if (x->u.size == 0 || y->u.size == 0)
    {
        radix_integer_zero(&res->u, radix);
        res->v = vz;
        res->N = Nz;
        return _padic_radix_finalize_assume_canonical(res, ctx);
    }

    xn = FLINT_ABS(x->u.size);
    yn = FLINT_ABS(y->u.size);
    target = _padic_radix_target_prec(vz, ctx);

    if (Nz == PADIC_RADIX_EXACT)
    {
        if (target == PADIC_RADIX_PREC_INF)
        {
            if (!_padic_radix_exact_digits_ok((xn + yn) * e, ctx))
                return GR_UNABLE;
            radix_integer_mul(&res->u, &x->u, &y->u, radix);
            res->N = PADIC_RADIX_EXACT;
        }
        else
        {
            slong digits = target - vz;
            slong maxD = (xn + yn) * e;
            slong minD = (xn + yn - 2) * e + 1;

            if (digits <= 0)
            {
                radix_integer_zero(&res->u, radix);
                res->N = target;
            }
            else if (maxD <= digits)
            {
                radix_integer_mul(&res->u, &x->u, &y->u, radix);
                res->N = PADIC_RADIX_EXACT;
            }
            else if (minD > digits)
            {
                _padic_radix_unit_mullow(res, x, y, radix_div_exp(digits + e - 1, radix), radix);
                res->N = target;
            }
            else
            {
                radix_integer_mul(&res->u, &x->u, &y->u, radix);
                res->N = PADIC_RADIX_EXACT;
            }
        }
    }
    else
    {
        slong Neff = (target == PADIC_RADIX_PREC_INF) ? _padic_radix_cap_prec(Nz, vz, ctx) : FLINT_MIN(Nz, target);

        if (Neff - vz <= 0)
        {
            radix_integer_zero(&res->u, radix);
            res->N = Neff;
        }
        else
        {
            int inrange = (vz <= PADIC_RADIX_ERR_MAX && vz >= -PADIC_RADIX_ERR_MAX
                           && Nz <= PADIC_RADIX_ERR_MAX);
            slong K = (inrange && Kcap < Neff) ? Kcap : Neff;

            if (K <= vz)
            {
                radix_integer_fit_limbs(&res->u, 1, radix)[0] = 1;
                res->u.size = 1;
                res->v = vz;
                res->N = Neff;
                return GR_SUCCESS;
            }

            _padic_radix_unit_mullow(res, x, y, radix_div_exp(K - vz + e - 1, radix), radix);
            res->v = vz;

            if (inrange)
            {
                res->N = Neff;
                if (!nomod)
                    _radix_integer_mod_digits_inplace(&res->u, K - vz, radix);
                return GR_SUCCESS;
            }

            res->N = Nz;
        }
    }

    res->v = vz;
    return _padic_radix_finalize_assume_canonical(res, ctx);
}

/*
    res = x * y, with truncating multiplication.

    The product unit u_x * u_y has valuation v_z = v_x + v_y and is needed only
    to the result's absolute precision; we therefore compute it with
    radix_integer_mullow_limbs to just the required number of limbs rather than
    forming the full product (cf. the reference above). The unit product of two
    canonical units is canonical, so no canonicalisation is needed.

    The error term N_z = min(v_y + N_x, v_x + N_y) (an exact factor contributing
    no error). When both inputs are exact the product may still be exact, so we
    detect that and avoid introducing an error radius unnecessarily, using the
    operand sizes:
      - if the product certainly fits the target ((xn+yn)*e <= digits), keep it
        exact (full product, no reduction);
      - if it certainly overflows ((xn+yn-1)*e + 1 > digits), it is inexact, so
        truncate to ceil(digits/e) limbs;
      - otherwise (a one-limb-wide borderline) form the full product and let the
        reducer's fits test decide.
    When an input is inexact the result is always inexact, and we truncate to the
    ceil(r/e) limbs that survive reduction (r = relative precision in digits).

    The final reduction/representability handling is shared with the reference
    via _padic_radix_finalize_assume_canonical.
*/
int
padic_radix_mul(padic_radix_t res, const padic_radix_t x, const padic_radix_t y, gr_ctx_t ctx)
{
    return _padic_radix_mul_capped(res, x, y, WORD_MAX, 0, ctx);
}

/* ------------------------------------------------------------------------- */
/*    Division                                                               */
/* ------------------------------------------------------------------------- */

/*
    res = a / b. The divisor b must be checked for being zero by the caller's
    contract: a definite zero gives GR_DOMAIN, an unknown-zero gives GR_UNABLE.
    For two exact inputs whose units divide over Z the quotient is returned
    exactly; otherwise the unit inverse is computed to the working precision.
*/
static int
_padic_radix_div(padic_radix_t res, const padic_radix_t a, const padic_radix_t b, gr_ctx_t ctx)
{
    radix_struct * radix = PADIC_RADIX_CTX_RADIX(ctx);
    slong e = radix->exp;
    slong va = a->v, vb = b->v, Na = a->N, Nb = b->N;
    slong prec_abs = PADIC_RADIX_CTX_PREC_ABS(ctx);
    slong prec_rel = PADIC_RADIX_CTX_PREC_REL(ctx);
    slong vq, rel_a, rel_b, rel, relcap, krel, nlimbs;
    radix_integer_t unit;
    truth_t bz;
    int invertible;

    bz = padic_radix_is_zero(b, ctx);
    if (bz == T_TRUE)
        return GR_DOMAIN;                       /* division by zero */
    if (bz == T_UNKNOWN)
        return GR_UNABLE;                        /* divisor not known to be nonzero */

    /* b is definitely nonzero, so its unit is a genuine p-adic unit */

    if (a->u.size == 0 && Na == PADIC_RADIX_EXACT)
        return padic_radix_zero(res, ctx);       /* exact 0 / b = exact 0 */

    vq = va - vb;

    /* exact / exact: return an exact quotient when the units divide over Z.

       TODO: depending on the operand sizes this exactness test may be better
       performed as a Hensel division with remainder (radix_divmod_bn with a
       non-NULL rem) rather than a Euclidean division. */
    if (Na == PADIC_RADIX_EXACT && Nb == PADIC_RADIX_EXACT)
    {
        radix_integer_t qint;
        int divides;

        radix_integer_init(qint, radix);
        divides = radix_integer_div(qint, &a->u, &b->u, radix);

        if (divides)
        {
            FLINT_SWAP(radix_integer_struct, res->u, *qint);
            res->v = vq;
            res->N = PADIC_RADIX_EXACT;
            radix_integer_clear(qint, radix);
            return _padic_radix_finalize_assume_canonical(res, ctx);
        }

        radix_integer_clear(qint, radix);
        /* fall through: infinite expansion, truncate to context precision */
    }

    /* relative precision of the quotient = min of the operand relative
       precisions, then capped by the context precisions */
    rel_a = (Na == PADIC_RADIX_EXACT) ? PADIC_RADIX_PREC_INF : (Na - va);
    rel_b = (Nb == PADIC_RADIX_EXACT) ? PADIC_RADIX_PREC_INF : (Nb - vb);
    rel = FLINT_MIN(rel_a, rel_b);

    relcap = PADIC_RADIX_PREC_INF;
    if (prec_rel != PADIC_RADIX_PREC_INF)
        relcap = prec_rel;
    if (prec_abs != PADIC_RADIX_PREC_INF)
    {
        slong c = prec_abs - vq;                 /* bounded: no overflow */
        relcap = (relcap == PADIC_RADIX_PREC_INF) ? c : FLINT_MIN(relcap, c);
    }

    /* (without context precision, at most size_limit digits) */
    if (relcap == PADIC_RADIX_PREC_INF && rel != PADIC_RADIX_PREC_INF && (ulong) rel > ctx->size_limit)
        rel = ctx->size_limit;

    if (rel == PADIC_RADIX_PREC_INF && relcap == PADIC_RADIX_PREC_INF)
        return GR_UNABLE;                        /* exact ring, infinite expansion */

    if (rel == PADIC_RADIX_PREC_INF)
        krel = relcap;
    else if (relcap == PADIC_RADIX_PREC_INF)
        krel = rel;
    else
        krel = FLINT_MIN(rel, relcap);

    if (krel > PADIC_RADIX_ERR_MAX)
        krel = PADIC_RADIX_ERR_MAX;              /* never compute more than we can store */

    if (krel <= 0)
    {
        /* quotient valuation is at/below the precision horizon */
        radix_integer_zero(&res->u, radix);
        res->v = 0;
        res->N = vq + krel;                      /* bounded; finalize clamps */
        return _padic_radix_finalize_assume_canonical(res, ctx);
    }

    nlimbs = radix_div_exp(krel + e - 1, radix);                 /* digits -> limbs (ceil) */

    radix_integer_init(unit, radix);

    if (a->u.size == 0)
    {
        /* a is an inexact zero: the quotient unit is zero. */
        radix_integer_zero(unit, radix);
        invertible = 1;                          /* b is a genuine unit */
    }
    else
    {
        slong an = FLINT_MIN(FLINT_ABS(a->u.size), nlimbs);
        slong bn = FLINT_MIN(FLINT_ABS(b->u.size), nlimbs);
        slong sgn = a->u.size ^ b->u.size;
        slong rn = nlimbs;
        nn_ptr ud = radix_integer_fit_limbs(unit, nlimbs, radix);

        invertible = radix_divmod_bn(ud, NULL, a->u.d, an, b->u.d, bn, nlimbs, radix);

        if (invertible)
        {
            MPN_NORM(ud, rn);
            unit->size = (sgn >= 0) ? rn : -rn;
        }
    }

    if (!invertible)
    {
        radix_integer_clear(unit, radix);
        return GR_UNABLE;                        /* defensive: should not happen */
    }

    FLINT_SWAP(radix_integer_struct, res->u, *unit);
    res->v = vq;
    res->N = vq + krel;                          /* bounded; finalize clamps + truncates */

    radix_integer_clear(unit, radix);

    return _padic_radix_finalize_assume_canonical(res, ctx);
}

int
padic_radix_div(padic_radix_t res, const padic_radix_t x, const padic_radix_t y, gr_ctx_t ctx)
{
    return _padic_radix_div(res, x, y, ctx);
}

int
padic_radix_inv(padic_radix_t res, const padic_radix_t x, gr_ctx_t ctx)
{
    radix_struct * radix = PADIC_RADIX_CTX_RADIX(ctx);
    slong e = radix->exp;
    slong vx = x->v, Nx = x->N;
    slong prec_abs = PADIC_RADIX_CTX_PREC_ABS(ctx);
    slong prec_rel = PADIC_RADIX_CTX_PREC_REL(ctx);
    slong vq, rel, relcap, krel, nlimbs;
    truth_t xz;
    int invertible;

    xz = padic_radix_is_zero(x, ctx);
    if (xz == T_TRUE)
        return GR_DOMAIN;                        /* 1 / 0 */
    if (xz == T_UNKNOWN)
        return GR_UNABLE;                        /* x not known to be nonzero */

    /* x is a genuine p-adic unit times p^vx, so 1/x has valuation -vx. */
    vq = -vx;

    /* The only units that invert exactly over Z are +-1; an exact such x gives
       an exact +-p^{-vx}. Any other unit inverts to an infinite expansion,
       handled by truncation to the context precision below. */
    if (Nx == PADIC_RADIX_EXACT && FLINT_ABS(x->u.size) == 1 && x->u.d[0] == 1)
    {
        radix_integer_fit_limbs(&res->u, 1, radix)[0] = 1;
        res->u.size = (x->u.size > 0) ? 1 : -1;
        res->v = vq;
        res->N = PADIC_RADIX_EXACT;
        return _padic_radix_finalize_assume_canonical(res, ctx);
    }

    /* relative precision of 1/x equals that of x, capped by the context
       precisions (mirrors _padic_radix_div with an exact unit numerator). */
    rel = (Nx == PADIC_RADIX_EXACT) ? PADIC_RADIX_PREC_INF : (Nx - vx);

    relcap = PADIC_RADIX_PREC_INF;
    if (prec_rel != PADIC_RADIX_PREC_INF)
        relcap = prec_rel;
    if (prec_abs != PADIC_RADIX_PREC_INF)
    {
        slong c = prec_abs - vq;                 /* bounded: no overflow */
        relcap = (relcap == PADIC_RADIX_PREC_INF) ? c : FLINT_MIN(relcap, c);
    }

    /* (without context precision, at most size_limit digits) */
    if (relcap == PADIC_RADIX_PREC_INF && rel != PADIC_RADIX_PREC_INF && (ulong) rel > ctx->size_limit)
        rel = ctx->size_limit;

    if (rel == PADIC_RADIX_PREC_INF && relcap == PADIC_RADIX_PREC_INF)
        return GR_UNABLE;                        /* exact ring, infinite expansion */

    if (rel == PADIC_RADIX_PREC_INF)
        krel = relcap;
    else if (relcap == PADIC_RADIX_PREC_INF)
        krel = rel;
    else
        krel = FLINT_MIN(rel, relcap);

    if (krel > PADIC_RADIX_ERR_MAX)
        krel = PADIC_RADIX_ERR_MAX;              /* never compute more than we can store */

    if (krel <= 0)
    {
        /* inverse valuation is at/below the precision horizon */
        radix_integer_zero(&res->u, radix);
        res->v = 0;
        res->N = vq + krel;                      /* bounded; finalize clamps */
        return _padic_radix_finalize_assume_canonical(res, ctx);
    }

    nlimbs = radix_div_exp(krel + e - 1, radix);                 /* digits -> limbs (ceil) */

    /* res->u = x_unit^{-1} (mod p^{e*nlimbs}); invmod_limbs preserves the sign
       of x, which is also the sign of 1/x. Handles res aliasing x internally. */
    invertible = radix_integer_invmod_limbs(&res->u, &x->u, nlimbs, radix);
    if (!invertible)
        flint_abort();  /* should not happen */

    res->v = vq;
    res->N = vq + krel;                          /* bounded; finalize clamps + truncates */

    return _padic_radix_finalize_assume_canonical(res, ctx);
}

/* TODO: optimize this */
/*
    If the exact unit u = x->u is a perfect square, store its terminating square
    root in res->u and return 1; otherwise return 0. On entry res->u holds one
    square root of u modulo B^M (M limbs, as produced by radix_integer_sqrtmod_
    limbs); the terminating root, if it exists, is that value or its negation
    modulo B^M (both are checked by squaring and comparing against the full u).
*/
static int
_padic_radix_exact_unit_sqrt(padic_radix_t res, const padic_radix_t x, slong M,
    const radix_t radix)
{
    ulong p = DIGIT_RADIX(radix);
    slong ncand = (p == 2) ? 4 : 2;     /* a unit square has 4 roots mod 2^n */
    nn_ptr root, c, sq, d;
    slong i, k, sqn, cn, rs, cmp_limbs;
    int found = 0, neg = 0;
    TMP_INIT;

    if (x->u.size <= 0)         /* a perfect-square exact unit is positive */
        return 0;

    TMP_START;
    root = TMP_ALLOC(M * sizeof(ulong));
    c    = TMP_ALLOC(M * sizeof(ulong));
    sq   = TMP_ALLOC(2 * M * sizeof(ulong));

    /* the canonical branch residue r, zero-extended to M limbs */
    rs = FLINT_MIN(FLINT_ABS(res->u.size), M);
    flint_mpn_copyi(root, res->u.d, rs);
    flint_mpn_zero(root + rs, M - rs);

    /* limbs of r that are reliable for the sign test: all but the top limb,
       whose top bit is the one the 2-adic root leaves undetermined. */
    cmp_limbs = (M > 1) ? (M - 1) : 1;

    /* candidates: root, -root, and (for p = 2) root + B^M/2, -root + B^M/2. The
       candidate that squares (as an exact integer) to x->u is the small positive
       integer root s; only s, not B^M - s, can match exactly. */
    for (k = 0; k < ncand && !found; k++)
    {
        if (k & 1)
            radix_neg(c, root, M, radix);
        else
            flint_mpn_copyi(c, root, M);

        if (k >= 2)             /* add B^M/2 = 2^{eM-1} in the top limb (p = 2) */
        {
            ulong h = UWORD(1) << (radix->exp - 1);
            ulong t = c[M - 1] + h;
            if (t >= LIMB_RADIX(radix))
                t -= LIMB_RADIX(radix);     /* mod B; carry out of B^M dropped */
            c[M - 1] = t;
        }

        radix_mulmid(sq, c, M, c, M, 0, 2 * M, radix);   /* c^2 (full product) */
        sqn = 2 * M;
        MPN_NORM(sq, sqn);

        if (sqn == x->u.size)
        {
            int eq = 1;
            for (i = 0; i < sqn; i++)
                if (sq[i] != x->u.d[i]) { eq = 0; break; }

            if (eq)
            {
                /* c == s, the small positive integer root. The two p-adic roots
                   are +s and -s; the canonical branch is whichever is congruent
                   to r (mod p^k). Store that signed integer so the exact root
                   agrees with the finite-precision result on every shared digit.
                   r == s (on the reliable limbs) selects +s, otherwise -s (whose
                   nonnegative residue B^M - s equals r). */
                neg = 0;
                for (i = 0; i < cmp_limbs; i++)
                    if (root[i] != c[i]) { neg = 1; break; }
                found = 1;
            }
        }
    }

    if (found)
    {
        cn = M;
        MPN_NORM(c, cn);
        if (cn == 0)
            cn = 1;                 /* the root of a unit is nonzero */
        d = radix_integer_fit_limbs(&res->u, cn, radix);
        flint_mpn_copyi(d, c, cn);
        res->u.size = neg ? -cn : cn;   /* signed, consistent with the branch */
    }

    TMP_END;
    return found;
}

/*
    square modulo p (modulo 8 for p = 2). The relative precision of the root
    equals that of x, capped by the context precisions, exactly as for inv/div.
    A definite zero gives the exact zero; an unknown (inexact) zero gives
    GR_UNABLE; an odd valuation or a non-residue unit gives GR_DOMAIN. The branch
    is pinned to a canonical representative, and when the input is exact and a
    perfect square the exact (terminating) root is returned regardless of the
    context precision.
*/
int
padic_radix_sqrt(padic_radix_t res, const padic_radix_t x, gr_ctx_t ctx)
{
    radix_struct * radix = PADIC_RADIX_CTX_RADIX(ctx);
    slong e = radix->exp;
    slong vx = x->v, Nx = x->N;
    slong prec_abs = PADIC_RADIX_CTX_PREC_ABS(ctx);
    slong prec_rel = PADIC_RADIX_CTX_PREC_REL(ctx);
    slong vq, rel, relcap, krel, nlimbs, loss, kin, ndet;
    truth_t xz;
    int issquare, input_exact;
    radix_integer_t unit;

    input_exact = (Nx == PADIC_RADIX_EXACT);

    xz = padic_radix_is_zero(x, ctx);
    if (xz == T_TRUE)
        return padic_radix_zero(res, ctx);       /* sqrt(0) = 0 (exact) */
    if (xz == T_UNKNOWN)
        return GR_UNABLE;                        /* x not known to be nonzero */

    /* A square has even valuation. */
    if (vx & WORD(1))
        return GR_DOMAIN;
    vq = vx / 2;

    /* x = + p^vx exactly (unit 1) has the exact root + p^{vq}. (Other exact
       units have an infinite root expansion, truncated to context precision
       below.) */
    if (Nx == PADIC_RADIX_EXACT && x->u.size == 1 && x->u.d[0] == 1)
    {
        radix_integer_fit_limbs(&res->u, 1, radix)[0] = 1;
        res->u.size = 1;
        res->v = vq;
        res->N = PADIC_RADIX_EXACT;
        return _padic_radix_finalize_assume_canonical(res, ctx);
    }

    /* Unit must be known mod 2^3 */
    if (DIGIT_RADIX(radix) == 2 && Nx != PADIC_RADIX_EXACT && Nx - vx < 3)
        return GR_UNABLE;

    /* Relative precision of the root. sqrt preserves relative precision for odd
       p, but the 2-adic square root loses one bit (d/dt t^2 = 2t has valuation 1),
       so an input known to relative precision 'rel' determines the root only to
       'rel - loss'. This loss applies to the input-limited term only: when the
       context's relative cap binds and the input still has bits to spare, the
       capped bits are all determined. Dropping the one unreliable bit is what
       makes sqrt a single precision-independent function (a lower-precision
       result is the truncation of a higher-precision one). */
    loss = (DIGIT_RADIX(radix) == 2) ? 1 : 0;

    rel = (Nx == PADIC_RADIX_EXACT) ? PADIC_RADIX_PREC_INF : (Nx - vx);
    if (rel != PADIC_RADIX_PREC_INF)
        rel = FLINT_MAX(rel - loss, 0);          /* input-limited root precision */

    relcap = PADIC_RADIX_PREC_INF;
    if (prec_rel != PADIC_RADIX_PREC_INF)
        relcap = prec_rel;
    if (prec_abs != PADIC_RADIX_PREC_INF)
    {
        slong c = prec_abs - vq;                 /* bounded: no overflow */
        relcap = (relcap == PADIC_RADIX_PREC_INF) ? c : FLINT_MIN(relcap, c);
    }

    /* (without context precision, at most size_limit digits) */
    if (relcap == PADIC_RADIX_PREC_INF && rel != PADIC_RADIX_PREC_INF && (ulong) rel > ctx->size_limit)
        rel = ctx->size_limit;

    if (rel == PADIC_RADIX_PREC_INF && relcap == PADIC_RADIX_PREC_INF)
        krel = PADIC_RADIX_PREC_INF;             /* exact ring, infinite expansion */
    else if (rel == PADIC_RADIX_PREC_INF)
        krel = relcap;
    else if (relcap == PADIC_RADIX_PREC_INF)
        krel = rel;
    else
        krel = FLINT_MIN(rel, relcap);

    if (krel != PADIC_RADIX_PREC_INF && krel > PADIC_RADIX_ERR_MAX)
        krel = PADIC_RADIX_ERR_MAX;

    /* Input precision to develop: one bit beyond the delivered root precision so
       that the unreliable top bit is computed for headroom and then truncated away
       by the finalizer, leaving every delivered bit pinned. */
    kin = (krel == PADIC_RADIX_PREC_INF) ? PADIC_RADIX_PREC_INF
        : (krel > PADIC_RADIX_ERR_MAX - loss ? PADIC_RADIX_ERR_MAX : krel + loss);

    /* Develop the unit root. If the input is exact the root may terminate (x a
       perfect square), and detecting that needs about half of u's limbs no matter
       what the context precision is -- this is what lets sqrt(4) = 2 exactly even
       in a finite-precision ring. Otherwise develop kin = krel + loss relative
       bits. Non-squareness is read from radix_integer_sqrtmod_limbs, which finds
       the root modulo p itself, so it is never recomputed here. */
    ndet = input_exact ? (FLINT_ABS(x->u.size) / 2 + 2) : 0;

    if (krel == PADIC_RADIX_PREC_INF)
        nlimbs = ndet;                              /* exact ring implies input exact */
    else if (kin <= 0)
        nlimbs = 1;
    else
        nlimbs = (kin + e - 1) / e;

    if (input_exact && ndet > nlimbs)
        nlimbs = ndet;                              /* enough to detect termination */

    /* Materialise the unit's nonnegative residue modulo p^{e*nlimbs} (this also
       resolves a signed exact unit), then take its square root. */
    radix_integer_init(unit, radix);
    radix_integer_mod_digits(unit, &x->u, e * nlimbs, radix);
    issquare = radix_integer_sqrtmod_limbs(&res->u, unit, nlimbs, radix);
    radix_integer_clear(unit, radix);

    if (!issquare)
        return GR_DOMAIN;                            /* unit is not a square mod p */

    /* Exact input + perfect-square unit: return the exact (terminating) root, with
       the sign matching the canonical branch so the exact result agrees with the
       finite-precision result on every shared digit. Works in any context. */
    if (input_exact &&
        _padic_radix_exact_unit_sqrt(res, x, (ndet > 0 ? ndet : nlimbs), radix))
    {
        res->v = vq;
        res->N = PADIC_RADIX_EXACT;
        return _padic_radix_finalize_assume_canonical(res, ctx);
    }

    if (krel == PADIC_RADIX_PREC_INF)
        return GR_UNABLE;       /* exact ring: an irrational root cannot be exact */

    if (krel <= 0)
    {
        /* root is zero to the tracked precision */
        radix_integer_zero(&res->u, radix);
        res->v = 0;
        res->N = vq + krel;                          /* bounded; finalize clamps */
        return _padic_radix_finalize_assume_canonical(res, ctx);
    }

    res->v = vq;
    res->N = vq + krel;                              /* bounded; finalize clamps + truncates */

    return _padic_radix_finalize_assume_canonical(res, ctx);
}

/*
    Whether x is a square. A square has even valuation and a unit that is a square
    in Z_p: for odd p this is decided by the single low digit (a quadratic residue
    modulo p, which always lifts), and that digit is always known; for p = 2 it
    needs the unit modulo 8, so a relative precision below 3 is undecided. A
    definite zero is a square; an inexact (unknown) zero is undecided.
*/
truth_t
padic_radix_is_square(const padic_radix_t x, gr_ctx_t ctx)
{
    radix_struct * radix = PADIC_RADIX_CTX_RADIX(ctx);
    ulong p = DIGIT_RADIX(radix);
    truth_t xz;

    xz = padic_radix_is_zero(x, ctx);
    if (xz == T_TRUE)
        return T_TRUE;                  /* 0 = 0^2 */
    if (xz == T_UNKNOWN)
        return T_UNKNOWN;               /* inexact zero: cannot tell */

    if (x->v & WORD(1))
        return T_FALSE;                 /* a square has even valuation */

    if (p == 2)
    {
        slong rel = (x->N == PADIC_RADIX_EXACT)
            ? PADIC_RADIX_PREC_INF : (x->N - x->v);
        if (rel != PADIC_RADIX_PREC_INF && rel < 3)
            return T_UNKNOWN;           /* need the unit modulo 8 */
    }

    ulong d = x->u.d[0];

    if (x->u.size < 0)
        d = p - d;

    if (p == 2)
        return ((d & 7) == 1) ? T_TRUE : T_FALSE;

    return (n_jacobi_unsigned(nmod_set_ui(d, radix->b), p) == 1) ? T_TRUE : T_FALSE;
}

/*
    Reciprocal square root x^{-1/2}. Computed directly from
    radix_integer_rsqrtmod_limbs (the same primitive that underlies sqrt), not by
    inverting sqrt. The unit's reciprocal root carries valuation -v/2; relative
    precision and the one-bit 2-adic loss are handled exactly as for sqrt. The
    branch is pinned to the reciprocal of the canonical square root, so
    rsqrt(x) * sqrt(x) == 1. GR_DOMAIN when x is not a square (odd valuation or a
    non-residue unit) or is a definite zero; GR_UNABLE when x is not known nonzero,
    or in the exact ring when the reciprocal root has an infinite expansion (every
    case except x = +p^v, whose reciprocal root is the exact + p^{-v/2}).
*/
int
padic_radix_rsqrt(padic_radix_t res, const padic_radix_t x, gr_ctx_t ctx)
{
    radix_struct * radix = PADIC_RADIX_CTX_RADIX(ctx);
    slong e = radix->exp;
    slong vx = x->v, Nx = x->N;
    slong prec_abs = PADIC_RADIX_CTX_PREC_ABS(ctx);
    slong prec_rel = PADIC_RADIX_CTX_PREC_REL(ctx);
    slong vq, vr, rel, relcap, krel, nlimbs, loss, kin;
    truth_t xz;
    int issquare;
    radix_integer_t unit;

    xz = padic_radix_is_zero(x, ctx);
    if (xz == T_TRUE)
        return GR_DOMAIN;                        /* 1/sqrt(0) is undefined */
    if (xz == T_UNKNOWN)
        return GR_UNABLE;                        /* x not known to be nonzero */

    if (vx & WORD(1))
        return GR_DOMAIN;                        /* a square has even valuation */
    vq = vx / 2;
    vr = -vq;                                    /* valuation of x^{-1/2} */

    /* x = + p^vx exactly (unit 1) has the exact reciprocal root + p^{-vq}. */
    if (Nx == PADIC_RADIX_EXACT && x->u.size == 1 && x->u.d[0] == 1)
    {
        radix_integer_fit_limbs(&res->u, 1, radix)[0] = 1;
        res->u.size = 1;
        res->v = vr;
        res->N = PADIC_RADIX_EXACT;
        return _padic_radix_finalize_assume_canonical(res, ctx);
    }

    loss = (DIGIT_RADIX(radix) == 2) ? 1 : 0;

    rel = (Nx == PADIC_RADIX_EXACT) ? PADIC_RADIX_PREC_INF : (Nx - vx);
    if (rel != PADIC_RADIX_PREC_INF)
        rel = FLINT_MAX(rel - loss, 0);

    relcap = PADIC_RADIX_PREC_INF;
    if (prec_rel != PADIC_RADIX_PREC_INF)
        relcap = prec_rel;
    if (prec_abs != PADIC_RADIX_PREC_INF)
    {
        slong c = prec_abs - vr;                 /* absolute cap relative to vr = -vq */
        relcap = (relcap == PADIC_RADIX_PREC_INF) ? c : FLINT_MIN(relcap, c);
    }

    /* (without context precision, at most size_limit digits) */
    if (relcap == PADIC_RADIX_PREC_INF && rel != PADIC_RADIX_PREC_INF && (ulong) rel > ctx->size_limit)
        rel = ctx->size_limit;

    if (rel == PADIC_RADIX_PREC_INF && relcap == PADIC_RADIX_PREC_INF)
        krel = PADIC_RADIX_PREC_INF;             /* exact ring */
    else if (rel == PADIC_RADIX_PREC_INF)
        krel = relcap;
    else if (relcap == PADIC_RADIX_PREC_INF)
        krel = rel;
    else
        krel = FLINT_MIN(rel, relcap);

    if (krel != PADIC_RADIX_PREC_INF && krel > PADIC_RADIX_ERR_MAX)
        krel = PADIC_RADIX_ERR_MAX;

    kin = (krel == PADIC_RADIX_PREC_INF) ? PADIC_RADIX_PREC_INF
        : (krel > PADIC_RADIX_ERR_MAX - loss ? PADIC_RADIX_ERR_MAX : krel + loss);

    /* Develop the reciprocal unit root. In the exact ring the result is infinite
       (the unit-1 case is already returned above), so one limb is enough to decide
       squareness; otherwise develop kin = krel + loss relative bits. Non-squareness
       is read directly from radix_integer_rsqrtmod_limbs. */
    if (krel == PADIC_RADIX_PREC_INF || kin <= 0)
        nlimbs = 1;
    else
        nlimbs = (kin + e - 1) / e;

    radix_integer_init(unit, radix);
    radix_integer_mod_digits(unit, &x->u, e * nlimbs, radix);
    issquare = radix_integer_rsqrtmod_limbs(&res->u, unit, nlimbs, radix);
    radix_integer_clear(unit, radix);

    if (!issquare)
        return GR_DOMAIN;                        /* unit is not a square mod p */

    if (krel == PADIC_RADIX_PREC_INF)
        return GR_UNABLE;       /* exact ring: an infinite reciprocal root */

    if (krel <= 0)
    {
        radix_integer_zero(&res->u, radix);
        res->v = 0;
        res->N = vr + krel;
        return _padic_radix_finalize_assume_canonical(res, ctx);
    }

    res->v = vr;
    res->N = vr + krel;
    return _padic_radix_finalize_assume_canonical(res, ctx);
}

/* ------------------------------------------------------------------------- */
/*    Predicates                                                             */
/* ------------------------------------------------------------------------- */

truth_t
padic_radix_is_zero(const padic_radix_t x, gr_ctx_t ctx)
{
    if (x->N == PADIC_RADIX_EXACT)
        return (x->u.size == 0) ? T_TRUE : T_FALSE;

    /* inexact: x = u p^v + O(p^N) */
    if (x->u.size == 0)
        return T_UNKNOWN;            /* 0 + O(p^N): might be zero or not */

    /* nonzero unit: its lowest (nonzero) digit sits at position v */
    return (x->N > x->v) ? T_FALSE : T_UNKNOWN;
}

truth_t
padic_radix_is_one(const padic_radix_t x, gr_ctx_t ctx)
{
    radix_struct * radix = PADIC_RADIX_CTX_RADIX(ctx);

    if (x->N == PADIC_RADIX_EXACT)
        return (x->v == 0 && radix_integer_is_one(&x->u, radix)) ? T_TRUE : T_FALSE;

    if (x->N <= 0)
        return T_UNKNOWN;            /* not even the units digit is known */

    /* one has valuation 0 and units digit 1 */
    if (x->u.size == 0)
        return T_UNKNOWN;            /* 0 + O(p^N): cannot rule one in or out beyond horizon */

    if (x->v != 0)
        return T_FALSE;              /* units digit (0) differs from 1, and N >= 1 so it is known */

    /* v == 0: compare the unit's residue modulo p^N to 1 (units are nonnegative) */
    return _radix_integer_residue_is_one_mod(&x->u, x->N, radix) ? T_UNKNOWN : T_FALSE;
}

truth_t
padic_radix_is_neg_one(const padic_radix_t x, gr_ctx_t ctx)
{
    radix_struct * radix = PADIC_RADIX_CTX_RADIX(ctx);

    if (x->N == PADIC_RADIX_EXACT)
        return (x->v == 0 && radix_integer_is_neg_one(&x->u, radix)) ? T_TRUE : T_FALSE;

    if (x->N <= 0)
        return T_UNKNOWN;            /* not even the units digit is known */

    /* -1 has valuation 0; its residue mod p^N is p^N - 1 (all digits p-1) */
    if (x->v != 0)
        return T_FALSE;              /* x has a known digit incompatible with -1 */

    /* v == 0: compare the unit's residue modulo p^N to p^N - 1 (units nonnegative) */
    return _radix_integer_residue_is_neg_one_mod(&x->u, x->N, radix) ? T_UNKNOWN : T_FALSE;
}

truth_t
padic_radix_equal(const padic_radix_t x, const padic_radix_t y, gr_ctx_t ctx)
{
    radix_struct * radix = PADIC_RADIX_CTX_RADIX(ctx);
    slong err, vx = x->v, vy = y->v;
    int xz, yz, eq;

    if (x->N == PADIC_RADIX_EXACT && y->N == PADIC_RADIX_EXACT)
    {
        if (vx != vy)
            return T_FALSE;
        return radix_integer_equal(&x->u, &y->u, radix) ? T_TRUE : T_FALSE;
    }

    /* at least one is inexact: they agree as ring elements modulo p^err */
    err = FLINT_MIN(x->N, y->N);   /* EXACT == WORD_MAX acts as +inf */

    if (err <= 0)
        return T_UNKNOWN;

    /* Decide x == y (mod p^err) structurally from the valuations and units. The
       lowest digit of a unit (either sign) sits at its valuation and is nonzero,
       so differing valuations disagree at min(vx, vy); equal valuations reduce to
       a unit congruence that _radix_units_congruent_mod settles in O(unit size),
       regardless of how large err or the valuations are. */
    xz = (x->u.size == 0);
    yz = (y->u.size == 0);

    if (xz && yz)
        eq = 1;
    else if (xz)
        eq = (vy >= err);                 /* y == 0 mod p^err iff its valuation reaches the horizon */
    else if (yz)
        eq = (vx >= err);
    else if (vx == vy)
        eq = _radix_units_congruent_mod(&x->u, &y->u, err - vx, radix);
    else
        eq = (FLINT_MIN(vx, vy) >= err);  /* differ at the lower valuation if it is within the horizon */

    return eq ? T_UNKNOWN : T_FALSE;
}

truth_t
padic_radix_is_invertible(const padic_radix_t x, gr_ctx_t ctx)
{
    /* the radix p-adic numbers form a field: x is invertible iff x is nonzero */
    return truth_not(padic_radix_is_zero(x, ctx));
}

/* ------------------------------------------------------------------------- */
/*    Randomisation                                                          */
/* ------------------------------------------------------------------------- */

int
padic_radix_randtest(padic_radix_t res, flint_rand_t state, gr_ctx_t ctx)
{
    radix_struct * radix = PADIC_RADIX_CTX_RADIX(ctx);
    int ctx_exact = (PADIC_RADIX_CTX_PREC_ABS(ctx) == PADIC_RADIX_PREC_INF
                  && PADIC_RADIX_CTX_PREC_REL(ctx) == PADIC_RADIX_PREC_INF);

    /* the limits mode is only meaningful for an inexact ring (it produces
       error terms and over-large valuations that an exact ring cannot hold) */
    int limits = (PADIC_RADIX_CTX_FLAGS(ctx) & PADIC_RADIX_TEST_LIMITS)
              && PADIC_RADIX_CTX_PREC_REL(ctx) != PADIC_RADIX_PREC_INF
              && (n_randint(state, 2) == 0);

    if (n_randint(state, 2))
        radix_integer_set_si(&res->u, n_randint(state, 5) - 2, radix);
    else
        radix_integer_randtest_limbs(&res->u, state, 1 + n_randint(state, 10), radix);

    if (!PADIC_RADIX_CTX_SIGNED(ctx))
        radix_integer_abs(&res->u, &res->u, radix);

    if (limits && (n_randint(state, 4) == 0))
    {
        /* valuation and precision at or just inside the representable range, so
           that subsequent operations push the intermediate exponents over the
           ERR_MAX boundary and exercise the clamping in _padic_radix_finalize */
        slong big = PADIC_RADIX_ERR_MAX - n_randint(state, 4);

        res->v = (n_randlimb(state) & 1) ? big : -big;

        switch (n_randint(state, 3))
        {
            case 0:  res->N = PADIC_RADIX_EXACT; break;
            case 1:  res->N = (n_randlimb(state) & 1) ? big : -big; break;
            default: res->N = res->v + n_randint(state, 6); break;
        }
    }
    else
    {
        res->v = n_randint(state, 4);
        res->N = PADIC_RADIX_EXACT;

        if (!ctx_exact && (n_randlimb(state) & 1))
            res->N = res->v + n_randint(state, 6);   /* random finite absolute precision */
    }

    if (_padic_radix_finalize(res, ctx) != GR_SUCCESS)
        padic_radix_zero(res, ctx);                  /* always return a valid element */

    return GR_SUCCESS;
}

/* ------------------------------------------------------------------------- */
/*    Output                                                                 */
/* ------------------------------------------------------------------------- */

int
padic_radix_write(gr_stream_t out, const padic_radix_t x, gr_ctx_t ctx)
{
    radix_struct * radix = PADIC_RADIX_CTX_RADIX(ctx);
    ulong p = GR_PADIC_RADIX_CTX(ctx)->p;
    int status = GR_SUCCESS;

    if (x->u.size == 0)
    {
        status |= gr_stream_write(out, "0");
    }
    else
    {
        slong size = x->u.size;
        char * str;

        if (PADIC_RADIX_CTX_DECIMAL(ctx))
            str = radix_get_str_decimal(NULL, x->u.d, FLINT_ABS(size), size < 0, radix);
        else
            str = radix_get_str_sum(NULL, x->u.d, FLINT_ABS(size), size < 0, 1, radix);

        if (x->v != 0)
            status |= gr_stream_write(out, "(");
        status |= gr_stream_write_free(out, str);
        if (x->v != 0)
            status |= gr_stream_write(out, ")");

        if (x->v != 0)
        {
            status |= gr_stream_write(out, " * ");
            status |= gr_stream_write_ui(out, p);
            status |= gr_stream_write(out, "^");
            status |= gr_stream_write_si(out, x->v);
        }
    }

    if (x->N != PADIC_RADIX_EXACT)
    {
        status |= gr_stream_write(out, " + O(");
        status |= gr_stream_write_ui(out, p);
        status |= gr_stream_write(out, "^");
        status |= gr_stream_write_si(out, x->N);
        status |= gr_stream_write(out, ")");
    }

    return status;
}

/* ------------------------------------------------------------------------- */
/*    Context                                                                */
/* ------------------------------------------------------------------------- */

void
padic_radix_ctx_clear(gr_ctx_t ctx)
{
    radix_clear(PADIC_RADIX_CTX_RADIX(ctx));
    flint_free(GR_PADIC_RADIX_CTX(ctx));
}

int
padic_radix_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    status |= gr_stream_write(out, "Radix ");
    status |= gr_stream_write_ui(out, GR_PADIC_RADIX_CTX(ctx)->p);
    status |= gr_stream_write(out, "-adic numbers");

    status |= gr_stream_write(out, " (rel prec ");
    if (PADIC_RADIX_CTX_PREC_REL(ctx) == PADIC_RADIX_PREC_INF)
        status |= gr_stream_write(out, "inf");
    else
        status |= gr_stream_write_si(out, PADIC_RADIX_CTX_PREC_REL(ctx));

    status |= gr_stream_write(out, ", abs prec ");
    if (PADIC_RADIX_CTX_PREC_ABS(ctx) == PADIC_RADIX_PREC_INF)
        status |= gr_stream_write(out, "inf");
    else
        status |= gr_stream_write_si(out, PADIC_RADIX_CTX_PREC_ABS(ctx));

    status |= gr_stream_write(out, ")");
    return status;
}

static truth_t
padic_radix_ctx_is_exact(gr_ctx_t ctx)
{
    return (PADIC_RADIX_CTX_PREC_ABS(ctx) == PADIC_RADIX_PREC_INF
         && PADIC_RADIX_CTX_PREC_REL(ctx) == PADIC_RADIX_PREC_INF) ? T_TRUE : T_FALSE;
}

int _padic_radix_methods_initialized = 0;

gr_static_method_table _padic_radix_methods;

gr_method_tab_input _padic_radix_methods_input[] =
{
    {GR_METHOD_CTX_CLEAR,       (gr_funcptr) padic_radix_ctx_clear},
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) padic_radix_ctx_write},
    {GR_METHOD_CTX_IS_RING,     (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FIELD,    (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_UNIQUE_FACTORIZATION_DOMAIN,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_RATIONAL_VECTOR_SPACE,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_EXACT,    (gr_funcptr) padic_radix_ctx_is_exact},
    {GR_METHOD_INIT,            (gr_funcptr) padic_radix_init},
    {GR_METHOD_CLEAR,           (gr_funcptr) padic_radix_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) padic_radix_swap},
    {GR_METHOD_SET_SHALLOW,     (gr_funcptr) padic_radix_set_shallow},
    {GR_METHOD_RANDTEST,        (gr_funcptr) padic_radix_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) padic_radix_write},
    {GR_METHOD_BIG_O_BASE_FMPZ, (gr_funcptr) padic_radix_big_o_base_fmpz},
    {GR_METHOD_ZERO,            (gr_funcptr) padic_radix_zero},
    {GR_METHOD_ONE,             (gr_funcptr) padic_radix_one},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) padic_radix_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) padic_radix_is_one},
    {GR_METHOD_IS_NEG_ONE,      (gr_funcptr) padic_radix_is_neg_one},
    {GR_METHOD_IS_INVERTIBLE,   (gr_funcptr) padic_radix_is_invertible},
    {GR_METHOD_EQUAL,           (gr_funcptr) padic_radix_equal},
    {GR_METHOD_SET,             (gr_funcptr) padic_radix_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) padic_radix_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) padic_radix_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) padic_radix_set_fmpz},
    {GR_METHOD_GET_FMPZ,        (gr_funcptr) padic_radix_get_fmpz},
    {GR_METHOD_NEG,             (gr_funcptr) padic_radix_neg},
    {GR_METHOD_ADD,             (gr_funcptr) padic_radix_add},
    {GR_METHOD_SUB,             (gr_funcptr) padic_radix_sub},
    {GR_METHOD_MUL,             (gr_funcptr) padic_radix_mul},
    {GR_METHOD_MUL_SI,          (gr_funcptr) padic_radix_mul_si},
    {GR_METHOD_MUL_UI,          (gr_funcptr) padic_radix_mul_ui},
    {GR_METHOD_MUL_FMPZ,        (gr_funcptr) padic_radix_mul_fmpz},
    {GR_METHOD_ADDMUL,          (gr_funcptr) padic_radix_addmul},
    {GR_METHOD_SUBMUL,          (gr_funcptr) padic_radix_submul},
    {GR_METHOD_ADDMUL_SI,       (gr_funcptr) padic_radix_addmul_si},
    {GR_METHOD_SUBMUL_SI,       (gr_funcptr) padic_radix_submul_si},
    {GR_METHOD_ADDMUL_UI,       (gr_funcptr) padic_radix_addmul_ui},
    {GR_METHOD_SUBMUL_UI,       (gr_funcptr) padic_radix_submul_ui},
    {GR_METHOD_ADDMUL_FMPZ,     (gr_funcptr) padic_radix_addmul_fmpz},
    {GR_METHOD_SUBMUL_FMPZ,     (gr_funcptr) padic_radix_submul_fmpz},
    {GR_METHOD_POW_UI,          (gr_funcptr) padic_radix_pow_ui},
    {GR_METHOD_POW_SI,          (gr_funcptr) padic_radix_pow_si},
    {GR_METHOD_POW_FMPZ,        (gr_funcptr) padic_radix_pow_fmpz},
    {GR_METHOD_INV,             (gr_funcptr) padic_radix_inv},
    {GR_METHOD_DIV,             (gr_funcptr) padic_radix_div},
    {GR_METHOD_SQRT,            (gr_funcptr) padic_radix_sqrt},
    {GR_METHOD_RSQRT,           (gr_funcptr) padic_radix_rsqrt},
    {GR_METHOD_IS_SQUARE,       (gr_funcptr) padic_radix_is_square},
    {GR_METHOD_EXP,             (gr_funcptr) padic_radix_exp},
    {GR_METHOD_LOG,             (gr_funcptr) padic_radix_log},
    {GR_METHOD_VEC_DOT,         (gr_funcptr) padic_radix_dot},
    {GR_METHOD_VEC_DOT_REV,     (gr_funcptr) padic_radix_dot_rev},
    {GR_METHOD_VEC_DOT_STRIDED, (gr_funcptr) padic_radix_dot_strided},
    {GR_METHOD_VEC_MUL_SCALAR,  (gr_funcptr) _padic_radix_vec_mul_scalar},
    {GR_METHOD_VEC_MUL_SCALAR_SI,   (gr_funcptr) _padic_radix_vec_mul_scalar_si},
    {GR_METHOD_VEC_MUL_SCALAR_UI,   (gr_funcptr) _padic_radix_vec_mul_scalar_ui},
    {GR_METHOD_VEC_MUL_SCALAR_FMPZ, (gr_funcptr) _padic_radix_vec_mul_scalar_fmpz},
    {GR_METHOD_VEC_ADDMUL_SCALAR,   (gr_funcptr) _padic_radix_vec_addmul_scalar},
    {GR_METHOD_VEC_SUBMUL_SCALAR,   (gr_funcptr) _padic_radix_vec_submul_scalar},
    {GR_METHOD_VEC_ADDMUL_SCALAR_SI,    (gr_funcptr) _padic_radix_vec_addmul_scalar_si},
    {GR_METHOD_VEC_SUBMUL_SCALAR_SI,    (gr_funcptr) _padic_radix_vec_submul_scalar_si},
    {GR_METHOD_VEC_ADDMUL_SCALAR_FMPZ,  (gr_funcptr) _padic_radix_vec_addmul_scalar_fmpz},
    /* todo: this ought to be automatic */
    /* todo: use different algorithm with exact precision */
    {GR_METHOD_POLY_MULLOW,     (gr_funcptr) _padic_radix_poly_mullow},
    {GR_METHOD_POLY_MULMID,     (gr_funcptr) _padic_radix_poly_mulmid},
    {GR_METHOD_MAT_MUL,         (gr_funcptr) padic_radix_mat_mul},
    {GR_METHOD_MAT_DET,         (gr_funcptr) padic_radix_mat_det},
    {GR_METHOD_MAT_FIND_NONZERO_PIVOT, (gr_funcptr) padic_radix_mat_find_nonzero_pivot},
    {GR_METHOD_MAT_NONSINGULAR_SOLVE, (gr_funcptr) padic_radix_mat_nonsingular_solve},
    {GR_METHOD_MAT_INV,         (gr_funcptr) padic_radix_mat_inv},
    {GR_METHOD_MAT_CHARPOLY,    (gr_funcptr) padic_radix_mat_charpoly},
    {0,                         (gr_funcptr) NULL},
};

int
gr_ctx_init_padic_radix(gr_ctx_t ctx, ulong p, slong prec_rel, slong prec_abs, int flags)
{
    padic_radix_ctx_struct * pctx;

    FLINT_ASSERT(n_is_prime(p));

    ctx->which_ring = GR_CTX_PADIC_RADIX;
    ctx->sizeof_elem = sizeof(padic_radix_struct);
    ctx->size_limit = WORD_MAX;

    GR_CTX_DATA_AS_PTR(ctx) = flint_malloc(sizeof(padic_radix_ctx_struct));
    pctx = GR_PADIC_RADIX_CTX(ctx);

    radix_init(&pctx->radix, p, 0);   /* limb radix p^e with e chosen automatically */
    pctx->p = p;

    if (prec_rel == PADIC_RADIX_PREC_INF)
        pctx->prec_rel = PADIC_RADIX_PREC_INF;
    else
        pctx->prec_rel = FLINT_MIN(FLINT_MAX(prec_rel, 0), PADIC_RADIX_ERR_MAX);

    if (prec_abs == PADIC_RADIX_PREC_INF)
        pctx->prec_abs = PADIC_RADIX_PREC_INF;
    else
        pctx->prec_abs = FLINT_MAX(FLINT_MIN(prec_abs, PADIC_RADIX_ERR_MAX),
                                   -PADIC_RADIX_ERR_MAX);

    pctx->flags = flags;

    ctx->methods = _padic_radix_methods;

    if (!_padic_radix_methods_initialized)
    {
        gr_method_tab_init(_padic_radix_methods, _padic_radix_methods_input);
        _padic_radix_methods_initialized = 1;
    }

    return GR_SUCCESS;
}

int
gr_ctx_init_padic_radix_randtest(gr_ctx_t ctx, flint_rand_t state, slong maxprec)
{
    ulong p = n_randtest_prime(state, 0);
    slong prec_abs, prec_rel;
    int flags = 0;

    if (n_randint(state, 2))
        flags |= PADIC_RADIX_SIGNED;

    switch (n_randint(state, 4))
    {
        case 0:   /* classic padic: absolute precision only */
            prec_abs = 1 + n_randint(state, maxprec);
            prec_rel = PADIC_RADIX_PREC_INF;
            break;
        case 1:   /* relative precision only */
            prec_abs = PADIC_RADIX_PREC_INF;
            prec_rel = 1 + n_randint(state, maxprec);
            break;
        case 2:   /* both */
            prec_abs = 1 + n_randint(state, maxprec);
            prec_rel = 1 + n_randint(state, maxprec);
            break;
        default:  /* exact ring */
            prec_abs = PADIC_RADIX_PREC_INF;
            prec_rel = PADIC_RADIX_PREC_INF;
            break;
    }

    if (n_randint(state, 2))
        flags |= PADIC_RADIX_TEST_LIMITS;

    if (n_randint(state, 4) == 0)
        flags |= PADIC_RADIX_FAST_LINALG;

    /* allow e.g. get_fmpz roundtrip to return unable for huge test values */
    if (flags & PADIC_RADIX_TEST_LIMITS)
        ctx->size_limit = (1 << 16);

    return gr_ctx_init_padic_radix(ctx, p, prec_rel, prec_abs, flags);
}

/* ------------------------------------------------------------------------- */
/*    Fused multiply-add, scalar and vector operations, powering            */
/* ------------------------------------------------------------------------- */

static inline void
_padic_radix_stack_init(padic_radix_struct * T, nn_ptr buf, slong alloc)
{
    T->u.d = buf;
    T->u.alloc = alloc;
    T->u.size = 0;
    T->v = 0;
    T->N = PADIC_RADIX_EXACT;
}

#define PADIC_RADIX_STACK_LIMBS 24

/* (x0 + x1 B) * p^d mod B^L for 0 <= d < L e, L in {1, 2}: one inline
   truncated product (d < e) or a limb shift and a one-limb product */
static inline void
_padic_radix_lsh_12(ulong * a0, ulong * a1, ulong x0, ulong x1, slong d, slong L,
    const radix_struct * radix)
{
    ulong t[2];

    if (d == 0)
    {
        *a0 = x0;
        *a1 = (L == 2) ? x1 : 0;
    }
    else if (d < radix->exp)
    {
        slong n = _padic_radix_mullow_12(t, x0, x1, radix->bpow[d], 0, L, radix->B);
        *a0 = t[0];
        *a1 = (L == 2 && n == 2) ? t[1] : 0;
    }
    else
    {
        _padic_radix_mullow_12(t, x0, 0, radix->bpow[d - radix->exp], 0, 1, radix->B);
        *a0 = 0;
        *a1 = t[0];
    }
}

/*
    Fused addmul for small results: an inline product of the units and one
    truncation, for positive units, an inexact product and a result of at
    most two limbs. The operand of higher valuation is first multiplied by
    p^(v_hi - v_lo) (an inline truncated product), and a cancelled lowest
    digit (v_t == v_res) needs a right digit shift.
*/
static inline int
_padic_radix_addmul_small(padic_radix_t res, const padic_radix_t x, const padic_radix_t y,
    int sub, gr_ctx_t ctx)
{
    const radix_struct * radix = PADIC_RADIX_CTX_RADIX(ctx);
    slong e = radix->exp;
    slong xs = x->u.size, ys = y->u.size, rs = res->u.size;
    slong vx, vy, vr, vt, Nx, Ny, Ntin, Nt, tt, tm, vmin, Nz, P, Ldig, L, rem, w;
    ulong B = radix->B.n, a0 = 0, a1 = 0, b0 = 0, b1 = 0, s0, s1 = 0;
    nn_ptr d;

    if (xs <= 0 || ys <= 0 || rs <= 0 || rs > 2)
        return 0;

    Nx = x->N;
    Ny = y->N;
    if (Nx == PADIC_RADIX_EXACT && Ny == PADIC_RADIX_EXACT)
        return 0;

    vx = x->v;
    vy = y->v;
    vr = res->v;
    vt = vx + vy;

    Ntin = PADIC_RADIX_EXACT;
    if (Nx != PADIC_RADIX_EXACT)
        Ntin = vy + Nx;
    if (Ny != PADIC_RADIX_EXACT)
        Ntin = FLINT_MIN(Ntin, vx + Ny);

    if (vt > PADIC_RADIX_ERR_MAX || vt < -PADIC_RADIX_ERR_MAX || Ntin > PADIC_RADIX_ERR_MAX
        || vr > PADIC_RADIX_ERR_MAX || vr < -PADIC_RADIX_ERR_MAX)
        return 0;

    tt = _padic_radix_target_prec(vt, ctx);
    Nt = (tt == PADIC_RADIX_PREC_INF) ? Ntin : FLINT_MIN(Ntin, tt);
    if (Nt <= vt)
        return 0;

    vmin = FLINT_MIN(vt, vr);
    Nz = FLINT_MIN(res->N, Nt);
    tm = _padic_radix_target_prec(vmin, ctx);
    P = (tm == PADIC_RADIX_PREC_INF) ? Nz : FLINT_MIN(Nz, tm);
    if (vmin >= P)
        return 0;

    Ldig = P - vmin;
    if (Ldig > 2 * e)
        return 0;
    L = (Ldig > e) ? 2 : 1;

    if (vr < P)
        _padic_radix_lsh_12(&a0, &a1, res->u.d[0], (rs > 1) ? res->u.d[1] : 0, vr - vmin, L, radix);

    if (vt < P)
    {
        ulong t[2];
        ulong x1 = (xs > 1) ? x->u.d[1] : 0, y1 = (ys > 1) ? y->u.d[1] : 0;
        slong m = radix_div_exp(P - vt + e - 1, radix);
        slong tn = _padic_radix_mullow_12(t, x->u.d[0], x1, y->u.d[0], y1, m, radix->B);
        _padic_radix_lsh_12(&b0, &b1, t[0], (tn == 2) ? t[1] : 0, vt - vmin, L, radix);
    }

    if (!sub)
    {
        ulong cy;
        s0 = a0 + b0;
        cy = (s0 < a0) | (s0 >= B);
        if (cy)
            s0 -= B;
        if (L == 2)
        {
            s1 = a1 + b1;
            if ((s1 < a1) | (s1 >= B))
                s1 -= B;
            s1 += cy;
            if (s1 == B)
                s1 = 0;
        }
    }
    else
    {
        int bw = (a0 < b0);
        s0 = a0 - b0;
        if (bw)
            s0 += B;
        if (L == 2)
        {
            s1 = a1 - b1;
            if (a1 < b1)
                s1 += B;
            if (bw)
                s1 = (s1 == 0) ? B - 1 : s1 - 1;
        }
    }

    rem = Ldig - (L - 1) * e;
    if (rem < e)
    {
        if (L == 1)
            s0 = n_rem_precomp(s0, radix->bpow[rem], radix->bpow_div + rem);
        else
            s1 = n_rem_precomp(s1, radix->bpow[rem], radix->bpow_div + rem);
    }

    if (s0 == 0 && s1 == 0)
    {
        if (P != Nz)
            return 0;
        res->u.size = 0;
        res->v = 0;
        res->N = P;
        return 1;
    }

    /* valuation of the sum (only nonzero when v_t == v_res) */
    if (s0 != 0)
        w = _padic_radix_divisible_p(s0, radix) ? _padic_radix_val1(s0, radix) : 0;
    else
        w = e + _padic_radix_val1(s1, radix);

    if (w != 0)
    {
        if (P != Nz)
            return 0;   /* defensive; cannot happen for valid inputs */

        /* S >>= w digits */
        if (w >= e)
        {
            s0 = n_div_precomp(s1, radix->bpow_div + (w - e));
            s1 = 0;
        }
        else
        {
            ulong q1 = n_div_precomp(s1, radix->bpow_div + w);
            ulong r1 = s1 - q1 * radix->bpow[w];
            s0 = n_div_precomp(s0, radix->bpow_div + w) + r1 * radix->bpow[e - w];
            s1 = q1;
        }
    }

    d = radix_integer_fit_limbs(&res->u, 2, radix);
    d[0] = s0;
    d[1] = s1;
    res->u.size = (s1 != 0) ? 2 : 1;
    res->v = vmin + w;
    res->N = P;
    return 1;
}

/* res = res +/- x y, identical to mul into a temporary then add/sub; see
   the comment below for the cap argument */
static int
_padic_radix_addmul_core(padic_radix_t res, const padic_radix_t x, const padic_radix_t y,
    int sub, nn_ptr buf, slong alloc, gr_ctx_t ctx)
{
    padic_radix_struct T;
    slong Kcap = WORD_MAX;
    int status;

    if (_padic_radix_addmul_small(res, x, y, sub, ctx))
        return GR_SUCCESS;

    if ((x->N != PADIC_RADIX_EXACT || y->N != PADIC_RADIX_EXACT)
        && x->u.size != 0 && y->u.size != 0)
    {
        slong vt = x->v + y->v;

        Kcap = res->N;
        if (res->u.size != 0)
        {
            slong t = _padic_radix_target_prec(FLINT_MIN(vt, res->v), ctx);
            Kcap = FLINT_MIN(Kcap, t);
        }
    }

    if (x->u.size < 0 || y->u.size < 0)
    {
        padic_radix_t H;
        padic_radix_init(H, ctx);
        status = _padic_radix_mul_capped(H, x, y, Kcap, 1, ctx);
        status |= _padic_radix_add_sub(res, res, H, sub, ctx);
        padic_radix_clear(H, ctx);
        return status;
    }

    _padic_radix_stack_init(&T, buf, alloc);
    status = _padic_radix_mul_capped(&T, x, y, Kcap, 1, ctx);
    status |= _padic_radix_add_sub(res, res, &T, sub, ctx);
    return status;
}

int
_padic_radix_addmul_sub(padic_radix_t res, const padic_radix_t x, const padic_radix_t y, int sub, gr_ctx_t ctx)
{
    slong need = FLINT_ABS(x->u.size) + FLINT_ABS(y->u.size) + 2;
    int status;

    if (need <= PADIC_RADIX_STACK_LIMBS)
    {
        ulong buf[PADIC_RADIX_STACK_LIMBS];
        status = _padic_radix_addmul_core(res, x, y, sub, buf, PADIC_RADIX_STACK_LIMBS, ctx);
    }
    else
    {
        nn_ptr buf;
        TMP_INIT;
        TMP_START;
        buf = TMP_ALLOC(need * sizeof(ulong));
        status = _padic_radix_addmul_core(res, x, y, sub, buf, need, ctx);
        TMP_END;
    }

    return status;
}

int padic_radix_addmul(padic_radix_t res, const padic_radix_t x, const padic_radix_t y, gr_ctx_t ctx)
{ return _padic_radix_addmul_sub(res, x, y, 0, ctx); }
int padic_radix_submul(padic_radix_t res, const padic_radix_t x, const padic_radix_t y, gr_ctx_t ctx)
{ return _padic_radix_addmul_sub(res, x, y, 1, ctx); }

/* ---- scalars: |a| in radix B, valuation stripped (digit shift), then the
   set_si rounding (_padic_radix_finish) ---- */
static int
_padic_radix_set_scalar(padic_radix_t C, ulong a, int neg, gr_ctx_t ctx)
{
    const radix_struct * radix = PADIC_RADIX_CTX_RADIX(ctx);
    ulong d0, d1, B = radix->B.n;
    slong w = 0, t, n;

    C->N = PADIC_RADIX_EXACT;

    if (a == 0)
    {
        C->u.size = 0;
        C->v = 0;
        return GR_SUCCESS;
    }

    if (a < B) { d0 = a; d1 = 0; }
    else       { d1 = a / B; d0 = a - d1 * B; }

    if (d0 == 0)
    {
        w = radix->exp;
        d0 = d1;
        d1 = 0;
    }

    t = _padic_radix_val1(d0, radix);
    if (t != 0)
    {
        ulong q1 = n_div_precomp(d1, radix->bpow_div + t);
        ulong r1 = d1 - q1 * radix->bpow[t];
        d0 = n_div_precomp(d0, radix->bpow_div + t) + r1 * radix->bpow[radix->exp - t];
        d1 = q1;
    }
    w += t;

    if (neg && !PADIC_RADIX_CTX_SIGNED(ctx))
    {
        slong target = _padic_radix_target_prec(w, ctx);
        slong need;

        if (target == PADIC_RADIX_PREC_INF)
            return GR_UNABLE;

        need = (slong) radix_div_exp(FLINT_MAX(target - w, 0) + radix->exp - 1, radix) + 2;
        if (need > C->u.alloc)
        {
            C->u.d = flint_malloc(need * sizeof(ulong));
            C->u.alloc = need;
        }
    }

    n = (d1 != 0) ? 2 : 1;
    C->u.d[0] = d0;
    C->u.d[1] = d1;
    C->u.size = neg ? -n : n;
    C->v = w;

    return _padic_radix_finish(C, ctx);
}

static int
_padic_radix_set_scalar_fmpz(padic_radix_t C, const fmpz_t a, gr_ctx_t ctx)
{
    if (!COEFF_IS_MPZ(*a))
    {
        slong s = *a;
        return _padic_radix_set_scalar(C, (s < 0) ? -(ulong) s : (ulong) s, s < 0, ctx);
    }
    else
    {
        padic_radix_t t;
        int status;
        padic_radix_init(t, ctx);
        status = padic_radix_set_fmpz(t, a, ctx);
        if (status == GR_SUCCESS)
        {
            C->u = t->u;
            C->v = t->v;
            C->N = t->N;
        }
        else
            padic_radix_clear(t, ctx);
        return status;
    }
}

#define PADIC_RADIX_SCALAR_DECL(C, Cbuf) \
    ulong Cbuf[2]; padic_radix_struct C[1]; _padic_radix_stack_init(C, Cbuf, 2)
#define PADIC_RADIX_SCALAR_CLEAR(C, Cbuf) \
    do { if (C->u.d != Cbuf) flint_free(C->u.d); } while (0)
#define PADIC_RADIX_SI(a) (a < 0) ? -(ulong) a : (ulong) a, a < 0

int
padic_radix_set_si(padic_radix_t res, slong a, gr_ctx_t ctx)
{
    PADIC_RADIX_SCALAR_DECL(C, Cbuf);
    int status = _padic_radix_set_scalar(C, PADIC_RADIX_SI(a), ctx);
    if (status == GR_SUCCESS)
    {
        radix_integer_set(&res->u, &C->u, PADIC_RADIX_CTX_RADIX(ctx));
        res->v = C->v;
        res->N = C->N;
    }
    PADIC_RADIX_SCALAR_CLEAR(C, Cbuf);
    return status;
}

int
padic_radix_set_ui(padic_radix_t res, ulong a, gr_ctx_t ctx)
{
    PADIC_RADIX_SCALAR_DECL(C, Cbuf);
    int status = _padic_radix_set_scalar(C, a, 0, ctx);
    if (status == GR_SUCCESS)
    {
        radix_integer_set(&res->u, &C->u, PADIC_RADIX_CTX_RADIX(ctx));
        res->v = C->v;
        res->N = C->N;
    }
    PADIC_RADIX_SCALAR_CLEAR(C, Cbuf);
    return status;
}

#define PADIC_RADIX_SCALAR_BINOP(name, ctype, SETCALL, BODY)                          \
int name(padic_radix_t res, const padic_radix_t x, ctype a, gr_ctx_t ctx)    \
{                                                                            \
    PADIC_RADIX_SCALAR_DECL(C, Cbuf);                                                 \
    int status = SETCALL;                                                    \
    if (status == GR_SUCCESS)                                                \
        status = BODY;                                                       \
    PADIC_RADIX_SCALAR_CLEAR(C, Cbuf);                                                \
    return status;                                                           \
}

PADIC_RADIX_SCALAR_BINOP(padic_radix_mul_si, slong, _padic_radix_set_scalar(C, PADIC_RADIX_SI(a), ctx), padic_radix_mul(res, x, C, ctx))
PADIC_RADIX_SCALAR_BINOP(padic_radix_mul_ui, ulong, _padic_radix_set_scalar(C, a, 0, ctx), padic_radix_mul(res, x, C, ctx))
PADIC_RADIX_SCALAR_BINOP(padic_radix_mul_fmpz, const fmpz_t, _padic_radix_set_scalar_fmpz(C, a, ctx), padic_radix_mul(res, x, C, ctx))
PADIC_RADIX_SCALAR_BINOP(padic_radix_addmul_si, slong, _padic_radix_set_scalar(C, PADIC_RADIX_SI(a), ctx), _padic_radix_addmul_sub(res, x, C, 0, ctx))
PADIC_RADIX_SCALAR_BINOP(padic_radix_submul_si, slong, _padic_radix_set_scalar(C, PADIC_RADIX_SI(a), ctx), _padic_radix_addmul_sub(res, x, C, 1, ctx))
PADIC_RADIX_SCALAR_BINOP(padic_radix_addmul_ui, ulong, _padic_radix_set_scalar(C, a, 0, ctx), _padic_radix_addmul_sub(res, x, C, 0, ctx))
PADIC_RADIX_SCALAR_BINOP(padic_radix_submul_ui, ulong, _padic_radix_set_scalar(C, a, 0, ctx), _padic_radix_addmul_sub(res, x, C, 1, ctx))
PADIC_RADIX_SCALAR_BINOP(padic_radix_addmul_fmpz, const fmpz_t, _padic_radix_set_scalar_fmpz(C, a, ctx), _padic_radix_addmul_sub(res, x, C, 0, ctx))
PADIC_RADIX_SCALAR_BINOP(padic_radix_submul_fmpz, const fmpz_t, _padic_radix_set_scalar_fmpz(C, a, ctx), _padic_radix_addmul_sub(res, x, C, 1, ctx))

/* ---- vectors ---- */
static int
_padic_radix_vec_addmul_submul_scalar(padic_radix_struct * vec1, const padic_radix_struct * vec2,
    slong len, const padic_radix_struct * c, int sub, gr_ctx_t ctx)
{
    ulong sbuf[PADIC_RADIX_STACK_LIMBS];
    nn_ptr buf = sbuf;
    slong alloc = PADIC_RADIX_STACK_LIMBS, i, maxn = 0, need;
    int status = GR_SUCCESS;

    if (len <= 0)
        return GR_SUCCESS;

    for (i = 0; i < len; i++)
        maxn = FLINT_MAX(maxn, FLINT_ABS(vec2[i].u.size));

    need = maxn + FLINT_ABS(c->u.size) + 2;
    if (need > alloc)
    {
        buf = flint_malloc(need * sizeof(ulong));
        alloc = need;
    }

    for (i = 0; i < len; i++)
    {
        need = FLINT_ABS(vec2[i].u.size) + FLINT_ABS(c->u.size) + 2;
        if (need > alloc)
        {
            if (buf != sbuf)
                flint_free(buf);
            buf = flint_malloc(need * sizeof(ulong));
            alloc = need;
        }
        status |= _padic_radix_addmul_core(vec1 + i, vec2 + i, c, sub, buf, alloc, ctx);
    }

    if (buf != sbuf)
        flint_free(buf);
    return status;
}

int _padic_radix_vec_addmul_scalar(padic_radix_struct * v1, const padic_radix_struct * v2, slong len, const padic_radix_t c, gr_ctx_t ctx)
{ return _padic_radix_vec_addmul_submul_scalar(v1, v2, len, c, 0, ctx); }
int _padic_radix_vec_submul_scalar(padic_radix_struct * v1, const padic_radix_struct * v2, slong len, const padic_radix_t c, gr_ctx_t ctx)
{ return _padic_radix_vec_addmul_submul_scalar(v1, v2, len, c, 1, ctx); }

int
_padic_radix_vec_mul_scalar(padic_radix_struct * v1, const padic_radix_struct * v2, slong len, const padic_radix_t c, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong i;
    for (i = 0; i < len; i++)
        status |= padic_radix_mul(v1 + i, v2 + i, c, ctx);
    return status;
}

#define PADIC_RADIX_VEC_SCALAR(name, ctype, SETCALL, BODY)                                  \
int name(padic_radix_struct * vec1, const padic_radix_struct * vec2, slong len,    \
    ctype a, gr_ctx_t ctx)                                                         \
{                                                                                  \
    PADIC_RADIX_SCALAR_DECL(C, Cbuf);                                                       \
    int status;                                                                    \
    if (len <= 0)                                                                  \
        return GR_SUCCESS;                                                         \
    status = SETCALL;                                                              \
    if (status == GR_SUCCESS)                                                      \
        status = BODY;                                                             \
    PADIC_RADIX_SCALAR_CLEAR(C, Cbuf);                                                      \
    return status;                                                                 \
}

PADIC_RADIX_VEC_SCALAR(_padic_radix_vec_addmul_scalar_si, slong, _padic_radix_set_scalar(C, PADIC_RADIX_SI(a), ctx), _padic_radix_vec_addmul_submul_scalar(vec1, vec2, len, C, 0, ctx))
PADIC_RADIX_VEC_SCALAR(_padic_radix_vec_submul_scalar_si, slong, _padic_radix_set_scalar(C, PADIC_RADIX_SI(a), ctx), _padic_radix_vec_addmul_submul_scalar(vec1, vec2, len, C, 1, ctx))
PADIC_RADIX_VEC_SCALAR(_padic_radix_vec_addmul_scalar_fmpz, const fmpz_t, _padic_radix_set_scalar_fmpz(C, a, ctx), _padic_radix_vec_addmul_submul_scalar(vec1, vec2, len, C, 0, ctx))
PADIC_RADIX_VEC_SCALAR(_padic_radix_vec_mul_scalar_si, slong, _padic_radix_set_scalar(C, PADIC_RADIX_SI(a), ctx), _padic_radix_vec_mul_scalar(vec1, vec2, len, C, ctx))
PADIC_RADIX_VEC_SCALAR(_padic_radix_vec_mul_scalar_ui, ulong, _padic_radix_set_scalar(C, a, 0, ctx), _padic_radix_vec_mul_scalar(vec1, vec2, len, C, ctx))
PADIC_RADIX_VEC_SCALAR(_padic_radix_vec_mul_scalar_fmpz, const fmpz_t, _padic_radix_set_scalar_fmpz(C, a, ctx), _padic_radix_vec_mul_scalar(vec1, vec2, len, C, ctx))

/* ---- powering on the unit (here the stored unit itself) ---- */

/* v_p(n) for n >= 1 */
static slong
_padic_radix_fmpz_val_p(const fmpz_t n, ulong p)
{
    if (!COEFF_IS_MPZ(*n))
    {
        ulong m = FLINT_ABS(*n), q;
        slong k = 0;
        while (m != 0 && m % p == 0)
        {
            q = m / p;
            m = q;
            k++;
        }
        return k;
    }
    else
    {
        fmpz_t t, pp;
        slong k;
        fmpz_init(t);
        fmpz_init_set_ui(pp, p);
        k = fmpz_remove(t, n, pp);
        fmpz_clear(t);
        fmpz_clear(pp);
        return k;
    }
}

/*
    x^n for n >= 2, x with a nonzero unit, rounded once:
        v = n v_x,   N = (n - 1) v_x + N_x + v_p(n)   (EXACT if x is exact).
    The extra v_p(n): for x = a + h with v(a) = v_x, v(h) >= N_x and
    r = N_x - v_x >= 1, every term C(n,k) a^(n-k) h^k (k >= 1) of
    x^n - a^n has valuation >= (n-1) v_x + N_x + v_p(n) + (k-1) r - v_p(k)
    (using v_p(C(n,k)) >= v_p(n) - v_p(k)), and v_p(k) <= k - 1 <= (k-1) r.
    Returns -1 if the exponents leave the representable range.
*/
static int
_padic_radix_pow_core(padic_radix_t res, const padic_radix_t x, const fmpz_t n, gr_ctx_t ctx)
{
    const radix_struct * radix = PADIC_RADIX_CTX_RADIX(ctx);
    slong e = radix->exp;
    slong vx = x->v, Nx = x->N, xn = FLINT_ABS(x->u.size);
    slong vout, Nout, target, Neff, r, L, un, Rn, i, nbits;
    int exact = (Nx == PADIC_RADIX_EXACT), neg, fits;
    nn_ptr u, R, T, rd;
    fmpz_t t;
    TMP_INIT;

    fmpz_init(t);
    fmpz_mul_si(t, n, vx);
    if (fmpz_cmp_si(t, PADIC_RADIX_ERR_MAX) > 0 || fmpz_cmp_si(t, -PADIC_RADIX_ERR_MAX) < 0)
    {
        fmpz_clear(t);
        return -1;
    }
    vout = fmpz_get_si(t);
    if (!exact)
    {
        fmpz_sub_si(t, t, vx);
        fmpz_add_si(t, t, Nx);
        fmpz_add_ui(t, t, _padic_radix_fmpz_val_p(n, GR_PADIC_RADIX_CTX(ctx)->p));
        if (fmpz_cmp_si(t, PADIC_RADIX_ERR_MAX) > 0 || fmpz_cmp_si(t, -PADIC_RADIX_ERR_MAX) < 0)
        {
            fmpz_clear(t);
            return -1;
        }
        Nout = fmpz_get_si(t);
    }
    else
        Nout = PADIC_RADIX_EXACT;
    fmpz_clear(t);

    target = _padic_radix_target_prec(vout, ctx);
    Neff = exact ? target : FLINT_MIN(Nout, target);
    if (!exact && target == PADIC_RADIX_PREC_INF)
        Neff = _padic_radix_cap_prec(Neff, vout, ctx);

    if (Neff == PADIC_RADIX_PREC_INF)
        return -1;

    if (Neff <= vout)
    {
        radix_integer_zero(&res->u, radix);
        res->v = 0;
        res->N = Neff;
        return _padic_radix_finalize_assume_canonical(res, ctx);
    }

    neg = (x->u.size < 0) && fmpz_is_odd(n);

    if (exact && xn == 1 && x->u.d[0] == 1)
    {
        radix_integer_fit_limbs(&res->u, 1, radix)[0] = 1;
        res->u.size = neg ? -1 : 1;
        res->v = vout;
        res->N = PADIC_RADIX_EXACT;
        return _padic_radix_finalize_assume_canonical(res, ctx);
    }

    r = Neff - vout;
    L = radix_div_exp(r + e - 1, radix);

    if (exact && fmpz_cmp_ui(n, L / (xn + 1)) < 0)
        L = (slong) fmpz_get_ui(n) * (xn + 1);

    TMP_START;
    u = TMP_ALLOC((L + 2) * sizeof(ulong));
    R = TMP_ALLOC((2 * L + 2) * sizeof(ulong));
    T = TMP_ALLOC((2 * L + 2) * sizeof(ulong));

    un = FLINT_MIN(xn, L);
    flint_mpn_copyi(u, x->u.d, un);
    fits = exact && (xn <= L);
    while (un > 0 && u[un - 1] == 0)
        un--;

    flint_mpn_copyi(R, u, un);
    Rn = un;
    nbits = fmpz_bits(n);

    for (i = nbits - 2; i >= 0; i--)
    {
        slong hi;

        hi = fits ? 2 * Rn : FLINT_MIN(2 * Rn, L);
        if (Rn <= 2 && hi <= 2)
            Rn = _padic_radix_mullow_12(T, R[0], Rn > 1 ? R[1] : 0, R[0], Rn > 1 ? R[1] : 0, hi, radix->B);
        else
        {
            radix_mulmid(T, R, Rn, R, Rn, 0, hi, radix);
            Rn = hi;
            while (Rn > 0 && T[Rn - 1] == 0)
                Rn--;
        }
        FLINT_SWAP(nn_ptr, R, T);
        if (fits && Rn > L)
        {
            fits = 0;
            Rn = L;
            while (Rn > 0 && R[Rn - 1] == 0)
                Rn--;
        }

        if (fmpz_tstbit(n, i))
        {
            hi = fits ? Rn + un : FLINT_MIN(Rn + un, L);
            if (Rn <= 2 && un <= 2 && hi <= 2)
                Rn = _padic_radix_mullow_12(T, R[0], Rn > 1 ? R[1] : 0, u[0], un > 1 ? u[1] : 0, hi, radix->B);
            else
            {
                radix_mulmid(T, R, Rn, u, un, 0, hi, radix);
                Rn = hi;
                while (Rn > 0 && T[Rn - 1] == 0)
                    Rn--;
            }
            FLINT_SWAP(nn_ptr, R, T);
            if (fits && Rn > L)
            {
                fits = 0;
                Rn = L;
                while (Rn > 0 && R[Rn - 1] == 0)
                    Rn--;
            }
        }
    }

    rd = radix_integer_fit_limbs(&res->u, Rn, radix);
    flint_mpn_copyi(rd, R, Rn);
    res->u.size = neg ? -Rn : Rn;
    res->v = vout;
    res->N = exact ? (fits ? PADIC_RADIX_EXACT : target) : Nout;

    TMP_END;
    return _padic_radix_finalize_assume_canonical(res, ctx);
}

static int
_padic_radix_pow_chain(padic_radix_t res, const padic_radix_t x, const fmpz_t n, gr_ctx_t ctx)
{
    padic_radix_t R, T;
    slong i;
    int status = GR_SUCCESS;

    padic_radix_init(R, ctx);
    padic_radix_init(T, ctx);
    radix_integer_set(&R->u, &x->u, PADIC_RADIX_CTX_RADIX(ctx));
    R->v = x->v;
    R->N = x->N;

    for (i = fmpz_bits(n) - 2; i >= 0 && status == GR_SUCCESS; i--)
    {
        status |= padic_radix_mul(T, R, R, ctx);
        if (status != GR_SUCCESS)
        {
            padic_radix_swap(R, T, ctx);
            break;
        }
        if (fmpz_tstbit(n, i))
            status |= padic_radix_mul(R, T, x, ctx);
        else
            padic_radix_swap(R, T, ctx);
    }

    padic_radix_swap(res, R, ctx);
    padic_radix_clear(R, ctx);
    padic_radix_clear(T, ctx);
    return status;
}

/*
    0 + O(p^N) to the n-th power: the ball p^N Z_p maps into p^(nN) Z_p
    (for any sign of N), so the result is 0 + O(p^(nN)). The generic chain
    gives O(p^N), which is too optimistic when N < 0 (the zero element is
    stored with v = 0, so products of two such elements claim min(N1, N2)
    instead of N1 + N2).
*/
static int
_padic_radix_pow_zero_inexact(padic_radix_t res, slong N, const fmpz_t n, gr_ctx_t ctx)
{
    fmpz_t t;
    slong M;

    fmpz_init(t);
    fmpz_mul_si(t, n, N);
    if (fmpz_cmp_si(t, -PADIC_RADIX_ERR_MAX) < 0)
    {
        fmpz_clear(t);
        return GR_UNABLE;
    }
    M = (fmpz_cmp_si(t, PADIC_RADIX_ERR_MAX) > 0) ? PADIC_RADIX_ERR_MAX : fmpz_get_si(t);
    fmpz_clear(t);

    radix_integer_zero(&res->u, PADIC_RADIX_CTX_RADIX(ctx));
    res->v = 0;
    res->N = M;
    return GR_SUCCESS;
}

/* an exact power computed by the multiplication chain in a context
   without finite precision (where nothing is truncated) has about n |u|
   limbs, which must be within the size limit of the context; in other
   contexts the products are truncated to the precision, and powers of +-1
   stay one limb (the valuation saturates) */
static int
_padic_radix_pow_chain_ok(const padic_radix_t x, const fmpz_t n, gr_ctx_t ctx)
{
    fmpz_t t;
    int ok;

    if (x->N != PADIC_RADIX_EXACT
        || PADIC_RADIX_CTX_PREC_REL(ctx) != PADIC_RADIX_PREC_INF
        || PADIC_RADIX_CTX_PREC_ABS(ctx) != PADIC_RADIX_PREC_INF
        || (FLINT_ABS(x->u.size) == 1 && x->u.d[0] == 1))
        return 1;

    fmpz_init(t);
    fmpz_mul_ui(t, n, FLINT_ABS(x->u.size));
    fmpz_mul_ui(t, t, PADIC_RADIX_CTX_RADIX(ctx)->exp);
    ok = fmpz_cmp_ui(t, ctx->size_limit) <= 0;
    fmpz_clear(t);
    return ok;
}

static int
_padic_radix_pow_fmpz_nonneg(padic_radix_t res, const padic_radix_t x, const fmpz_t n, gr_ctx_t ctx)
{
    int status;

    if (fmpz_is_zero(n))
        return padic_radix_one(res, ctx);

    if (fmpz_is_one(n))
        return padic_radix_set(res, x, ctx);

    if (x->u.size == 0)
    {
        if (x->N == PADIC_RADIX_EXACT)
            return padic_radix_zero(res, ctx);
        return _padic_radix_pow_zero_inexact(res, x->N, n, ctx);
    }

    if (res == x)
    {
        padic_radix_t t;
        padic_radix_init(t, ctx);
        status = _padic_radix_pow_core(t, x, n, ctx);
        if (status == -1)
            status = _padic_radix_pow_chain_ok(x, n, ctx)
                ? _padic_radix_pow_chain(t, x, n, ctx) : GR_UNABLE;
        padic_radix_swap(res, t, ctx);
        padic_radix_clear(t, ctx);
        return status;
    }

    status = _padic_radix_pow_core(res, x, n, ctx);
    if (status == -1)
    {
        if (!_padic_radix_pow_chain_ok(x, n, ctx))
            return GR_UNABLE;
        status = _padic_radix_pow_chain(res, x, n, ctx);
    }
    return status;
}

int
padic_radix_pow_fmpz(padic_radix_t res, const padic_radix_t x, const fmpz_t n, gr_ctx_t ctx)
{
    if (fmpz_sgn(n) >= 0)
        return _padic_radix_pow_fmpz_nonneg(res, x, n, ctx);
    else
    {
        fmpz_t m;
        int status;
        fmpz_init(m);
        fmpz_neg(m, n);
        status = padic_radix_inv(res, x, ctx);
        if (status == GR_SUCCESS)
            status = _padic_radix_pow_fmpz_nonneg(res, res, m, ctx);
        fmpz_clear(m);
        return status;
    }
}

int
padic_radix_pow_ui(padic_radix_t res, const padic_radix_t x, ulong n, gr_ctx_t ctx)
{
    fmpz_t m;
    int status;
    fmpz_init_set_ui(m, n);
    status = _padic_radix_pow_fmpz_nonneg(res, x, m, ctx);
    fmpz_clear(m);
    return status;
}

int
padic_radix_pow_si(padic_radix_t res, const padic_radix_t x, slong n, gr_ctx_t ctx)
{
    int status;
    if (n >= 0)
        return padic_radix_pow_ui(res, x, n, ctx);
    status = padic_radix_inv(res, x, ctx);
    if (status == GR_SUCCESS && n != -1)
        status = padic_radix_pow_ui(res, res, -(ulong) n, ctx);
    return status;
}
