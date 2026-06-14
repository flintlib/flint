/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "radix.h"
#include "fmpz.h"
#include "gr.h"

/* ------------------------------------------------------------------------- */
/*    Internal helpers                                                       */
/* ------------------------------------------------------------------------- */

/* Normalise a limb array length. */
static slong
_norm(nn_srcptr d, slong n)
{
    while (n > 0 && d[n - 1] == 0)
        n--;
    return n;
}

/*
    res = (nonnegative residue of x) mod p^k, as a radix_integer with digit
    range [0, k). Built from a limb-level modulo plus a partial reduction of
    the top limb. Handles negative x (the residue is the nonnegative
    representative in [0, p^k)).
*/
static void
radix_integer_mod_digits_nonneg(radix_integer_t res, const radix_integer_t x,
    slong k, const radix_t radix)
{
    slong e = radix->exp;
    slong limbs, rem, tot;

    if (k <= 0)
    {
        radix_integer_zero(res, radix);
        return;
    }

    limbs = k / e;
    rem = k - limbs * e;       /* 0 <= rem < e */
    tot = limbs + (rem != 0);

    radix_integer_mod_limbs(res, x, tot, radix);   /* nonnegative, size <= tot */

    if (rem != 0 && res->size == tot)
    {
        /* reduce the top limb (index `limbs`) modulo p^rem */
        ulong top = res->d[limbs];
        ulong m = radix->bpow[rem];
        ulong t = n_rem_precomp(top, m, radix->bpow_div + rem);

        res->d[limbs] = t;
        res->size = _norm(res->d, tot);
    }
}

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
    slong e = radix->exp;
    slong un = FLINT_ABS(x->size);
    slong limbs = k / e;
    slong rem = k - limbs * e;

    if (rem == 0)
        return un <= limbs;
    if (un <= limbs)
        return 1;
    if (un == limbs + 1)
        return x->d[limbs] < radix->bpow[rem];
    return 0;
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
    radix_integer_mod_digits_nonneg(rx, ux, k, radix);
    radix_integer_mod_digits_nonneg(ry, uy, k, radix);
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
_radix_padic_canonicalise(radix_padic_t x, gr_ctx_t ctx)
{
    radix_struct * radix = RADIX_PADIC_CTX_RADIX(ctx);
    slong w;

    if (x->u.size == 0)
    {
        x->v = 0;
        return;
    }

    w = radix_integer_valuation_digits(&x->u, radix);

    if (w != 0)
    {
        radix_integer_rshift_digits(&x->u, &x->u, w, radix);
        x->v += w;
    }
}

/*
    Target absolute precision for a value of valuation v under the context
    default precisions: min(prec_abs, v + prec_rel), with PREC_INF acting as
    +infinity.
*/
static slong
_radix_padic_target_prec(slong v, gr_ctx_t ctx)
{
    slong a = RADIX_PADIC_CTX_PREC_ABS(ctx);
    slong r = RADIX_PADIC_CTX_PREC_REL(ctx);
    slong t;

    if (r == RADIX_PADIC_PREC_INF)
    {
        t = RADIX_PADIC_PREC_INF;
    }
    else
    {
        /* r is finite and clamped to [0, ERR_MAX]; v is bounded in practice */
        t = v + r;
        if (t < 0 || t > RADIX_PADIC_ERR_MAX)  /* overflow / very large -> treat as inf-ish */
            t = (v >= 0) ? RADIX_PADIC_PREC_INF : t;
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
_radix_padic_reduce(radix_padic_t x, gr_ctx_t ctx)
{
    radix_struct * radix = RADIX_PADIC_CTX_RADIX(ctx);
    slong target, Neff, r;

    if (x->u.size == 0)
    {
        x->v = 0;
        if (x->N != RADIX_PADIC_EXACT)
            x->N = FLINT_MIN(x->N, RADIX_PADIC_ERR_MAX);
        return;
    }

    target = _radix_padic_target_prec(x->v, ctx);
    Neff = (x->N == RADIX_PADIC_EXACT) ? target : FLINT_MIN(x->N, target);

    if (Neff == RADIX_PADIC_PREC_INF)
        return;   /* stays exact */

    r = Neff - x->v;

    if (r <= 0)
    {
        radix_integer_zero(&x->u, radix);
        x->v = 0;
        x->N = FLINT_MIN(Neff, RADIX_PADIC_ERR_MAX);
        return;
    }

    /* exact value that fits: keep it exact (allows signed units) */
    if (x->N == RADIX_PADIC_EXACT && _radix_integer_fits_digits(&x->u, r, radix))
        return;

    /* otherwise the result is inexact: store the nonnegative residue */
    radix_integer_mod_digits_nonneg(&x->u, &x->u, r, radix);
    if (x->u.size == 0)
        x->v = 0;
    x->N = FLINT_MIN(Neff, RADIX_PADIC_ERR_MAX);
}

/*
    Finish a conversion: a freshly built (possibly signed) canonical element
    with N == EXACT. Handles the unsigned-mode representability of negative
    values and then imposes the context precision.
*/
static int
_radix_padic_finish(radix_padic_t res, gr_ctx_t ctx)
{
    if (res->u.size < 0 && !RADIX_PADIC_CTX_SIGNED(ctx))
    {
        slong target = _radix_padic_target_prec(res->v, ctx);

        if (target == RADIX_PADIC_PREC_INF)
            return GR_UNABLE;   /* no finite residue: cannot represent exactly */

        res->N = target;        /* force inexact; reduce takes the nonneg residue */
    }

    _radix_padic_reduce(res, ctx);
    return GR_SUCCESS;
}

/* ------------------------------------------------------------------------- */
/*    Memory management                                                      */
/* ------------------------------------------------------------------------- */

void
radix_padic_init(radix_padic_t res, gr_ctx_t ctx)
{
    radix_integer_init(&res->u, RADIX_PADIC_CTX_RADIX(ctx));
    res->v = 0;
    res->N = RADIX_PADIC_EXACT;
}

void
radix_padic_clear(radix_padic_t res, gr_ctx_t ctx)
{
    radix_integer_clear(&res->u, RADIX_PADIC_CTX_RADIX(ctx));
}

void
radix_padic_swap(radix_padic_t x, radix_padic_t y, gr_ctx_t ctx)
{
    FLINT_SWAP(radix_padic_struct, *x, *y);
}

void
radix_padic_set_shallow(radix_padic_t res, const radix_padic_t x, gr_ctx_t ctx)
{
    *res = *x;
}

/* ------------------------------------------------------------------------- */
/*    Error term access                                                      */
/* ------------------------------------------------------------------------- */

slong
radix_padic_get_error(const radix_padic_t x, gr_ctx_t ctx)
{
    return x->N;
}

truth_t
radix_padic_is_exact(const radix_padic_t x, gr_ctx_t ctx)
{
    return (x->N == RADIX_PADIC_EXACT) ? T_TRUE : T_FALSE;
}

/* ------------------------------------------------------------------------- */
/*    Assignments                                                            */
/* ------------------------------------------------------------------------- */

int
radix_padic_zero(radix_padic_t res, gr_ctx_t ctx)
{
    radix_integer_zero(&res->u, RADIX_PADIC_CTX_RADIX(ctx));
    res->v = 0;
    res->N = RADIX_PADIC_EXACT;
    return GR_SUCCESS;
}

int
radix_padic_one(radix_padic_t res, gr_ctx_t ctx)
{
    radix_integer_one(&res->u, RADIX_PADIC_CTX_RADIX(ctx));
    res->v = 0;
    res->N = RADIX_PADIC_EXACT;
    _radix_padic_reduce(res, ctx);   /* becomes 0 + O(p^N) if target precision <= 0 */
    return GR_SUCCESS;
}

int
radix_padic_set(radix_padic_t res, const radix_padic_t x, gr_ctx_t ctx)
{
    if (res != x)
    {
        radix_integer_set(&res->u, &x->u, RADIX_PADIC_CTX_RADIX(ctx));
        res->v = x->v;
        res->N = x->N;
    }

    _radix_padic_reduce(res, ctx);
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
_radix_padic_exact_set_finish(radix_padic_t res, gr_ctx_t ctx)
{
    _radix_padic_canonicalise(res, ctx);

    if (res->u.size != 0
        && (res->v > RADIX_PADIC_ERR_MAX || res->v < -RADIX_PADIC_ERR_MAX))
        return GR_UNABLE;

    return GR_SUCCESS;
}

int
radix_padic_exact_set_fmpz(radix_padic_t res, const fmpz_t c, gr_ctx_t ctx)
{
    radix_integer_set_fmpz(&res->u, c, RADIX_PADIC_CTX_RADIX(ctx));
    res->v = 0;
    res->N = RADIX_PADIC_EXACT;
    return _radix_padic_exact_set_finish(res, ctx);
}

int
radix_padic_exact_set_ui(radix_padic_t res, ulong c, gr_ctx_t ctx)
{
    radix_integer_set_ui(&res->u, c, RADIX_PADIC_CTX_RADIX(ctx));
    res->v = 0;
    res->N = RADIX_PADIC_EXACT;
    return _radix_padic_exact_set_finish(res, ctx);
}

int
radix_padic_exact_set_si(radix_padic_t res, slong c, gr_ctx_t ctx)
{
    radix_integer_set_si(&res->u, c, RADIX_PADIC_CTX_RADIX(ctx));
    res->v = 0;
    res->N = RADIX_PADIC_EXACT;
    return _radix_padic_exact_set_finish(res, ctx);
}

int
radix_padic_set_fmpz(radix_padic_t res, const fmpz_t c, gr_ctx_t ctx)
{
    int status = radix_padic_exact_set_fmpz(res, c, ctx);
    if (status != GR_SUCCESS)
        return status;
    return _radix_padic_finish(res, ctx);        /* round to context precision */
}

int
radix_padic_set_ui(radix_padic_t res, ulong c, gr_ctx_t ctx)
{
    int status = radix_padic_exact_set_ui(res, c, ctx);
    if (status != GR_SUCCESS)
        return status;
    return _radix_padic_finish(res, ctx);
}

int
radix_padic_set_si(radix_padic_t res, slong c, gr_ctx_t ctx)
{
    int status = radix_padic_exact_set_si(res, c, ctx);
    if (status != GR_SUCCESS)
        return status;
    return _radix_padic_finish(res, ctx);
}

int
radix_padic_get_fmpz(fmpz_t res, const radix_padic_t x, gr_ctx_t ctx)
{
    radix_struct * radix = RADIX_PADIC_CTX_RADIX(ctx);

    if (x->N != RADIX_PADIC_EXACT)
        return GR_UNABLE;       /* not a definite integer */

    if (x->v < 0)
        return GR_DOMAIN;       /* genuine p-adic non-integer */

    if (x->u.size == 0)
    {
        fmpz_zero(res);
        return GR_SUCCESS;
    }

    if (x->v > ctx->size_limit)
        return GR_UNABLE;

    radix_integer_get_fmpz(res, &x->u, radix);

    if (x->v != 0)
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_ui_pow_ui(t, GR_RADIX_PADIC_CTX(ctx)->p, x->v);
        fmpz_mul(res, res, t);
        fmpz_clear(t);
    }

    return GR_SUCCESS;
}

/* ------------------------------------------------------------------------- */
/*    Arithmetic                                                             */
/* ------------------------------------------------------------------------- */

int
radix_padic_neg(radix_padic_t res, const radix_padic_t x, gr_ctx_t ctx)
{
    radix_struct * radix = RADIX_PADIC_CTX_RADIX(ctx);

    if (res != x)
        radix_padic_set(res, x, ctx);

    radix_integer_neg(&res->u, &res->u, radix);

    if (res->N == RADIX_PADIC_EXACT)
        return _radix_padic_finish(res, ctx);

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
            radix_integer_mod_digits_nonneg(&res->u, &res->u, r, radix);
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
    the canonicalisation scan that _radix_padic_finalize performs.
*/
static int
_radix_padic_finalize_assume_canonical(radix_padic_t res, gr_ctx_t ctx)
{
    radix_struct * radix = RADIX_PADIC_CTX_RADIX(ctx);

    /*
        Enforce the representability invariant |v| <= ERR_MAX and (N == EXACT
        or |N| <= ERR_MAX). Keeping both exponents within ERR_MAX = WORD_MAX/4
        guarantees that the sums formed by add/sub/mul (v_x + v_y, v_y + N_x,
        ...) never overflow a signed word.
    */
    if (res->N != RADIX_PADIC_EXACT)
    {
        if (res->N < -RADIX_PADIC_ERR_MAX)
            return GR_UNABLE;                       /* precision too low to represent */
        if (res->N > RADIX_PADIC_ERR_MAX)
            res->N = RADIX_PADIC_ERR_MAX;           /* clamp down: always sound */
    }

    if (res->u.size != 0)
    {
        if (res->v > RADIX_PADIC_ERR_MAX)
        {
            /* u * p^v with v > ERR_MAX is == 0 modulo p^ERR_MAX */
            slong newN = (res->N == RADIX_PADIC_EXACT)
                       ? RADIX_PADIC_ERR_MAX
                       : FLINT_MIN(res->N, RADIX_PADIC_ERR_MAX);
            radix_integer_zero(&res->u, radix);
            res->v = 0;
            res->N = newN;                          /* now finite */
        }
        else if (res->v < -RADIX_PADIC_ERR_MAX)
        {
            return GR_UNABLE;                       /* valuation too small to represent */
        }
    }

    if (res->N == RADIX_PADIC_EXACT)
        return _radix_padic_finish(res, ctx);

    _radix_padic_reduce(res, ctx);
    return GR_SUCCESS;
}

static int
_radix_padic_finalize(radix_padic_t res, gr_ctx_t ctx)
{
    _radix_padic_canonicalise(res, ctx);
    return _radix_padic_finalize_assume_canonical(res, ctx);
}

/* res = x + y  (sub != 0 selects x - y). */
static int
_radix_padic_add_sub(radix_padic_t res, const radix_padic_t x,
    const radix_padic_t y, int sub, gr_ctx_t ctx)
{
    radix_struct * radix = RADIX_PADIC_CTX_RADIX(ctx);
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
    target = _radix_padic_target_prec(vmin, ctx);
    if (Nz == RADIX_PADIC_EXACT)
        P = target;
    else if (target == RADIX_PADIC_PREC_INF)
        P = Nz;
    else
        P = FLINT_MIN(Nz, target);

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
        return _radix_padic_finalize_assume_canonical(res, ctx);
    }

    /* whole result is below the precision horizon: it is 0 + O(p^P) */
    if (P != RADIX_PADIC_PREC_INF && vmin >= P)
    {
        radix_integer_zero(&res->u, radix);
        res->v = 0;
        res->N = P;                 /* known only modulo p^P (P <= Nz) */
        return _radix_padic_finalize_assume_canonical(res, ctx);   /* zero is canonical */
    }

    /* Equal valuations */
    if (vx == vy)
    {
        if (sub)
            radix_integer_sub(&res->u, &x->u, &y->u, radix);
        else
            radix_integer_add(&res->u, &x->u, &y->u, radix);
        res->v = vmin;
        res->N = Nz;
        return _radix_padic_finalize(res, ctx);
    }

    radix_integer_init(T, radix);

    Nresult = Nz;

    if (P != RADIX_PADIC_PREC_INF && vmax >= P)
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
    return _radix_padic_finalize_assume_canonical(res, ctx);
}

int
radix_padic_add(radix_padic_t res, const radix_padic_t x, const radix_padic_t y, gr_ctx_t ctx)
{
    return _radix_padic_add_sub(res, x, y, 0, ctx);
}

int
radix_padic_sub(radix_padic_t res, const radix_padic_t x, const radix_padic_t y, gr_ctx_t ctx)
{
    return _radix_padic_add_sub(res, x, y, 1, ctx);
}

/*
    Reference (unoptimized) multiplication: forms the full product of the units
    and lets the finalizer reduce/truncate. Kept as a correctness oracle for the
    optimized radix_padic_mul below (the test suite cross-checks the two). The
    unit product of two canonical units is itself canonical, so this uses the
    canonical-assuming finalizer.
*/
int
_radix_padic_mul_reference(radix_padic_t res, const radix_padic_t x, const radix_padic_t y, gr_ctx_t ctx)
{
    radix_struct * radix = RADIX_PADIC_CTX_RADIX(ctx);
    slong vx = x->v, vy = y->v, Nx = x->N, Ny = y->N;
    slong vz, Nz;
    radix_integer_t U;

    /* an exact zero annihilates exactly, even against an inexact factor */
    if ((x->u.size == 0 && Nx == RADIX_PADIC_EXACT)
     || (y->u.size == 0 && Ny == RADIX_PADIC_EXACT))
        return radix_padic_zero(res, ctx);

    vz = vx + vy;

    Nz = RADIX_PADIC_EXACT;
    if (Nx != RADIX_PADIC_EXACT)
        Nz = FLINT_MIN(Nz, vy + Nx);
    if (Ny != RADIX_PADIC_EXACT)
        Nz = FLINT_MIN(Nz, vx + Ny);

    radix_integer_init(U, radix);
    radix_integer_mul(U, &x->u, &y->u, radix);

    FLINT_SWAP(radix_integer_struct, res->u, *U);
    res->v = vz;
    res->N = Nz;

    radix_integer_clear(U, radix);

    return _radix_padic_finalize_assume_canonical(res, ctx);
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
    via _radix_padic_finalize_assume_canonical.
*/
int
radix_padic_mul(radix_padic_t res, const radix_padic_t x, const radix_padic_t y, gr_ctx_t ctx)
{
    radix_struct * radix = RADIX_PADIC_CTX_RADIX(ctx);
    slong e = radix->exp;
    slong vx = x->v, vy = y->v, Nx = x->N, Ny = y->N;
    slong vz, Nz, xn, yn, target;

    //return _radix_padic_mul_reference(res, x, y, ctx);

    /* an exact zero annihilates exactly */
    if ((x->u.size == 0 && Nx == RADIX_PADIC_EXACT)
     || (y->u.size == 0 && Ny == RADIX_PADIC_EXACT))
        return radix_padic_zero(res, ctx);

    vz = vx + vy;

    Nz = RADIX_PADIC_EXACT;
    if (Nx != RADIX_PADIC_EXACT)
        Nz = FLINT_MIN(Nz, vy + Nx);
    if (Ny != RADIX_PADIC_EXACT)
        Nz = FLINT_MIN(Nz, vx + Ny);

    /* inexact zero operand: the product unit is zero */
    if (x->u.size == 0 || y->u.size == 0)
    {
        radix_integer_zero(&res->u, radix);
        res->v = vz;
        res->N = Nz;
        return _radix_padic_finalize_assume_canonical(res, ctx);
    }

    xn = FLINT_ABS(x->u.size);
    yn = FLINT_ABS(y->u.size);
    target = _radix_padic_target_prec(vz, ctx);

    /*
        The product is written straight into res->u. All scalar operand fields
        (vx, vy, Nx, Ny, xn, yn) are captured above, and radix_integer_mul /
        radix_integer_mullow_limbs handle res aliasing x or y, so no temporary
        is required.
    */
    if (Nz == RADIX_PADIC_EXACT)
    {
        /* both inputs exact */
        if (target == RADIX_PADIC_PREC_INF)
        {
            radix_integer_mul(&res->u, &x->u, &y->u, radix);   /* exact ring: full product */
            res->N = RADIX_PADIC_EXACT;
        }
        else
        {
            slong digits = target - vz;                    /* exact digits permitted */
            slong maxD = (xn + yn) * e;                    /* upper bound on product digits */
            slong minD = (xn + yn - 2) * e + 1;            /* lower bound on product digits */

            if (digits <= 0)
            {
                radix_integer_zero(&res->u, radix);        /* below the horizon: 0 + O(p^target) */
                res->N = target;
            }
            else if (maxD <= digits)
            {
                radix_integer_mul(&res->u, &x->u, &y->u, radix); /* certainly fits: exact */
                res->N = RADIX_PADIC_EXACT;
            }
            else if (minD > digits)
            {
                slong n = (digits + e - 1) / e;            /* certainly inexact: truncate */
                radix_integer_mullow_limbs(&res->u, &x->u, &y->u, n, radix);
                res->N = target;
            }
            else
            {
                radix_integer_mul(&res->u, &x->u, &y->u, radix); /* borderline: reducer decides */
                res->N = RADIX_PADIC_EXACT;
            }
        }
    }
    else
    {
        /* at least one inexact input: the product is inexact */
        slong Neff = (target == RADIX_PADIC_PREC_INF) ? Nz : FLINT_MIN(Nz, target);
        slong r = Neff - vz;                               /* surviving digits */

        if (r <= 0)
        {
            radix_integer_zero(&res->u, radix);
            res->N = Neff;
        }
        else
        {
            slong n = (r + e - 1) / e;
            radix_integer_mullow_limbs(&res->u, &x->u, &y->u, n, radix);
            res->N = Nz;                                   /* reducer recomputes Neff & truncates */
        }
    }

    res->v = vz;

    return _radix_padic_finalize_assume_canonical(res, ctx);
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
_radix_padic_div(radix_padic_t res, const radix_padic_t a, const radix_padic_t b, gr_ctx_t ctx)
{
    radix_struct * radix = RADIX_PADIC_CTX_RADIX(ctx);
    slong e = radix->exp;
    slong va = a->v, vb = b->v, Na = a->N, Nb = b->N;
    slong prec_abs = RADIX_PADIC_CTX_PREC_ABS(ctx);
    slong prec_rel = RADIX_PADIC_CTX_PREC_REL(ctx);
    slong vq, rel_a, rel_b, rel, relcap, krel, nlimbs;
    radix_integer_t unit;
    truth_t bz;
    int invertible;

    bz = radix_padic_is_zero(b, ctx);
    if (bz == T_TRUE)
        return GR_DOMAIN;                       /* division by zero */
    if (bz == T_UNKNOWN)
        return GR_UNABLE;                        /* divisor not known to be nonzero */

    /* b is definitely nonzero, so its unit is a genuine p-adic unit */

    if (a->u.size == 0 && Na == RADIX_PADIC_EXACT)
        return radix_padic_zero(res, ctx);       /* exact 0 / b = exact 0 */

    vq = va - vb;

    /* exact / exact: return an exact quotient when the units divide over Z.

       TODO: depending on the operand sizes this exactness test may be better
       performed as a Hensel division with remainder (radix_divmod_bn with a
       non-NULL rem) rather than a Euclidean division. */
    if (Na == RADIX_PADIC_EXACT && Nb == RADIX_PADIC_EXACT)
    {
        radix_integer_t qint;
        int divides;

        radix_integer_init(qint, radix);
        divides = radix_integer_div(qint, &a->u, &b->u, radix);

        if (divides)
        {
            FLINT_SWAP(radix_integer_struct, res->u, *qint);
            res->v = vq;
            res->N = RADIX_PADIC_EXACT;
            radix_integer_clear(qint, radix);
            return _radix_padic_finalize_assume_canonical(res, ctx);
        }

        radix_integer_clear(qint, radix);
        /* fall through: infinite expansion, truncate to context precision */
    }

    /* relative precision of the quotient = min of the operand relative
       precisions, then capped by the context precisions */
    rel_a = (Na == RADIX_PADIC_EXACT) ? RADIX_PADIC_PREC_INF : (Na - va);
    rel_b = (Nb == RADIX_PADIC_EXACT) ? RADIX_PADIC_PREC_INF : (Nb - vb);
    rel = FLINT_MIN(rel_a, rel_b);

    relcap = RADIX_PADIC_PREC_INF;
    if (prec_rel != RADIX_PADIC_PREC_INF)
        relcap = prec_rel;
    if (prec_abs != RADIX_PADIC_PREC_INF)
    {
        slong c = prec_abs - vq;                 /* bounded: no overflow */
        relcap = (relcap == RADIX_PADIC_PREC_INF) ? c : FLINT_MIN(relcap, c);
    }

    if (rel == RADIX_PADIC_PREC_INF && relcap == RADIX_PADIC_PREC_INF)
        return GR_UNABLE;                        /* exact ring, infinite expansion */

    if (rel == RADIX_PADIC_PREC_INF)
        krel = relcap;
    else if (relcap == RADIX_PADIC_PREC_INF)
        krel = rel;
    else
        krel = FLINT_MIN(rel, relcap);

    if (krel > RADIX_PADIC_ERR_MAX)
        krel = RADIX_PADIC_ERR_MAX;              /* never compute more than we can store */

    if (krel <= 0)
    {
        /* quotient valuation is at/below the precision horizon */
        radix_integer_zero(&res->u, radix);
        res->v = 0;
        res->N = vq + krel;                      /* bounded; finalize clamps */
        return _radix_padic_finalize_assume_canonical(res, ctx);
    }

    nlimbs = (krel + e - 1) / e;                 /* digits -> limbs (ceil) */

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

    return _radix_padic_finalize_assume_canonical(res, ctx);
}

int
radix_padic_div(radix_padic_t res, const radix_padic_t x, const radix_padic_t y, gr_ctx_t ctx)
{
    return _radix_padic_div(res, x, y, ctx);
}

int
radix_padic_inv(radix_padic_t res, const radix_padic_t x, gr_ctx_t ctx)
{
    radix_struct * radix = RADIX_PADIC_CTX_RADIX(ctx);
    slong e = radix->exp;
    slong vx = x->v, Nx = x->N;
    slong prec_abs = RADIX_PADIC_CTX_PREC_ABS(ctx);
    slong prec_rel = RADIX_PADIC_CTX_PREC_REL(ctx);
    slong vq, rel, relcap, krel, nlimbs;
    truth_t xz;
    int invertible;

    xz = radix_padic_is_zero(x, ctx);
    if (xz == T_TRUE)
        return GR_DOMAIN;                        /* 1 / 0 */
    if (xz == T_UNKNOWN)
        return GR_UNABLE;                        /* x not known to be nonzero */

    /* x is a genuine p-adic unit times p^vx, so 1/x has valuation -vx. */
    vq = -vx;

    /* The only units that invert exactly over Z are +-1; an exact such x gives
       an exact +-p^{-vx}. Any other unit inverts to an infinite expansion,
       handled by truncation to the context precision below. */
    if (Nx == RADIX_PADIC_EXACT && FLINT_ABS(x->u.size) == 1 && x->u.d[0] == 1)
    {
        radix_integer_fit_limbs(&res->u, 1, radix)[0] = 1;
        res->u.size = (x->u.size > 0) ? 1 : -1;
        res->v = vq;
        res->N = RADIX_PADIC_EXACT;
        return _radix_padic_finalize_assume_canonical(res, ctx);
    }

    /* relative precision of 1/x equals that of x, capped by the context
       precisions (mirrors _radix_padic_div with an exact unit numerator). */
    rel = (Nx == RADIX_PADIC_EXACT) ? RADIX_PADIC_PREC_INF : (Nx - vx);

    relcap = RADIX_PADIC_PREC_INF;
    if (prec_rel != RADIX_PADIC_PREC_INF)
        relcap = prec_rel;
    if (prec_abs != RADIX_PADIC_PREC_INF)
    {
        slong c = prec_abs - vq;                 /* bounded: no overflow */
        relcap = (relcap == RADIX_PADIC_PREC_INF) ? c : FLINT_MIN(relcap, c);
    }

    if (rel == RADIX_PADIC_PREC_INF && relcap == RADIX_PADIC_PREC_INF)
        return GR_UNABLE;                        /* exact ring, infinite expansion */

    if (rel == RADIX_PADIC_PREC_INF)
        krel = relcap;
    else if (relcap == RADIX_PADIC_PREC_INF)
        krel = rel;
    else
        krel = FLINT_MIN(rel, relcap);

    if (krel > RADIX_PADIC_ERR_MAX)
        krel = RADIX_PADIC_ERR_MAX;              /* never compute more than we can store */

    if (krel <= 0)
    {
        /* inverse valuation is at/below the precision horizon */
        radix_integer_zero(&res->u, radix);
        res->v = 0;
        res->N = vq + krel;                      /* bounded; finalize clamps */
        return _radix_padic_finalize_assume_canonical(res, ctx);
    }

    nlimbs = (krel + e - 1) / e;                 /* digits -> limbs (ceil) */

    /* res->u = x_unit^{-1} (mod p^{e*nlimbs}); invmod_limbs preserves the sign
       of x, which is also the sign of 1/x. Handles res aliasing x internally. */
    invertible = radix_integer_invmod_limbs(&res->u, &x->u, nlimbs, radix);
    if (!invertible)
        flint_abort();  /* should not happen */

    res->v = vq;
    res->N = vq + krel;                          /* bounded; finalize clamps + truncates */

    return _radix_padic_finalize_assume_canonical(res, ctx);
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
_radix_padic_exact_unit_sqrt(radix_padic_t res, const radix_padic_t x, slong M,
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
radix_padic_sqrt(radix_padic_t res, const radix_padic_t x, gr_ctx_t ctx)
{
    radix_struct * radix = RADIX_PADIC_CTX_RADIX(ctx);
    slong e = radix->exp;
    slong vx = x->v, Nx = x->N;
    slong prec_abs = RADIX_PADIC_CTX_PREC_ABS(ctx);
    slong prec_rel = RADIX_PADIC_CTX_PREC_REL(ctx);
    slong vq, rel, relcap, krel, nlimbs, loss, kin, ndet;
    truth_t xz;
    int issquare, input_exact;
    radix_integer_t unit;

    input_exact = (Nx == RADIX_PADIC_EXACT);

    xz = radix_padic_is_zero(x, ctx);
    if (xz == T_TRUE)
        return radix_padic_zero(res, ctx);       /* sqrt(0) = 0 (exact) */
    if (xz == T_UNKNOWN)
        return GR_UNABLE;                        /* x not known to be nonzero */

    /* A square has even valuation. */
    if (vx & WORD(1))
        return GR_DOMAIN;
    vq = vx / 2;

    /* x = + p^vx exactly (unit 1) has the exact root + p^{vq}. (Other exact
       units have an infinite root expansion, truncated to context precision
       below.) */
    if (Nx == RADIX_PADIC_EXACT && x->u.size == 1 && x->u.d[0] == 1)
    {
        radix_integer_fit_limbs(&res->u, 1, radix)[0] = 1;
        res->u.size = 1;
        res->v = vq;
        res->N = RADIX_PADIC_EXACT;
        return _radix_padic_finalize_assume_canonical(res, ctx);
    }

    /* Relative precision of the root. sqrt preserves relative precision for odd
       p, but the 2-adic square root loses one bit (d/dt t^2 = 2t has valuation 1),
       so an input known to relative precision 'rel' determines the root only to
       'rel - loss'. This loss applies to the input-limited term only: when the
       context's relative cap binds and the input still has bits to spare, the
       capped bits are all determined. Dropping the one unreliable bit is what
       makes sqrt a single precision-independent function (a lower-precision
       result is the truncation of a higher-precision one). */
    loss = (DIGIT_RADIX(radix) == 2) ? 1 : 0;

    rel = (Nx == RADIX_PADIC_EXACT) ? RADIX_PADIC_PREC_INF : (Nx - vx);
    if (rel != RADIX_PADIC_PREC_INF)
        rel = FLINT_MAX(rel - loss, 0);          /* input-limited root precision */

    relcap = RADIX_PADIC_PREC_INF;
    if (prec_rel != RADIX_PADIC_PREC_INF)
        relcap = prec_rel;
    if (prec_abs != RADIX_PADIC_PREC_INF)
    {
        slong c = prec_abs - vq;                 /* bounded: no overflow */
        relcap = (relcap == RADIX_PADIC_PREC_INF) ? c : FLINT_MIN(relcap, c);
    }

    if (rel == RADIX_PADIC_PREC_INF && relcap == RADIX_PADIC_PREC_INF)
        krel = RADIX_PADIC_PREC_INF;             /* exact ring, infinite expansion */
    else if (rel == RADIX_PADIC_PREC_INF)
        krel = relcap;
    else if (relcap == RADIX_PADIC_PREC_INF)
        krel = rel;
    else
        krel = FLINT_MIN(rel, relcap);

    if (krel != RADIX_PADIC_PREC_INF && krel > RADIX_PADIC_ERR_MAX)
        krel = RADIX_PADIC_ERR_MAX;

    /* Input precision to develop: one bit beyond the delivered root precision so
       that the unreliable top bit is computed for headroom and then truncated away
       by the finalizer, leaving every delivered bit pinned. */
    kin = (krel == RADIX_PADIC_PREC_INF) ? RADIX_PADIC_PREC_INF
        : (krel > RADIX_PADIC_ERR_MAX - loss ? RADIX_PADIC_ERR_MAX : krel + loss);

    /* Develop the unit root. If the input is exact the root may terminate (x a
       perfect square), and detecting that needs about half of u's limbs no matter
       what the context precision is -- this is what lets sqrt(4) = 2 exactly even
       in a finite-precision ring. Otherwise develop kin = krel + loss relative
       bits. Non-squareness is read from radix_integer_sqrtmod_limbs, which finds
       the root modulo p itself, so it is never recomputed here. */
    ndet = input_exact ? (FLINT_ABS(x->u.size) / 2 + 2) : 0;

    if (krel == RADIX_PADIC_PREC_INF)
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
    radix_integer_mod_digits_nonneg(unit, &x->u, e * nlimbs, radix);
    issquare = radix_integer_sqrtmod_limbs(&res->u, unit, nlimbs, radix);
    radix_integer_clear(unit, radix);

    if (!issquare)
        return GR_DOMAIN;                            /* unit is not a square mod p */

    /* Exact input + perfect-square unit: return the exact (terminating) root, with
       the sign matching the canonical branch so the exact result agrees with the
       finite-precision result on every shared digit. Works in any context. */
    if (input_exact &&
        _radix_padic_exact_unit_sqrt(res, x, (ndet > 0 ? ndet : nlimbs), radix))
    {
        res->v = vq;
        res->N = RADIX_PADIC_EXACT;
        return _radix_padic_finalize_assume_canonical(res, ctx);
    }

    if (krel == RADIX_PADIC_PREC_INF)
        return GR_UNABLE;       /* exact ring: an irrational root cannot be exact */

    if (krel <= 0)
    {
        /* root is zero to the tracked precision */
        radix_integer_zero(&res->u, radix);
        res->v = 0;
        res->N = vq + krel;                          /* bounded; finalize clamps */
        return _radix_padic_finalize_assume_canonical(res, ctx);
    }

    res->v = vq;
    res->N = vq + krel;                              /* bounded; finalize clamps + truncates */

    return _radix_padic_finalize_assume_canonical(res, ctx);
}

/*
    Whether x is a square. A square has even valuation and a unit that is a square
    in Z_p: for odd p this is decided by the single low digit (a quadratic residue
    modulo p, which always lifts), and that digit is always known; for p = 2 it
    needs the unit modulo 8, so a relative precision below 3 is undecided. A
    definite zero is a square; an inexact (unknown) zero is undecided.
*/
truth_t
radix_padic_is_square(const radix_padic_t x, gr_ctx_t ctx)
{
    radix_struct * radix = RADIX_PADIC_CTX_RADIX(ctx);
    ulong p = DIGIT_RADIX(radix);
    truth_t xz;

    xz = radix_padic_is_zero(x, ctx);
    if (xz == T_TRUE)
        return T_TRUE;                  /* 0 = 0^2 */
    if (xz == T_UNKNOWN)
        return T_UNKNOWN;               /* inexact zero: cannot tell */

    if (x->v & WORD(1))
        return T_FALSE;                 /* a square has even valuation */

    if (p == 2)
    {
        slong rel = (x->N == RADIX_PADIC_EXACT)
            ? RADIX_PADIC_PREC_INF : (x->N - x->v);
        if (rel != RADIX_PADIC_PREC_INF && rel < 3)
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
radix_padic_rsqrt(radix_padic_t res, const radix_padic_t x, gr_ctx_t ctx)
{
    radix_struct * radix = RADIX_PADIC_CTX_RADIX(ctx);
    slong e = radix->exp;
    slong vx = x->v, Nx = x->N;
    slong prec_abs = RADIX_PADIC_CTX_PREC_ABS(ctx);
    slong prec_rel = RADIX_PADIC_CTX_PREC_REL(ctx);
    slong vq, vr, rel, relcap, krel, nlimbs, loss, kin;
    truth_t xz;
    int issquare;
    radix_integer_t unit;

    xz = radix_padic_is_zero(x, ctx);
    if (xz == T_TRUE)
        return GR_DOMAIN;                        /* 1/sqrt(0) is undefined */
    if (xz == T_UNKNOWN)
        return GR_UNABLE;                        /* x not known to be nonzero */

    if (vx & WORD(1))
        return GR_DOMAIN;                        /* a square has even valuation */
    vq = vx / 2;
    vr = -vq;                                    /* valuation of x^{-1/2} */

    /* x = + p^vx exactly (unit 1) has the exact reciprocal root + p^{-vq}. */
    if (Nx == RADIX_PADIC_EXACT && x->u.size == 1 && x->u.d[0] == 1)
    {
        radix_integer_fit_limbs(&res->u, 1, radix)[0] = 1;
        res->u.size = 1;
        res->v = vr;
        res->N = RADIX_PADIC_EXACT;
        return _radix_padic_finalize_assume_canonical(res, ctx);
    }

    loss = (DIGIT_RADIX(radix) == 2) ? 1 : 0;

    rel = (Nx == RADIX_PADIC_EXACT) ? RADIX_PADIC_PREC_INF : (Nx - vx);
    if (rel != RADIX_PADIC_PREC_INF)
        rel = FLINT_MAX(rel - loss, 0);

    relcap = RADIX_PADIC_PREC_INF;
    if (prec_rel != RADIX_PADIC_PREC_INF)
        relcap = prec_rel;
    if (prec_abs != RADIX_PADIC_PREC_INF)
    {
        slong c = prec_abs - vr;                 /* absolute cap relative to vr = -vq */
        relcap = (relcap == RADIX_PADIC_PREC_INF) ? c : FLINT_MIN(relcap, c);
    }

    if (rel == RADIX_PADIC_PREC_INF && relcap == RADIX_PADIC_PREC_INF)
        krel = RADIX_PADIC_PREC_INF;             /* exact ring */
    else if (rel == RADIX_PADIC_PREC_INF)
        krel = relcap;
    else if (relcap == RADIX_PADIC_PREC_INF)
        krel = rel;
    else
        krel = FLINT_MIN(rel, relcap);

    if (krel != RADIX_PADIC_PREC_INF && krel > RADIX_PADIC_ERR_MAX)
        krel = RADIX_PADIC_ERR_MAX;

    kin = (krel == RADIX_PADIC_PREC_INF) ? RADIX_PADIC_PREC_INF
        : (krel > RADIX_PADIC_ERR_MAX - loss ? RADIX_PADIC_ERR_MAX : krel + loss);

    /* Develop the reciprocal unit root. In the exact ring the result is infinite
       (the unit-1 case is already returned above), so one limb is enough to decide
       squareness; otherwise develop kin = krel + loss relative bits. Non-squareness
       is read directly from radix_integer_rsqrtmod_limbs. */
    if (krel == RADIX_PADIC_PREC_INF || kin <= 0)
        nlimbs = 1;
    else
        nlimbs = (kin + e - 1) / e;

    radix_integer_init(unit, radix);
    radix_integer_mod_digits_nonneg(unit, &x->u, e * nlimbs, radix);
    issquare = radix_integer_rsqrtmod_limbs(&res->u, unit, nlimbs, radix);
    radix_integer_clear(unit, radix);

    if (!issquare)
        return GR_DOMAIN;                        /* unit is not a square mod p */

    if (krel == RADIX_PADIC_PREC_INF)
        return GR_UNABLE;       /* exact ring: an infinite reciprocal root */

    if (krel <= 0)
    {
        radix_integer_zero(&res->u, radix);
        res->v = 0;
        res->N = vr + krel;
        return _radix_padic_finalize_assume_canonical(res, ctx);
    }

    res->v = vr;
    res->N = vr + krel;
    return _radix_padic_finalize_assume_canonical(res, ctx);
}

/* ------------------------------------------------------------------------- */
/*    Predicates                                                             */
/* ------------------------------------------------------------------------- */

truth_t
radix_padic_is_zero(const radix_padic_t x, gr_ctx_t ctx)
{
    if (x->N == RADIX_PADIC_EXACT)
        return (x->u.size == 0) ? T_TRUE : T_FALSE;

    /* inexact: x = u p^v + O(p^N) */
    if (x->u.size == 0)
        return T_UNKNOWN;            /* 0 + O(p^N): might be zero or not */

    /* nonzero unit: its lowest (nonzero) digit sits at position v */
    return (x->N > x->v) ? T_FALSE : T_UNKNOWN;
}

truth_t
radix_padic_is_one(const radix_padic_t x, gr_ctx_t ctx)
{
    radix_struct * radix = RADIX_PADIC_CTX_RADIX(ctx);

    if (x->N == RADIX_PADIC_EXACT)
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
radix_padic_is_neg_one(const radix_padic_t x, gr_ctx_t ctx)
{
    radix_struct * radix = RADIX_PADIC_CTX_RADIX(ctx);

    if (x->N == RADIX_PADIC_EXACT)
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
radix_padic_equal(const radix_padic_t x, const radix_padic_t y, gr_ctx_t ctx)
{
    radix_struct * radix = RADIX_PADIC_CTX_RADIX(ctx);
    slong err, vx = x->v, vy = y->v;
    int xz, yz, eq;

    if (x->N == RADIX_PADIC_EXACT && y->N == RADIX_PADIC_EXACT)
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
radix_padic_is_invertible(const radix_padic_t x, gr_ctx_t ctx)
{
    /* the radix p-adic numbers form a field: x is invertible iff x is nonzero */
    return truth_not(radix_padic_is_zero(x, ctx));
}

/* ------------------------------------------------------------------------- */
/*    Randomisation                                                          */
/* ------------------------------------------------------------------------- */

int
radix_padic_randtest(radix_padic_t res, flint_rand_t state, gr_ctx_t ctx)
{
    radix_struct * radix = RADIX_PADIC_CTX_RADIX(ctx);
    int ctx_exact = (RADIX_PADIC_CTX_PREC_ABS(ctx) == RADIX_PADIC_PREC_INF
                  && RADIX_PADIC_CTX_PREC_REL(ctx) == RADIX_PADIC_PREC_INF);

    /* the limits mode is only meaningful for an inexact ring (it produces
       error terms and over-large valuations that an exact ring cannot hold) */
    int limits = (RADIX_PADIC_CTX_FLAGS(ctx) & RADIX_PADIC_TEST_LIMITS)
              && RADIX_PADIC_CTX_PREC_REL(ctx) != RADIX_PADIC_PREC_INF
              && (n_randint(state, 2) == 0);

    if (n_randint(state, 2))
        radix_integer_set_si(&res->u, n_randint(state, 5) - 2, radix);
    else
        radix_integer_randtest_limbs(&res->u, state, 1 + n_randint(state, 10), radix);

    if (!RADIX_PADIC_CTX_SIGNED(ctx))
        radix_integer_abs(&res->u, &res->u, radix);

    if (limits && (n_randint(state, 4) == 0))
    {
        /* valuation and precision at or just inside the representable range, so
           that subsequent operations push the intermediate exponents over the
           ERR_MAX boundary and exercise the clamping in _radix_padic_finalize */
        slong big = RADIX_PADIC_ERR_MAX - n_randint(state, 4);

        res->v = (n_randlimb(state) & 1) ? big : -big;

        switch (n_randint(state, 3))
        {
            case 0:  res->N = RADIX_PADIC_EXACT; break;
            case 1:  res->N = (n_randlimb(state) & 1) ? big : -big; break;
            default: res->N = res->v + n_randint(state, 6); break;
        }
    }
    else
    {
        res->v = n_randint(state, 4);
        res->N = RADIX_PADIC_EXACT;

        if (!ctx_exact && (n_randlimb(state) & 1))
            res->N = res->v + n_randint(state, 6);   /* random finite absolute precision */
    }

    if (_radix_padic_finalize(res, ctx) != GR_SUCCESS)
        radix_padic_zero(res, ctx);                  /* always return a valid element */

    return GR_SUCCESS;
}

/* ------------------------------------------------------------------------- */
/*    Output                                                                 */
/* ------------------------------------------------------------------------- */

int
radix_padic_write(gr_stream_t out, const radix_padic_t x, gr_ctx_t ctx)
{
    radix_struct * radix = RADIX_PADIC_CTX_RADIX(ctx);
    ulong p = GR_RADIX_PADIC_CTX(ctx)->p;
    int status = GR_SUCCESS;

    if (x->u.size == 0)
    {
        status |= gr_stream_write(out, "0");
    }
    else
    {
        slong size = x->u.size;
        char * str;

        if (RADIX_PADIC_CTX_DECIMAL(ctx))
            str = radix_get_str_decimal(NULL, x->u.d, FLINT_ABS(size), size < 0, radix);
        else
            str = radix_get_str_sum(NULL, x->u.d, FLINT_ABS(size), size < 0, 1, radix);

        status |= gr_stream_write(out, "(");
        status |= gr_stream_write_free(out, str);
        status |= gr_stream_write(out, ")");

        if (x->v != 0)
        {
            status |= gr_stream_write(out, " * ");
            status |= gr_stream_write_ui(out, p);
            status |= gr_stream_write(out, "^");
            status |= gr_stream_write_si(out, x->v);
        }
    }

    if (x->N != RADIX_PADIC_EXACT)
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
radix_padic_ctx_clear(gr_ctx_t ctx)
{
    radix_clear(RADIX_PADIC_CTX_RADIX(ctx));
    flint_free(GR_RADIX_PADIC_CTX(ctx));
}

int
radix_padic_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    status |= gr_stream_write(out, "Radix ");
    status |= gr_stream_write_ui(out, GR_RADIX_PADIC_CTX(ctx)->p);
    status |= gr_stream_write(out, "-adic numbers");

    status |= gr_stream_write(out, " (rel prec ");
    if (RADIX_PADIC_CTX_PREC_REL(ctx) == RADIX_PADIC_PREC_INF)
        status |= gr_stream_write(out, "inf");
    else
        status |= gr_stream_write_si(out, RADIX_PADIC_CTX_PREC_REL(ctx));

    status |= gr_stream_write(out, ", abs prec ");
    if (RADIX_PADIC_CTX_PREC_ABS(ctx) == RADIX_PADIC_PREC_INF)
        status |= gr_stream_write(out, "inf");
    else
        status |= gr_stream_write_si(out, RADIX_PADIC_CTX_PREC_ABS(ctx));

    status |= gr_stream_write(out, ")");
    return status;
}

static truth_t
radix_padic_ctx_is_exact(gr_ctx_t ctx)
{
    return (RADIX_PADIC_CTX_PREC_ABS(ctx) == RADIX_PADIC_PREC_INF
         && RADIX_PADIC_CTX_PREC_REL(ctx) == RADIX_PADIC_PREC_INF) ? T_TRUE : T_FALSE;
}

int _radix_padic_methods_initialized = 0;

gr_static_method_table _radix_padic_methods;

gr_method_tab_input _radix_padic_methods_input[] =
{
    {GR_METHOD_CTX_CLEAR,       (gr_funcptr) radix_padic_ctx_clear},
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) radix_padic_ctx_write},
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
    {GR_METHOD_CTX_IS_EXACT,    (gr_funcptr) radix_padic_ctx_is_exact},
    {GR_METHOD_INIT,            (gr_funcptr) radix_padic_init},
    {GR_METHOD_CLEAR,           (gr_funcptr) radix_padic_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) radix_padic_swap},
    {GR_METHOD_SET_SHALLOW,     (gr_funcptr) radix_padic_set_shallow},
    {GR_METHOD_RANDTEST,        (gr_funcptr) radix_padic_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) radix_padic_write},
    {GR_METHOD_ZERO,            (gr_funcptr) radix_padic_zero},
    {GR_METHOD_ONE,             (gr_funcptr) radix_padic_one},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) radix_padic_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) radix_padic_is_one},
    {GR_METHOD_IS_NEG_ONE,      (gr_funcptr) radix_padic_is_neg_one},
    {GR_METHOD_IS_INVERTIBLE,   (gr_funcptr) radix_padic_is_invertible},
    {GR_METHOD_EQUAL,           (gr_funcptr) radix_padic_equal},
    {GR_METHOD_SET,             (gr_funcptr) radix_padic_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) radix_padic_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) radix_padic_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) radix_padic_set_fmpz},
    {GR_METHOD_GET_FMPZ,        (gr_funcptr) radix_padic_get_fmpz},
    {GR_METHOD_NEG,             (gr_funcptr) radix_padic_neg},
    {GR_METHOD_ADD,             (gr_funcptr) radix_padic_add},
    {GR_METHOD_SUB,             (gr_funcptr) radix_padic_sub},
    {GR_METHOD_MUL,             (gr_funcptr) radix_padic_mul},
    {GR_METHOD_INV,             (gr_funcptr) radix_padic_inv},
    {GR_METHOD_DIV,             (gr_funcptr) radix_padic_div},
    {GR_METHOD_SQRT,            (gr_funcptr) radix_padic_sqrt},
    {GR_METHOD_RSQRT,           (gr_funcptr) radix_padic_rsqrt},
    {GR_METHOD_IS_SQUARE,       (gr_funcptr) radix_padic_is_square},
    {0,                         (gr_funcptr) NULL},
};

int
gr_ctx_init_radix_padic(gr_ctx_t ctx, ulong p, slong prec_rel, slong prec_abs, int flags)
{
    radix_padic_ctx_struct * pctx;

    FLINT_ASSERT(n_is_prime(p));

    ctx->which_ring = GR_CTX_RADIX_PADIC;
    ctx->sizeof_elem = sizeof(radix_padic_struct);
    ctx->size_limit = WORD_MAX;

    GR_CTX_DATA_AS_PTR(ctx) = flint_malloc(sizeof(radix_padic_ctx_struct));
    pctx = GR_RADIX_PADIC_CTX(ctx);

    radix_init(&pctx->radix, p, 0);   /* limb radix p^e with e chosen automatically */
    pctx->p = p;

    if (prec_rel == RADIX_PADIC_PREC_INF)
        pctx->prec_rel = RADIX_PADIC_PREC_INF;
    else
        pctx->prec_rel = FLINT_MIN(FLINT_MAX(prec_rel, 0), RADIX_PADIC_ERR_MAX);

    if (prec_abs == RADIX_PADIC_PREC_INF)
        pctx->prec_abs = RADIX_PADIC_PREC_INF;
    else
        pctx->prec_abs = FLINT_MAX(FLINT_MIN(prec_abs, RADIX_PADIC_ERR_MAX),
                                   -RADIX_PADIC_ERR_MAX);

    pctx->flags = flags;

    ctx->methods = _radix_padic_methods;

    if (!_radix_padic_methods_initialized)
    {
        gr_method_tab_init(_radix_padic_methods, _radix_padic_methods_input);
        _radix_padic_methods_initialized = 1;
    }

    return GR_SUCCESS;
}

