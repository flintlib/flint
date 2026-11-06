/*
    Copyright (C) 2018 arbguest
    Copyright (C) 2018, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <float.h>
#include <math.h>
#include "arb.h"
#include "arb_mat.h"
#include "gr.h"
#include "gr_mat.h"
#include "nfloat.h"
#include "dfloat.h"
#include "arb_mat/impl.h"

/* The approximate LU decomposition and solving functions of arb_mat are
   the generic algorithms of gr_mat over a floating-point ring holding the
   midpoints: dfloat or nfloat when they are faster and the exponents
   fit, arf otherwise or when those fail (see impl.h). */

/* ---- the policy and the conversions, shared with acb_mat ---- */

/* The size cutoffs for dfloat come from timings of the inverse against
   nfloat at the next multiple of 64 bits, where nfloat wins for small
   matrices with two or more doubles. */

#define NFLOAT_MAX_PREC (NFLOAT_MAX_LIMBS * FLINT_BITS)

int
_arb_mat_approx_backend(slong n, slong prec, slong range, int FLINT_UNUSED(cplx), int op)
{
    int N;

    if (n < ARB_MAT_APPROX_MIN_N || range < 0)
        return ARB_MAT_APPROX_ARF;

    if (dfloat_is_supported())
    {
        for (N = 1; N <= DFLOAT_MAX_N; N++)
            if (prec + ARB_MAT_APPROX_DFLOAT_PAD(N) <= 53 * N)
                break;

        if (N <= DFLOAT_MAX_N && range <= ARB_MAT_APPROX_DFLOAT_RANGE(N)
            && (op == ARB_MAT_APPROX_OP_EIG || N == 1 || n >= 8 * N))
            return N;
    }

    if (prec <= NFLOAT_MAX_PREC)
        return ARB_MAT_APPROX_NFLOAT;

    return ARB_MAT_APPROX_ARF;
}

int
_arb_mat_approx_backend_next(int kind, slong prec)
{
    if (kind > 0 && prec <= NFLOAT_MAX_PREC)
        return ARB_MAT_APPROX_NFLOAT;

    if (kind == ARB_MAT_APPROX_ARF)
        return ARB_MAT_APPROX_NONE;

    return ARB_MAT_APPROX_ARF;
}

int
_arb_mat_approx_ctx_init(gr_ctx_t ctx, int kind, slong prec, int cplx)
{
    if (kind > 0)
        return gr_ctx_init_dfloat(ctx, kind, cplx ? DFLOAT_COMPLEX : 0);

    if (kind == ARB_MAT_APPROX_ARF)
    {
        /* as dfloat and nfloat do for small prec, compute with at least
           53 bits: at a handful of bits the iterative algorithms (the
           QR iteration) give meaningless results */
        prec = FLINT_MAX(prec, 53);

        if (cplx)
            gr_ctx_init_complex_float_acf(ctx, prec);
        else
            gr_ctx_init_real_float_arf(ctx, prec);
        return GR_SUCCESS;
    }

    prec = FLINT_BITS * ((prec + FLINT_BITS - 1) / FLINT_BITS);

    return cplx ? nfloat_complex_ctx_init(ctx, prec, 0)
                : nfloat_ctx_init(ctx, prec, 0);
}

int
_arb_mat_approx_exp_update(slong * emin, slong * emax, const arf_t x)
{
    slong e;

    if (arf_is_zero(x))
        return 1;

    if (arf_is_special(x) || COEFF_IS_MPZ(*ARF_EXPREF(x)))
        return 0;

    e = *ARF_EXPREF(x);

    if (e > ARB_MAT_APPROX_MAX_EXP || e < -ARB_MAT_APPROX_MAX_EXP)
        return 0;

    *emin = FLINT_MIN(*emin, e);
    *emax = FLINT_MAX(*emax, e);
    return 1;
}

slong
_arb_mat_approx_range(slong emin, slong emax)
{
    return (emax < emin) ? 0 : emax - emin;
}

int
_arb_mat_approx_set_arf(gr_ptr res, const arf_t x, slong scale, int kind, gr_ctx_t ctx)
{
    arf_struct t;

    if (kind == ARB_MAT_APPROX_ARF)
    {
        arf_set(res, x);
        return GR_SUCCESS;
    }

    /* a shallow copy with the exponent changed; the exponent is small
       (checked by _arb_mat_approx_exp_update) */
    t = *x;
    if (scale != 0 && !arf_is_special(x))
        *ARF_EXPREF(&t) += scale;

    if (kind > 0)
    {
        _dfloat_set_arf((double *) res, kind, NULL, &t);
        return (fabs(((double *) res)[0]) <= DBL_MAX) ? GR_SUCCESS : GR_UNABLE;
    }
    else
    {
        return nfloat_set_arf(res, &t, ctx);
    }
}

int
_arb_mat_approx_is_finite(gr_srcptr x, int kind, gr_ctx_t FLINT_UNUSED(ctx))
{
    int i;

    if (kind > 0)
    {
        for (i = 0; i < kind; i++)
            if (!(fabs(((const double *) x)[i]) <= DBL_MAX))
                return 0;
        return 1;
    }
    else if (kind == ARB_MAT_APPROX_NFLOAT)
    {
        return !NFLOAT_IS_SPECIAL(x) || NFLOAT_IS_ZERO(x);
    }
    else
    {
        /* arf has no exponent limits; nonfinite values are passed on */
        return 1;
    }
}

int
_arb_mat_approx_get_arf(arf_t res, gr_srcptr x, slong scale, int kind, gr_ctx_t ctx)
{
    if (kind == ARB_MAT_APPROX_ARF)
    {
        /* moved out: the backend element is discarded afterwards */
        arf_swap(res, (arf_ptr) x);
        arf_mul_2exp_si(res, res, scale);
        return 1;
    }

    if (!_arb_mat_approx_is_finite(x, kind, ctx))
        return 0;

    if (kind > 0)
        _dfloat_get_arf(res, (const double *) x, kind);
    else if (nfloat_get_arf(res, x, ctx) != GR_SUCCESS)
        return 0;

    arf_mul_2exp_si(res, res, scale);
    return 1;
}

gr_ptr
_arb_mat_approx_imag_ptr(gr_ptr x, int kind, gr_ctx_t ctx)
{
    if (kind > 0)
        return (double *) x + kind;
    else if (kind == ARB_MAT_APPROX_NFLOAT)
        return NFLOAT_COMPLEX_IM(x, ctx);
    else
        return (arf_ptr) x + 1;     /* acf_imagref */
}

/* ---- real matrices ---- */

static int
_exp_range(slong * emin, slong * emax, const arb_mat_t A, int part)
{
    slong i, j;

    *emin = WORD_MAX;
    *emax = WORD_MIN;

    for (i = 0; i < arb_mat_nrows(A); i++)
        for (j = 0; j < arb_mat_ncols(A); j++)
            if (_arb_mat_approx_in_part(i, j, part))
                if (!_arb_mat_approx_exp_update(emin, emax,
                        arb_midref(arb_mat_entry(A, i, j))))
                    return 0;

    return 1;
}

/* the scaling exponent bringing the largest entry into [1/2, 1) */
static slong
_scale(slong emax)
{
    return (emax == WORD_MIN) ? 0 : -emax;
}

static int
_get_gr(gr_mat_t res, const arb_mat_t A, slong scale, int part, int kind, gr_ctx_t ctx)
{
    slong i, j, sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    for (i = 0; i < arb_mat_nrows(A) && status == GR_SUCCESS; i++)
        for (j = 0; j < arb_mat_ncols(A); j++)
            if (_arb_mat_approx_in_part(i, j, part))
                status |= _arb_mat_approx_set_arf(GR_MAT_ENTRY(res, i, j, sz),
                    arb_midref(arb_mat_entry(A, i, j)), scale, kind, ctx);

    return status;
}

/* An input matrix for the ring: the converted midpoints of the entries
   in part1 (scaled by 2^scale1) and part2 (scaled by 2^scale2; -1 for
   none), zero elsewhere. For arf, which needs no conversion, the entries
   are shallow copies of the midpoints, so that the matrix must only be
   read and must be freed with _input_clear. */
static int
_input(gr_mat_t res, const arb_mat_t A, int part1, slong scale1, int part2, slong scale2, int kind, gr_ctx_t ctx)
{
    slong i, j, r = arb_mat_nrows(A), c = arb_mat_ncols(A);
    int status = GR_SUCCESS;

    if (kind == ARB_MAT_APPROX_ARF)
    {
        res->entries = flint_calloc(r * c + 1, sizeof(arf_struct));
        res->r = r;
        res->c = c;
        res->stride = c;

        for (i = 0; i < r; i++)
            for (j = 0; j < c; j++)
                if (_arb_mat_approx_in_part(i, j, part1)
                    || (part2 != -1 && _arb_mat_approx_in_part(i, j, part2)))
                    ((arf_ptr) res->entries)[i * c + j] = *arb_midref(arb_mat_entry(A, i, j));

        return GR_SUCCESS;
    }

    gr_mat_init(res, r, c, ctx);
    status |= _get_gr(res, A, scale1, part1, kind, ctx);
    if (part2 != -1)
        status |= _get_gr(res, A, scale2, part2, kind, ctx);
    return status;
}

static void
_input_clear(gr_mat_t res, int kind, gr_ctx_t ctx)
{
    if (kind == ARB_MAT_APPROX_ARF)
        flint_free(res->entries);
    else
        gr_mat_clear(res, ctx);
}

static int
_is_finite(const gr_mat_t A, int part, int kind, gr_ctx_t ctx)
{
    slong i, j, sz = ctx->sizeof_elem;

    if (kind == ARB_MAT_APPROX_ARF)
        return 1;

    for (i = 0; i < A->r; i++)
        for (j = 0; j < A->c; j++)
            if (_arb_mat_approx_in_part(i, j, part))
                if (!_arb_mat_approx_is_finite(GR_MAT_ENTRY(A, i, j, sz), kind, ctx))
                    return 0;

    return 1;
}

/* Write the midpoints of res (with zero radii), scaling the entries in
   the given part by 2^scale1 and the others by 2^scale2. Only called once
   the result is known to be finite. */
static void
_set_gr(arb_mat_t res, const gr_mat_t A, int part, slong scale1, slong scale2, int kind, gr_ctx_t ctx)
{
    slong i, j, sz = ctx->sizeof_elem;
    arb_ptr t;

    for (i = 0; i < A->r; i++)
    {
        for (j = 0; j < A->c; j++)
        {
            t = arb_mat_entry(res, i, j);
            _arb_mat_approx_get_arf(arb_midref(t), GR_MAT_ENTRY(A, i, j, sz),
                _arb_mat_approx_in_part(i, j, part) ? scale1 : scale2, kind, ctx);
            mag_zero(arb_radref(t));
        }
    }
}

static void
_set_nan(arb_mat_t X)
{
    slong i, j;

    for (i = 0; i < arb_mat_nrows(X); i++)
        for (j = 0; j < arb_mat_ncols(X); j++)
        {
            arf_nan(arb_midref(arb_mat_entry(X, i, j)));
            mag_zero(arb_radref(arb_mat_entry(X, i, j)));
        }
}

int
arb_mat_approx_lu(slong * P, arb_mat_t LU, const arb_mat_t A, slong prec)
{
    slong m = arb_mat_nrows(A), n = arb_mat_ncols(A), emin, emax, rank, s;
    int kind, status, result = 0, done = 0;
    gr_ctx_t ctx;
    gr_mat_t GA, T;

    if (m == 0 || n == 0)
        return 1;

    if (_exp_range(&emin, &emax, A, ARB_MAT_APPROX_ALL) && emax != WORD_MIN)
    {
        kind = _arb_mat_approx_backend(FLINT_MIN(m, n), prec,
            _arb_mat_approx_range(emin, emax), 0, ARB_MAT_APPROX_OP_LU);
        s = _scale(emax);
    }
    else
    {
        kind = ARB_MAT_APPROX_ARF;
        s = 0;
    }

    for ( ; !done && kind != ARB_MAT_APPROX_NONE; kind = _arb_mat_approx_backend_next(kind, prec))
    {
        if (kind == ARB_MAT_APPROX_ARF)
            s = 0;

        if (_arb_mat_approx_ctx_init(ctx, kind, prec, 0) != GR_SUCCESS)
            continue;

        status = _input(GA, A, ARB_MAT_APPROX_ALL, s, -1, 0, kind, ctx);
        gr_mat_init(T, m, n, ctx);
        if (status == GR_SUCCESS)
            status = gr_mat_lu(&rank, P, T, GA, m == n, ctx);

        /* U is scaled, L is not */
        if (status == GR_SUCCESS && rank == FLINT_MIN(m, n)
            && _is_finite(T, ARB_MAT_APPROX_ALL, kind, ctx))
        {
            _set_gr(LU, T, ARB_MAT_APPROX_UPPER, -s, 0, kind, ctx);
            result = 1;
            done = 1;
        }
        else if (kind == ARB_MAT_APPROX_ARF)
        {
            done = 1;
        }

        _input_clear(GA, kind, ctx);
        gr_mat_clear(T, ctx);
        gr_ctx_clear(ctx);
    }

    return result;
}

/* Solve T X = B with T of the given part (ARB_MAT_APPROX_ALL for a
   general matrix, or LU data for lu_precomp, or triangular), scaled
   as follows: for a unit triangular T (and the L part of LU data) the
   off-diagonal entries are not scaled, since T contains the ones on the
   diagonal; otherwise T (or U) is scaled to bring its largest entry
   into [1/2, 1). B is scaled likewise. Returns 0 if the system cannot
   be solved even with arf (a singular matrix, or a zero on the diagonal
   of a triangular one); then X is not written for a general system and
   set to NaN otherwise. */

#define SOLVE_GENERAL 0
#define SOLVE_LU_PRECOMP 1
#define SOLVE_TRIL 2
#define SOLVE_TRIU 3

static int
_solve(arb_mat_t X, const slong * perm, const arb_mat_t T,
    const arb_mat_t B, int which, int unit, slong prec)
{
    slong n = arb_mat_nrows(T), m = arb_mat_ncols(B);
    slong emin, emax, uemin, uemax, bemin, bemax, range, sT, sB;
    int kind, status, op, upart, lpart, ok, result = 0, done = 0;
    gr_ctx_t ctx;
    gr_mat_t GT, GB, GX;

    if (n == 0 || m == 0)
        return 1;

    op = (which == SOLVE_GENERAL) ? ARB_MAT_APPROX_OP_LU : ARB_MAT_APPROX_OP_TRI;

    /* lpart: unscaled part (unit), upart: scaled part */
    if (which == SOLVE_GENERAL)
    {
        lpart = -1;
        upart = ARB_MAT_APPROX_ALL;
    }
    else if (which == SOLVE_LU_PRECOMP)
    {
        lpart = ARB_MAT_APPROX_STRICT_LOWER;
        upart = ARB_MAT_APPROX_UPPER;
    }
    else if (unit)
    {
        lpart = (which == SOLVE_TRIL) ? ARB_MAT_APPROX_STRICT_LOWER : ARB_MAT_APPROX_STRICT_UPPER;
        upart = -1;
    }
    else
    {
        lpart = -1;
        upart = (which == SOLVE_TRIL) ? ARB_MAT_APPROX_LOWER : ARB_MAT_APPROX_UPPER;
    }

    range = 0;
    sT = sB = 0;
    ok = !(op == ARB_MAT_APPROX_OP_TRI && m < ARB_MAT_APPROX_TRI_MIN_RHS);

    if (ok && lpart != -1)
    {
        ok = _exp_range(&emin, &emax, T, lpart);
        /* include the unit diagonal (exponent 1) */
        if (ok)
            range = FLINT_MAX(range, _arb_mat_approx_range(FLINT_MIN(emin, 1), FLINT_MAX(emax, 1)));
    }

    if (ok && upart != -1)
    {
        /* a zero matrix (singular) is left to arf */
        ok = _exp_range(&uemin, &uemax, T, upart) && uemax != WORD_MIN;
        if (ok)
        {
            range = FLINT_MAX(range, _arb_mat_approx_range(uemin, uemax));
            sT = _scale(uemax);
        }
    }

    if (ok)
    {
        ok = _exp_range(&bemin, &bemax, B, ARB_MAT_APPROX_ALL);
        if (ok)
        {
            range = FLINT_MAX(range, _arb_mat_approx_range(bemin, bemax));
            sB = _scale(bemax);
        }
    }

    kind = ok ? _arb_mat_approx_backend(n, prec, range, 0, op) : ARB_MAT_APPROX_ARF;

    for ( ; !done && kind != ARB_MAT_APPROX_NONE; kind = _arb_mat_approx_backend_next(kind, prec))
    {
        if (kind == ARB_MAT_APPROX_ARF)
            sT = sB = 0;

        if (_arb_mat_approx_ctx_init(ctx, kind, prec, 0) != GR_SUCCESS)
            continue;

        if (lpart != -1)
            status = _input(GT, T, lpart, 0, upart, sT, kind, ctx);
        else
            status = _input(GT, T, upart, sT, -1, 0, kind, ctx);
        status |= _input(GB, B, ARB_MAT_APPROX_ALL, sB, -1, 0, kind, ctx);
        gr_mat_init(GX, n, m, ctx);

        if (status == GR_SUCCESS)
        {
            if (which == SOLVE_GENERAL)
                status = gr_mat_nonsingular_solve(GX, GT, GB, ctx);
            else if (which == SOLVE_LU_PRECOMP)
                status = gr_mat_nonsingular_solve_lu_precomp(GX, perm, GT, GB, ctx);
            else if (which == SOLVE_TRIL)
                status = gr_mat_nonsingular_solve_tril(GX, GT, GB, unit, ctx);
            else
                status = gr_mat_nonsingular_solve_triu(GX, GT, GB, unit, ctx);
        }

        /* the inputs may be shallow copies of entries of X */
        _input_clear(GT, kind, ctx);
        _input_clear(GB, kind, ctx);

        /* T' = 2^sT T, B' = 2^sB B, so X = 2^(sT - sB) X' */
        if (status == GR_SUCCESS && _is_finite(GX, ARB_MAT_APPROX_ALL, kind, ctx))
        {
            _set_gr(X, GX, ARB_MAT_APPROX_ALL, sT - sB, 0, kind, ctx);
            result = 1;
            done = 1;
        }
        else if (kind == ARB_MAT_APPROX_ARF)
        {
            if (which != SOLVE_GENERAL)
                _set_nan(X);
            done = 1;
        }

        gr_mat_clear(GX, ctx);
        gr_ctx_clear(ctx);
    }

    return result;
}

int
arb_mat_approx_solve(arb_mat_t X, const arb_mat_t A, const arb_mat_t B, slong prec)
{
    return _solve(X, NULL, A, B, SOLVE_GENERAL, 0, prec);
}

void
arb_mat_approx_solve_lu_precomp(arb_mat_t X, const slong * perm, const arb_mat_t A, const arb_mat_t B, slong prec)
{
    _solve(X, perm, A, B, SOLVE_LU_PRECOMP, 1, prec);
}

void
arb_mat_approx_solve_tril(arb_mat_t X, const arb_mat_t L, const arb_mat_t B, int unit, slong prec)
{
    _solve(X, NULL, L, B, SOLVE_TRIL, unit, prec);
}

void
arb_mat_approx_solve_triu(arb_mat_t X, const arb_mat_t U, const arb_mat_t B, int unit, slong prec)
{
    _solve(X, NULL, U, B, SOLVE_TRIU, unit, prec);
}
