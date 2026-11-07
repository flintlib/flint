/*
    Copyright (C) 2013 Timo Hartmann
    Copyright (C) 2018 arbguest
    Copyright (C) 2018, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_mat.h"
#include "acf.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_mat.h"
#include "arb_mat/impl.h"

/* The approximate LU decomposition, solving and eigenvalue functions of
   acb_mat: the generic algorithms of gr_mat over dfloat, nfloat_complex
   or acf; see arb_mat/impl.h and arb_mat/approx_gr.c. */

static int
_entry_is_finite(gr_srcptr x, int kind, gr_ctx_t ctx)
{
    return _arb_mat_approx_is_finite(x, kind, ctx) &&
        _arb_mat_approx_is_finite(_arb_mat_approx_imag_ptr((gr_ptr) x, kind, ctx), kind, ctx);
}

/* write the midpoint (with zero radius) of an entry known to be finite */
static void
_set_entry(acb_t res, gr_srcptr x, slong scale, int kind, gr_ctx_t ctx)
{
    _arb_mat_approx_get_arf(arb_midref(acb_realref(res)), x, scale, kind, ctx);
    _arb_mat_approx_get_arf(arb_midref(acb_imagref(res)),
        _arb_mat_approx_imag_ptr((gr_ptr) x, kind, ctx), scale, kind, ctx);
    mag_zero(arb_radref(acb_realref(res)));
    mag_zero(arb_radref(acb_imagref(res)));
}

static int
_exp_range(slong * emin, slong * emax, const acb_mat_t A, int part)
{
    slong i, j;

    *emin = WORD_MAX;
    *emax = WORD_MIN;

    for (i = 0; i < acb_mat_nrows(A); i++)
        for (j = 0; j < acb_mat_ncols(A); j++)
            if (_arb_mat_approx_in_part(i, j, part))
                if (!_arb_mat_approx_exp_update(emin, emax,
                        arb_midref(acb_realref(acb_mat_entry(A, i, j))))
                    || !_arb_mat_approx_exp_update(emin, emax,
                        arb_midref(acb_imagref(acb_mat_entry(A, i, j)))))
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
_get_gr(gr_mat_t res, const acb_mat_t A, slong scale, int part, int kind, gr_ctx_t ctx)
{
    slong i, j, sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;
    gr_ptr e;

    for (i = 0; i < acb_mat_nrows(A) && status == GR_SUCCESS; i++)
    {
        for (j = 0; j < acb_mat_ncols(A); j++)
        {
            if (_arb_mat_approx_in_part(i, j, part))
            {
                e = GR_MAT_ENTRY(res, i, j, sz);
                status |= _arb_mat_approx_set_arf(e,
                    arb_midref(acb_realref(acb_mat_entry(A, i, j))), scale, kind, ctx);
                status |= _arb_mat_approx_set_arf(_arb_mat_approx_imag_ptr(e, kind, ctx),
                    arb_midref(acb_imagref(acb_mat_entry(A, i, j))), scale, kind, ctx);
            }
        }
    }

    return status;
}

/* An input matrix for the ring, as in arb_mat/approx_gr.c: for acf,
   shallow copies of the midpoints, to be freed with _input_clear. */
static int
_input(gr_mat_t res, const acb_mat_t A, int part1, slong scale1, int part2, slong scale2, int kind, gr_ctx_t ctx)
{
    slong i, j, r = acb_mat_nrows(A), c = acb_mat_ncols(A);
    int status = GR_SUCCESS;
    acf_struct * e;

    if (kind == ARB_MAT_APPROX_ARF)
    {
        res->entries = flint_calloc(r * c + 1, sizeof(acf_struct));
        res->r = r;
        res->c = c;
        res->stride = c;

        for (i = 0; i < r; i++)
        {
            for (j = 0; j < c; j++)
            {
                if (_arb_mat_approx_in_part(i, j, part1)
                    || (part2 != -1 && _arb_mat_approx_in_part(i, j, part2)))
                {
                    e = ((acf_struct *) res->entries) + i * c + j;
                    *acf_realref(e) = *arb_midref(acb_realref(acb_mat_entry(A, i, j)));
                    *acf_imagref(e) = *arb_midref(acb_imagref(acb_mat_entry(A, i, j)));
                }
            }
        }

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
                if (!_entry_is_finite(GR_MAT_ENTRY(A, i, j, sz), kind, ctx))
                    return 0;

    return 1;
}

/* Write the midpoints of res (with zero radii), scaling the entries in
   the given part by 2^scale1 and the others by 2^scale2. Only called once
   the result is known to be finite. */
static void
_set_gr(acb_mat_t res, const gr_mat_t A, int part, slong scale1, slong scale2, int kind, gr_ctx_t ctx)
{
    slong i, j, sz = ctx->sizeof_elem;

    for (i = 0; i < A->r; i++)
        for (j = 0; j < A->c; j++)
            _set_entry(acb_mat_entry(res, i, j), GR_MAT_ENTRY(A, i, j, sz),
                _arb_mat_approx_in_part(i, j, part) ? scale1 : scale2, kind, ctx);
}

static void
_set_nan(acb_mat_t X)
{
    slong i, j;

    for (i = 0; i < acb_mat_nrows(X); i++)
        for (j = 0; j < acb_mat_ncols(X); j++)
        {
            acb_indeterminate(acb_mat_entry(X, i, j));
            mag_zero(arb_radref(acb_realref(acb_mat_entry(X, i, j))));
            mag_zero(arb_radref(acb_imagref(acb_mat_entry(X, i, j))));
        }
}

int
acb_mat_approx_lu(slong * P, acb_mat_t LU, const acb_mat_t A, slong prec)
{
    slong m = acb_mat_nrows(A), n = acb_mat_ncols(A), emin, emax, rank, s;
    int kind, status, result = 0, done = 0;
    gr_ctx_t ctx;
    gr_mat_t GA, T;

    if (m == 0 || n == 0)
        return 1;

    if (_exp_range(&emin, &emax, A, ARB_MAT_APPROX_ALL) && emax != WORD_MIN)
    {
        kind = _arb_mat_approx_backend(FLINT_MIN(m, n), prec,
            _arb_mat_approx_range(emin, emax), 1, ARB_MAT_APPROX_OP_LU);
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

        if (_arb_mat_approx_ctx_init(ctx, kind, prec, 1) != GR_SUCCESS)
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

/* As in arb_mat/approx_gr.c. */

#define SOLVE_GENERAL 0
#define SOLVE_LU_PRECOMP 1
#define SOLVE_TRIL 2
#define SOLVE_TRIU 3

static int
_solve(acb_mat_t X, const slong * perm, const acb_mat_t T,
    const acb_mat_t B, int which, int unit, slong prec)
{
    slong n = acb_mat_nrows(T), m = acb_mat_ncols(B);
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
        /* a zero matrix (singular) is left to acf */
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

    kind = ok ? _arb_mat_approx_backend(n, prec, range, 1, op) : ARB_MAT_APPROX_ARF;

    for ( ; !done && kind != ARB_MAT_APPROX_NONE; kind = _arb_mat_approx_backend_next(kind, prec))
    {
        if (kind == ARB_MAT_APPROX_ARF)
            sT = sB = 0;

        if (_arb_mat_approx_ctx_init(ctx, kind, prec, 1) != GR_SUCCESS)
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
acb_mat_approx_solve(acb_mat_t X, const acb_mat_t A, const acb_mat_t B, slong prec)
{
    return _solve(X, NULL, A, B, SOLVE_GENERAL, 0, prec);
}

void
acb_mat_approx_solve_lu_precomp(acb_mat_t X, const slong * perm, const acb_mat_t A, const acb_mat_t B, slong prec)
{
    _solve(X, perm, A, B, SOLVE_LU_PRECOMP, 1, prec);
}

void
acb_mat_approx_solve_tril(acb_mat_t X, const acb_mat_t L, const acb_mat_t B, int unit, slong prec)
{
    _solve(X, NULL, L, B, SOLVE_TRIL, unit, prec);
}

void
acb_mat_approx_solve_triu(acb_mat_t X, const acb_mat_t U, const acb_mat_t B, int unit, slong prec)
{
    _solve(X, NULL, U, B, SOLVE_TRIU, unit, prec);
}

/* Eigenvalues: gr_mat_eig_qr on the scaled midpoints. With dfloat and
   nfloat, any failure (including no convergence) moves on to the next
   type; with acf, which comes last, the approximations are written even
   if the iteration has not converged, and the return value is 0. */
int
acb_mat_approx_eig_qr(acb_ptr E, acb_mat_t L, acb_mat_t R, const acb_mat_t A, const mag_t tol, slong maxiter, slong prec)
{
    slong n = acb_mat_nrows(A), emin, emax, temin, temax, i, sz, s;
    int kind, status, converged, result = 0, done = 0;
    gr_ctx_t ctx;
    gr_mat_t GA, GL, GR;
    gr_ptr GE, eps;
    arf_t t;

    if (n != acb_mat_ncols(A))
        flint_throw(FLINT_ERROR, "acb_mat_approx_eig_qr: the matrix must be square\n");

    if (n == 0)
        return 1;

    /* The same defaults at prec bits with every type, rather than at the
       precision of the type; the deflation tolerance is relative to the
       norm, so it is unaffected by the scaling. */
    if (maxiter <= 0)
        maxiter = 14 * n + prec / 10;

    arf_init(t);

    if (tol == NULL)
    {
        arf_one(t);
        arf_mul_2exp_si(t, t, -prec);
        arf_div_ui(t, t, 100 * n, MAG_BITS, ARF_RND_DOWN);
    }
    else
    {
        arf_set_mag(t, tol);
    }

    temin = WORD_MAX;
    temax = WORD_MIN;

    if (_exp_range(&emin, &emax, A, ARB_MAT_APPROX_ALL) && emax != WORD_MIN
        && _arb_mat_approx_exp_update(&temin, &temax, t))
    {
        kind = _arb_mat_approx_backend(n, prec,
            _arb_mat_approx_range(emin, emax), 1, ARB_MAT_APPROX_OP_EIG);
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

        if (_arb_mat_approx_ctx_init(ctx, kind, prec, 1) != GR_SUCCESS)
            continue;

        sz = ctx->sizeof_elem;
        GE = gr_heap_init_vec(n, ctx);
        eps = gr_heap_init(ctx);
        gr_mat_init(GL, (L != NULL) ? n : 0, (L != NULL) ? n : 0, ctx);
        gr_mat_init(GR, (R != NULL) ? n : 0, (R != NULL) ? n : 0, ctx);

        status = _input(GA, A, ARB_MAT_APPROX_ALL, s, -1, 0, kind, ctx);
        status |= _arb_mat_approx_set_arf(eps, t, 0, kind, ctx);

        if (status == GR_SUCCESS)
            status = _gr_mat_eig_qr(GE, (L != NULL) ? GL : NULL,
                (R != NULL) ? GR : NULL, GA, eps, maxiter,
                (kind == ARB_MAT_APPROX_ARF) ? &converged : NULL, ctx);

        /* the input may be shallow copies of entries of L or R */
        _input_clear(GA, kind, ctx);

        if (status == GR_SUCCESS && kind != ARB_MAT_APPROX_ARF)
        {
            converged = 1;
            for (i = 0; i < n && status == GR_SUCCESS; i++)
                if (!_entry_is_finite(GR_ENTRY(GE, i, sz), kind, ctx))
                    status = GR_UNABLE;
            if ((L != NULL && !_is_finite(GL, ARB_MAT_APPROX_ALL, kind, ctx))
                || (R != NULL && !_is_finite(GR, ARB_MAT_APPROX_ALL, kind, ctx)))
                status = GR_UNABLE;
        }

        if (status == GR_SUCCESS)
        {
            /* the eigenvalues of 2^s A are those of A times 2^s; the
               eigenvectors are the same */
            for (i = 0; i < n; i++)
                _set_entry(E + i, GR_ENTRY(GE, i, sz), -s, kind, ctx);
            if (L != NULL)
                _set_gr(L, GL, ARB_MAT_APPROX_ALL, 0, 0, kind, ctx);
            if (R != NULL)
                _set_gr(R, GR, ARB_MAT_APPROX_ALL, 0, 0, kind, ctx);
            result = converged;
            done = 1;
        }
        else if (kind == ARB_MAT_APPROX_ARF)
        {
            /* an arithmetic failure (nonfinite input) */
            for (i = 0; i < n; i++)
            {
                acb_indeterminate(E + i);
                mag_zero(arb_radref(acb_realref(E + i)));
                mag_zero(arb_radref(acb_imagref(E + i)));
            }
            if (L != NULL)
                _set_nan(L);
            if (R != NULL)
                _set_nan(R);
            done = 1;
        }

        gr_heap_clear_vec(GE, n, ctx);
        gr_heap_clear(eps, ctx);
        gr_mat_clear(GL, ctx);
        gr_mat_clear(GR, ctx);
        gr_ctx_clear(ctx);
    }

    arf_clear(t);
    return result;
}
