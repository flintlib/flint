/*
    Copyright (C) 2020 Daniel Schultz
    This file is part of FLINT.
    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_mat.h"
#include "nmod_vec.h"
#include "thread_support.h"

#if HAVE_BLAS

#include "cblas.h"

/*
    b (and hence the signed remainder) definitely fits into an int, which
    helps the conversion to a double
*/
static void _lift_vec(double * a, ulong * b, slong len, ulong n)
{
    slong i;
    for (i = 0; i < len; i++)
        a[i] = (int)(b[i] - (n & FLINT_SIGN_EXT(n/2 - b[i])));
}

static void _lift_vec_sp(float * a, ulong * b, slong len, ulong n)
{
    slong i;
    for (i = 0; i < len; i++)
        a[i] = (int)(b[i] - (n & FLINT_SIGN_EXT(n/2 - b[i])));
}

static void _lift_vec_crt(double * a, ulong * b, slong len, nmod_t ctx)
{
    slong i;
    for (i = 0; i < len; i++)
    {
        mp_limb_t bn;
        NMOD_RED(bn, b[i], ctx);
        a[i] = (int)(bn - (ctx.n & FLINT_SIGN_EXT(ctx.n/2 - bn)));
    }
}


static int _nmod_mat_mul_blas_crt(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
{
    slong i, j, pi, pj;
    slong m = A->r;
    slong k = A->c;
    slong n = B->c;
    double * dC, * dA, * dB;
    nmod_t ctx = C->mod;
    ulong s, t, hi, lo, reshi, reslo, shifts[16];
    slong crtnum;
    nmod_t crtmod[12]; /* not nec prime */
    mp_limb_t pmodinv[12][12];
    mp_limb_t q[12], u[12], v[12];
    fmpz_t prodmod, maxentry;

    fmpz_init_set_ui(maxentry, k);
    fmpz_mul_ui(maxentry, maxentry, ctx.n - 1);
    fmpz_mul_ui(maxentry, maxentry, ctx.n - 1);

    t = n_sqrt((UWORD(1)<<55)/k - 1);
    fmpz_init_set_ui(prodmod, t);
    nmod_init(crtmod + 0, t);
    crtnum = 1;

    do {
        t = crtmod[crtnum - 1].n;
        do {
            if (crtnum > 11 || t < 100)
            {
                fmpz_clear(maxentry);
                fmpz_clear(prodmod);
                return 0;
            }
            t--;
        } while (n_gcd(fmpz_fdiv_ui(prodmod, t), t) != 1);
        fmpz_mul_ui(prodmod, prodmod, t);
        nmod_init(crtmod + crtnum, t);
        crtnum++;
    } while (fmpz_cmp(prodmod, maxentry) <= 0);

    /* note that if k is sufficiently big, i.e. k >= 4, then crtnum >= 3 */

    /* arange the crt moduli in increasing order */
    for (pi = 0; pi < crtnum/2; pi++)
    {
        nmod_t tmp = crtmod[pi];
        crtmod[pi] = crtmod[crtnum - 1 - pi];
        crtmod[crtnum - 1 - pi] = tmp;
    }

    dA = flint_malloc(m*k*sizeof(double));
    dB = flint_malloc(k*n*sizeof(double));
    dC = flint_calloc(crtnum*m*n, sizeof(double));

    for (pi = 0; pi < crtnum; pi++)
    {
        for (i = 0; i < m; i++)
            _lift_vec_crt(dA + i*k, A->rows[i], k, crtmod[pi]);

        for (i = 0; i < k; i++)
            _lift_vec_crt(dB + i*n, B->rows[i], n, crtmod[pi]);

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
                                       1.0, dA, k, dB, n, 0.0, dC + pi*m*n, n);
    }

    /*
        set p_i = crtmod[i].n
        for finding u given its image u[i] mod p_i, first solve for the v[i]:
            u = v[0] + v[1]*p_0 + v[2]*p_0*p_1 + ... + v[crtnum-1]*p_0*...*p_{crtnum-1}
        then evaluate this dot product modulo ctx.n to find u mod ctx.n
    */

    for (pi = 0; pi < crtnum; pi++)
    {
        t = 1;
        s = 1;
        for (pj = pi - 1; pj >= 0; pj--)
        {
            t = nmod_mul(t, crtmod[pj].n, crtmod[pi]);
            FLINT_ASSERT(crtmod[pj].n < ctx.n);
            s = nmod_mul(s, crtmod[pj].n, ctx);
            pmodinv[pi][pj] = nmod_neg(nmod_inv(t, crtmod[pi]), crtmod[pi]);
        }
        q[pi] = s; /* q[i] = p_0 * ... * p[i-1] mod ctx.n */

        /* for double -> nmod conversion */
        shifts[pi] = ((UWORD(1)<<54)/crtmod[pi].n)*crtmod[pi].n;
    }

    for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
    {
        for (pi = 0; pi < crtnum; pi++)
        {
            slong a = (slong) dC[i*n + j + pi*m*n];
            mp_limb_t b = (a < 0) ? a + shifts[pi] : a;
            NMOD_RED(u[pi], b, crtmod[pi]);
        }

        reslo = u[0];
        reshi = 0;
        for (pi = 1; pi < crtnum; pi++)
        {
            FLINT_ASSERT(u[pi] < crtmod[pi].n);
            FLINT_ASSERT(u[0] < crtmod[pi].n);
            t = pmodinv[pi][0]*nmod_sub(u[0], u[pi], crtmod[pi]);
            for (pj = 1; pj < pi; pj++)
                t += pmodinv[pi][pj]*v[pj];
            NMOD_RED(v[pi], t, crtmod[pi]);
            umul_ppmm(hi, lo, v[pi], q[pi]);
            add_ssaaaa(reshi, reslo, reshi, reslo, hi, lo);
        }

        if (reshi < ctx.n)
            NMOD_RED2(C->rows[i][j], reshi, reslo, ctx);
        else
            NMOD2_RED2(C->rows[i][j], reshi, reslo, ctx);
    }

    flint_free(dA);
    flint_free(dB);
    flint_free(dC);

    fmpz_clear(maxentry);
    fmpz_clear(prodmod);

    return 1;
}


static int _nmod_mat_mul_blas_sp(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
{
    slong i, j;
    slong m = A->r;
    slong k = A->c;
    slong n = B->c;
    float * dC, * dA, * dB;
    ulong shift;
    nmod_t ctx = C->mod;

    dA = flint_malloc(m*k*sizeof(float));
    dB = flint_malloc(k*n*sizeof(float));
    dC = flint_calloc(m*n, sizeof(float));

    for (i = 0; i < m; i++)
        _lift_vec_sp(dA + i*k, A->rows[i], k, ctx.n);

    for (i = 0; i < k; i++)
        _lift_vec_sp(dB + i*n, B->rows[i], n, ctx.n);

    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
                                                1.0, dA, k, dB, n, 0.0, dC, n);

    /*
        the shift for negative a must satisfy
            ctx.n divides shift, and
            2^24 <= shift <= 2^64
    */
    shift = ((UWORD(1)<<30)/ctx.n)*ctx.n;
    for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
    {
        slong a = (slong) dC[i*n + j];
        mp_limb_t b = (a < 0) ? a + shift : a;
        NMOD_RED(C->rows[i][j], b, ctx);
    }

    flint_free(dA);
    flint_free(dB);
    flint_free(dC);

    return 1;
}


int nmod_mat_mul_blas(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
{
    slong i, j;
    slong m = A->r;
    slong k = A->c;
    slong n = B->c;
    double * dC, * dA, * dB;
    ulong hi, lo, shift;
    nmod_t ctx = C->mod;

    if (m < 1 || k < 1 || n < 1 || m > INT_MAX || k > INT_MAX || n > INT_MAX)
        return 0;

    /*
        Each nmod has |_lift()| <= floor(mod.n/2) in the signed representation.
        Want k*floor(mod.n/2)^2 < 2^53 **assuming** no BLAS rounding errors,
        otherwise, must resort to modular methods.
    */

    if (FLINT_BITS != 64)
        return 0;

    umul_ppmm(hi, lo, ctx.n/2, ctx.n/2);
    if (hi != 0)
        return _nmod_mat_mul_blas_crt(C, A, B);

    umul_ppmm(hi, lo, lo, k);
    if (hi != 0 || lo >= (UWORD(1) << 53))
        return _nmod_mat_mul_blas_crt(C, A, B);

    if (lo < UWORD(1) << 24)
        return _nmod_mat_mul_blas_sp(C, A, B);

    dA = flint_malloc(m*k*sizeof(double));
    dB = flint_malloc(k*n*sizeof(double));
    dC = flint_calloc(m*n, sizeof(double));

    for (i = 0; i < m; i++)
        _lift_vec(dA + i*k, A->rows[i], k, ctx.n);

    for (i = 0; i < k; i++)
        _lift_vec(dB + i*n, B->rows[i], n, ctx.n);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,
                                                1.0, dA, k, dB, n, 0.0, dC, n);

    /*
        the shift for negative a must satisfy
            ctx.n divides shift, and
            2^53 <= shift <= 2^64
    */
    shift = ((UWORD(1)<<54)/ctx.n)*ctx.n;
    for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
    {
        slong a = (slong) dC[i*n + j];
        mp_limb_t b = (a < 0) ? a + shift : a;
        NMOD_RED(C->rows[i][j], b, ctx);
    }

    flint_free(dA);
    flint_free(dB);
    flint_free(dC);

    return 1;
}

#else

int nmod_mat_mul_blas(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
{
    return 0;
}

#endif
