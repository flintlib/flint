/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpzi.h"
#include "gr.h"
#include "fft_small.h"
#include "ulong_extras.h"

/* Middle products of Gaussian integer polynomials with the coefficient
   arithmetic done pointwise on fft_small transforms: the same rolling
   schoolbook as _fmpz_poly_mulmid_classical_fft_small, with two
   transforms per coefficient (real and imaginary parts) and four
   pointwise operations per complex product,
       re += ra rb,  re -= ia ib,  im += ra ib,  im += ia rb.
   Squaring shares transforms, pairs symmetric terms through doubled
   accumulators, and computes diagonal squares with three
   multiplications ((r + ii)^2 = r^2 - i^2 + 2ri i, the cross term
   riding the doubled accumulator). Returns 1 on success, 0 when the
   transformed representation is unavailable or unprofitable. */

#if FLINT_HAVE_FFT_SMALL

static nn_srcptr
_coeff_limbs(const fmpz * f, ulong * scratch, slong * n, int * sgn)
{
    fmpz c = *f;

    if (!COEFF_IS_MPZ(c))
    {
        *scratch = FLINT_ABS(c);
        *n = 1;
        *sgn = c < 0;
        return scratch;
    }
    else
    {
        mpz_srcptr m = COEFF_TO_PTR(c);

        *n = FLINT_ABS(m->_mp_size);
        *sgn = m->_mp_size < 0;
        return m->_mp_d;
    }
}

static int
_set_coeff(gr_ptr elem, const fmpz * f, gr_ctx_t ctx)
{
    ulong scratch;
    slong n;
    int sgn;
    nn_srcptr p = _coeff_limbs(f, &scratch, &n, &sgn);

    return gr_transformed_mpn_set(elem, p, n, sgn, ctx);
}

/* the accumulated element, converted out directly into an fmpz */
/* converts the accumulator out destructively -- the caller re-initializes
   it (borrowed storage, so that is free) before the next output */
static int
_get_coeff(fmpz * f, gr_ptr elem, gr_ctx_t ctx)
{
    slong w = gr_transformed_mpn_get_limbs(ctx, elem);
    slong zl;
    int sg, status;
    mpz_ptr m;

    if (w <= 0)
        w = 1;

    m = _fmpz_promote(f);
    status = gr_transformed_mpn_get_destructive(FLINT_MPZ_REALLOC(m, w), w,
                                                &zl, &sg, elem, ctx);
    if (status != GR_SUCCESS)
    {
        _fmpz_demote(f);
        fmpz_zero(f);
        return status;
    }

    m->_mp_size = sg ? -zl : zl;
    _fmpz_demote_val(f);
    return GR_SUCCESS;
}


int
_fmpzi_poly_mulmid_classical_fft_small(fmpzi_struct * res,
    const fmpzi_struct * poly1, slong len1,
    const fmpzi_struct * poly2, slong len2,
    slong nlo, slong nhi)
{
    gr_ctx_t ctx;
    gr_ptr E, B, X, accR, accI, accR2, accI2;
    double * accd[4];
    slong i, j, k, sz, abits, bbits;
    slong * xidx;
    int share = (poly1 == poly2);
    int status = GR_SUCCESS;
    int ok;

    if (len2 < 1 || len1 < 1 || nlo >= nhi)
        return 0;

    if (len1 < len2)
    {
        const fmpzi_struct * tp = poly1;
        slong tl = len1;
        poly1 = poly2;
        len1 = len2;
        poly2 = tp;
        len2 = tl;
    }
    k = len2;

    /* bounds over both components; the subtraction in the real part
       makes the representation signed regardless of input signs */
    {
        slong m, t;
        abits = 0;
        for (i = 0; i < len1; i++)
        {
            t = fmpz_bits(fmpzi_realref(poly1 + i));
            m = fmpz_bits(fmpzi_imagref(poly1 + i));
            abits = FLINT_MAX(abits, FLINT_MAX(t, m));
        }
        if (share)
            bbits = abits;
        else
        {
            bbits = 0;
            for (i = 0; i < len2; i++)
            {
                t = fmpz_bits(fmpzi_realref(poly2 + i));
                m = fmpz_bits(fmpzi_imagref(poly2 + i));
                bbits = FLINT_MAX(bbits, FLINT_MAX(t, m));
            }
        }
    }
    if (abits == 0 || bbits == 0)
    {
        for (i = 0; i < nhi - nlo; i++)
            fmpzi_zero(res + i);
        return 1;
    }

    /* Decline if too sparse; to be moved to the algorithm selection in
       a future _fmpzi_poly_mulmid when this  function has been integrated.
    {
        slong wa = 0, wb = 0;
        double density, need;
        for (i = 0; i < len1; i++)
            wa += fmpz_bits(fmpzi_realref(poly1 + i))
                + fmpz_bits(fmpzi_imagref(poly1 + i));
        if (share)
            wb = wa;
        else
            for (i = 0; i < len2; i++)
                wb += fmpz_bits(fmpzi_realref(poly2 + i))
                    + fmpz_bits(fmpzi_imagref(poly2 + i));
        density = (double) (wa + wb)
                  / (2.0 * ((double) len1 * abits
                            + (double) len2 * bbits));
        need = (k <= 2) ? 0.90 : (k <= 4) ? 0.62
             : (k <= 8) ? 0.50 : 0.30;
        if (density < need)
            return 0;
    } */

    /* each output component accumulates at most 2 k product
       magnitudes of one product's size */
    if (gr_ctx_init_transformed_mpn(ctx,
            abits + bbits, 2 * k, 1, 4 * k + 4) != GR_SUCCESS)
        return 0;

    sz = ctx->sizeof_elem;
    /* 2k kept transforms (re, im interleaved), 2k rolling slots, and
       four accumulators on borrowed slabs */
    E = flint_malloc((4 * k + 4) * sz);
    B = E;
    X = GR_ENTRY(E, 2 * k, sz);
    accR = GR_ENTRY(E, 4 * k, sz);
    accI = GR_ENTRY(E, 4 * k + 1, sz);
    accR2 = GR_ENTRY(E, 4 * k + 2, sz);
    accI2 = GR_ENTRY(E, 4 * k + 3, sz);
    for (i = 0; i < 4 * k; i++)
        gr_init(GR_ENTRY(E, i, sz), ctx);
    {
        ulong slab = n_round_up(gr_transformed_mpn_sizeof_data(ctx),
                                FLINT_FFT_SMALL_ALIGNMENT);
        for (i = 0; i < 4; i++)
            accd[i] = flint_aligned_alloc(FLINT_FFT_SMALL_ALIGNMENT,
                                          slab);
        gr_transformed_mpn_init_borrowed(accR, accd[0], ctx);
        gr_transformed_mpn_init_borrowed(accI, accd[1], ctx);
        gr_transformed_mpn_init_borrowed(accR2, accd[2], ctx);
        gr_transformed_mpn_init_borrowed(accI2, accd[3], ctx);
    }
    xidx = flint_malloc(k * sizeof(slong));
    for (i = 0; i < k; i++)
        xidx[i] = -1;

    for (j = 0; j < k && status == GR_SUCCESS; j++)
    {
        const fmpzi_struct * c = (share ? poly1 : poly2) + j;
        status |= _set_coeff(GR_ENTRY(B, 2 * j, sz),
                             fmpzi_realref(c), ctx);
        status |= _set_coeff(GR_ENTRY(B, 2 * j + 1, sz),
                             fmpzi_imagref(c), ctx);
    }

#define A_PAIR(t, outr, outi) \
    do { \
        if (share && (t) < k) \
        { \
            (outr) = GR_ENTRY(B, 2 * (t), sz); \
            (outi) = GR_ENTRY(B, 2 * (t) + 1, sz); \
        } \
        else \
        { \
            slong _s = share ? ((t) - k) % k : (t) % k; \
            gr_ptr _er = GR_ENTRY(X, 2 * _s, sz); \
            gr_ptr _ei = GR_ENTRY(X, 2 * _s + 1, sz); \
            if (xidx[_s] != (t)) \
            { \
                status |= _set_coeff(_er, \
                    fmpzi_realref(poly1 + (t)), ctx); \
                status |= _set_coeff(_ei, \
                    fmpzi_imagref(poly1 + (t)), ctx); \
                xidx[_s] = (t); \
            } \
            (outr) = _er; \
            (outi) = _ei; \
        } \
    } while (0)

    /* one complex product accumulated into the (re, im) target pair */
#define CPROD(tr, ti, trused, tiused, ar, ai, br, bi) \
    do { \
        status |= (trused) ? gr_addmul(tr, ar, br, ctx) \
                           : gr_mul(tr, ar, br, ctx); \
        (trused) = 1; \
        status |= gr_submul(tr, ai, bi, ctx); \
        status |= (tiused) ? gr_addmul(ti, ar, bi, ctx) \
                           : gr_mul(ti, ar, bi, ctx); \
        (tiused) = 1; \
        status |= gr_addmul(ti, ai, br, ctx); \
    } while (0)

    for (i = nlo; i < nhi && status == GR_SUCCESS; i++)
    {
        slong jlo = FLINT_MAX(0, i - len1 + 1);
        slong jhi = FLINT_MIN(i, k - 1);
        int r1u = 0, i1u = 0, r2u = 0, i2u = 0;
        gr_ptr ar, ai, br, bi;

        if (!share)
        {
            for (j = jlo; j <= jhi; j++)
            {
                A_PAIR(i - j, ar, ai);
                br = GR_ENTRY(B, 2 * j, sz);
                bi = GR_ENTRY(B, 2 * j + 1, sz);
                if (status != GR_SUCCESS)
                    break;
                CPROD(accR, accI, r1u, i1u, ar, ai, br, bi);
            }
        }
        else
        {
            for (j = jlo; j <= jhi; j++)
            {
                slong jp = i - j;
                int paired = (jp >= jlo && jp <= jhi);

                if (paired && j > jp)
                    continue;

                A_PAIR(i - j, ar, ai);
                A_PAIR(j, br, bi);
                if (status != GR_SUCCESS)
                    break;

                if (paired && j < jp)
                    CPROD(accR2, accI2, r2u, i2u, ar, ai, br, bi);
                else if (paired)
                {
                    /* the diagonal square: r^2 - i^2 into the single
                       accumulator, r i into the doubled one */
                    status |= r1u ? gr_addmul(accR, ar, ar, ctx)
                                  : gr_mul(accR, ar, ar, ctx);
                    r1u = 1;
                    status |= gr_submul(accR, ai, ai, ctx);
                    status |= i2u ? gr_addmul(accI2, ar, ai, ctx)
                                  : gr_mul(accI2, ar, ai, ctx);
                    i2u = 1;
                }
                else
                    CPROD(accR, accI, r1u, i1u, ar, ai, br, bi);
            }
        }

        if (status == GR_SUCCESS && r2u)
        {
            if (!r1u)
            {
                status |= gr_set(accR, accR2, ctx);
                status |= gr_add(accR, accR, accR2, ctx);
            }
            else
            {
                status |= gr_add(accR, accR, accR2, ctx);
                status |= gr_add(accR, accR, accR2, ctx);
            }
            r1u = 1;
        }
        if (status == GR_SUCCESS && i2u)
        {
            if (!i1u)
            {
                status |= gr_set(accI, accI2, ctx);
                status |= gr_add(accI, accI, accI2, ctx);
            }
            else
            {
                status |= gr_add(accI, accI, accI2, ctx);
                status |= gr_add(accI, accI, accI2, ctx);
            }
            i1u = 1;
        }

        if (status == GR_SUCCESS)
        {
            if (r1u)
                status |= _get_coeff(fmpzi_realref(res + i - nlo),
                                     accR, ctx);
            else
                fmpz_zero(fmpzi_realref(res + i - nlo));
            if (i1u)
                status |= _get_coeff(fmpzi_imagref(res + i - nlo),
                                     accI, ctx);
            else
                fmpz_zero(fmpzi_imagref(res + i - nlo));
            /* borrowed storage: re-initialization is free */
            gr_transformed_mpn_init_borrowed(accR, accd[0], ctx);
            gr_transformed_mpn_init_borrowed(accI, accd[1], ctx);
            gr_transformed_mpn_init_borrowed(accR2, accd[2], ctx);
            gr_transformed_mpn_init_borrowed(accI2, accd[3], ctx);
        }
    }

#undef A_PAIR
#undef CPROD

    ok = (status == GR_SUCCESS);
    for (i = 0; i < 4 * k; i++)
        gr_clear(GR_ENTRY(E, i, sz), ctx);
    flint_free(E);
    for (i = 0; i < 4; i++)
        flint_free(accd[i]);
    flint_free(xidx);
    gr_ctx_clear(ctx);
    return ok;
}

#else

int
_fmpzi_poly_mulmid_classical_fft_small(fmpzi_struct * res,
    const fmpzi_struct * poly1, slong len1,
    const fmpzi_struct * poly2, slong len2,
    slong nlo, slong nhi)
{
    return 0;
}

#endif
