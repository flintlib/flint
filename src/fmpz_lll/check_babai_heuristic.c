/*
    Copyright (C) 2009, 2010 William Hart
    Copyright (C) 2009, 2010 Andy Novocin
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <float.h>
#include "gmpcompat.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "fmpz_lll.h"
#ifdef GM
#undef GM
#endif
#define GM ((fl->rt == Z_BASIS) ? A->exactSP : B)

#include "gr.h"
#include "gr_vec.h"
#include "gr_mat.h"
#include "nfloat.h"

static int _gr_cmp_d(gr_srcptr x, double y, gr_ctx_t ctx)
{
    gr_ptr t;
    int status;
    int cmp;
    GR_TMP_INIT(t, ctx);
    status = gr_set_d(t, y, ctx);
    status |= gr_cmp(&cmp, x, t, ctx);
    GR_TMP_CLEAR(t, ctx);
    GR_MUST_SUCCEED(status);
    return cmp;
}

static int _gr_sgn(gr_srcptr x, gr_ctx_t ctx)
{
    gr_ptr t;
    int sgn;
    GR_TMP_INIT(t, ctx);
    GR_MUST_SUCCEED(gr_cmp(&sgn, x, t, ctx));
    GR_TMP_CLEAR(t, ctx);
    return sgn;
}

static int _gr_vec_set_fmpz_vec(gr_srcptr res, const fmpz * vec, slong len, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong i, sz = ctx->sizeof_elem;

    for (i = 0; i < len; i++)
        status |= gr_set_fmpz(GR_ENTRY(res, i, sz), vec + i, ctx);

    return status;
}

static int _gr_vec_norm2(gr_ptr res, gr_srcptr vec, slong len, gr_ctx_t ctx)
{
    /* todo */
    return _gr_vec_dot(res, NULL, 0, vec, vec, len, ctx);
}


/* XXX: dubious use of DBL_MIN */

int
fmpz_lll_check_babai_heuristic(int kappa, fmpz_mat_t B, fmpz_mat_t U,
                               gr_mat_t mu, gr_mat_t r, gr_ptr s,
                               gr_mat_t appB, fmpz_gram_t A, int a, int zeros,
                               int kappamax, int n, gr_ptr tmp, gr_ptr rtmp,
                               gr_ctx_t ctx, const fmpz_lll_t fl)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

#define ENTRY(mat, ii, jj) GR_MAT_ENTRY(mat, ii, jj, sz)
#define ROW(mat, ii) GR_MAT_ENTRY(mat, ii, 0, sz)

    if (fl->rt == Z_BASIS && fl->gt == APPROX)
    {
        int i, j, k, test, aa;
        slong max_expo = WORD_MAX;
        fmpz_t ztmp;
        double halfplus, onedothalfplus;
        ulong loops;

        fmpz_init(ztmp);

        aa = (a > zeros) ? a : zeros + 1;

        halfplus = (fl->eta + 0.5) / 2;
        onedothalfplus = 1.0 + halfplus;

        loops = 0;

        do
        {
            test = 0;

            /* ************************************** */
            /* Step2: compute the GSO for stage kappa */
            /* ************************************** */

            for (j = aa; j < kappa; j++)
            {
                if (_gr_cmp_d(ENTRY(A->appSP2, kappa, j), DBL_MIN, ctx) == 0)
                {
                    status |= _gr_vec_dot(ENTRY(A->appSP2, kappa, j), NULL, 0, ROW(appB, kappa), ROW(appB, j), n, ctx);

#if 0
                    /* If a heuristic told us that some cancellation probably happened,
                       recompute the scalar product at full precision */
                    _fmpz_vec_dot(ztmp, fmpz_mat_row(B, kappa), fmpz_mat_row(B, j), n);
                    status |= gr_set_fmpz(ENTRY(A->appSP2, kappa, j), ztmp, ctx);
#endif
                }

                /* we write to rtmp instead of ENTRY(r, kappa, j) directly
                   to avoid possible aliasing issues */
                status |= _gr_vec_dot(rtmp, ENTRY(A->appSP2, kappa, j), 1, ENTRY(mu, j, zeros + 1), ENTRY(r, kappa, zeros + 1), j - 1 - zeros, ctx);
                status |= gr_set(ENTRY(r, kappa, j), rtmp, ctx);
                status |= gr_div(ENTRY(mu, kappa, j), ENTRY(r, kappa, j), ENTRY(r, j, j), ctx);
            }

            if (loops >= 20)
            {
                slong new_max_expo = WORD_MIN;
                for (j = 0; j < kappa; j++)
                {
                    slong expo2;
                    double FLINT_SET_BUT_UNUSED(mant);
                    status |= gr_get_d_2exp_si(&mant, &expo2, ENTRY(mu, kappa, j), ctx);
                    new_max_expo = FLINT_MAX(new_max_expo, expo2);
                }
                if (new_max_expo > max_expo - SIZE_RED_FAILURE_THRESH)
                {
                    fmpz_clear(ztmp);
                    return -1;
                }
                max_expo = new_max_expo;
            }

            /* **************************** */
            /* Step3--5: compute the X_j's  */
            /* **************************** */

            for (j = kappa - 1; j > zeros; j--)
            {
                /* test of the relaxed size-reduction condition */
                status |= gr_abs(tmp, ENTRY(mu, kappa, j), ctx);

                if (_gr_cmp_d(tmp, halfplus, ctx) > 0)
                {
                    test = 1;

                    /* we consider separately the cases X = +-1 */
                    if (_gr_cmp_d(tmp, onedothalfplus, ctx) <= 0)
                    {
                        int sgn = _gr_sgn(ENTRY(mu, kappa, j), ctx);
                        if (sgn >= 0)   /* in this case, X is 1 */
                        {
                            status |= _gr_vec_sub(ENTRY(mu, kappa, zeros + 1), ENTRY(mu, kappa, zeros + 1), ENTRY(mu, j, zeros + 1), j - (zeros + 1), ctx);
                            _fmpz_vec_sub(fmpz_mat_row(B, kappa), fmpz_mat_row(B, kappa), fmpz_mat_row(B, j), n);
                            if (U != NULL)
                                _fmpz_vec_sub(fmpz_mat_row(U, kappa), fmpz_mat_row(U, kappa), fmpz_mat_row(U, j), U->c);
                        }
                        else    /* otherwise X is -1 */
                        {
                            status |= _gr_vec_add(ENTRY(mu, kappa, zeros + 1), ENTRY(mu, kappa, zeros + 1), ENTRY(mu, j, zeros + 1), j - (zeros + 1), ctx);
                            _fmpz_vec_add(fmpz_mat_row(B, kappa), fmpz_mat_row(B, kappa), fmpz_mat_row(B, j), n);
                            if (U != NULL)
                                _fmpz_vec_add(fmpz_mat_row(U, kappa), fmpz_mat_row(U, kappa), fmpz_mat_row(U, j), U->c);
                        }
                    }
                    else        /* we must have |X| >= 2 */
                    {
                        status |= gr_set(tmp, ENTRY(mu, kappa, j), ctx);
                        status |= gr_set_d(rtmp, 0.5, ctx);

                        if (_gr_sgn(tmp, ctx) < 0)
                        {
                            status |= gr_sub(tmp, tmp, rtmp, ctx);
                            status |= gr_ceil(tmp, tmp, ctx);
                        }
                        else
                        {
                            status |= gr_add(tmp, tmp, rtmp, ctx);
                            status |= gr_floor(tmp, tmp, ctx);
                        }

                        /* TODO: consider converting to integer if high precision? */
                        /* TODO: or make nfloat check for small multiplier */
                        status |= _gr_vec_submul_scalar(ENTRY(mu, kappa, zeros + 1), ENTRY(mu, j, zeros + 1), j - (zeros + 1), tmp, ctx);
                        status |= gr_get_fmpz(ztmp, tmp, ctx);

                        _fmpz_vec_scalar_submul_fmpz(fmpz_mat_row(B, kappa), fmpz_mat_row(B, j), n, ztmp);
                        if (U != NULL)
                            _fmpz_vec_scalar_submul_fmpz(fmpz_mat_row(U, kappa), fmpz_mat_row(U, j), U->c, ztmp);
                    }
                }
            }

            if (test)           /* Anything happened? */
            {
                status |= _gr_vec_set_fmpz_vec(ROW(appB, kappa), fmpz_mat_row(B, kappa), n, ctx);
                aa = zeros + 1;

                for (i = zeros + 1; i <= kappa; i++)
                    status |= gr_set_d(ENTRY(A->appSP2, kappa, i), DBL_MIN, ctx);

                for (i = kappa + 1; i <= kappamax; i++)
                    status |= gr_set_d(ENTRY(A->appSP2, i, kappa), DBL_MIN, ctx);
            }
            loops++;
        } while (test);

        if (_gr_cmp_d(ENTRY(A->appSP2, kappa, kappa), DBL_MIN, ctx) == 0)
            status |= _gr_vec_norm2(ENTRY(A->appSP2, kappa, kappa), ROW(appB, kappa), n, ctx);

        status |= gr_set(GR_ENTRY(s, zeros + 1, sz), ENTRY(A->appSP2, kappa, kappa), ctx);

        for (k = zeros + 1; k < kappa - 1; k++)
        {
            status |= gr_mul(tmp, ENTRY(mu, kappa, k), ENTRY(r, kappa, k), ctx);
            status |= gr_sub(GR_ENTRY(s, k + 1, sz), GR_ENTRY(s, k, sz), tmp, ctx);
        }

        fmpz_clear(ztmp);
    }
    else
    {
        int i, j, k, test, aa;
        slong max_expo = WORD_MAX;
        fmpz_t ztmp;
        double halfplus, onedothalfplus;
        ulong loops;

        fmpz_init(ztmp);

        aa = (a > zeros) ? a : zeros + 1;

        halfplus = (fl->eta + 0.5) / 2;
        onedothalfplus = 1.0 + halfplus;

        loops = 0;

        do
        {
            fmpz *x;

            test = 0;

            /* ************************************** */
            /* Step2: compute the GSO for stage kappa */
            /* ************************************** */

            for (j = aa; j < kappa; j++)
            {
                if (j > zeros + 2)
                {
                    status |= gr_set_fmpz(rtmp, fmpz_mat_entry(GM, kappa, j), ctx);
                    status |= _gr_vec_dot(rtmp, rtmp, 1, ENTRY(mu, j, zeros + 1), ENTRY(r, kappa, zeros + 1), (j - 1) - (zeros + 1), ctx);
                    status |= gr_mul(tmp, ENTRY(mu, j, j - 1), ENTRY(r, kappa, j - 1), ctx);
                    status |= gr_sub(ENTRY(r, kappa, j), rtmp, tmp, ctx);
                }
                else if (j == zeros + 2)
                {
                    status |= gr_mul(tmp, ENTRY(mu, j, zeros + 1), ENTRY(r, kappa, zeros + 1), ctx);
                    status |= gr_set_fmpz(ENTRY(r, kappa, j), fmpz_mat_entry(GM, kappa, j), ctx);
                    status |= gr_sub(ENTRY(r, kappa, j), ENTRY(r, kappa, j), tmp, ctx);
                }
                else
                    status |= gr_set_fmpz(ENTRY(r, kappa, j), fmpz_mat_entry(GM, kappa, j), ctx);

                status |= gr_div(ENTRY(mu, kappa, j), ENTRY(r, kappa, j), ENTRY(r, j, j), ctx);
            }

            if (loops >= 20)
            {
                slong new_max_expo = WORD_MIN;
                for (j = 0; j < kappa; j++)
                {
                    double mant;
                    slong expo2;
                    status |= gr_get_d_2exp_si(&mant, &expo2, ENTRY(mu, kappa, j), ctx);
                    new_max_expo = FLINT_MAX(new_max_expo, expo2);
                }
                if (new_max_expo > max_expo - SIZE_RED_FAILURE_THRESH)
                {
                    fmpz_clear(ztmp);
                    return -1;
                }
                max_expo = new_max_expo;
            }

            /* **************************** */
            /* Step3--5: compute the X_j's  */
            /* **************************** */

            x = _fmpz_vec_init(kappa - 1 - zeros);
            for (j = kappa - 1; j > zeros; j--)
            {
                /* test of the relaxed size-reduction condition */
                status |= gr_abs(tmp, ENTRY(mu, kappa, j), ctx);

                if (_gr_cmp_d(tmp, halfplus, ctx) > 0)
                {
                    test = 1;

                    /* we consider separately the cases X = +-1 */
                    if (_gr_cmp_d(tmp, onedothalfplus, ctx) <= 0)
                    {
                        int sgn = _gr_sgn(ENTRY(mu, kappa, j), ctx);
                        if (sgn >= 0)   /* in this case, X is 1 */
                        {
                            fmpz_one(x + j);

                            status |= _gr_vec_sub(ENTRY(mu, kappa, zeros + 1),
                                ENTRY(mu, kappa, zeros + 1), ENTRY(mu, j, zeros + 1), j - (zeros + 1), ctx);

                            if (fl->rt == Z_BASIS && B != NULL)
                            {
                                _fmpz_vec_sub(fmpz_mat_row(B, kappa), fmpz_mat_row(B, kappa),
                                              fmpz_mat_row(B, j), n);
                            }
                            if (U != NULL)
                            {
                                _fmpz_vec_sub(fmpz_mat_row(U, kappa),
                                              fmpz_mat_row(U, kappa), fmpz_mat_row(U, j),
                                              U->c);
                            }
                        }
                        else    /* otherwise X is -1 */
                        {
                            fmpz_set_si(x + j, -WORD(1));

                            status |= _gr_vec_add(ENTRY(mu, kappa, zeros + 1),
                                ENTRY(mu, kappa, zeros + 1), ENTRY(mu, j, zeros + 1), j - (zeros + 1), ctx);

                            if (fl->rt == Z_BASIS && B != NULL)
                            {
                                _fmpz_vec_add(fmpz_mat_row(B, kappa),
                                              fmpz_mat_row(B, kappa), fmpz_mat_row(B, j), n);
                            }
                            if (U != NULL)
                            {
                                _fmpz_vec_add(fmpz_mat_row(U, kappa),
                                              fmpz_mat_row(U, kappa), fmpz_mat_row(U, j),
                                              U->c);
                            }
                        }
                    }
                    else        /* we must have |X| >= 2 */
                    {
                        status |= gr_set(tmp, ENTRY(mu, kappa, j), ctx);
                        status |= gr_set_d(rtmp, 0.5, ctx);

                        if (_gr_sgn(tmp, ctx) < 0)
                        {
                            status |= gr_sub(tmp, tmp, rtmp, ctx);
                            status |= gr_ceil(tmp, tmp, ctx);
                        }
                        else
                        {
                            status |= gr_add(tmp, tmp, rtmp, ctx);
                            status |= gr_floor(tmp, tmp, ctx);
                        }

                        /* TODO: consider converting to integer if high precision? */
                        /* TODO: or make nfloat check for small multiplier */
                        status |= _gr_vec_submul_scalar(ENTRY(mu, kappa, zeros + 1), ENTRY(mu, j, zeros + 1), j - (zeros + 1), tmp, ctx);
                        status |= gr_get_fmpz(x + j, tmp, ctx);

                        if (fl->rt == Z_BASIS && B != NULL)
                        {
                            _fmpz_vec_scalar_submul_fmpz(fmpz_mat_row(B, kappa),
                                                         fmpz_mat_row(B, j), n,
                                                         x + j);
                        }
                        if (U != NULL)
                        {
                            _fmpz_vec_scalar_submul_fmpz(fmpz_mat_row(U, kappa),
                                                         fmpz_mat_row(U, j),
                                                         U->c, x + j);
                        }
                    }
                }
            }

            if (test)           /* Anything happened? */
            {
                aa = zeros + 1;

                for (j = zeros + 1; j < kappa; j++)
                {
                    fmpz_pow_ui(ztmp, x + j, 2);
                    fmpz_addmul(fmpz_mat_entry(GM, kappa, kappa),
                                ztmp, fmpz_mat_entry(GM, j, j));

                    fmpz_mul(ztmp, x + j, fmpz_mat_entry(GM, kappa, j));
                    fmpz_mul_2exp(ztmp, ztmp, 1);
                    fmpz_sub(fmpz_mat_entry(GM, kappa, kappa),
                             fmpz_mat_entry(GM, kappa, kappa), ztmp);

                    for (i = zeros + 1; i < j; i++)
                    {
                        fmpz_mul(ztmp, x + i, x + j);
                        fmpz_mul(ztmp, ztmp, fmpz_mat_entry(GM, j, i));
                        fmpz_mul_2exp(ztmp, ztmp, 1);
                        fmpz_add(fmpz_mat_entry(GM, kappa, kappa),
                                 fmpz_mat_entry(GM, kappa, kappa), ztmp);
                    }
                }

                for (i = zeros + 1; i < kappa; i++)
                {
                    for (j = zeros + 1; j <= i; j++)
                        fmpz_submul(fmpz_mat_entry(GM, kappa, i),
                                    x + j, fmpz_mat_entry(GM, i, j));
                    for (j = i + 1; j < kappa; j++)
                        fmpz_submul(fmpz_mat_entry(GM, kappa, i),
                                    x + j, fmpz_mat_entry(GM, j, i));
                }

                for (i = kappa + 1; i < GM->r; i++)
                {
                    for (j = zeros + 1; j < kappa; j++)
                        fmpz_submul(fmpz_mat_entry(GM, i, kappa),
                                    x + j, fmpz_mat_entry(GM, i, j));
                }
            }

            _fmpz_vec_clear(x, kappa - 1 - zeros);
            loops++;
        } while (test);

        status |= gr_set_fmpz(GR_ENTRY(s, zeros + 1, sz), fmpz_mat_entry(GM, kappa, kappa), ctx);

        for (k = zeros + 1; k < kappa - 1; k++)
        {
            status |= gr_mul(tmp, ENTRY(mu, kappa, k), ENTRY(r, kappa, k), ctx);
            status |= gr_sub(GR_ENTRY(s, k + 1, sz), GR_ENTRY(s, k, sz), tmp, ctx);
        }

        fmpz_clear(ztmp);
    }

    /* Debugging */
    GR_MUST_SUCCEED(status);

    if (status != GR_SUCCESS)
        return -1;

    return 0;
}

#undef GM
