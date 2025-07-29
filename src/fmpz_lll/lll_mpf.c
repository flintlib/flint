/*
    Copyright (C) 2009, 2010 William Hart
    Copyright (C) 2009, 2010 Andy Novocin
    Copyright (C) 2014 Abhinav Baid
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <float.h>
#include <gmp.h>
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "fmpz_lll.h"
#include "double_extras.h"


#include "gr.h"
#include "gr_vec.h"
#include "gr_mat.h"
#include "nfloat.h"

static void
fmpz_mat_move_row(fmpz_mat_t A, slong i, slong j)
{
    gr_ctx_t ctx;
    gr_ctx_init_fmpz(ctx);
    GR_MUST_SUCCEED(gr_mat_move_row((gr_mat_struct *) A, i, j, ctx));
}

static int _gr_cmp(gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
{
    int sgn, status;
    status = gr_cmp(&sgn, x, y, ctx);
    GR_MUST_SUCCEED(status);
    return sgn;
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

static int _gr_vec_norm2(gr_ptr res, gr_srcptr vec, slong len, gr_ctx_t ctx)
{
    /* todo */
    return _gr_vec_dot(res, NULL, 0, vec, vec, len, ctx);
}

#ifdef GM
#undef GM
#endif
#define GM ((fl->rt == Z_BASIS) ? A->exactSP : B)

int fmpz_lll_mpf2_with_removal(fmpz_mat_t B, fmpz_mat_t U, flint_bitcnt_t prec, const fmpz_t gs_B, const fmpz_lll_t fl)
{
    int newd = 0;
    int ok = 1;
    fmpz_t rii;
    int with_removal = (gs_B != NULL);

    gr_ctx_t ctx;
    int status = GR_SUCCESS;

#if 0
    gr_ctx_init_mpf(ctx, prec);
#else
    if (nfloat_ctx_init(ctx, prec, 0) != GR_SUCCESS)
        gr_ctx_init_real_float_arf(ctx, prec);
#endif

    slong sz = ctx->sizeof_elem;

#define ENTRY(mat, ii, jj) GR_MAT_ENTRY(mat, ii, jj, sz)
#define ROW(mat, ii) GR_MAT_ENTRY(mat, ii, 0, sz)

    if (fl->rt == Z_BASIS && fl->gt == APPROX)
    {
        int kappa, kappa2, d, n, i, j, zeros, kappamax;
        gr_mat_t mu, r, appB;
        fmpz_gram_t A;
        gr_ptr s, appSPtmp;
        gr_ptr ctt, tmp, rtmp;
        int *alpha;

        n = B->c;
        d = B->r;

        GR_TMP_INIT3(ctt, tmp, rtmp, ctx);

        status |= gr_set_d(ctt, (fl->delta + 1) / 2, ctx);

        alpha = (int *) flint_malloc(d * sizeof(int));

        gr_mat_init(mu, d, d, ctx);
        gr_mat_init(r, d, d, ctx);
        gr_mat_init(appB, d, n, ctx);
        gr_mat_init(A->appSP2, d, d, ctx);

        if (U != NULL)
        {
            if (U->r != d)
            {
                flint_throw(FLINT_ERROR, "(fmpz_lll_mpf*): Incompatible dimensions of capturing matrix.\n");
            }
        }

        GR_TMP_INIT_VEC(s, d, ctx);
        GR_TMP_INIT_VEC(appSPtmp, d, ctx);

        for (i = 0; i < d; i++)
            for (j = 0; j < d; j++)
                status |= gr_set_d(ENTRY(A->appSP2, i, j), DBL_MIN, ctx);

        /* ************************** */
        /* Step1: Initialization Step */
        /* ************************** */

        status |= gr_mat_set_fmpz_mat(appB, B, ctx);

        /* ********************************* */
        /* Step2: Initializing the main loop */
        /* ********************************* */

        kappamax = 0;
        i = 0;

        do
        {
            status |= _gr_vec_norm2(ENTRY(A->appSP2, i, i), ROW(appB, i), n, ctx);
        } while ((_gr_sgn(ENTRY(A->appSP2, i, i), ctx) == 0)
                 && (++i < d));

        zeros = i - 1;          /* all vectors B[i] with i <= zeros are zero vectors */
        kappa = i + 1;
        kappamax = kappa;

        if (zeros < d - 1)
        {
            status |= gr_set(ENTRY(r, i, i), ENTRY(A->appSP2, i, i), ctx);
        }

        for (i = zeros + 1; i < d; i++)
            alpha[i] = 0;

        while (kappa < d)
        {
            int babai_ok = 0;

            if (kappa > kappamax)
                kappamax = kappa;

            /* ********************************** */
            /* Step3: Call to the Babai algorithm */
            /* ********************************** */
            babai_ok =
                fmpz_lll_check_babai_heuristic(kappa, B, U, mu, r, s, appB,
                                               A, alpha[kappa], zeros,
                                               kappamax, n, tmp, rtmp,
                                               ctx, fl);

            if (babai_ok == -1)
            {
                flint_free(alpha);
                GR_TMP_CLEAR3(ctt, tmp, rtmp, ctx);
                gr_mat_clear(mu, ctx);
                gr_mat_clear(r, ctx);
                gr_mat_clear(appB, ctx);
                gr_mat_clear(A->appSP2, ctx);
                GR_TMP_CLEAR_VEC(s, d, ctx);
                GR_TMP_CLEAR_VEC(appSPtmp, d, ctx);
                /* Need to switch to mpf / arb */
                return -1;
            }

            /* ************************************ */
            /* Step4: Success of Lovasz's condition */
            /* ************************************ */

            status |= gr_mul(tmp, ENTRY(r, kappa - 1, kappa - 1), ctt, ctx);

            if (_gr_cmp(tmp, GR_ENTRY(s, kappa - 1, sz), ctx) <= 0)
            {
                alpha[kappa] = kappa;
                status |= gr_mul(tmp, ENTRY(mu, kappa, kappa - 1), ENTRY(r, kappa, kappa - 1), ctx);
                status |= gr_sub(ENTRY(r, kappa, kappa), GR_ENTRY(s, kappa - 1, sz), tmp, ctx);
                kappa++;
            }
            else
            {

                /* ******************************************* */
                /* Step5: Find the right insertion index kappa */
                /* kappa2 remains the initial kappa            */
                /* ******************************************* */

                kappa2 = kappa;

                if (with_removal)
                {
                    if (kappa == d - 1 && gs_B != NULL)
                    {
                        fmpz_init(rii);
                        status |= gr_mul(tmp, ENTRY(mu, kappa, kappa - 1), ENTRY(r, kappa, kappa - 1), ctx);
                        status |= gr_mul_2exp_si(tmp, tmp, 1, ctx);
                        status |= gr_sub(tmp, GR_ENTRY(s, kappa - 1, sz), tmp, ctx);
                        status |= gr_trunc(tmp, tmp, ctx);
                        status |= gr_get_fmpz(rii, tmp, ctx); /* using a heuristic lower bound on the final GS norm */
                        if (fmpz_cmp(rii, gs_B) > 0)
                        {
                            d--;
                        }
                        fmpz_clear(rii);
                        if (kappa >= d)
                        {
                            break;
                        }
                    }
                }

                do
                {
                    kappa--;
                    if (kappa > zeros + 1)
                    {
                        status |= gr_mul(tmp, ENTRY(r, kappa - 1, kappa - 1), ctt, ctx);
                    }
                } while ((kappa >= zeros + 2)
                         && (_gr_cmp(GR_ENTRY(s, kappa - 1, sz), tmp, ctx) <= 0));

                for (i = kappa; i < kappa2; i++)
                    if (kappa <= alpha[i])
                        alpha[i] = kappa;

                for (i = kappa2; i > kappa; i--)
                    alpha[i] = alpha[i - 1];

                for (i = kappa2 + 1; i <= kappamax; i++)
                    if (kappa < alpha[i])
                        alpha[i] = kappa;

                alpha[kappa] = kappa;

                /* ****************************** */
                /* Step6: Update the mu's and r's */
                /* ****************************** */

                status |= gr_mat_move_row(mu, kappa2, kappa, ctx);
                status |= gr_mat_move_row(r, kappa2, kappa, ctx);
                status |= gr_set(ENTRY(r, kappa, kappa), GR_ENTRY(s, kappa, sz), ctx);

                /* ************************ */
                /* Step7: Update B and appB */
                /* ************************ */

                fmpz_mat_move_row(B, kappa2, kappa);

                if (U != NULL)
                    fmpz_mat_move_row(U, kappa2, kappa);

                status |= gr_mat_move_row(appB, kappa2, kappa, ctx);

                /* *************************** */
                /* Step8: Update appSP: tricky */
                /* *************************** */

                status |= _gr_vec_set(appSPtmp, ENTRY(A->appSP2, kappa2, 0), kappa2 + 1, ctx);

                for (i = kappa2 + 1; i <= kappamax; i++)
                    status |= gr_set(GR_ENTRY(appSPtmp, i, sz), ENTRY(A->appSP2, i, kappa2), ctx);

                for (i = kappa2; i > kappa; i--)
                {
                    for (j = 0; j < kappa; j++)
                        status |= gr_set(ENTRY(A->appSP2, i, j), ENTRY(A->appSP2, i - 1, j), ctx);

                    status |= gr_set(ENTRY(A->appSP2, i, kappa), GR_ENTRY(appSPtmp, i - 1, sz), ctx);

                    for (j = kappa + 1; j <= i; j++)
                        status |= gr_set(ENTRY(A->appSP2, i, j), ENTRY(A->appSP2, i - 1, j - 1), ctx);

                    for (j = kappa2 + 1; j <= kappamax; j++)
                        status |= gr_set(ENTRY(A->appSP2, j, i), ENTRY(A->appSP2, j, i - 1), ctx);
                }

                for (i = 0; i < kappa; i++)
                    status |= gr_set(ENTRY(A->appSP2, kappa, i), GR_ENTRY(appSPtmp, i, sz), ctx);

                status |= gr_set(ENTRY(A->appSP2, kappa, kappa), GR_ENTRY(appSPtmp, kappa2, sz), ctx);

                for (i = kappa2 + 1; i <= kappamax; i++)
                    status |= gr_set(ENTRY(A->appSP2, i, kappa), GR_ENTRY(appSPtmp, i, sz), ctx);

                if (_gr_sgn(ENTRY(r, kappa, kappa), ctx) <= 0)
                {
                    zeros++;
                    kappa++;
                    status |= _gr_vec_norm2(ENTRY(A->appSP2, kappa, kappa), ROW(appB, kappa), n, ctx);
                    status |= gr_set(ENTRY(r, kappa, kappa), ENTRY(A->appSP2, kappa, kappa), ctx);
                }

                kappa++;
            }
        }

        if (with_removal)
        {
            if (gs_B != NULL)
            {
                newd = d;
                fmpz_init(rii);
                for (i = d - 1; (i >= 0) && (ok > 0); i--)
                {
                    status |= gr_mul_2exp_si(tmp, ENTRY(r, i, i), -1, ctx);
                    /* rii is the G-S length of ith vector divided by 2 */
                    status |= gr_trunc(tmp, tmp, ctx);
                    status |= gr_get_fmpz(rii, tmp, ctx);

                    if ((ok = fmpz_cmp(rii, gs_B)) > 0)
                    {
                        newd--;
                    }
                }
                fmpz_clear(rii);
            }
        }

        flint_free(alpha);
        GR_TMP_CLEAR3(ctt, tmp, rtmp, ctx);
        gr_mat_clear(mu, ctx);
        gr_mat_clear(r, ctx);
        gr_mat_clear(appB, ctx);
        gr_mat_clear(A->appSP2, ctx);
        GR_TMP_CLEAR_VEC(s, d, ctx);
        GR_TMP_CLEAR_VEC(appSPtmp, d, ctx);
    }
    else
    {
        int kappa, kappa2, d, n, i, j, zeros, kappamax, update_b = 1;
        gr_mat_t mu, r;
        fmpz_gram_t A;
        gr_ptr s;
        gr_ptr ctt, tmp, rtmp;
        int *alpha;

        n = B->c;
        d = B->r;

        GR_TMP_INIT3(ctt, tmp, rtmp, ctx);

        status |= gr_set_d(ctt, (fl->delta + 1) / 2, ctx);

        alpha = (int *) flint_malloc(d * sizeof(int));

        gr_mat_init(mu, d, d, ctx);
        gr_mat_init(r, d, d, ctx);
        if (fl->rt == Z_BASIS)
        {
            fmpz_mat_init(A->exactSP, d, d);
        }

        GR_TMP_INIT_VEC(s, d, ctx);

        if (U != NULL)
        {
            if (U->r != d)
            {
                flint_throw(FLINT_ERROR, "(fmpz_lll_mpf*): Incompatible dimensions of capturing matrix.\n");
            }
            else if (U->c == d && n > d && fmpz_mat_is_one(U))
            {
                update_b = 0;
            }
        }

        if (fl->rt == Z_BASIS)
        {
            fmpz_mat_gram(A->exactSP, B);
        }

        /* ********************************* */
        /* Step2: Initializing the main loop */
        /* ********************************* */

        kappamax = 0;
        i = 0;

        do
            ;
        while ((fmpz_cmp_ui(fmpz_mat_entry(GM, i, i), 0) <= 0) && (++i < d));

        zeros = i - 1;          /* all vectors B[i] with i <= zeros are zero vectors */
        kappa = i + 1;
        kappamax = kappa;

        if (zeros < d - 1)
        {
            status |= gr_set_fmpz(ENTRY(r, i, i), fmpz_mat_entry(GM, i, i), ctx);
        }

        for (i = zeros + 1; i < d; i++)
            alpha[i] = 0;

        while (kappa < d)
        {
            int babai_ok = 0;

            if (kappa > kappamax)
                kappamax = kappa;

            /* ********************************** */
            /* Step3: Call to the Babai algorithm */
            /* ********************************** */
            babai_ok =
                fmpz_lll_check_babai_heuristic(kappa, (update_b ? B : NULL), U,
                                               mu, r, s, NULL, A, alpha[kappa],
                                               zeros, kappamax, n, tmp, rtmp,
                                               ctx, fl);

            if (babai_ok == -1)
            {
                flint_free(alpha);
                GR_TMP_CLEAR3(ctt, tmp, rtmp, ctx);
                gr_mat_clear(mu, ctx);
                gr_mat_clear(r, ctx);
                if (fl->rt == Z_BASIS)
                {
                    fmpz_mat_clear(A->exactSP);
                }
                GR_TMP_CLEAR_VEC(s, d, ctx);
                /* Need to switch to mpf / arb */
                return -1;
            }

            /* ************************************ */
            /* Step4: Success of Lovasz's condition */
            /* ************************************ */

            status |= gr_mul(tmp, ENTRY(r, kappa - 1, kappa - 1), ctt, ctx);

            if (_gr_cmp(tmp, GR_ENTRY(s, kappa - 1, sz), ctx) <= 0)
            {
                alpha[kappa] = kappa;
                status |= gr_mul(tmp, ENTRY(mu, kappa, kappa - 1), ENTRY(r, kappa, kappa - 1), ctx);
                status |= gr_sub(ENTRY(r, kappa, kappa), GR_ENTRY(s, kappa - 1, sz), tmp, ctx);
                kappa++;
            }
            else
            {

                /* ******************************************* */
                /* Step5: Find the right insertion index kappa */
                /* kappa2 remains the initial kappa            */
                /* ******************************************* */

                kappa2 = kappa;

                if (with_removal)
                {
                    if (kappa == d - 1 && gs_B != NULL)
                    {
                        fmpz_init(rii);
                        status |= gr_mul(tmp, ENTRY(mu, kappa, kappa - 1), ENTRY(r, kappa, kappa - 1), ctx);
                        status |= gr_mul_2exp_si(tmp, tmp, 1, ctx);
                        status |= gr_sub(tmp, GR_ENTRY(s, kappa - 1, sz), tmp, ctx);
                        status |= gr_trunc(tmp, tmp, ctx);
                        status |= gr_get_fmpz(rii, tmp, ctx); /* using a heuristic lower bound on the final GS norm */
                        if (fmpz_cmp(rii, gs_B) > 0)
                        {
                            d--;
                        }
                        fmpz_clear(rii);
                        if (kappa >= d)
                        {
                            break;
                        }
                    }
                }

                do
                {
                    kappa--;
                    if (kappa > zeros + 1)
                    {
                        status |= gr_mul(tmp, ENTRY(r, kappa - 1, kappa - 1), ctt, ctx);
                    }
                } while ((kappa >= zeros + 2)
                         && (_gr_cmp(GR_ENTRY(s, kappa - 1, sz), tmp, ctx) <= 0));

                for (i = kappa; i < kappa2; i++)
                    if (kappa <= alpha[i])
                        alpha[i] = kappa;

                for (i = kappa2; i > kappa; i--)
                    alpha[i] = alpha[i - 1];

                for (i = kappa2 + 1; i <= kappamax; i++)
                    if (kappa < alpha[i])
                        alpha[i] = kappa;

                alpha[kappa] = kappa;

                /* ****************************** */
                /* Step6: Update the mu's and r's */
                /* ****************************** */

                status |= gr_mat_move_row(mu, kappa2, kappa, ctx);
                status |= gr_mat_move_row(r, kappa2, kappa, ctx);
                status |= gr_set(ENTRY(r, kappa, kappa), GR_ENTRY(s, kappa, sz), ctx);

                /* *************** */
                /* Step7: Update B */
                /* *************** */

                if (fl->rt == Z_BASIS && update_b)
                    fmpz_mat_move_row(B, kappa2, kappa);

                if (U != NULL)
                    fmpz_mat_move_row(U, kappa2, kappa);

                /* ********************* */
                /* Step8: Update exactSP */
                /* ********************* */

                for (j = kappa2; j > kappa; j--)
                {
                    for (i = kappa2; i < d; i++)
                        fmpz_swap(fmpz_mat_entry(GM, i, j),
                                  fmpz_mat_entry(GM, i, j - 1));
                    for (i = 0; i < kappa; i++)
                        fmpz_swap(fmpz_mat_entry(GM, j, i),
                                  fmpz_mat_entry(GM, j - 1, i));
                }
                for (j = kappa2; j > kappa; j--)
                {
                    for (i = j; i > kappa; i--)
                        fmpz_swap(fmpz_mat_entry(GM, j, i),
                                  fmpz_mat_entry(GM, j - 1, i - 1));
                }
                for (j = 0; 2 * j < kappa2 - kappa; j++)
                    fmpz_swap(fmpz_mat_entry(GM, kappa + j, kappa),
                              fmpz_mat_entry(GM, kappa2 - j, kappa));

                if (_gr_sgn(ENTRY(r, kappa, kappa), ctx) <= 0)
                {
                    zeros++;
                    kappa++;
                    status |= gr_set_fmpz(ENTRY(r, kappa, kappa), fmpz_mat_entry(GM, kappa, kappa), ctx);
                }

                kappa++;
            }
        }

        if (fl->rt == GRAM)
        {
            for (i = 0; i < B->r - 1; i++)
            {
                for (j = i + 1; j < B->r; j++)
                {
                    fmpz_set(fmpz_mat_entry(B, i, j), fmpz_mat_entry(B, j, i));
                }
            }
        }
        else if (!update_b)
        {
            fmpz_mat_mul(B, U, B);
        }


        if (with_removal)
        {
            if (gs_B != NULL)
            {
                newd = d;
                fmpz_init(rii);
                for (i = d - 1; (i >= 0) && (ok > 0); i--)
                {
                    status |= gr_mul_2exp_si(tmp, ENTRY(r, i, i), -1, ctx);
                    status |= gr_trunc(tmp, tmp, ctx);
                    /* rii is the G-S length of ith vector divided by 2 */
                    status |= gr_get_fmpz(rii, tmp, ctx);
                    if ((ok = fmpz_cmp(rii, gs_B)) > 0)
                    {
                        newd--;
                    }
                }
                fmpz_clear(rii);
            }
        }

        flint_free(alpha);
        GR_TMP_CLEAR3(ctt, tmp, rtmp, ctx);
        gr_mat_clear(mu, ctx);
        gr_mat_clear(r, ctx);
        if (fl->rt == Z_BASIS)
        {
            fmpz_mat_clear(A->exactSP);
        }
        GR_TMP_CLEAR_VEC(s, d, ctx);
    }

    /* Debugging */
    GR_MUST_SUCCEED(status);

    return newd;
}

int fmpz_lll_mpf2(fmpz_mat_t B, fmpz_mat_t U, flint_bitcnt_t prec, const fmpz_lll_t fl)
{
    return fmpz_lll_mpf2_with_removal(B, U, prec, NULL, fl);
}

int
fmpz_lll_mpf(fmpz_mat_t B, fmpz_mat_t U, const fmpz_lll_t fl)
{
    flint_bitcnt_t prec = 0;
    int result, num_loops = 0;

    do
    {
        /* todo: why not always double prec? */
        if (num_loops < 20)
            prec += 64;
        else
            prec *= 2;

        result = fmpz_lll_mpf2(B, U, prec, fl);

        num_loops++;
    } while (((result == -1) || (!fmpz_lll_is_reduced(B, fl, prec)))
             && (prec < UWORD_MAX));

    return result;
}

int
fmpz_lll_mpf_with_removal(fmpz_mat_t B, fmpz_mat_t U, const fmpz_t gs_B,
                          const fmpz_lll_t fl)
{
    flint_bitcnt_t prec = 0;
    int result, num_loops = 0;

    do
    {
        /* todo: why not always double prec? */
        if (num_loops < 20)
            prec += 64;
        else
            prec *= 2;

        result = fmpz_lll_mpf2_with_removal(B, U, prec, gs_B, fl);

        num_loops++;
    } while (((result == -1)
              ||
              (!fmpz_lll_is_reduced_with_removal(B, fl, gs_B, result, prec)))
             && (prec < UWORD_MAX));

    return result;
}

