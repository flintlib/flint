/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr.h"
#include "gr_vec.h"
#include "gr_mat.h"
#include "gr_poly.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"

/* Some helper methods to compute s += a * c in-place for segments of vectors,
   allowing different lengths for s and a and updating the length of s. */

static int
_gr_poly_addmul_1(slong * new_slen, gr_ptr s, slong slen, gr_srcptr a, slong alen, gr_srcptr c, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    if (gr_is_one(c, ctx) == T_TRUE)
    {
        status |= _gr_vec_add(s, s, a, FLINT_MIN(slen, alen), ctx);
        if (alen > slen)
            status |= _gr_vec_set(GR_ENTRY(s, slen, sz), GR_ENTRY(a, slen, sz), alen - slen, ctx);
    }
    else if (gr_is_neg_one(c, ctx) == T_TRUE)
    {
        status |= _gr_vec_sub(s, s, a, FLINT_MIN(slen, alen), ctx);
        if (alen > slen)
            status |= _gr_vec_neg(GR_ENTRY(s, slen, sz), GR_ENTRY(a, slen, sz), alen - slen, ctx);
    }
    else
    {
        status |= _gr_vec_addmul_scalar(s, a, FLINT_MIN(slen, alen), c, ctx);
        if (alen > slen)
            status |= _gr_vec_mul_scalar(GR_ENTRY(s, slen, sz), GR_ENTRY(a, slen, sz), alen - slen, c, ctx);
    }

    *new_slen = FLINT_MAX(slen, alen);
    return status;
}

static int
_gr_poly_addmul_1_shift(slong * new_slen, nn_ptr s, slong slen, slong shift, nn_srcptr a, slong alen, gr_srcptr c, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    if (shift > slen)
    {
        status |= _gr_vec_zero(GR_ENTRY(s, slen, sz), shift - slen, ctx);
        slen = shift;
    }

    if (alen > 0)
        status |= _gr_poly_addmul_1(new_slen, GR_ENTRY(s, shift, sz), slen - shift, a, alen, c, ctx);

    *new_slen = FLINT_MAX(slen, shift + alen);
    return status;
}

static int
_gr_poly_addmul_1_fmpz(slong * new_slen, gr_ptr s, slong slen, gr_srcptr a, slong alen, const fmpz_t c, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    if (fmpz_is_one(c))
    {
        status |= _gr_vec_add(s, s, a, FLINT_MIN(slen, alen), ctx);
        if (alen > slen)
            status |= _gr_vec_set(GR_ENTRY(s, slen, sz), GR_ENTRY(a, slen, sz), alen - slen, ctx);
    }
    else if (fmpz_equal_si(c, -1))
    {
        status |= _gr_vec_sub(s, s, a, FLINT_MIN(slen, alen), ctx);
        if (alen > slen)
            status |= _gr_vec_neg(GR_ENTRY(s, slen, sz), GR_ENTRY(a, slen, sz), alen - slen, ctx);
    }
    else
    {
        status |= _gr_vec_addmul_scalar_fmpz(s, a, FLINT_MIN(slen, alen), c, ctx);
        if (alen > slen)
            status |= _gr_vec_mul_scalar_fmpz(GR_ENTRY(s, slen, sz), GR_ENTRY(a, slen, sz), alen - slen, c, ctx);
    }

    *new_slen = FLINT_MAX(slen, alen);
    return status;
}

static int
_gr_poly_addmul_1_fmpz_shift(slong * new_slen, nn_ptr s, slong slen, slong shift, nn_srcptr a, slong alen, const fmpz_t c, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    if (shift > slen)
    {
        status |= _gr_vec_zero(GR_ENTRY(s, slen, sz), shift - slen, ctx);
        slen = shift;
    }

    if (alen > 0)
        status |= _gr_poly_addmul_1_fmpz(new_slen, GR_ENTRY(s, shift, sz), slen - shift, a, alen, c, ctx);

    *new_slen = FLINT_MAX(slen, shift + alen);
    return status;
}

/* Version using gr_mat interpolation matrix, suitable e.g. over Z/nZ.
   To do: this uses ctx itself for the interpolation matrix, but in some cases
   it would be more appropriate to use either a subring or an extension ring. */
static int
_gr_poly_mullow_toom_serial_gr_mat(gr_ptr res, gr_srcptr A, slong len1, gr_srcptr B, slong len2, slong trunc, slong r, slong m, slong mlen1, slong mlen2, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    slong i, j;
    gr_ptr T, U, V;
    slong lenR, lenT, lenU, lenV;
    int squaring;

    gr_mat_t M, Minv;

    len1 = FLINT_MIN(len1, trunc);
    len2 = FLINT_MIN(len2, trunc);

    squaring = (A == B) && (len1 == len2);

    gr_mat_init(M, r, r, ctx);
    gr_mat_init(Minv, r, r, ctx);

    for (i = 0; i < r; i++)
    {
        if (i == 0)
        {
            /* First evaluation point is 0 */
            status |= gr_one(gr_mat_entry_ptr(M, i, 0, ctx), ctx);
            status |= _gr_vec_zero(gr_mat_entry_ptr(M, i, 1, ctx), r - 1, ctx);
        }
        else if (i == r - 1)
        {
            /* Last evaluation point is inf */
            status |= _gr_vec_zero(gr_mat_entry_ptr(M, i, 0, ctx), r - 1, ctx);
            status |= gr_one(gr_mat_entry_ptr(M, i, r - 1, ctx), ctx);
        }
        else
        {
            /* Middle evaluation points are 1,-1,2,-2,3,-3,... */
            gr_ptr C = gr_mat_entry_ptr(M, r - 1, 0, ctx);  /* scratch space */
            status |= gr_set_si(C, (i % 2) ? (i + 1) / 2 : -((i + 1) / 2), ctx);
            status |= _gr_vec_set_powers(gr_mat_entry_ptr(M, i, 0, ctx), C, r, ctx);
        }
    }

    status = gr_mat_inv(Minv, M, ctx);
    if (status != GR_SUCCESS)
    {
        status = GR_UNABLE;
        goto cleanup;
    }

    T = gr_heap_init_vec(m, ctx);
    if (!squaring)
        U = gr_heap_init_vec(m, ctx);
    V = gr_heap_init_vec((2 * m - 1), ctx);

    lenR = 0;

    for (i = 0; i < r; i++)
    {
        gr_srcptr TT, UU;

        if (i == 0)  /* Evaluate at 0 */
        {
            lenT = FLINT_MIN(m, len1);
            lenU = FLINT_MIN(m, len2);
            TT = A;
            UU = B;
        }
        else if (i == r - 1)  /* Evaluate at inf */
        {
            if (mlen1 + mlen2 - 1 >= r)
            {
                lenT = FLINT_MIN(m, len1 - (mlen1 - 1) * m);
                lenU = FLINT_MIN(m, len2 - (mlen2 - 1) * m);

                TT = GR_ENTRY(A, (mlen1 - 1) * m, sz);
                UU = GR_ENTRY(B, (mlen2 - 1) * m, sz);
            }
            else
            {
                lenT = lenU = 0;
                continue;
            }
        }
        else
        {
            /* M * A */
            lenT = 0;
            for (j = 0; j * m < len1; j++)
            {
                if (gr_is_zero(gr_mat_entry_srcptr(M, i, j, ctx), ctx) == T_TRUE)
                    continue;

                status |= _gr_poly_addmul_1(&lenT, T, lenT, GR_ENTRY(A, j * m, sz),
                    FLINT_MIN(m, len1 - j * m), gr_mat_entry_srcptr(M, i, j, ctx), ctx);
            }

            /* M * B */
            if (!squaring)
            {
                lenU = 0;
                for (j = 0; j * m < len2; j++)
                {
                    if (gr_is_zero(gr_mat_entry_srcptr(M, i, j, ctx), ctx) == T_TRUE)
                        continue;

                    status |= _gr_poly_addmul_1(&lenU, U, lenU, GR_ENTRY(B, j * m, sz),
                        FLINT_MIN(m, len2 - j * m), gr_mat_entry_srcptr(M, i, j, ctx), ctx);
                }
            }

            TT = T;
            UU = U;
        }

        if (squaring)
        {
            UU = TT;
            lenU = lenT;
        }

        lenV = lenT + lenU - 1;
        status |= _gr_poly_mul(V, TT, lenT, UU, lenU, ctx);

        /* Minv^T * C */
        for (j = 0; j < r && j * m < trunc; j++)
        {
            if (gr_is_zero(gr_mat_entry_srcptr(Minv, j, i, ctx), ctx) == T_TRUE)
                continue;

            status |= _gr_poly_addmul_1_shift(&lenR, res, lenR, j * m, V,
                FLINT_MIN(lenV, trunc - j * m), gr_mat_entry_srcptr(Minv, j, i, ctx), ctx);
        }
    }

    gr_heap_clear_vec(T, m, ctx);
    if (!squaring)
        gr_heap_clear_vec(U, m, ctx);
    gr_heap_clear_vec(V, (2 * m - 1), ctx);

cleanup:
    gr_mat_clear(M, ctx);
    gr_mat_clear(Minv, ctx);

    return status;
}

/* Version using an fmpz_mat interpolation matrix with a common denominator,
   suitable over integral domains in characteristic 0. Todo: consider using
   reciprocal evaluation points 1/2, etc.; this should give smaller
   denominators. */

static int
_gr_poly_mullow_toom_serial_fmpz_mat(gr_ptr res, gr_srcptr A, slong len1, gr_srcptr B, slong len2, slong trunc, slong r, slong m, slong mlen1, slong mlen2, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    slong i, j;
    gr_ptr T, U, V;
    slong lenR, lenT, lenU, lenV;
    int squaring;

    fmpz_mat_t M, Minv;
    fmpz_t den;

    squaring = (A == B) && (len1 == len2);

    fmpz_mat_init(M, r, r);
    fmpz_mat_init(Minv, r, r);
    fmpz_init(den);

    for (i = 0; i < r; i++)
    {
        if (i == 0)
        {
            /* First evaluation point is 0 */
            fmpz_one(fmpz_mat_entry(M, i, 0));
            _fmpz_vec_zero(fmpz_mat_entry(M, i, 1), r - 1);
        }
        else if (i == r - 1)
        {
            /* Last evaluation point is inf */
            _fmpz_vec_zero(fmpz_mat_entry(M, i, 0), r - 1);
            fmpz_one(fmpz_mat_entry(M, i, r - 1));
        }
        else
        {
            /* Middle evaluation points are 1,-1,2,-2,3,-3,... */
            slong C = (i % 2) ? (i + 1) / 2 : -((i + 1) / 2);

            for (j = 0; j < r; j++)
            {
                if (j == 0)
                    fmpz_one(fmpz_mat_entry(M, i, j));
                else if (j == 1)
                    fmpz_set_si(fmpz_mat_entry(M, i, j), C);
                else
                    fmpz_mul_si(fmpz_mat_entry(M, i, j), fmpz_mat_entry(M, i, j - 1), C);
            }
        }
    }

    if (!fmpz_mat_inv(Minv, den, M))
        flint_throw(FLINT_ERROR, "_gr_poly_mullow_toom_serial: interpolation matrix not invertible over Z");

    T = gr_heap_init_vec(m, ctx);
    if (!squaring)
        U = gr_heap_init_vec(m, ctx);
    else
        U = NULL;
    V = gr_heap_init_vec((2 * m - 1), ctx);

    lenR = 0;

    for (i = 0; i < r; i++)
    {
        gr_srcptr TT, UU;

        if (i == 0)  /* Evaluate at 0 */
        {
            lenT = FLINT_MIN(m, len1);
            lenU = FLINT_MIN(m, len2);
            TT = A;
            UU = B;
        }
        else if (i == r - 1)  /* Evaluate at inf */
        {
            if (mlen1 + mlen2 - 1 >= r)
            {
                lenT = FLINT_MIN(m, len1 - (mlen1 - 1) * m);
                lenU = FLINT_MIN(m, len2 - (mlen2 - 1) * m);

                TT = GR_ENTRY(A, (mlen1 - 1) * m, sz);
                UU = GR_ENTRY(B, (mlen2 - 1) * m, sz);
            }
            else
            {
                lenT = lenU = 0;
                continue;
            }
        }
        else
        {
            /* M * A */
            lenT = 0;
            for (j = 0; j * m < len1; j++)
            {
                if (fmpz_is_zero(fmpz_mat_entry(M, i, j)))
                    continue;

                status |= _gr_poly_addmul_1_fmpz(&lenT, T, lenT, GR_ENTRY(A, j * m, sz),
                    FLINT_MIN(m, len1 - j * m), fmpz_mat_entry(M, i, j), ctx);
            }

            /* M * B */
            if (!squaring)
            {
                lenU = 0;
                for (j = 0; j * m < len2; j++)
                {
                    if (fmpz_is_zero(fmpz_mat_entry(M, i, j)))
                        continue;

                    status |= _gr_poly_addmul_1_fmpz(&lenU, U, lenU, GR_ENTRY(B, j * m, sz),
                        FLINT_MIN(m, len2 - j * m), fmpz_mat_entry(M, i, j), ctx);
                }
            }

            TT = T;
            UU = U;
        }

        if (squaring)
        {
            UU = TT;
            lenU = lenT;
        }

        lenV = lenT + lenU - 1;
        status |= _gr_poly_mul(V, TT, lenT, UU, lenU, ctx);

        /* Minv^T * C */
        for (j = 0; j < r && j * m < trunc; j++)
        {
            if (fmpz_is_zero(fmpz_mat_entry(Minv, j, i)))
                continue;

            status |= _gr_poly_addmul_1_fmpz_shift(&lenR, res, lenR, j * m, V,
                FLINT_MIN(lenV, trunc - j * m), fmpz_mat_entry(Minv, j, i), ctx);
        }
    }

    if (!fmpz_is_one(den))
        status |= _gr_vec_divexact_scalar_fmpz(res, res, trunc, den, ctx);

    gr_heap_clear_vec(T, m, ctx);
    if (!squaring)
        gr_heap_clear_vec(U, m, ctx);
    gr_heap_clear_vec(V, (2 * m - 1), ctx);

    fmpz_mat_clear(M);
    fmpz_mat_clear(Minv);
    fmpz_clear(den);

    return status;
}

static ulong _cdiv(ulong a, ulong b)
{
    return (a + b - 1) / b;
}


int
_gr_poly_mullow_toom_serial(gr_ptr res, gr_srcptr A, slong len1, gr_srcptr B, slong len2, slong trunc, slong num_points, slong splitting_degree, gr_ctx_t ctx)
{
    slong r, m, m1, m2, mlen1, mlen2;

    len1 = FLINT_MIN(len1, trunc);
    len2 = FLINT_MIN(len2, trunc);

    if ((num_points == 0) == (splitting_degree == 0))
    {
        flint_throw(FLINT_ERROR, "_gr_poly_mullow_toom_serial: exactly one of num_points and splitting_degree must be given");
    }

    if (num_points != 0)
    {
        r = num_points;
        r = FLINT_MIN(r, len1 + len2 - 1);
        r = FLINT_MAX(r, 1);

        /* Viewing A and B as polynomials in x^m of length mlen1 and mlen2,
           select m along with mlen1, mlen2 such that:

           mlen1 + mlen2 - 1 <= r
           mlen1 * m >= len1
           mlen2 * m >= len2
        */

        /* choose mlen1 ~= r * len1/(len1+len2) */
        mlen1 = _cdiv(len1 * r, len1 + len2);
        mlen2 = r + 1 - mlen1;

        m1 = _cdiv(len1, mlen1);
        m2 = _cdiv(len2, mlen2);
        m = FLINT_MAX(m1, m2);
        mlen1 = _cdiv(len1, m);
        mlen2 = _cdiv(len2, m);
    }
    else
    {
        m = splitting_degree;
        mlen1 = _cdiv(len1, m);
        mlen2 = _cdiv(len2, m);
        r = mlen1 + mlen2 - 1;
    }

    // flint_printf("r = %wd   len1 = %wd  len2 = %wd  m = %wd  mlen1 = %wd  mlen2 = %wd\n", r, len1, len2, m, mlen1, mlen2);

    if (gr_ctx_is_finite_characteristic(ctx) == T_TRUE)
        return _gr_poly_mullow_toom_serial_gr_mat(res, A, len1, B, len2, trunc, r, m, mlen1, mlen2, ctx);
    else
        return _gr_poly_mullow_toom_serial_fmpz_mat(res, A, len1, B, len2, trunc, r, m, mlen1, mlen2, ctx);
}

int
gr_poly_mullow_toom_serial(gr_poly_t res, const gr_poly_t poly1,
                                            const gr_poly_t poly2,
                                                slong n, slong num_points, slong splitting_degree, gr_ctx_t ctx)
{
    slong len_out;
    int status;

    if (poly1->length == 0 || poly2->length == 0 || n == 0)
        return gr_poly_zero(res, ctx);

    len_out = poly1->length + poly2->length - 1;
    n = FLINT_MIN(n, len_out);

    if (res == poly1 || res == poly2)
    {
        gr_poly_t t;
        gr_poly_init2(t, n, ctx);
        status = _gr_poly_mullow_toom_serial(t->coeffs, poly1->coeffs, poly1->length, poly2->coeffs, poly2->length, n, num_points, splitting_degree, ctx);
        gr_poly_swap(res, t, ctx);
        gr_poly_clear(t, ctx);
    }
    else
    {
        gr_poly_fit_length(res, n, ctx);
        status = _gr_poly_mullow_toom_serial(res->coeffs, poly1->coeffs, poly1->length, poly2->coeffs, poly2->length, n, num_points, splitting_degree, ctx);
    }

    _gr_poly_set_length(res, n, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}

