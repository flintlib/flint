/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "thread_support.h"
#include "longlong.h"
#include "mpn_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "gr.h"
#include "gr_mat.h"

static void
_si_vec_sub(slong * s, const slong * x, slong n)
{
    slong i;

    for (i = 0; i < n; i++)
        s[i] -= x[i];
}

static void
_si_vec_add(slong * s, const slong * x, slong n)
{
    slong i;

    for (i = 0; i < n; i++)
        s[i] += x[i];
}

#define SET_1 \
        t = a; if (a < 0) { negative = !negative; t = -t; } \
        p[0] = t; \

#define MULTIPLY_1(prod) \
        a = (prod); t = a; if (a < 0) { negative = !negative; t = -t; } \
        cy = p[plen] = mpn_mul_1(p, p, plen, t); \
        plen += (cy != 0); \

/* {res, rn} += (-1)^negative * (x[0] * ... * x[n-1]) assuming that |x[i]| < 2^sbits*/
static void
_nn_addproduct_1(nn_ptr res, slong rn, const slong * x, slong n, slong sbits, int negative)
{
    slong i, plen = 1, a = 0;
    ulong p[FLINT_BITS], t, cy;

    if (sbits < FLINT_BITS / 8)
    {
        switch (n % 8)
        {
            case 1: a = x[n - 1]; n -= 1; break;
            case 2: a = x[n - 1] * x[n - 2]; n -= 2; break;
            case 3: a = x[n - 1] * x[n - 2] * x[n - 3]; n -= 3; break;
            case 4: a = x[n - 1] * x[n - 2] * x[n - 3] * x[n - 4]; n -= 4; break;
            case 5: a = x[n - 1] * x[n - 2] * x[n - 3] * x[n - 4] * x[n - 5]; n -= 5; break;
            case 6: a = x[n - 1] * x[n - 2] * x[n - 3] * x[n - 4] * x[n - 5] * x[n - 6]; n -= 6; break;
            case 7: a = x[n - 1] * x[n - 2] * x[n - 3] * x[n - 4] * x[n - 5] * x[n - 6] * x[n - 7]; n -= 7; break;
            case 0: a = x[n - 1] * x[n - 2] * x[n - 3] * x[n - 4] * x[n - 5] * x[n - 6] * x[n - 7] * x[n - 8]; n -= 8; break;
        }

        SET_1

        for (i = 0; i < n; i += 8)
        {
            MULTIPLY_1(x[i] * x[i + 1] * x[i + 2] * x[i + 3] * x[i + 4] * x[i + 5] * x[i + 6] * x[i + 7])
        }
    }
    else if (sbits < FLINT_BITS / 4)
    {
        switch (n % 4)
        {
            case 1: a = x[n - 1]; n -= 1; break;
            case 2: a = x[n - 1] * x[n - 2]; n -= 2; break;
            case 3: a = x[n - 1] * x[n - 2] * x[n - 3]; n -= 3; break;
            case 0: a = x[n - 1] * x[n - 2] * x[n - 3] * x[n - 4]; n -= 4; break;
        }

        SET_1

        for (i = 0; i < n; i += 4)
        {
            MULTIPLY_1(x[i] * x[i + 1] * x[i + 2] * x[i + 3])
        }
    }
    else if (sbits < FLINT_BITS / 2)
    {
        switch (n % 2)
        {
            case 1: a = x[n - 1]; n -= 1; break;
            case 0: a = x[n - 1] * x[n - 2]; n -= 2; break;
        }

        SET_1

        for (i = 0; i < n; i += 2)
        {
            MULTIPLY_1(x[i] * x[i + 1])
        }
    }
    else
    {
        a = x[n - 1]; n -= 1;
        SET_1

        for (i = 0; i < n; i++)
        {
            MULTIPLY_1(x[i])
        }
    }

    if (negative)
        mpn_sub(res, res, rn, p, plen);
    else
        mpn_add(res, res, rn, p, plen);
}

static void
fmpz_mat_permanent_glynn1(fmpz_t res, const fmpz_mat_t A, slong sbits, slong rn)
{
    slong n = A->r;
    slong i, j, jprev, d;
    slong *s, *A2;
    nn_ptr rr;
    TMP_INIT;

    TMP_START;

    s = TMP_ALLOC(sizeof(slong) * n);
    A2 = TMP_ALLOC(sizeof(slong) * n * n);
    rr = TMP_ALLOC(sizeof(ulong) * rn);
    flint_mpn_zero(rr, rn);

    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            A2[j * n + i] = *fmpz_mat_entry(A, i, j) << 1;

    for (i = 0; i < n; i++)
    {
        s[i] = 0;
        for (j = 0; j < n; j++)
            s[i] += *fmpz_mat_entry(A, i, j);
    }

    _nn_addproduct_1(rr, rn, s, n, sbits, 0);

    jprev = 0;
    for (i = 1; i < (WORD(1) << (n - 1)); i++)
    {
        j = i ^ (i >> 1);

        if (j > jprev)
        {
            d = FLINT_BIT_COUNT(j - jprev) - 1;
            _si_vec_sub(s, A2 + d * n, n);
        }
        else
        {
            d = FLINT_BIT_COUNT(jprev - j) - 1;
            _si_vec_add(s, A2 + d * n, n);
        }

        _nn_addproduct_1(rr, rn, s, n, sbits, i & 1);
        jprev = j;
    }

    fmpz_set_signed_ui_array(res, rr, rn);
    fmpz_tdiv_q_2exp(res, res, n - 1);

    TMP_END;
}

typedef struct
{
    const fmpz_mat_struct * A;
    slong * A2;
    fmpz * rp;
    slong istart;
    slong istop;
    slong sbits;
    slong rn;
}
glynn_args_t;

static void
glynn_worker1(slong chunk, void * args1)
{
    glynn_args_t * args = (glynn_args_t *) args1;
    const fmpz_mat_struct * A = (const fmpz_mat_struct *) args[chunk].A;
    const slong * A2 = args[chunk].A2;
    slong n = A->r;
    slong istart = args[chunk].istart, istop = args[chunk].istop;
    fmpz * rp;
    slong i, j, jprev, k, l, d;
    slong * s;
    nn_ptr rr;

    slong sbits = args->sbits;
    slong rn = args->rn;

    rp = args[chunk].rp;

    s = flint_malloc(sizeof(slong) * n);
    rr = flint_calloc(rn, sizeof(ulong));

    for (i = istart; i < istop; i++)
    {
        if (i == 0)
        {
            for (j = 0; j < n; j++)
            {
                s[j] = 0;
                for (k = 0; k < n; k++)
                    s[j] += *fmpz_mat_entry(A, j, k);
            }

            _nn_addproduct_1(rr, rn, s, n, sbits, 0);
        }
        else
        {
            j = i ^ (i >> 1);

            if (i == istart)
            {
                for (k = 0; k < n; k++)
                {
                    s[k] = 0;

                    for (l = 0; l < n; l++)
                    {
                        if (j & (WORD(1) << l))
                            s[k] -= *fmpz_mat_entry(A, k, l);
                        else
                            s[k] += *fmpz_mat_entry(A, k, l);
                    }
                }
            }
            else
            {
                jprev = (i - 1) ^ ((i - 1) >> 1);

                if (j > jprev)
                {
                    d = FLINT_BIT_COUNT(j - jprev) - 1;
                    _si_vec_sub(s, A2 + d * n, n);
                }
                else
                {
                    d = FLINT_BIT_COUNT(jprev - j) - 1;
                    _si_vec_add(s, A2 + d * n, n);
                }
            }

            _nn_addproduct_1(rr, rn, s, n, sbits, i & 1);
        }
    }

    fmpz_set_signed_ui_array(rp, rr, rn);

    flint_free(s);
    flint_free(rr);
}

static void
fmpz_mat_permanent_glynn1_threaded(fmpz_t res, const fmpz_mat_t A, slong sbits, slong rn)
{
    slong n = A->r;
    slong * A2;
    slong chunk, chunks;
    slong i, j, start, stop;
    fmpz * r;
    glynn_args_t * args;

    chunks = flint_get_num_available_threads();

    A2 = flint_malloc(n * n * sizeof(slong));
    r = _fmpz_vec_init(chunks);

    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            A2[j * n + i] = *fmpz_mat_entry(A, i, j) << 1;

    start = 0;
    stop = WORD(1) << (n - 1);

    args = flint_malloc(sizeof(glynn_args_t) * chunks);

    for (chunk = 0; chunk < chunks; chunk++)
    {
        slong istart, istop;

        istart = chunk * (stop - start + (chunks - 1)) / chunks;
        istop  = (chunk + 1) * (stop - start + (chunks - 1)) / chunks;
        istop = FLINT_MIN(istop, stop);

        args[chunk].A = A;
        args[chunk].A2 = A2;
        args[chunk].rp = r + chunk;
        args[chunk].istart = istart;
        args[chunk].istop = istop;
        args[chunk].sbits = sbits;
        args[chunk].rn = rn;
    }

    flint_parallel_do(glynn_worker1, args, chunks, chunks, FLINT_PARALLEL_UNIFORM);

    _fmpz_vec_sum(res, r, chunks);
    fmpz_tdiv_q_2exp(res, res, n - 1);

    flint_free(args);

    flint_free(A2);
    _fmpz_vec_clear(r, chunks);
}

static int
_fmpz_mat_permanent_generic(fmpz_t res, const fmpz_mat_t A)
{
    slong n = A->r;
    gr_ctx_t ctx;
    int status;

    gr_ctx_init_fmpz(ctx);

    if (n <= 4)
        status = gr_mat_permanent_cofactor(res, (const gr_mat_struct *) A, ctx);
    else if (n <= 10 || flint_get_num_available_threads() == 1)
        status = gr_mat_permanent_glynn(res, (const gr_mat_struct *) A, ctx);
    else
        status = gr_mat_permanent_glynn_threaded(res, (const gr_mat_struct *) A, ctx);

    return (status == GR_SUCCESS);
}

int
fmpz_mat_permanent(fmpz_t res, const fmpz_mat_t A)
{
    slong n = A->r;
    slong Abits, sbits, pbits, rbits, rn;

    if (n <= 2)
    {
        if (n == 0)
            fmpz_one(res);
        else if (n == 1)
            fmpz_set(res, fmpz_mat_entry(A, 0, 0));
        else if (n == 2)
            fmpz_fmma(res, fmpz_mat_entry(A, 0, 0), fmpz_mat_entry(A, 1, 1),
                           fmpz_mat_entry(A, 0, 1), fmpz_mat_entry(A, 1, 0));
        return 1;
    }

    if (n == 3 || n > FLINT_BITS - 2)
        return _fmpz_mat_permanent_generic(res, A);

    /* Todo: for non-uniform matrices, work with rowwise or columnwise
             bounds which may be tighter than matrixwise bounds. */
    Abits = fmpz_mat_max_bits(A);
    Abits = FLINT_ABS(Abits);

    /* Bound bits in sum of n rows or columns with alternating signs */
    sbits = Abits + FLINT_BIT_COUNT(n);

    /* A must have small coefficients; max(2,n)*A must fit in slong */
    if (Abits > SMALL_FMPZ_BITCOUNT_MAX || sbits > FLINT_BITS - 1)
        return _fmpz_mat_permanent_generic(res, A);

    /* Bound bits in product of n sums */
    pbits = n * sbits;
    /* Bound bits in sum of 2^(n-1) products, + 1 bit for signs */
    rbits = pbits + (n - 1) + 1;
    /* Number of limbs */
    rn = (rbits + FLINT_BITS - 1) / FLINT_BITS;

    /* Todo: cutoff could depend on Abits */
    if (n <= 10 || flint_get_num_available_threads() == 1)
        fmpz_mat_permanent_glynn1(res, A, sbits, rn);
    else
        fmpz_mat_permanent_glynn1_threaded(res, A, sbits, rn);

    return 1;
}

