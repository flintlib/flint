/*
    Copyright (C) 2010,2011,2018 Fredrik Johansson
    Copyright (C) 2016 Aaditya Thakkar
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

static void _dot1(fmpz_t z, fmpz * a, slong * b, slong len)
{
    slong i;
    slong s0 = 0;

    for (i = 0; i < len; i++)
        s0 += a[i]*b[i];

    fmpz_set_si(z, s0);
}

static void _dot2(fmpz_t z, fmpz * a, slong * b, slong len)
{
    slong i;
    ulong p1, p0, s1, s0;

    s1 = s0 = 0;
    for (i = 0; i < len; i++)
    {
        smul_ppmm(p1, p0, a[i], b[i]);
        add_ssaaaa(s1, s0, s1, s0, p1, p0);
    }

    fmpz_set_signed_uiui(z, s1, s0);
}

static void _dot3(fmpz_t z, fmpz * a, slong * b, slong len)
{
    slong i;
    ulong p1, p0, s2, s1, s0;

    s2 = s1 = s0 = 0;
    for (i = 0; i < len; i++)
    {
        smul_ppmm(p1, p0, a[i], b[i]);
        add_sssaaaaaa(s2, s1, s0, s2, s1, s0, FLINT_SIGN_EXT(p1), p1, p0);
    }

    fmpz_set_signed_uiuiui(z, s2, s1, s0);
}

static void _dot_add1(ulong * s, fmpz * a, slong * b, slong len)
{
    slong i;
    slong s0 = s[0];

    for (i = 0; i < len; i++)
        s0 += a[i]*b[i];

    s[0] = s0;
}

static void _dot_add2(ulong * s, fmpz * a, slong * b, slong len)
{
    slong i;
    ulong p1, p0, s1, s0, t1, t0;

    FLINT_ASSERT(len > 0);

    s0 = s[0];
    s1 = s[1];

    if (len & 1)
    {
        smul_ppmm(t1, t0, a[0], b[0]);
        a++;
        b++;
    }
    else
    {
        t1 = t0 = 0;
    }

    for (i = 0; i < len/2; i++)
    {
        smul_ppmm(p1, p0, a[2*i+0], b[2*i+0]);
        add_ssaaaa(s1, s0, s1, s0, p1, p0);
        smul_ppmm(p1, p0, a[2*i+1], b[2*i+1]);
        add_ssaaaa(t1, t0, t1, t0, p1, p0);
    }

    add_ssaaaa(s1, s0, s1, s0, t1, t0);

    s[0] = s0;
    s[1] = s1;
}

static void _dot_add3(ulong * s, fmpz * a, slong * b, slong len)
{
    slong i;
    ulong p1, p0, s2, s1, s0, t2, t1, t0;

    FLINT_ASSERT(len > 0);

    s0 = s[0];
    s1 = s[1];
    s2 = s[2];

    if (len & 1)
    {
        smul_ppmm(t1, t0, a[0], b[0]);
        t2 = FLINT_SIGN_EXT(t1);
        a++;
        b++;
    }
    else
    {
        t2 = t1 = t0 = 0;
    }

    for (i = 0; i < len/2; i++)
    {
        smul_ppmm(p1, p0, a[2*i+0], b[2*i+0]);
        add_sssaaaaaa(s2, s1, s0, s2, s1, s0, FLINT_SIGN_EXT(p1), p1, p0);
        smul_ppmm(p1, p0, a[2*i+1], b[2*i+1]);
        add_sssaaaaaa(t2, t1, t0, t2, t1, t0, FLINT_SIGN_EXT(p1), p1, p0);
    }

    add_sssaaaaaa(s2, s1, s0, s2, s1, s0, t2, t1, t0);

    s[0] = s0;
    s[1] = s1;
    s[2] = s2;
}


typedef struct {
    slong Astartrow;
    slong Astoprow;
    slong Bstartcol;
    slong Bstopcol;
    slong n;
    slong k;
    slong m;
    slong k_blk_sz;
    slong m_blk_sz;
    fmpz ** Crows;
    fmpz ** Arows;
    fmpz ** Brows;
    slong * BL;
    int words;
} _worker_arg;

static void _tr_worker(void * varg)
{
    _worker_arg * arg = (_worker_arg *) varg;
    slong Bstartcol = arg->Bstartcol;
    slong Bstopcol = arg->Bstopcol;
    slong k = arg->k;
    slong n = arg->n;
    slong k_blk_sz = arg->k_blk_sz;
    fmpz ** Brows = arg->Brows;
    slong * BL = arg->BL;
    slong i, iq, ir, j;

    iq = ir = 0;
    for (i = 0; i < k; i++)
    {
        for (j = Bstartcol; j < Bstopcol; j++)
            BL[iq*n*k_blk_sz + j*k_blk_sz + ir] = Brows[i][j];

        /* (iq, ir) = divrem(i, k_blk_sz) */
        ir++;
        if (ir >= k_blk_sz)
        {
            iq++;
            ir = 0;
        }
    }
}

static void _mul_worker(void * varg)
{
    _worker_arg * arg = (_worker_arg *) varg;
    slong Astartrow = arg->Astartrow;
    slong Astoprow = arg->Astoprow;
    slong n = arg->n;
    slong k = arg->k;
    slong m_blk_sz = arg->m_blk_sz;
    slong k_blk_sz = arg->k_blk_sz;
    fmpz ** Crows = arg->Crows;
    fmpz ** Arows = arg->Arows;
    slong * BL = arg->BL;
    slong * TA;
    ulong * TC;
    slong h, hh, i, ii, j;
    int words = arg->words;
    TMP_INIT;

    if (k <= k_blk_sz)
    {
        /* no blocking overhead: the B matrix is fully transposed in BL */
        if (words == 1)
        {
            for (h = Astartrow; h < Astoprow; h++)
                for (j = 0; j < n; j++)
                    _dot1(&Crows[h][j], Arows[h], &BL[j*k_blk_sz], k);
        }
        else if (words == 2)
        {
            for (h = Astartrow; h < Astoprow; h++)
                for (j = 0; j < n; j++)
                    _dot2(&Crows[h][j], Arows[h], &BL[j*k_blk_sz], k);
        }
        else
        {
            for (h = Astartrow; h < Astoprow; h++)
                for (j = 0; j < n; j++)
                    _dot3(&Crows[h][j], Arows[h], &BL[j*k_blk_sz], k);
        }

        return;
    }

    TMP_START;

    TA = TMP_ARRAY_ALLOC(m_blk_sz*k_blk_sz, slong);
    TC = TMP_ARRAY_ALLOC(n*m_blk_sz*words, ulong);

    for (h = Astartrow; h < Astoprow; h += m_blk_sz)
    {
        slong hstop = FLINT_MIN(Astoprow - h, m_blk_sz);

        /* TC is a compressed block for the answer C[h:h+hhstop-1, all] */
        for (j = 0; j < n*hstop*words; j++)
            TC[j] = 0;

        for (i = 0; i < k; i += k_blk_sz)
        {
            slong istop = FLINT_MIN(k_blk_sz, k - i);

            /* get a compressed copy of A[h:h+hhstop, i:i+iistop] into TA */
            for (hh = 0; hh < hstop; hh++)
            for (ii = 0; ii < istop; ii++)
                TA[hh*k_blk_sz + ii] = Arows[h + hh][i + ii];

            /* addmul into answer block */
            if (words == 1)
            {
                for (j = 0; j < n; j++)
                for (hh = 0; hh < hstop; hh++)
                    _dot_add1(&TC[1*(hh + hstop*j)], &TA[hh*k_blk_sz],
                                                 &BL[i*n + j*k_blk_sz], istop);
            }
            else if (words == 2)
            {
                for (j = 0; j < n; j++)
                for (hh = 0; hh < hstop; hh++)
                    _dot_add2(&TC[2*(hh + hstop*j)], &TA[hh*k_blk_sz],
                                                 &BL[i*n + j*k_blk_sz], istop);
            }
            else
            {
                for (j = 0; j < n; j++)
                for (hh = 0; hh < hstop; hh++)
                    _dot_add3(&TC[3*(hh + hstop*j)], &TA[hh*k_blk_sz],
                                                 &BL[i*n + j*k_blk_sz], istop);
            }
        }

        /* copy out answer for C[h:h+hhstop-1, all] */
        if (words == 1)
        {
            for (j = 0; j < n; j++)
            for (hh = 0; hh < hstop; hh++)
                fmpz_set_si(&Crows[h + hh][j], (slong)TC[hh + hstop*j]);
        }
        else if (words == 2)
        {
            for (j = 0; j < n; j++)
            for (hh = 0; hh < hstop; hh++)
                fmpz_set_signed_uiui(&Crows[h + hh][j],
                                     TC[2*(hh + hstop*j) + 1],
                                     TC[2*(hh + hstop*j) + 0]);
        }
        else
        {
            for (j = 0; j < n; j++)
            for (hh = 0; hh < hstop; hh++)
                fmpz_set_signed_uiuiui(&Crows[h + hh][j],
                                       TC[3*(hh + hstop*j) + 2],
                                       TC[3*(hh + hstop*j) + 1],
                                       TC[3*(hh + hstop*j) + 0]);
        }
    }

    TMP_END;
    return;
}


void _fmpz_mat_mul_small_internal(
    fmpz_mat_t C,
    const fmpz_mat_t A,
    const fmpz_mat_t B,
    flint_bitcnt_t Cbits)
{
    slong i;
    slong m = fmpz_mat_nrows(A);
    slong k = fmpz_mat_nrows(B);
    slong n = fmpz_mat_ncols(B);
    _worker_arg mainarg;
    thread_pool_handle * handles;
    slong num_workers;
    _worker_arg * args;
    slong limit;
    slong k_blk_ct;
    slong k_blk_sz;
    slong m_blk_sz;
    TMP_INIT;

    FLINT_ASSERT(m > 0 && k > 0 && n > 0);

    TMP_START;

    /* _dot1 moves through the data quickly, so it is worth blocking */

    /* block size in the m direction */
    m_blk_sz = 16;

    /* choose a block size in the k direction */
    if (k <= 128)
    {
        k_blk_sz = k;
        k_blk_ct = 1;
    }
    else
    {
        k_blk_sz = 128;
        k_blk_ct = (k + k_blk_sz - 1)/k_blk_sz;
    }

    mainarg.m_blk_sz = m_blk_sz;
    mainarg.k_blk_sz = k_blk_sz;
    mainarg.Astartrow = 0;
    mainarg.Astoprow = m;
    mainarg.Bstartcol = 0;
    mainarg.Bstopcol = n;

    mainarg.k = k;
    mainarg.m = m;
    mainarg.n = n;
    mainarg.Crows = C->rows;
    mainarg.Arows = A->rows;
    mainarg.Brows = B->rows;
    mainarg.BL = TMP_ARRAY_ALLOC(n*k_blk_ct*k_blk_sz, slong);

    if (Cbits <= SMALL_FMPZ_BITCOUNT_MAX)
        mainarg.words = 1;
    else if (Cbits <= 2*FLINT_BITS - 1)
        mainarg.words = 2;
    else
        mainarg.words = 3;

    /* limit on number of threads */
    limit = FLINT_MAX(k, n);
    limit = FLINT_MIN(limit, m);
    limit = limit < 32 ? 0 : (limit - 32)/16;

    if (limit < 2)
    {
use_one_thread:

        _tr_worker(&mainarg);
        _mul_worker(&mainarg);

        TMP_END;
        return;
    }

    num_workers = flint_request_threads(&handles, limit);
    if (num_workers < 1)
    {
        flint_give_back_threads(handles, num_workers);
        goto use_one_thread;
    }

    args = FLINT_ARRAY_ALLOC(num_workers, _worker_arg);

    for (i = 0; i < num_workers; i++)
    {
        args[i].m_blk_sz = mainarg.m_blk_sz;
        args[i].k_blk_sz = mainarg.k_blk_sz;
        args[i].Astartrow = (i + 0)*m/(num_workers + 1);
        args[i].Astoprow = (i + 1)*m/(num_workers + 1);
        args[i].Bstartcol = (i + 0)*n/(num_workers + 1);
        args[i].Bstopcol = (i + 1)*n/(num_workers + 1);
        args[i].k = mainarg.k;
        args[i].m = mainarg.m;
        args[i].n = mainarg.n;
        args[i].Crows = mainarg.Crows;
        args[i].Arows = mainarg.Arows;
        args[i].Brows = mainarg.Brows;
        args[i].BL = mainarg.BL;
        args[i].words = mainarg.words;
    }

    i = num_workers;
    mainarg.Astartrow = (i + 0)*m/(num_workers + 1);
    mainarg.Astoprow = (i + 1)*m/(num_workers + 1);
    mainarg.Bstartcol = (i + 0)*n/(num_workers + 1);
    mainarg.Bstopcol = (i + 1)*n/(num_workers + 1);

    for (i = 0; i < num_workers; i++)
        thread_pool_wake(global_thread_pool, handles[i], 0, _tr_worker, &args[i]);
    _tr_worker(&mainarg);
    for (i = 0; i < num_workers; i++)
        thread_pool_wait(global_thread_pool, handles[i]);

    for (i = 0; i < num_workers; i++)
        thread_pool_wake(global_thread_pool, handles[i], 0, _mul_worker, &args[i]);
    _mul_worker(&mainarg);
    for (i = 0; i < num_workers; i++)
        thread_pool_wait(global_thread_pool, handles[i]);

    flint_give_back_threads(handles, num_workers);
    flint_free(args);

    TMP_END;
}

void _fmpz_mat_mul_small(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
{
    slong Abits, Bbits;
    flint_bitcnt_t Cbits;

    if (fmpz_mat_is_empty(A) || fmpz_mat_is_empty(B))
    {
        fmpz_mat_zero(C);
        return;
    }

    Abits = fmpz_mat_max_bits(A);
    Bbits = fmpz_mat_max_bits(B);
    Abits = FLINT_ABS(Abits);
    Bbits = FLINT_ABS(Bbits);

    Cbits = Abits + Bbits + FLINT_BIT_COUNT(A->c);

    _fmpz_mat_mul_small_internal(C, A, B, Cbits);
}

