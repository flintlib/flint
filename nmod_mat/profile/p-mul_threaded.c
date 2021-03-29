/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_mat.h"
#include "ulong_extras.h"
#include "thread_support.h"
#include "profiler.h"

#define MAC12(h, l, a, b)                   \
do {                                        \
    mp_limb_t pp0 = a*b;                    \
    add_ssaaaa(h, l, h, l, UWORD(0), pp0);  \
} while (0)

#define MAC22(h, l, a, b)               \
do {                                    \
    mp_limb_t pp1, pp0;                 \
    umul_ppmm(pp1, pp0, a, b);          \
    add_ssaaaa(h, l, h, l, pp1, pp0);   \
} while (0)

#define MAC23(h, m, l, a, b)                                \
do {                                                        \
    mp_limb_t pp1, pp0;                                     \
    umul_ppmm(pp1, pp0, a, b);                              \
    add_sssaaaaaa(h, m, l, h, m, l, UWORD(0), pp1, pp0);    \
} while (0)


mp_limb_t _my_dot_add(
    mp_limb_t r,
    const mp_limb_t * x,
    const mp_limb_t * y,
    slong len,
    nmod_t mod,
    int nlimbs)
{
    slong i;
    mp_limb_t a2, a1, a0, b2, b1, b0, c1, c0, d1, d0;

    a0 = r;
    a2 = a1 = 0;
    b2 = b1 = b0 = 0;
    c1 = c0 = 0;
    d1 = d0 = 0;

    if (nlimbs == 2)
    {
        if (mod.norm < FLINT_BITS/2)
        {
            for (i = 0; i + 4 <= len; i += 4)
            {
                MAC22(a1, a0, x[i+0], y[i+0]);
                MAC22(b1, b0, x[i+1], y[i+1]);
                MAC22(c1, c0, x[i+2], y[i+2]);
                MAC22(d1, d0, x[i+3], y[i+3]);
            }
        }
        else
        {
            for (i = 0; i + 4 <= len; i += 4)
            {
                MAC12(a1, a0, x[i+0], y[i+0]);
                MAC12(b1, b0, x[i+1], y[i+1]);
                MAC12(c1, c0, x[i+2], y[i+2]);
                MAC12(d1, d0, x[i+3], y[i+3]);
            }
        }

        add_ssaaaa(a1, a0, a1, a0, b1, b0);
        add_ssaaaa(c1, c0, c1, c0, d1, d0);
        add_ssaaaa(a1, a0, a1, a0, c1, c0);

        for ( ; i < len; i++)
            MAC22(a1, a0, x[i], y[i]);

        NMOD2_RED2(r, a1, a0, mod);
        return r;
    }
    else if (nlimbs == 3)
    {
        for (i = 0; i + 2 <= len; i += 2)
        {
            MAC23(a2, a1, a0, x[i+0], y[i+0]);
            MAC23(b2, b1, b0, x[i+1], y[i+1]);
        }

        add_sssaaaaaa(a2, a1, a0, a2, a1, a0, b2, b1, b0);

        for ( ; i < len; i++)
            MAC23(a2, a1, a0, x[i], y[i]);

        NMOD_RED3(r, a2, a1, a0, mod); 
        return r;
    }
    else
    {
        for (i = 0; i + 2 <= len; i += 2)
        {
            a0 += x[i+0]*y[i+0];
            b0 += x[i+1]*y[i+1];
        }

        a0 += b0;

        for ( ; i < len; i++)
            a0 += x[i] * y[i];

        NMOD_RED(r, a0, mod);
        return r;
    }
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
    mp_limb_t ** Crows;
    mp_limb_t ** Arows;
    mp_limb_t ** Brows;
    mp_limb_t * BL;
    nmod_t mod;
} _worker_arg;

static void _red_worker(void * varg)
{
    _worker_arg * arg = (_worker_arg *) varg;
    slong Bstartcol = arg->Bstartcol;
    slong Bstopcol = arg->Bstopcol;
    slong k = arg->k;
    slong n = arg->n;
    slong k_blk_sz = arg->k_blk_sz;
    mp_limb_t ** Brows = arg->Brows;
    mp_limb_t * BL = arg->BL;
    slong i, iq, ir, j;

    for (j = Bstartcol; j < Bstopcol; j++)
    {
        for (i = 0; i < k; i++)
        {
            iq = i/k_blk_sz;
            ir = i%k_blk_sz;
            BL[iq*n*k_blk_sz + j*k_blk_sz + ir] = Brows[i][j];
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
    slong m = arg->m;
    slong m_blk_sz = arg->m_blk_sz;
    slong k_blk_sz = arg->k_blk_sz;
    mp_limb_t ** Crows = arg->Crows;
    mp_limb_t ** Arows = arg->Arows;
    mp_limb_t * BL = arg->BL;
    nmod_t mod = arg->mod;
    mp_limb_t * TA, * TC;
    slong h, hh, i, ii, j;
    int nlimbs;
    TMP_INIT;

    /* +1 for feed-in */
    nlimbs = _nmod_vec_dot_bound_limbs(k_blk_sz + 1, mod);

    TMP_START;

    if (k <= k_blk_sz)
    {
        TA = TMP_ARRAY_ALLOC(2*m_blk_sz*k, mp_limb_t);

        for (h = Astartrow; h < Astoprow; h += m_blk_sz)
        {
            slong hstop = FLINT_MIN(m - h, m_blk_sz);

            for (hh = 0; hh < hstop; hh++)
            for (i = 0; i < k; i++)
                TA[hh*k + i] = Arows[h + hh][i];

            for (j = 0; j < n; j++)
            for (hh = 0; hh < hstop; hh++)
                Crows[h + hh][j] = _my_dot_add(0, &TA[hh*k],
                                        &BL[0*n + j*k_blk_sz], k, mod, nlimbs);        
        }

        TMP_END;
        return;
    }

    TA = TMP_ARRAY_ALLOC(m_blk_sz*k_blk_sz, mp_limb_t);
    TC = TMP_ARRAY_ALLOC(n*m_blk_sz, mp_limb_t);

    for (h = Astartrow; h < Astoprow; h += m_blk_sz)
    {
        slong hstop = FLINT_MIN(m - h, m_blk_sz);

        /* TC is a compressed block for the answer C[h:h+hhstop-1, all] */
        for (j = 0; j < n; j++)
        for (hh = 0; hh < hstop; hh++)
            TC[hh + m_blk_sz*j] = 0;

        for (i = 0; i < k; i += k_blk_sz)
        {
            slong istop = FLINT_MIN(k_blk_sz, k - i);

            /* get a compressed copy of A[h:h+hhstop, i:i+iistop] into TA */
            for (hh = 0; hh < hstop; hh++)
            for (ii = 0; ii < istop; ii++)
                TA[hh*k_blk_sz + ii] = Arows[h + hh][i + ii];

            /* addmul into answer block */
            for (j = 0; j < n; j++)
            for (hh = 0; hh < hstop; hh++)
                TC[hh + m_blk_sz*j] = _my_dot_add(TC[hh + m_blk_sz*j],
                                                  &TA[hh*k_blk_sz],
                                                  &BL[i*n + j*k_blk_sz],
                                                  istop, mod, nlimbs);
        }

        /* copy out answer for C[h:h+hhstop-1, all] */
        for (j = 0; j < n; j++)
        for (hh = 0; hh < hstop; hh++)
            Crows[h + hh][j] = TC[hh + m_blk_sz*j];
    }

    TMP_END;
    return;
}


void my_nmod_mat_mul(
    nmod_mat_t C,
    const nmod_mat_t A,
    const nmod_mat_t B)
{
    slong i;
    slong m = nmod_mat_nrows(A);
    slong k = nmod_mat_nrows(B);
    slong n = nmod_mat_ncols(B);
    _worker_arg mainarg;
    thread_pool_handle * handles;
    slong num_workers;
    _worker_arg * args;
    slong limit;
    slong k_blk_ct;
    slong k_blk_sz;
    slong m_blk_sz;
    TMP_INIT;

    if (m < 1 || k < 1 || n < 1)
    {
        nmod_mat_zero(C);
        return;
    }

    TMP_START;

    /* block size in the m direction */
    m_blk_sz = 16;

    /* choose a block size in the k direction */
    if (k <= 256)
    {
        k_blk_sz = k;
    }
    else
    {
        /*
            k_blk_sz = 256 is fine here, but we have some fun trying to
            maximize the relative size of the last block
        */
        slong l;
        double best_ratio = 0;

        k_blk_sz = 256-8*8;

        for (l = k_blk_sz; l < 256+4*8; l += 8)
        {
            double ratio = (double)(((k+l-1)%(l) + 1))/(double)(l);
            if (ratio > best_ratio)
            {
                k_blk_sz = l;
                best_ratio = ratio;
                if (best_ratio > 0.75)
                    break;
            }
        }
    }

    k_blk_ct = (k + k_blk_sz - 1)/k_blk_sz;

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
    mainarg.BL = TMP_ARRAY_ALLOC(n*k_blk_ct*k_blk_sz, mp_limb_t);
    mainarg.mod = A->mod;

    /* limit on number of threads */
    limit = FLINT_MAX(k, n);
    limit = FLINT_MIN(limit, m);
    limit = limit < 32 ? 0 : (limit - 32)/16;

    if (limit < 2)
    {
use_one_thread:

        _red_worker(&mainarg);
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
        args[i].mod = mainarg.mod;
    }

    i = num_workers;
    mainarg.Astartrow = (i + 0)*m/(num_workers + 1);
    mainarg.Astoprow = (i + 1)*m/(num_workers + 1);
    mainarg.Bstartcol = (i + 0)*n/(num_workers + 1);
    mainarg.Bstopcol = (i + 1)*n/(num_workers + 1);

    for (i = 0; i < num_workers; i++)
        thread_pool_wake(global_thread_pool, handles[i], 0, _red_worker, &args[i]);
    _red_worker(&mainarg);
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




void nmod_mat_mul_check(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
{
    slong i, j, k;
    mp_limb_t s0, s1, s2;
    mp_limb_t t0, t1;

    for (i = 0; i < A->r; i++)
    for (j = 0; j < B->c; j++)
    {
        s0 = s1 = s2 = UWORD(0);

        for (k = 0; k < A->c; k++)
        {
            umul_ppmm(t1, t0, A->rows[i][k], B->rows[k][j]);
            add_sssaaaaaa(s2, s1, s0, s2, s1, s0, 0, t1, t0);
        }

        NMOD_RED(s2, s2, C->mod);
        NMOD_RED3(s0, s2, s1, s0, C->mod);
        C->rows[i][j] = s0;
    }
}

int
main(void)
{
    slong m, k, n, i;
    FLINT_TEST_INIT(state);

    flint_printf("p-mul_threaded:\n");
    fflush(stdout);

    flint_set_num_threads(8);

    for (m = 4000; m >= 500; m = m*5/6)
    for (k = m; k >= m; k = k/2)
    for (n = m; n >= m; n = n/2)
    {
        mp_limb_t p;
        slong old_time, new_time;
        nmod_mat_t A, B, C, D;
        timeit_t timer;

        p = n_nextprime(UWORD(1)<<50, 1);

        nmod_mat_init(A, m, k, p);
        nmod_mat_init(B, k, n, p);
        nmod_mat_init(C, m, n, p);
        nmod_mat_init(D, m, n, p);

        nmod_mat_randfull(A, state);
        nmod_mat_randfull(B, state);
        nmod_mat_randfull(C, state);
        nmod_mat_randfull(D, state);

        flint_printf("------ dim %wd x %wd x %wd, %wd threads -------\n",
                     m, k, n, flint_get_num_threads());

        for (i = 0; i < 2; i++)
        {
            timeit_start(timer);
            my_nmod_mat_mul(D, A, B);
            timeit_stop(timer);
            new_time = timer->wall;


            timeit_start(timer);
            nmod_mat_mul_classical_threaded(C, A, B);
            timeit_stop(timer);
            old_time = timer->wall;

            flint_printf("old: %wd  new: %wd  ratio: %f\n",
                         old_time, new_time, (double)old_time/new_time);

            if (!nmod_mat_equal(C, D))
            {
                flint_printf("oops C != D\n");
                flint_abort();
            }
        }

        nmod_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(C);
        nmod_mat_clear(D);
    }

    return 0;
}

