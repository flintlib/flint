/*
    Copyright (C) 2020 Rudolph
    Copyright (C) 2020 D.H.J. Polymath

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"
#include "pthread.h"

slong platt_get_smk_index(slong B, const fmpz_t j, slong prec);
void get_smk_points(slong * res, slong A, slong B);

void _platt_smk(acb_ptr table, acb_ptr startvec, acb_ptr stopvec,
        const slong * smk_points, const arb_t t0, slong A, slong B,
        const fmpz_t jstart, const fmpz_t jstop, slong mstart, slong mstop,
        slong K, slong prec);

void _acb_dirichlet_platt_multieval(arb_ptr out, acb_srcptr S_table,
        const arb_t t0, slong A, slong B, const arb_t h, const fmpz_t J,
        slong K, slong sigma, slong prec);

typedef struct
{
    acb_ptr S;
    acb_ptr startvec;
    acb_ptr stopvec;
    const slong * smk_points;
    arb_srcptr t0;
    slong A;
    slong B;
    slong K;
    fmpz_t jstart;
    fmpz_t jstop;
    slong mstart;
    slong mstop;
    slong prec;
}
platt_smk_arg_t;

static void *
_platt_smk_thread(void * arg_ptr)
{
    platt_smk_arg_t *p = (platt_smk_arg_t *) arg_ptr;
    _platt_smk(p->S, p->startvec, p->stopvec, p->smk_points, p->t0, p->A, p->B,
               p->jstart, p->jstop, p->mstart, p->mstop, p->K, p->prec);
    flint_cleanup();
    return NULL;
}


void
acb_dirichlet_platt_multieval_threaded(arb_ptr out, const fmpz_t T, slong A,
        slong B, const arb_t h, const fmpz_t J, slong K,
        slong sigma, slong prec)
{
    slong i, num_threads, N;
    fmpz * smk_points;
    pthread_t * threads;
    platt_smk_arg_t * args;
    acb_ptr S;
    arb_t t0;
    fmpz_t threadtasks;

    num_threads = flint_get_num_threads();
    if (num_threads < 1)
    {
        flint_printf("no threads available\n");
        flint_abort();
    }

    N = A*B;
    fmpz_init(threadtasks);
    threads = flint_malloc(sizeof(pthread_t) * num_threads);
    args = flint_malloc(sizeof(platt_smk_arg_t) * num_threads);
    fmpz_add_si(threadtasks, J, num_threads - 1);
    fmpz_tdiv_q_ui(threadtasks, threadtasks, (ulong) num_threads);
    smk_points = _fmpz_vec_init(N);
    arb_init(t0);

    get_smk_points(smk_points, A, B);
    arb_set_fmpz(t0, T);

    S =  _acb_vec_init(K*N);
    for (i = 0; i < num_threads; i++)
    {
        args[i].S = S;
        args[i].startvec = _acb_vec_init(K);
        args[i].stopvec = _acb_vec_init(K);
        args[i].smk_points = smk_points;
        args[i].t0 = t0;
        args[i].A = A;
        args[i].B = B;
        args[i].K = K;
        args[i].prec = prec;
        fmpz_init(args[i].jstart);
        fmpz_init(args[i].jstop);
        fmpz_mul_si(args[i].jstart, threadtasks, i);
        fmpz_add_ui(args[i].jstart, args[i].jstart, 1);
        fmpz_mul_si(args[i].jstop, threadtasks, i + 1);
        args[i].mstart = platt_get_smk_index(B, args[i].jstart, prec);
        args[i].mstop = platt_get_smk_index(B, args[i].jstop, prec);
    }
    fmpz_set(args[num_threads-1].jstop, J);
    args[num_threads-1].mstop = platt_get_smk_index(B, J, prec);

    for (i = 0; i < num_threads; i++)
    {
        pthread_create(&threads[i], NULL, _platt_smk_thread, &args[i]);
    }

    for (i = 0; i < num_threads; i++)
    {
        pthread_join(threads[i], NULL);
    }

    for (i = 0; i < num_threads; i++)
    {
        slong k;
        for (k = 0; k < K; k++)
        {
            acb_ptr z;
            z = S + N*k + args[i].mstart;
            acb_add(z, z, args[i].startvec + k, prec);
            z = S + N*k + args[i].mstop;
            acb_add(z, z, args[i].stopvec + k, prec);
        }
        _acb_vec_clear(args[i].startvec, K);
        _acb_vec_clear(args[i].stopvec, K);
        fmpz_clear(args[i].jstart);
        fmpz_clear(args[i].jstop);
    }

    _acb_dirichlet_platt_multieval(out, S, t0, A, B, h, J, K, sigma, prec);

    arb_clear(t0);
    _acb_vec_clear(S, K*N);
    _fmpz_vec_clear(smk_points, N);

    flint_free(args);
    flint_free(threads);
}
