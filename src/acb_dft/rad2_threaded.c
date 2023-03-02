/*
    Copyright (C) 2016 Pascal Molin
    Copyright (C) 2020 D.H.J. Polymath

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dft.h"
#include "pthread.h"

void acb_dft_rad2_reorder(acb_ptr v, slong n);

typedef struct
{
    acb_ptr v;
    acb_ptr vend;
    slong k;
    slong l;
    slong jstart;
    slong jend;
    acb_srcptr w;
    slong prec;
}
acb_dft_rad2_arg_t;

void *
_acb_dft_rad2_thread(void * arg_ptr)
{
    acb_dft_rad2_arg_t arg = *((acb_dft_rad2_arg_t *) arg_ptr);
    slong j, rstart, pstep;
    acb_ptr p, r;
    acb_t tmp;

    acb_init(tmp);
    rstart = arg.jstart / arg.l;
    pstep = 2 * arg.k;

    for (p = arg.v; p < arg.vend; p += pstep)
    {
        for (r = p + rstart, j = arg.jstart; j < arg.jend; j += arg.l, r++)
        {
            acb_mul(tmp, r + arg.k, arg.w + j, arg.prec);
            acb_sub(r + arg.k, r + 0, tmp, arg.prec);
            acb_add(r + 0, r + 0, tmp, arg.prec);
        }
    }

    acb_clear(tmp);
    flint_cleanup();
    return NULL;
}

void
acb_dft_rad2_precomp_inplace_threaded(acb_ptr v, const acb_dft_rad2_t rad2, slong prec)
{
    slong num_threads;
    pthread_t * threads;
    acb_dft_rad2_arg_t * args;

    slong t, logt, logk, logl;
    slong n = rad2->n;
    slong nz = rad2->nz; /* always n/2 ? */
    slong logn = rad2->e;

    num_threads = FLINT_MIN(flint_get_num_threads(), nz);
    for (logt = 0; WORD(1) << (logt + 1) <= num_threads; logt++);
    t = WORD(1) << logt;

    threads = flint_malloc(sizeof(pthread_t) * t);
    args = flint_malloc(sizeof(acb_dft_rad2_arg_t) * t);

    acb_dft_rad2_reorder(v, n);

    for (logk = 0, logl = logn - 1; logk < logn; logk += 1, logl -= 1)
    {
        slong i, j, p;
        slong logpstep = logk + 1 + FLINT_MAX(0, logl - logt);
        slong logjstep = logl + FLINT_MIN(logk, logn - 1 - logt);
        slong pstep = WORD(1) << logpstep;
        slong jstep = WORD(1) << logjstep;
        i = 0;
        for (p = 0; p < n; p += pstep)
        {
            for (j = 0; j < nz ; j += jstep)
            {
                args[i].v = v + p;
                args[i].vend = v + p + pstep;
                args[i].jstart = j;
                args[i].jend = j + jstep;
                args[i].k = WORD(1) << logk;
                args[i].l = WORD(1) << logl;
                args[i].w = rad2->z;
                args[i].prec = prec;
                pthread_create(&threads[i], NULL, _acb_dft_rad2_thread, &args[i]);
                i++;
            }
        }
        if (i != t)
        {
            flint_printf("threaded dft error: unequal i=%wd t=%wd\n", i, t);
            flint_abort();
        }
        for (i = 0; i < t; i++)
        {
            pthread_join(threads[i], NULL);
        }
    }

    flint_free(threads);
    flint_free(args);
}

void
acb_dft_rad2_inplace_threaded(acb_ptr v, int e, slong prec)
{
    acb_dft_rad2_t rad2;
    acb_dft_rad2_init(rad2, e, prec);
    acb_dft_rad2_precomp_inplace_threaded(v, rad2, prec);
    acb_dft_rad2_clear(rad2);
}
