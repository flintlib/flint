/*
    Copyright (C) 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_support.h"
#include "ulong_extras.h"
#include "nmod_poly.h"

typedef struct
{
   volatile slong * j;
   slong k;
   slong n;
   slong glen;
   slong ginvlen;
   mp_srcptr g;
   mp_srcptr ginv;
   mp_ptr * res;
   nmod_t mod;
#if FLINT_USES_PTHREAD
   pthread_mutex_t * mutex;
#endif
} powers_preinv_arg_t;

void
_nmod_poly_powers_mod_preinv_worker(void * arg_ptr)
{
    powers_preinv_arg_t arg = *((powers_preinv_arg_t *) arg_ptr);
    slong i, j, k = arg.k, n = arg.n;
    slong glen = arg.glen, ginvlen = arg.ginvlen;
    mp_ptr * res = arg.res;
    mp_srcptr g = arg.g, ginv = arg.ginv;
    const nmod_t mod = arg.mod;

    while (1)
    {
#if FLINT_USES_PTHREAD
        pthread_mutex_lock(arg.mutex);
#endif
	j = *arg.j + k;
        *arg.j = j;
#if FLINT_USES_PTHREAD
        pthread_mutex_unlock(arg.mutex);
#endif

        if (j >= n)
            return;

        if (glen == 2) /* special case, constant polynomials */
        {
            for (i = j + 1; i < j + k && i < n; i++)
                res[i][0] = n_mulmod2_preinv(res[j][0], res[i - j][0],
				                              mod.n, mod.ninv);
        } else
        {
            for (i = j + 1; i < j + k && i < n; i++)
                _nmod_poly_mulmod_preinv(res[i], res[j],
                  glen - 1, res[i - j], glen - 1, g, glen, ginv, ginvlen, mod);
        }
    }
}

/*
    compute f^0, f^1, ..., f^(n-1) mod g, where g has length glen and f is
    reduced mod g and has length flen (possibly zero spaced)
    assumes res is an array of n arrays each with space for at least glen - 1
    coefficients and that flen > 0
    {ginv, ginvlen} must be set to the power series inverse of the reverse of g
*/
void
_nmod_poly_powers_mod_preinv_threaded_pool(mp_ptr * res, mp_srcptr f,
		 slong flen, slong n, mp_srcptr g, slong glen,
           mp_srcptr ginv, slong ginvlen, const nmod_t mod,
	                       thread_pool_handle * threads, slong num_threads)
{
    slong i, k, shared_j = 0;
    powers_preinv_arg_t * args;
#if FLINT_USES_PTHREAD
    pthread_mutex_t mutex;
#endif

    if (n == 0)
        return;

    if (n == 1)
    {
        if (glen > 1)
            res[0][0] = 1;

        if (glen > 2)
            flint_mpn_zero(res[0] + 1, glen - 2);

        return;
    }

    k = n_sqrt(n);

    /* compute baby steps */

    _nmod_poly_powers_mod_preinv_naive(res, f, flen, k + 1,
		                                  g, glen, ginv, ginvlen, mod);

    /* compute giant steps */

    /* f^(k*i) = f^(k*(i - 1))*f^k */
    if (glen == 2) /* special case, constant polys */
    {
        for (i = 2*k; i < n; i += k)
            res[i][0] = n_mulmod2_preinv(res[i - k][0], res[k][0],
                                                              mod.n, mod.ninv);
    } else
    {
        for (i = 2*k; i < n; i += k)
            _nmod_poly_mulmod_preinv(res[i], res[i - k], glen - 1,
			        res[k], glen - 1, g, glen, ginv, ginvlen, mod);
    }

    args = (powers_preinv_arg_t *)
                 flint_malloc(sizeof(powers_preinv_arg_t) * (num_threads + 1));

    for (i = 0; i < num_threads + 1; i++)
    {
        args[i].j       = &shared_j;
        args[i].k       = k;
        args[i].n       = n;
        args[i].glen    = glen;
        args[i].ginvlen = ginvlen;
        args[i].g       = g;
        args[i].ginv    = ginv;
        args[i].res     = res;
	args[i].mod     = mod;
#if FLINT_USES_PTHREAD
        args[i].mutex   = &mutex;
#endif
    }

#if FLINT_USES_PTHREAD
    pthread_mutex_init(&mutex, NULL);
#endif

    for (i = 0; i < num_threads; i++)
        thread_pool_wake(global_thread_pool, threads[i], 0,
                _nmod_poly_powers_mod_preinv_worker, &args[i]);

    _nmod_poly_powers_mod_preinv_worker(&args[num_threads]);

    for (i = 0; i < num_threads; i++)
        thread_pool_wait(global_thread_pool, threads[i]);

#if FLINT_USES_PTHREAD
    pthread_mutex_destroy(&mutex);
#endif

    flint_free(args);
}

void
_nmod_poly_powers_mod_preinv_threaded(mp_ptr * res, mp_srcptr f,
                 slong flen, slong n, mp_srcptr g, slong glen,
            mp_srcptr ginv, slong ginvlen, const nmod_t mod)
{
    thread_pool_handle * threads;
    slong num_threads = flint_request_threads(&threads, flint_get_num_threads());

   _nmod_poly_powers_mod_preinv_threaded_pool(res, f, flen, n,
                             g, glen, ginv, ginvlen, mod, threads, num_threads);

    flint_give_back_threads(threads, num_threads);
}

void
nmod_poly_powers_mod_bsgs(nmod_poly_struct * res,
		             const nmod_poly_t f, slong n, const nmod_poly_t g)
{
    slong i;

    nmod_poly_t ginv;
    mp_ptr * res_arr;

    if (nmod_poly_length(g) == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (nmod_poly_powers_mod_naive). Divide by zero.\n");
    }

    if (nmod_poly_length(f) == 0 || nmod_poly_length(g) == 1)
    {
        if (n > 0)
           nmod_poly_one(res + 0);

        for (i = 1; i < n; i++)
           nmod_poly_zero(res + i);

        return;
    }

    if (nmod_poly_length(f) >= nmod_poly_length(g))
    {
        nmod_poly_t q, r;

        nmod_poly_init_mod(q, f->mod);
        nmod_poly_init_mod(r, f->mod);

        nmod_poly_divrem(q, r, f, g);
        nmod_poly_powers_mod_bsgs(res, r, n, g);

        nmod_poly_clear(q);
        nmod_poly_clear(r);

        return;
    }

    res_arr = (mp_ptr *) flint_malloc(n*sizeof(mp_ptr));
    nmod_poly_init_mod(ginv, g->mod);

    for (i = 0; i < n; i++)
    {
       nmod_poly_fit_length(res + i, nmod_poly_length(g) - 1);
       res_arr[i] = res[i].coeffs;
       _nmod_poly_set_length(res + i, nmod_poly_length(g) - 1);
    }

    nmod_poly_reverse(ginv, g, nmod_poly_length(g));
    nmod_poly_inv_series(ginv, ginv, nmod_poly_length(g));

    _nmod_poly_powers_mod_preinv_threaded(res_arr, f->coeffs, f->length, n,
                     g->coeffs, g->length, ginv->coeffs, ginv->length, g->mod);

    for (i = 0; i < n; i++)
       _nmod_poly_normalise(res + i);

    nmod_poly_clear(ginv);
    flint_free(res_arr);
}
