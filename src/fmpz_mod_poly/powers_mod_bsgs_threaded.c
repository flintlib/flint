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
#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"

typedef struct
{
   volatile slong * j;
   slong k;
   slong n;
   slong glen;
   slong ginvlen;
   const fmpz * g;
   const fmpz * ginv;
   fmpz ** res;
   const fmpz_mod_ctx_struct * ctx;
#if FLINT_USES_PTHREAD
   pthread_mutex_t * mutex;
#endif
} fmpz_powers_preinv_arg_t;

void
_fmpz_mod_poly_powers_mod_preinv_worker(void * arg_ptr)
{
    fmpz_powers_preinv_arg_t arg = *((fmpz_powers_preinv_arg_t *) arg_ptr);
    slong i, j, k = arg.k, n = arg.n;
    slong glen = arg.glen, ginvlen = arg.ginvlen;
    fmpz ** res = arg.res;
    const fmpz * g = arg.g, * ginv = arg.ginv;
    const fmpz_mod_ctx_struct * ctx = arg.ctx;

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
            {
                fmpz_mod_mul(res[i] + 0, res[j] + 0, res[i - j] + 0, ctx);
            }
        } else
        {
            for (i = j + 1; i < j + k && i < n; i++)
                _fmpz_mod_poly_mulmod_preinv(res[i], res[j],
                  glen - 1, res[i - j], glen - 1, g, glen, ginv, ginvlen, ctx);
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
_fmpz_mod_poly_powers_mod_preinv_threaded_pool(fmpz ** res, const fmpz * f,
		 slong flen, slong n, const fmpz * g, slong glen,
           const fmpz * ginv, slong ginvlen, const fmpz_mod_ctx_t ctx,
	                       thread_pool_handle * threads, slong num_threads)
{
    slong i, k, shared_j = 0;
    fmpz_powers_preinv_arg_t * args;
#if FLINT_USES_PTHREAD
    pthread_mutex_t mutex;
#endif

    if (n == 0)
        return;

    if (n == 1)
    {
        if (glen > 1)
            fmpz_set_ui(res[0] + 0, 1);

        if (glen > 2)
        {

            for (i = 1; i < glen - 1; i++)
                 fmpz_zero(res[0] + i);
        }

        return;
    }

    k = n_sqrt(n);

    /* compute baby steps */

    _fmpz_mod_poly_powers_mod_preinv_naive(res, f, flen, k + 1,
		                                  g, glen, ginv, ginvlen, ctx);

    /* compute giant steps */

    /* f^(k*i) = f^(k*(i - 1))*f^k */
    if (glen == 2) /* special case, constant polys */
    {
        for (i = 2*k; i < n; i += k)
        {
            fmpz_mod_mul(res[i] + 0, res[i - k] + 0, res[k] + 0, ctx);
        }
    }
    else
    {
        for (i = 2*k; i < n; i += k)
            _fmpz_mod_poly_mulmod_preinv(res[i], res[i - k], glen - 1,
			        res[k], glen - 1, g, glen, ginv, ginvlen, ctx);
    }

    args = (fmpz_powers_preinv_arg_t *)
                 flint_malloc(sizeof(fmpz_powers_preinv_arg_t) * (num_threads + 1));

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
        args[i].ctx     = ctx;
#if FLINT_USES_PTHREAD
        args[i].mutex   = &mutex;
#endif
    }

#if FLINT_USES_PTHREAD
    pthread_mutex_init(&mutex, NULL);
#endif

    for (i = 0; i < num_threads; i++)
        thread_pool_wake(global_thread_pool, threads[i], 0,
                _fmpz_mod_poly_powers_mod_preinv_worker, &args[i]);

    _fmpz_mod_poly_powers_mod_preinv_worker(&args[num_threads]);

    for (i = 0; i < num_threads; i++)
        thread_pool_wait(global_thread_pool, threads[i]);

#if FLINT_USES_PTHREAD
    pthread_mutex_destroy(&mutex);
#endif

    flint_free(args);
}

void
fmpz_mod_poly_powers_mod_bsgs(fmpz_mod_poly_struct * res,
    	     const fmpz_mod_poly_t f, slong n, const fmpz_mod_poly_t g,
                                                      const fmpz_mod_ctx_t ctx)
{
    slong i;

    fmpz_mod_poly_t ginv;
    fmpz ** res_arr;
    thread_pool_handle * threads;
    slong num_threads;

    if (fmpz_mod_poly_length(g, ctx) == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_mod_poly_powers_mod_bsgs). Divide by zero.\n");
    }

    if (fmpz_mod_poly_length(f, ctx) == 0 || fmpz_mod_poly_length(g, ctx) == 1)
    {
        if (n > 0)
           fmpz_mod_poly_one(res + 0, ctx);

        for (i = 1; i < n; i++)
           fmpz_mod_poly_zero(res + i, ctx);

        return;
    }

    if (fmpz_mod_poly_length(f, ctx) >= fmpz_mod_poly_length(g, ctx))
    {
        fmpz_mod_poly_t q, r;

        fmpz_mod_poly_init(q, ctx);
        fmpz_mod_poly_init(r, ctx);

        fmpz_mod_poly_divrem(q, r, f, g, ctx);
        fmpz_mod_poly_powers_mod_naive(res, r, n, g, ctx);

        fmpz_mod_poly_clear(q, ctx);
        fmpz_mod_poly_clear(r, ctx);

        return;
    }

    res_arr = (fmpz **) flint_malloc(n*sizeof(fmpz *));
    fmpz_mod_poly_init(ginv, ctx);

    for (i = 0; i < n; i++)
    {
       fmpz_mod_poly_fit_length(res + i, fmpz_mod_poly_length(g, ctx) - 1, ctx);
       res_arr[i] = res[i].coeffs;
       _fmpz_mod_poly_set_length(res + i, fmpz_mod_poly_length(g, ctx) - 1);
    }

    fmpz_mod_poly_reverse(ginv, g, fmpz_mod_poly_length(g, ctx), ctx);
    fmpz_mod_poly_inv_series(ginv, ginv, fmpz_mod_poly_length(g, ctx), ctx);

    num_threads = flint_request_threads(&threads, flint_get_num_threads());

    _fmpz_mod_poly_powers_mod_preinv_threaded_pool(res_arr, f->coeffs,
    	       f->length, n, g->coeffs, g->length, ginv->coeffs, ginv->length,
                              ctx, threads, num_threads);

    flint_give_back_threads(threads, num_threads);

    for (i = 0; i < n; i++)
       _fmpz_mod_poly_normalise(res + i);

    fmpz_mod_poly_clear(ginv, ctx);
    flint_free(res_arr);
}
