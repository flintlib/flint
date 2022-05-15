/*
    Copyright (C) 2020 Daniel Schultz
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef THREAD_SUPPORT_H
#define THREAD_SUPPORT_H

#include "flint.h"
#include "thread_pool.h"

#ifdef __cplusplus
 extern "C" {
#endif

#define FLINT_DEFAULT_THREAD_LIMIT 99999

FLINT_DLL slong flint_request_threads(thread_pool_handle ** handles,
                                                           slong thread_limit);

FLINT_DLL void flint_give_back_threads(thread_pool_handle * handles,
                                                            slong num_handles);

#define FLINT_PARALLEL_UNIFORM 1
#define FLINT_PARALLEL_STRIDED 2
#define FLINT_PARALLEL_DYNAMIC 4
#define FLINT_PARALLEL_BSPLIT_LEFT_INPLACE 8
#define FLINT_PARALLEL_VERBOSE 512

typedef void (* do_func_t)(slong i, void * args);

FLINT_DLL void flint_parallel_do(do_func_t f, void * args, slong n, int thread_limit, int flags);

typedef void (* bsplit_merge_func_t)(void *, void *, void *, void *);
typedef void (* bsplit_basecase_func_t)(void *, slong, slong, void *);
typedef void (* bsplit_init_func_t)(void *, void *);
typedef void (* bsplit_clear_func_t)(void *, void *);

FLINT_DLL void flint_parallel_binary_splitting(void * res, bsplit_basecase_func_t basecase, bsplit_merge_func_t merge,
    size_t sizeof_res, bsplit_init_func_t init, bsplit_clear_func_t clear, void * args, slong a, slong b, slong basecase_cutoff, int thread_limit, int flags);

#ifdef __cplusplus
}
#endif

#endif
