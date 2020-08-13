/*
    Copyright (C) 2020 Daniel Schultz

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

#ifdef __cplusplus
}
#endif

#endif
