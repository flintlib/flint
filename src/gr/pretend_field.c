/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr.h"

#if FLINT_USES_PTHREAD
#include <pthread.h>

static pthread_mutex_t _gr_zero_divisor_mutex = PTHREAD_MUTEX_INITIALIZER;

void _gr_ctx_zero_divisor_lock(void) { pthread_mutex_lock(&_gr_zero_divisor_mutex); }
void _gr_ctx_zero_divisor_unlock(void) { pthread_mutex_unlock(&_gr_zero_divisor_mutex); }

#else

void _gr_ctx_zero_divisor_lock(void) { }
void _gr_ctx_zero_divisor_unlock(void) { }

#endif
