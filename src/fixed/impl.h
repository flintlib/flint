/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Library-internal view of the cached reduction tables.

   These thread-local objects are shared across the translation units
   of the module (the dispatch files, the reduction, and the
   specialized per-size implementations) but are deliberately NOT
   declared in fixed.h: Windows DLLs cannot export thread-local data,
   so external consumers -- the test suite -- go through the
   _fixed_exp_logs_entry / _fixed_atans_entry accessors instead.

   Each entry occupies _fixed_{exp_logs,atans}_n limbs: the value
   limbs with one guard limb below them.  Consumers wanting the top n
   limbs of entry i read tab + i * stride + (stride - n). */

#ifndef FIXED_IMPL_H
#define FIXED_IMPL_H

#include "flint.h"

#ifdef __cplusplus
extern "C" {
#endif

extern FLINT_TLS_PREFIX nn_ptr _fixed_exp_logs;
extern FLINT_TLS_PREFIX slong _fixed_exp_logs_n;
extern FLINT_TLS_PREFIX slong _fixed_exp_logs_r;

extern FLINT_TLS_PREFIX nn_ptr _fixed_atans;
extern FLINT_TLS_PREFIX slong _fixed_atans_n;
extern FLINT_TLS_PREFIX slong _fixed_atans_r;

#ifdef __cplusplus
}
#endif

#endif
