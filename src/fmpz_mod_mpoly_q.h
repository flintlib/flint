/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_MPOLY_Q_H
#define FMPZ_MPOLY_Q_H

#ifdef FMPZ_MPOLY_Q_INLINES_C
#define FMPZ_MPOLY_Q_INLINE
#else
#define FMPZ_MPOLY_Q_INLINE static inline
#endif

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mpoly.h"
#include "acb_types.h"

#ifdef __cplusplus
extern "C" {
#endif

#define fmpz_mod_mpoly_q_numref(x) (&((x)->num))
#define fmpz_mod_mpoly_q_denref(x) (&((x)->den))

/* Memory management */

void fmpz_mod_mpoly_q_init(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_q_clear(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_ctx_t ctx);



#ifdef __cplusplus
}
#endif

#endif
