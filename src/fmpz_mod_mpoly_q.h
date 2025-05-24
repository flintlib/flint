/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_MOD_MPOLY_Q_H
#define FMPZ_MOD_MPOLY_Q_H

#ifdef FMPZ_MOD_MPOLY_Q_INLINES_C
#define FMPZ_MOD_MPOLY_Q_INLINE
#else
#define FMPZ_MOD_MPOLY_Q_INLINE static inline
#endif

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mod_mpoly.h"
#include "acb_types.h"

#ifdef __cplusplus
extern "C" {
#endif

#define fmpz_mod_mpoly_q_numref(x) (&((x)->num))
#define fmpz_mod_mpoly_q_denref(x) (&((x)->den))

/* Memory management */

void fmpz_mod_mpoly_q_init(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_ctx_t ctx);

void fmpz_mod_mpoly_q_clear(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_ctx_t ctx);

/* Properties */


FMPZ_MOD_MPOLY_Q_INLINE int
fmpz_mod_mpoly_q_is_zero(const fmpz_mod_mpoly_q_t x, const fmpz_mod_mpoly_ctx_t ctx)
{
    return fmpz_mod_mpoly_is_zero(fmpz_mod_mpoly_q_numref(x), ctx);
}

FMPZ_MOD_MPOLY_Q_INLINE int
fmpz_mpoly_q_is_one(const fmpz_mod_mpoly_q_t x, const fmpz_mod_mpoly_ctx_t ctx)
{
    return fmpz_mod_mpoly_is_one(fmpz_mod_mpoly_q_numref(x), ctx) &&
           fmpz_mod_mpoly_is_one(fmpz_mod_mpoly_q_denref(x), ctx);
}

FMPZ_MOD_MPOLY_Q_INLINE int
fmpz_mod_mpoly_q_is_fmpz(const fmpz_mod_mpoly_q_t x, const fmpz_mod_mpoly_ctx_t ctx)
{
    return fmpz_mod_mpoly_is_fmpz(fmpz_mod_mpoly_q_numref(x), ctx) &&
           fmpz_mod_mpoly_is_one(fmpz_mod_mpoly_q_denref(x), ctx);
}

FMPZ_MOD_MPOLY_Q_INLINE int
fmpz_mod_mpoly_q_is_fmpq(const fmpz_mod_mpoly_q_t x, const fmpz_mod_mpoly_ctx_t ctx)
{
    return fmpz_mod_mpoly_is_fmpz(fmpz_mod_mpoly_q_numref(x), ctx) &&
           fmpz_mod_mpoly_is_fmpz(fmpz_mod_mpoly_q_denref(x), ctx);
}

void fmpz_mod_mpoly_q_used_vars(int * used, const fmpz_mod_mpoly_q_t f, const fmpz_mod_mpoly_ctx_t ctx);
void fmpz_mod_mpoly_q_used_vars_num(int * used, const fmpz_mod_mpoly_q_t f, const fmpz_mod_mpoly_ctx_t ctx);
void fmpz_mod_mpoly_q_used_vars_den(int * used, const fmpz_mod_mpoly_q_t f, const fmpz_mod_mpoly_ctx_t ctx);

/* Special values */

FMPZ_MOD_MPOLY_Q_INLINE void
fmpz_mod_mpoly_q_zero(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_mod_mpoly_zero(fmpz_mod_mpoly_q_numref(res), ctx);
    fmpz_mod_mpoly_one(fmpz_mod_mpoly_q_denref(res), ctx);
}

FMPZ_MOD_MPOLY_Q_INLINE void
fmpz_mod_mpoly_q_one(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_mod_mpoly_one(fmpz_mod_mpoly_q_numref(res), ctx);
    fmpz_mod_mpoly_one(fmpz_mod_mpoly_q_denref(res), ctx);
}

FMPZ_MOD_MPOLY_Q_INLINE void
fmpz_mod_mpoly_q_gen(fmpz_mod_mpoly_q_t res, slong i, const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_mod_mpoly_gen(fmpz_mod_mpoly_q_numref(res), i, ctx);
    fmpz_mod_mpoly_one(fmpz_mod_mpoly_q_denref(res), ctx);
}

/* Input and output */

void fmpz_mod_mpoly_q_print_pretty(const fmpz_mod_mpoly_q_t f, const char ** x, const fmpz_mod_mpoly_ctx_t ctx);
char * fmpz_mod_mpoly_q_get_str_pretty(const fmpz_mod_mpoly_q_t f, const char ** vars, const fmpz_mod_mpoly_ctx_t ctx);
int fmpz_mod_mpoly_q_set_str_pretty(fmpz_mod_mpoly_q_t res, const char * s, const char ** vars, fmpz_mod_mpoly_ctx_t ctx);


#ifdef __cplusplus
}
#endif

#endif
