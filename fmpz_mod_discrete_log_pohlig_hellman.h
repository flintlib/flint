/*
    Copyright (C) 2017 - 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_MOD_DISCRETE_LOG_POHLIG_HELLMAN_H
#define FMPZ_MOD_DISCRETE_LOG_POHLIG_HELLMAN_H

#ifdef FMPZ_MOD_DISCRETE_LOG_POHLIG_HELLMAN_INLINES_C
#define FMPZ_MOD_DISCRETE_LOG_POHLIG_HELLMAN_INLINE FLINT_DLL
#else
#define FMPZ_MOD_DISCRETE_LOG_POHLIG_HELLMAN_INLINE static __inline__
#endif

#include "flint.h"

#ifdef __cplusplus
 extern "C" {
#endif

typedef struct {
    fmpz_t gammapow;
    ulong cm;
} fmpz_mod_discrete_log_pohlig_hellman_table_entry_struct;

typedef struct {
    slong exp;
    ulong prime;
    fmpz_t gamma;
    fmpz_t gammainv;
    fmpz_t startingbeta;
    fmpz_t co;
    fmpz_t startinge;
    fmpz_t idem;
    ulong cbound;
    ulong dbound;
    fmpz_mod_discrete_log_pohlig_hellman_table_entry_struct * table; /* length cbound */
} fmpz_mod_discrete_log_pohlig_hellman_entry_struct;

typedef struct {
    fmpz_mod_ctx_t fpctx;
    fmpz_t pm1;      /* p - 1 */
    fmpz_t alpha;    /* p.r. of p */
    fmpz_t alphainv;
    slong num_factors;  /* factors of p - 1*/
    fmpz_mod_discrete_log_pohlig_hellman_entry_struct * entries;
} fmpz_mod_discrete_log_pohlig_hellman_struct;
typedef fmpz_mod_discrete_log_pohlig_hellman_struct fmpz_mod_discrete_log_pohlig_hellman_t[1];

FLINT_DLL void fmpz_mod_discrete_log_pohlig_hellman_init(
                    fmpz_mod_discrete_log_pohlig_hellman_t L);

FLINT_DLL void fmpz_mod_discrete_log_pohlig_hellman_clear(
                    fmpz_mod_discrete_log_pohlig_hellman_t L);

FLINT_DLL double fmpz_mod_discrete_log_pohlig_hellman_precompute_prime(
                    fmpz_mod_discrete_log_pohlig_hellman_t L,
                    const fmpz_t p);

FLINT_DLL void fmpz_mod_discrete_log_pohlig_hellman_run(
                    fmpz_t x,
                    const fmpz_mod_discrete_log_pohlig_hellman_t L,
                    const fmpz_t y);

FMPZ_MOD_DISCRETE_LOG_POHLIG_HELLMAN_INLINE
const fmpz * fmpz_mod_discrete_log_pohlig_hellman_primitive_root(
                    fmpz_mod_discrete_log_pohlig_hellman_t L)
{
    return L->alpha;
}

#ifdef __cplusplus
}
#endif

#endif
