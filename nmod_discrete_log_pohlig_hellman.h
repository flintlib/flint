/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef NMOD_DISCRETE_LOG_POHLIG_HELLMAN_H
#define NMOD_DISCRETE_LOG_POHLIG_HELLMAN_H

#ifdef NMOD_DISCRETE_LOG_POHLIG_HELLMAN_INLINES_C
#define NMOD_DISCRETE_LOG_POHLIG_HELLMAN_INLINE FLINT_DLL
#else
#define NMOD_DISCRETE_LOG_POHLIG_HELLMAN_INLINE static __inline__
#endif

#include "flint.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    ulong gammapow;
    ulong cm;
} nmod_discrete_log_pohlig_hellman_table_entry_struct;

typedef struct {
    slong exp;
    ulong prime;
    ulong gamma;
    ulong gammainv;
    ulong startingbeta;
    ulong co;
    ulong startinge;
    ulong idem;
    ulong cbound;
    ulong dbound;
    nmod_discrete_log_pohlig_hellman_table_entry_struct * table; /* length cbound */
} nmod_discrete_log_pohlig_hellman_entry_struct;

typedef struct {
    nmod_t mod;         /* p is mod.n */
    ulong alpha;    /* p.r. of p */
    ulong alphainv;
    slong num_factors;  /* factors of p - 1*/
    nmod_discrete_log_pohlig_hellman_entry_struct * entries;
} nmod_discrete_log_pohlig_hellman_struct;

typedef nmod_discrete_log_pohlig_hellman_struct nmod_discrete_log_pohlig_hellman_t[1];

FLINT_DLL void nmod_discrete_log_pohlig_hellman_init(
                nmod_discrete_log_pohlig_hellman_t L);

FLINT_DLL void nmod_discrete_log_pohlig_hellman_clear(
                nmod_discrete_log_pohlig_hellman_t L);

FLINT_DLL double nmod_discrete_log_pohlig_hellman_precompute_prime(
                nmod_discrete_log_pohlig_hellman_t L,
                ulong p);

FLINT_DLL ulong nmod_discrete_log_pohlig_hellman_run(
                const nmod_discrete_log_pohlig_hellman_t L,
                ulong y);

NMOD_DISCRETE_LOG_POHLIG_HELLMAN_INLINE
ulong nmod_discrete_log_pohlig_hellman_primitive_root(
        const nmod_discrete_log_pohlig_hellman_t L)
{
    return L->alpha;
}

#ifdef __cplusplus
}
#endif

#endif
