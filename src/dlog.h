/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef DLOG_H
#define DLOG_H

#ifdef DLOG_INLINES_C
#define DLOG_INLINE
#else
#define DLOG_INLINE static inline
#endif

#include "ulong_extras.h"

#ifdef __cplusplus
extern "C" {
#endif

enum
{
    DLOG_MODPE, DLOG_CRT, DLOG_POWER, DLOG_BSGS, DLOG_TABLE, DLOG_23
};

typedef struct dlog_precomp_struct dlog_precomp_struct;
typedef struct dlog_precomp_struct * dlog_precomp_ptr;

/* log in (1+pZ/p^eZ), e small: use recursion formulas
 * could use padic log instead but exponent is small
 * for ulongs */
typedef struct
{
    ulong inv1p;         /* 1 / (1 + p) */
    ulong invloga1;      /* 1 / log(a^(p-1),1+p) */
}
dlog_1modpe_struct;

typedef dlog_1modpe_struct dlog_1modpe_t[1];

/* log in (Z/p^eZ)^* */
typedef struct
{
    ulong p;
    ulong e;
    ulong pe1;                   /* p^(e-1) */
    ulong inva;
    nmod_t pe;
    dlog_precomp_struct * modp;
    dlog_1modpe_t modpe;
}
dlog_modpe_struct;

typedef dlog_modpe_struct dlog_modpe_t[1];

/* all logs precomputed in (Z/modZ)^ast */
typedef struct
{
    ulong mod;
    ulong * table;
}
dlog_table_struct;

typedef dlog_table_struct dlog_table_t[1];

/* bsgs table, already in flint */

typedef struct apow {
    ulong k;
    ulong ak;
} apow_t;

int apow_cmp(const apow_t * x, const apow_t * y);

typedef struct {
    nmod_t mod;
    ulong m;
    ulong am;
    ulong g;
    apow_t * table;
} dlog_bsgs_struct;

typedef dlog_bsgs_struct dlog_bsgs_t[1];
/* typedef bsgs_t dlog_bsgs_t; */

/* Pollard rho */
typedef struct {
    ulong a;
    nmod_t n;
    nmod_t mod;
    int nisprime;
} dlog_rho_struct;

typedef dlog_rho_struct dlog_rho_t[1];

/* CRT decomposition (Pohlig-Hellman) */
typedef struct
{
    nmod_t mod;
    nmod_t n;
    ulong num;
    ulong * expo;
    ulong * crt_coeffs;
    dlog_precomp_ptr pre;
}
dlog_crt_struct;

typedef dlog_crt_struct dlog_crt_t[1];

/* dlog when generator has prime power order */
typedef struct
{
    nmod_t mod;
    ulong p;
    ulong e;
    ulong * apk;
    dlog_precomp_struct * pre;
}
dlog_power_struct;

typedef dlog_power_struct dlog_power_t[1];

typedef ulong dlog_order23_t[1];

/* generic decomposition */
/*typedef */
struct
dlog_precomp_struct
{
    int type;
    ulong cost;
    union
    {
        dlog_table_t table;
        dlog_bsgs_t bsgs;
        dlog_crt_t crt;
        dlog_power_t power;
        dlog_modpe_t modpe;
        dlog_order23_t order23;
    } t;
};

typedef dlog_precomp_struct dlog_precomp_t[1];

void dlog_precomp_modpe_init(dlog_precomp_t pre, ulong a, ulong p, ulong e, ulong pe, ulong num);
void dlog_precomp_small_init(dlog_precomp_t pre, ulong a, ulong mod, ulong n, ulong num);
void dlog_precomp_n_init(dlog_precomp_t pre, ulong a, ulong mod, ulong n, ulong num);
void dlog_precomp_p_init(dlog_precomp_t pre, ulong a, ulong mod, ulong p, ulong num);
void dlog_precomp_pe_init(dlog_precomp_t pre, ulong a, ulong mod, ulong p, ulong e, ulong pe, ulong num);
void dlog_precomp_clear(dlog_precomp_t pre);

ulong dlog_precomp(const dlog_precomp_t pre, ulong b);

ulong dlog_order23_init(dlog_order23_t t, ulong a);
ulong dlog_table_init(dlog_table_t t, ulong a, ulong mod);
ulong dlog_crt_init(dlog_crt_t t, ulong a, ulong mod, ulong n, ulong num);
ulong dlog_power_init(dlog_power_t t, ulong a, ulong mod, ulong p, ulong e, ulong num);
ulong dlog_modpe_init(dlog_modpe_t t, ulong a, ulong p, ulong e, ulong pe, ulong num);
ulong dlog_bsgs_init(dlog_bsgs_t t, ulong a, ulong mod, ulong n, ulong m);
void dlog_1modpe_init(dlog_1modpe_t t, ulong a1, ulong p, ulong e, nmod_t pe);
void dlog_rho_init(dlog_rho_t t, ulong a, ulong mod, ulong n);
/*#define dlog_bsgs_init(t, a, n, m) bsgs_table_init(t, a, n, m)*/

ulong dlog_once(ulong b, ulong a, const nmod_t mod, ulong n);

DLOG_INLINE void
dlog_order23_clear(dlog_order23_t t)
{
    return;
}

DLOG_INLINE void
dlog_table_clear(dlog_table_t t)
{
    flint_free(t->table);
}

void dlog_crt_clear(dlog_crt_t t);

DLOG_INLINE void
dlog_power_clear(dlog_power_t t)
{
    flint_free(t->apk);
    dlog_precomp_clear(t->pre);
    flint_free(t->pre);
}

DLOG_INLINE void
dlog_modpe_clear(dlog_modpe_t t)
{
    dlog_precomp_clear(t->modp);
    flint_free(t->modp);
}

DLOG_INLINE void
dlog_bsgs_clear(dlog_bsgs_t t)
{
    flint_free(t->table);
}

DLOG_INLINE void
dlog_rho_clear(dlog_rho_t t)
{
    return;
}

DLOG_INLINE ulong
dlog_bsgs_size(ulong n, ulong num)
{
    if (2 * num < n)
        return (1 + n_sqrt(n)) * (1 + n_sqrt(num));
    else
        return n;
}

/*#define dlog_bsgs_clear(t) bsgs_table_clear(t)*/

ulong dlog_order23(const dlog_order23_t t, ulong b);
ulong dlog_table(const dlog_table_t t, ulong b);
ulong dlog_crt(const dlog_crt_t t, ulong b);
ulong dlog_power(const dlog_power_t t, ulong b);
ulong dlog_modpe(const dlog_modpe_t t, ulong b);
ulong dlog_bsgs(const dlog_bsgs_t t, ulong b);
ulong dlog_rho(const dlog_rho_t t, ulong b);
ulong dlog_1modpe_1modp(ulong b1, ulong p, ulong e, ulong inv1p, nmod_t pe);
ulong dlog_1modpe(const dlog_1modpe_t t, ulong b1, ulong p, ulong e, nmod_t pe);
ulong dlog_mod2e_1mod4(ulong b1, ulong e, ulong inva, nmod_t pe);
ulong dlog_mod2e(const dlog_modpe_t t, ulong b);
/*#define dlog_bsgs(t, b) n_discrete_log_bsgs_table(t, b)*/

#define DLOG_SMALL_LIM 50
#define DLOG_TABLE_LIM 50
#define DLOG_TABLE_P_LIM 50
#define DLOG_TABLE_MODPE_LIM 50
#define DLOG_TABLE_PE_LIM 50
#define DLOG_TABLE_N_LIM 50
#define DLOG_BSGS_LIM  500

#define DLOG_LOOP_MAX_FACTOR 6
#define DLOG_G_SMALL 0
#define DLOG_G_BIG 1
void dlog_n_factor_group(n_factor_t * fac, ulong bound);

#define DLOG_NOT_FOUND UWORD_MAX
#define DLOG_NONE UWORD_MAX

ulong dlog_vec_pindex_factorgcd(ulong * v, ulong nv, ulong p, nmod_t mod, ulong a, ulong na, ulong loga, ulong logm1, nmod_t order, int maxtry);
void dlog_vec_fill(ulong * v, ulong nv, ulong x);
void dlog_vec_set_not_found(ulong *v, ulong nv, nmod_t mod);
void dlog_vec_loop(ulong * v, ulong nv, ulong a, ulong va, nmod_t mod, ulong na, nmod_t order);
void dlog_vec_loop_add(ulong * v, ulong nv, ulong a, ulong va, nmod_t mod, ulong na, nmod_t order);
void dlog_vec_eratos_add(ulong * v, ulong nv, ulong a, ulong va, nmod_t mod, ulong na, nmod_t order);
void dlog_vec_eratos(ulong * v, ulong nv, ulong a, ulong va, nmod_t mod, ulong na, nmod_t order);
void dlog_vec_sieve_add(ulong * v, ulong nv, ulong a, ulong va, nmod_t mod, ulong na, nmod_t order);
void dlog_vec_sieve(ulong * v, ulong nv, ulong a, ulong va, nmod_t mod, ulong na, nmod_t order);
void dlog_vec_add(ulong * v, ulong nv, ulong a, ulong va, nmod_t mod, ulong na, nmod_t order);
void dlog_vec(ulong * v, ulong nv, ulong a, ulong va, nmod_t mod, ulong na, nmod_t order);


void dlog_vec_sieve_precomp(ulong *v, ulong nv, dlog_precomp_t pre,  ulong a, ulong va, nmod_t mod, ulong na, nmod_t order);
void dlog_vec_sieve_add_precomp(ulong *v, ulong nv, dlog_precomp_t pre, ulong a, ulong va, nmod_t mod, ulong na, nmod_t order);
void dlog_vec_add_precomp(ulong *v, ulong nv, dlog_precomp_t pre, ulong a, ulong va, nmod_t mod, ulong na, nmod_t order);

#ifdef __cplusplus
}
#endif

#endif
