/*
    Copyright (C) 2017 - 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_MOD_H
#define FMPZ_MOD_H

#ifdef FMPZ_MOD_INLINES_C
#define FMPZ_MOD_INLINE
#else
#define FMPZ_MOD_INLINE static inline
#endif

#include "fmpz_mod_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* all of the data we need to do arithmetic mod n ****************************/

/*
    Currently operations are special cased according to
        if      n < 2^FLINT_BITS     -> add1, sub1, mul1 using nmod
        else if n = 2^FLINT_BITS     -> add2s, sub2s, mul2s
        else if n < 2^(2*FLINT_BITS) -> add2, sub2, mul2
        else                         -> addN, subN, mulN

    A special case for the multiplication for 3-word shows no signs of
    diminishing returns, but it is not implemented currently.
*/

void fmpz_mod_ctx_init(fmpz_mod_ctx_t ctx, const fmpz_t n);
void fmpz_mod_ctx_init_ui(fmpz_mod_ctx_t ctx, ulong n);
void fmpz_mod_ctx_init_rand_bits(fmpz_mod_ctx_t ctx, flint_rand_t state, flint_bitcnt_t max_bits);
void fmpz_mod_ctx_init_rand_bits_prime(fmpz_mod_ctx_t ctx, flint_rand_t state, flint_bitcnt_t max_bits);

void fmpz_mod_ctx_clear(fmpz_mod_ctx_t ctx);

FMPZ_MOD_INLINE const fmpz * fmpz_mod_ctx_modulus(const fmpz_mod_ctx_t ctx)
{
    return ctx->n;
}

void fmpz_mod_ctx_set_modulus(fmpz_mod_ctx_t ctx, const fmpz_t n);
void fmpz_mod_ctx_set_modulus_ui(fmpz_mod_ctx_t ctx, ulong n);

int fmpz_mod_is_canonical(const fmpz_t a, const fmpz_mod_ctx_t ctx);

void fmpz_mod_assert_canonical(const fmpz_t a, const fmpz_mod_ctx_t ctx);

int fmpz_mod_is_one(const fmpz_t a, const fmpz_mod_ctx_t ctx);

int fmpz_mod_equal_fmpz(const fmpz_t a, const fmpz_t b, const fmpz_mod_ctx_t ctx);
int fmpz_mod_equal_si(const fmpz_t a, slong b, const fmpz_mod_ctx_t ctx);

void fmpz_mod_set_fmpz(fmpz_t a, const fmpz_t b, const fmpz_mod_ctx_t ctx);
void fmpz_mod_set_ui(fmpz_t a, ulong b, const fmpz_mod_ctx_t ctx);
void fmpz_mod_set_si(fmpz_t a, slong b, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_add1(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx);
void _fmpz_mod_add2s(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx);
void _fmpz_mod_add2(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx);
void _fmpz_mod_addN(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx);

FMPZ_MOD_INLINE
void fmpz_mod_add(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx)
{
    (ctx->add_fxn)(a, b, c, ctx);
}

void fmpz_mod_add_fmpz(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx);
void fmpz_mod_add_ui(fmpz_t a, const fmpz_t b, ulong c, const fmpz_mod_ctx_t ctx);
void fmpz_mod_add_si(fmpz_t a, const fmpz_t b, slong c, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_sub1(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx);
void _fmpz_mod_sub2s(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx);
void _fmpz_mod_sub2(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx);
void _fmpz_mod_subN(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx);

FMPZ_MOD_INLINE
void fmpz_mod_sub(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx)
{
    (ctx->sub_fxn)(a, b, c, ctx);
}

void fmpz_mod_sub_fmpz(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx);
void fmpz_mod_sub_ui(fmpz_t a, const fmpz_t b, ulong c, const fmpz_mod_ctx_t ctx);
void fmpz_mod_sub_si(fmpz_t a, const fmpz_t b, slong c, const fmpz_mod_ctx_t ctx);
void fmpz_mod_fmpz_sub(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx);
void fmpz_mod_ui_sub(fmpz_t a, ulong b, const fmpz_t c, const fmpz_mod_ctx_t ctx);
void fmpz_mod_si_sub(fmpz_t a, slong b, const fmpz_t c, const fmpz_mod_ctx_t ctx);

void fmpz_mod_neg(fmpz_t a, const fmpz_t b, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_mul1(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx);
void _fmpz_mod_mul2s(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx);
void _fmpz_mod_mul2(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx);
void _fmpz_mod_mulN(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx);

FMPZ_MOD_INLINE
void fmpz_mod_mul(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx)
{
    (ctx->mul_fxn)(a, b, c, ctx);
}

void fmpz_mod_mul_fmpz(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx);
void fmpz_mod_mul_ui(fmpz_t a, const fmpz_t b, ulong c, const fmpz_mod_ctx_t ctx);
void fmpz_mod_mul_si(fmpz_t a, const fmpz_t b, slong c, const fmpz_mod_ctx_t ctx);

void fmpz_mod_addmul(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_t d, const fmpz_mod_ctx_t ctx);

int fmpz_mod_is_invertible(const fmpz_t a, const fmpz_mod_ctx_t ctx);

void fmpz_mod_inv(fmpz_t a, const fmpz_t b, const fmpz_mod_ctx_t ctx);

int fmpz_mod_divides(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx);

void fmpz_mod_pow_ui(fmpz_t a, const fmpz_t b, ulong pow, const fmpz_mod_ctx_t ctx);
int fmpz_mod_pow_fmpz(fmpz_t a, const fmpz_t b, const fmpz_t pow, const fmpz_mod_ctx_t ctx);

void fmpz_mod_rand(fmpz_t a, flint_rand_t state, const fmpz_mod_ctx_t ctx);
void fmpz_mod_rand_not_zero(fmpz_t a, flint_rand_t state, const fmpz_mod_ctx_t ctx);

/* discrete logs a la Pohlig - Hellman ***************************************/

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

void fmpz_mod_discrete_log_pohlig_hellman_init(fmpz_mod_discrete_log_pohlig_hellman_t L);

void fmpz_mod_discrete_log_pohlig_hellman_clear(fmpz_mod_discrete_log_pohlig_hellman_t L);

double fmpz_mod_discrete_log_pohlig_hellman_precompute_prime(fmpz_mod_discrete_log_pohlig_hellman_t L, const fmpz_t p);

void fmpz_mod_discrete_log_pohlig_hellman_run(fmpz_t x, const fmpz_mod_discrete_log_pohlig_hellman_t L, const fmpz_t y);

FMPZ_MOD_INLINE
const fmpz * fmpz_mod_discrete_log_pohlig_hellman_primitive_root(fmpz_mod_discrete_log_pohlig_hellman_t L)
{
    return L->alpha;
}

int fmpz_next_smooth_prime(fmpz_t a, const fmpz_t b);

/* Declare dead functions *****************************************************/

#define fmpz_mod_ctx_get_modulus_mpz_read_only _Pragma("GCC error \"'fmpz_mod_ctx_get_modulus_mpz_read_only' is deprecated. Use 'fmpz_mod_ctx_modulus' instead.\"")

#ifdef __cplusplus
}
#endif

#endif
