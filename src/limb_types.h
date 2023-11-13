/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef LIMB_TYPES_H
#define LIMB_TYPES_H

#include "flint.h"

#ifdef __cplusplus
extern "C" {
#endif

#define FLINT_MAX_FACTORS_IN_LIMB 15

typedef struct
{
   int num;
   int exp[FLINT_MAX_FACTORS_IN_LIMB];
   ulong p[FLINT_MAX_FACTORS_IN_LIMB];
}
n_factor_t;

typedef struct
{
    slong small_i;
    slong small_num;
    unsigned int * small_primes;

    ulong sieve_a;
    ulong sieve_b;
    slong sieve_i;
    slong sieve_num;
    char * sieve;
}
n_primes_struct;

typedef n_primes_struct n_primes_t[1];

/* arrays of ulong */
typedef struct
{
    mp_limb_t * coeffs;
    slong alloc;
    slong length;
} n_poly_struct;

typedef n_poly_struct n_poly_t[1];
typedef n_poly_struct n_fq_poly_struct;
typedef n_poly_t n_fq_poly_t;

/* arrays of arrays of ulong */
typedef struct
{
    n_poly_struct * coeffs;
    slong alloc;
    slong length;
} n_bpoly_struct;

typedef n_bpoly_struct n_bpoly_t[1];
typedef n_bpoly_struct n_fq_bpoly_struct;
typedef n_bpoly_t n_fq_bpoly_t;

/* sparse arrays of ulong */
typedef struct
{
    ulong * exps;
    mp_limb_t * coeffs;
    slong length;
    slong alloc;
} n_polyu_struct;

typedef n_polyu_struct n_polyu_t[1];
typedef n_polyu_struct n_fq_polyu_struct;
typedef n_polyu_t n_fq_polyu_t;

/*
    sparse arrays of arrays of ulong
    n_polyu1n => one exponent is in the exps[i]
    n_polyu2n => two exponents are packed into the exps[i]
    ...
*/
typedef struct
{
    n_poly_struct * coeffs;
    ulong * exps;
    slong length;
    slong alloc;
} n_polyun_struct;

typedef n_polyun_struct n_polyun_t[1];
typedef n_polyun_struct n_fq_polyun_struct;
typedef n_polyun_t n_fq_polyun_t;

/* arrays of arrays of arrays of ulong */
typedef struct
{
    n_bpoly_struct * coeffs;
    slong alloc;
    slong length;
} n_tpoly_struct;

typedef n_tpoly_struct n_tpoly_t[1];
typedef n_tpoly_struct n_fq_tpoly_struct;
typedef n_tpoly_t n_fq_tpoly_t;

/* n_poly stack */
typedef struct
{
    n_poly_struct ** array;
    slong alloc;
    slong top;
} n_poly_stack_struct;

typedef n_poly_stack_struct n_poly_stack_t[1];

/* n_bpoly stack */
typedef struct
{
    n_bpoly_struct ** array;
    slong alloc;
    slong top;
} n_bpoly_stack_struct;

typedef n_bpoly_stack_struct n_bpoly_stack_t[1];

typedef struct
{
    n_poly_stack_t poly_stack;
    n_bpoly_stack_t bpoly_stack;
} n_poly_bpoly_stack_struct;

typedef n_poly_bpoly_stack_struct n_poly_bpoly_stack_t[1];

typedef struct
{
    mp_limb_t * M;
    mp_limb_t * T;
    mp_limb_t * Q;
    mp_limb_t * array;
    slong alloc;
    slong d;
    slong radix;
    mp_limb_t w;
} nmod_eval_interp_struct;

typedef nmod_eval_interp_struct nmod_eval_interp_t[1];

#ifdef __cplusplus
}
#endif

#endif /* LIMB_TYPES_H */
