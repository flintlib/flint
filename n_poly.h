/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef N_POLY_H
#define N_POLY_H

#ifdef N_POLY_INLINES_C
#define N_POLY_INLINE FLINT_DLL
#else
#define N_POLY_INLINE static __inline__
#endif

#undef ulong
#define ulong ulongxx /* interferes with system includes */
#include <stdio.h>
#undef ulong
#include <gmp.h>
#define ulong mp_limb_t

#include "nmod_mpoly.h"

#ifdef __cplusplus
 extern "C" {
#endif

typedef struct
{
    mp_limb_t * coeffs;
    slong alloc;
    slong length;
} n_poly_struct;

typedef n_poly_struct n_poly_t[1];


typedef struct
{
    n_poly_struct * coeffs;
    slong alloc;
    slong length;
} n_bpoly_struct;

typedef n_bpoly_struct n_bpoly_t[1];


typedef struct
{
    n_bpoly_struct * coeffs;
    slong alloc;
    slong length;
} n_tpoly_struct;

typedef n_tpoly_struct n_tpoly_t[1];


typedef struct
{
    ulong exp;
    mp_limb_t coeff;
} n_polyu_term_struct;

typedef struct
{
    ulong * exps;
    mp_limb_t * coeffs;
    slong length;
    slong alloc;
} n_polyu_struct;

typedef n_polyu_struct n_polyu_t[1];


typedef struct
{
    ulong exp;
    n_poly_t coeff;
} n_polyun_term_struct;

typedef struct
{
    n_polyun_term_struct * terms;
    slong length;
    slong alloc;
} n_polyun_struct;

typedef n_polyun_struct n_polyun_t[1];


typedef struct
{
    n_poly_struct ** array;
    slong alloc;
    slong top;
} n_poly_stack_struct;

typedef n_poly_stack_struct n_poly_stack_t[1];


typedef struct
{
    n_bpoly_struct ** array;
    slong alloc;
    slong top;
} n_bpoly_stack_struct;

typedef n_bpoly_stack_struct n_bpoly_stack_t[1];


typedef struct {
    n_poly_stack_t poly_stack;
    n_bpoly_stack_t bpoly_stack;
} n_poly_bpoly_stack_struct;

typedef n_poly_bpoly_stack_struct n_poly_bpoly_stack_t[1];

/*****************************************************************************/

N_POLY_INLINE
mp_limb_t nmod_addmul(mp_limb_t a, mp_limb_t b, mp_limb_t c, nmod_t mod)
{
    NMOD_ADDMUL(a, b, c, mod);
    return a;
}


N_POLY_INLINE
ulong pack_exp2(ulong e0, ulong e1)
{
    return (e0 << (1*(FLINT_BITS/2))) +
           (e1 << (0*(FLINT_BITS/2)));
}

N_POLY_INLINE
ulong pack_exp3(ulong e0, ulong e1, ulong e2)
{
    return (e0 << (2*(FLINT_BITS/3))) +
           (e1 << (1*(FLINT_BITS/3))) +
           (e2 << (0*(FLINT_BITS/3)));
}

N_POLY_INLINE
ulong extract_exp(ulong e, int idx, int nvars)
{
    return (e >> (idx*(FLINT_BITS/nvars))) &
            ((-UWORD(1)) >> (FLINT_BITS - FLINT_BITS/nvars));
}

/*****************************************************************************/

N_POLY_INLINE
void n_poly_init(n_poly_t A)
{
    A->length = 0;
    A->alloc = 0;
    A->coeffs = NULL;
}

N_POLY_INLINE
void n_poly_init2(n_poly_t A, slong alloc)
{
    A->length = 0;
    A->alloc = alloc;
    A->coeffs = NULL;
    if (alloc > 0)
        A->coeffs = (mp_limb_t *) flint_malloc(alloc*sizeof(mp_limb_t));
}

N_POLY_INLINE
void n_poly_clear(n_poly_t A)
{
    FLINT_ASSERT(A->alloc != 0 || A->coeffs == NULL);
    if (A->alloc > 0)
        flint_free(A->coeffs);
}

FLINT_DLL void n_poly_realloc(n_poly_t A, slong len);

FLINT_DLL void n_poly_print_pretty(const n_poly_t A, const char * x);

N_POLY_INLINE
void n_poly_fit_length(n_poly_t A, slong len)
{
    if (len > A->alloc)
        n_poly_realloc(A, len);
}


N_POLY_INLINE
void nmod_poly_mock(nmod_poly_t a, const n_poly_t b, nmod_t mod)
{
    a->coeffs = b->coeffs;
    a->length = b->length;
    a->alloc = b->alloc;
    a->mod = mod;
}

N_POLY_INLINE
void n_poly_mock(n_poly_t a, const nmod_poly_t b)
{
    a->coeffs = b->coeffs;
    a->length = b->length;
    a->alloc = b->alloc;
}


N_POLY_INLINE
void n_poly_set(n_poly_t A, const n_poly_t B)
{
    n_poly_fit_length(A, B->length);
    flint_mpn_copyi(A->coeffs, B->coeffs, B->length);
    A->length = B->length;
}


N_POLY_INLINE
void n_poly_swap(n_poly_t A, n_poly_t B)
{
    n_poly_struct t = *B;
    *B = *A;
    *A = t;
}

N_POLY_INLINE
void _n_poly_normalise(n_poly_t A)
{
    while (A->length > 0 && A->coeffs[A->length - 1] == 0)
        A->length--;
}

N_POLY_INLINE
slong n_poly_degree(const n_poly_t A)
{
    FLINT_ASSERT(A->length >= 0);
    return A->length - 1;
}

N_POLY_INLINE
int n_poly_is_one(const n_poly_t A)
{
    return A->length == 1 && A->coeffs[0] == 1;
}

N_POLY_INLINE
void n_poly_one(n_poly_t A)
{
    n_poly_fit_length(A, 1);
    A->length = 1;
    A->coeffs[0] = 1;
}

N_POLY_INLINE
void n_poly_set_ui(n_poly_t A, mp_limb_t c)
{
    n_poly_fit_length(A, 1);
    A->coeffs[0] = c;
    A->length = (c != 0);
}

N_POLY_INLINE
int n_poly_is_zero(const n_poly_t poly)
{
    return poly->length == 0;
}

N_POLY_INLINE
void n_poly_zero(n_poly_t res)
{
    res->length = 0;
}

N_POLY_INLINE
int n_poly_equal(const n_poly_t a, const n_poly_t b)
{
    if (a->length != b->length)
        return 0;

    if (a != b)
    {
        if (!_nmod_vec_equal(a->coeffs, b->coeffs, a->length))
            return 0;
    }

    return 1;
}

/*****************************************************************************/

FLINT_DLL int n_poly_mod_is_canonical(const n_poly_t A, nmod_t mod);

N_POLY_INLINE
void n_poly_mod_make_monic(n_poly_t A, const n_poly_t B, nmod_t mod)
{
    FLINT_ASSERT(B->length > 0);
    n_poly_fit_length(A, B->length);
    A->length = B->length;
    _nmod_poly_make_monic(A->coeffs, B->coeffs, B->length, mod);
}

N_POLY_INLINE
void n_poly_mod_taylor_shift(n_poly_t g, mp_limb_t c, nmod_t mod)
{
    _nmod_poly_taylor_shift(g->coeffs, c, g->length, mod);
}

N_POLY_INLINE
ulong n_poly_get_coeff(const n_poly_t poly, slong j)
{
    return (j >= poly->length) ? 0 : poly->coeffs[j];
}

N_POLY_INLINE
void n_poly_set_coeff_nonzero(n_poly_t A, slong j, ulong c)
{
    FLINT_ASSERT(c != 0);
    if (j >= A->length)
    {
        n_poly_fit_length(A, j + 1);
        flint_mpn_zero(A->coeffs + A->length, j - A->length);
        A->length = j + 1;
    }
    A->coeffs[j] = c;
}

FLINT_DLL void n_poly_set_coeff(n_poly_t A, slong e, ulong c);

FLINT_DLL void n_poly_mod_set_coeff_ui(n_poly_t A, slong j, ulong c, nmod_t mod);

N_POLY_INLINE
void n_poly_set_nmod_poly(n_poly_t a, const nmod_poly_t b)
{
    n_poly_fit_length(a, b->length);
    flint_mpn_copyi(a->coeffs, b->coeffs, b->length);
    a->length = b->length;
}

N_POLY_INLINE
void nmod_poly_set_n_poly(nmod_poly_t a, const n_poly_t b)
{
    nmod_poly_fit_length(a, b->length);
    flint_mpn_copyi(a->coeffs, b->coeffs, b->length);
    a->length = b->length;
}

N_POLY_INLINE
void n_poly_shift_left(n_poly_t A, const n_poly_t B, slong k)
{
    n_poly_fit_length(A, B->length + k);
    _nmod_poly_shift_left(A->coeffs, B->coeffs, B->length, k);
    A->length = B->length + k;
}

N_POLY_INLINE
void n_poly_shift_right(n_poly_t res, const n_poly_t poly, slong k)
{
    if (k >= poly->length)
    {
        res->length = 0;
    }
    else
    {
        const slong len = poly->length - k;
        n_poly_fit_length(res, len);
        _nmod_poly_shift_right(res->coeffs, poly->coeffs, len, k);
        res->length = len;
    }
}

N_POLY_INLINE
void n_poly_truncate(n_poly_t poly, slong len)
{
    if (poly->length > len)
    {
        poly->length = len;
        _n_poly_normalise(poly);
    }
}

N_POLY_INLINE
void _n_poly_mod_scalar_mul_nmod(n_poly_t A, const n_poly_t B, mp_limb_t c,
                                                                    nmod_t mod)
{
    FLINT_ASSERT(B->length <= B->alloc);
    n_poly_fit_length(A, B->length);
    _nmod_vec_scalar_mul_nmod(A->coeffs, B->coeffs, B->length, c, mod);
    A->length = B->length;
}

N_POLY_INLINE
void _n_poly_mod_scalar_mul_nmod_inplace(n_poly_t A, mp_limb_t c, nmod_t mod)
{
    _nmod_vec_scalar_mul_nmod(A->coeffs, A->coeffs, A->length, c, mod);
}

FLINT_DLL void n_poly_mod_scalar_mul_ui(n_poly_t A, const n_poly_t B,
                                                      mp_limb_t c, nmod_t ctx);

FLINT_DLL mp_limb_t n_poly_mod_eval_step2(n_poly_t Acur, const n_poly_t Ainc,
                                                                   nmod_t mod);

N_POLY_INLINE
mp_limb_t n_poly_mod_evaluate_nmod(const n_poly_t A, mp_limb_t c, nmod_t mod)
{
    return _nmod_poly_evaluate_nmod(A->coeffs, A->length, c, mod);
}

N_POLY_INLINE
void n_poly_mod_neg(n_poly_t A, const n_poly_t B, nmod_t mod)
{
    n_poly_fit_length(A, B->length);
    _nmod_vec_neg(A->coeffs, B->coeffs, B->length, mod);
    A->length = B->length;
}

N_POLY_INLINE
void n_poly_mod_add(n_poly_t A, const n_poly_t B, const n_poly_t C, nmod_t mod)
{
    slong Alen = FLINT_MAX(B->length, C->length);
    n_poly_fit_length(A, Alen);
    _nmod_poly_add(A->coeffs, B->coeffs, B->length, C->coeffs, C->length, mod);
    A->length = Alen;
    _n_poly_normalise(A);
}

FLINT_DLL void n_poly_mod_add_ui(n_poly_t res, const n_poly_t poly, ulong c, nmod_t ctx);


N_POLY_INLINE
void n_poly_mod_sub(n_poly_t A, const n_poly_t B, const n_poly_t C, nmod_t mod)
{
    slong Alen = FLINT_MAX(B->length, C->length);
    n_poly_fit_length(A, Alen);
    _nmod_poly_sub(A->coeffs, B->coeffs, B->length, C->coeffs, C->length, mod);
    A->length = Alen;
    _n_poly_normalise(A);
}

N_POLY_INLINE
void n_poly_mod_product_roots_nmod_vec(n_poly_t A, mp_srcptr r, slong n, nmod_t mod)
{
    n_poly_fit_length(A, n + 1);
    A->length = n + 1;
    _nmod_poly_product_roots_nmod_vec(A->coeffs, r, n, mod);
}

FLINT_DLL mp_limb_t n_poly_mod_div_root(n_poly_t Q, 
                                    const n_poly_t A, mp_limb_t c, nmod_t ctx);

FLINT_DLL void n_poly_mod_shift_left_scalar_addmul(n_poly_t A, slong k,
                                                      mp_limb_t c, nmod_t mod);

FLINT_DLL void n_poly_mod_addmul_linear(n_poly_t A, const n_poly_t B,
                     const n_poly_t C, mp_limb_t d1, mp_limb_t d0, nmod_t mod);

FLINT_DLL void n_poly_mod_eval2_pow(mp_limb_t * vp, mp_limb_t * vm,
                              const n_poly_t P, n_poly_t alphapow, nmod_t mod);

FLINT_DLL mp_limb_t n_poly_mod_div_root(n_poly_t Q, 
                                    const n_poly_t A, mp_limb_t c, nmod_t ctx);

N_POLY_INLINE
void _n_poly_mod_mul(n_poly_t A, const n_poly_t B, const n_poly_t C, nmod_t mod)
{
    slong Blen = B->length;
    slong Clen = C->length;
    slong Alen = Blen + Clen - 1;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(A != C);

    if (Clen <= 0 || Blen <= 0)
    {
        A->length = 0;
        return;
    }

    n_poly_fit_length(A, Alen);
    A->length = Alen;

    if (Blen >= Clen)
        _nmod_poly_mul(A->coeffs, B->coeffs, Blen, C->coeffs, Clen, mod);
    else
        _nmod_poly_mul(A->coeffs, C->coeffs, Clen, B->coeffs, Blen, mod);
}

N_POLY_INLINE
void _n_poly_mod_div(n_poly_t Q, const n_poly_t A, const n_poly_t B, nmod_t mod)
{
    const slong lenA = A->length, lenB = B->length;
    FLINT_ASSERT(lenB > 0);
    FLINT_ASSERT(Q != A && Q != B);
    if (lenA < lenB)
    {
        n_poly_zero(Q);
        return;
    }
    n_poly_fit_length(Q, lenA - lenB + 1);
    _nmod_poly_div(Q->coeffs, A->coeffs, lenA, B->coeffs, lenB, mod);
    Q->length = lenA - lenB + 1;
}

N_POLY_INLINE
void _n_poly_mod_rem(n_poly_t R, const n_poly_t A, const n_poly_t B, nmod_t mod)
{
    const slong lenA = A->length, lenB = B->length;
    FLINT_ASSERT(R != A && R != B);
    FLINT_ASSERT(lenB > 0);
    if (lenA < lenB)
    {
        n_poly_set(R, A);
        return;
    }
    n_poly_fit_length(R, lenB - 1);
    _nmod_poly_rem(R->coeffs, A->coeffs, lenA, B->coeffs, lenB, mod);
    R->length = lenB - 1;
    _n_poly_normalise(R);
}

N_POLY_INLINE
void _n_poly_mod_divrem(n_poly_t Q, n_poly_t R, const n_poly_t A,
                                                  const n_poly_t B, nmod_t mod)
{
    const slong lenA = A->length, lenB = B->length;

    FLINT_ASSERT(lenB > 0);
    FLINT_ASSERT(Q != A && Q != B);
    FLINT_ASSERT(R != A && R != B);

    if (lenA < lenB)
    {
        n_poly_set(R, A);
        n_poly_zero(Q);
        return;
    }

    n_poly_fit_length(Q, lenA - lenB + 1);
    n_poly_fit_length(R, lenB - 1);
    _nmod_poly_divrem(Q->coeffs, R->coeffs, A->coeffs, lenA, B->coeffs, lenB, mod);
    Q->length = lenA - lenB + 1;
    R->length = lenB - 1;
    _n_poly_normalise(R);
}

FLINT_DLL ulong n_poly_mod_remove(n_poly_t f, const n_poly_t p, nmod_t ctx);

FLINT_DLL void n_poly_mod_pow(n_poly_t res, const n_poly_t poly, ulong e,
                                                                   nmod_t ctx);

FLINT_DLL void n_poly_mod_mul(n_poly_t A, const n_poly_t B, const n_poly_t C,
                                                                   nmod_t mod);

FLINT_DLL void n_poly_mod_mullow(n_poly_t A, const n_poly_t B,
                                        const n_poly_t C, slong n, nmod_t mod);

FLINT_DLL void n_poly_mod_div(n_poly_t Q, const n_poly_t A, const n_poly_t B,
                                                                   nmod_t mod);

FLINT_DLL void n_poly_mod_rem(n_poly_t R, const n_poly_t A, const n_poly_t B,
                                                                   nmod_t mod);

FLINT_DLL void n_poly_mod_divrem(n_poly_t Q, n_poly_t R, const n_poly_t A,
                                                 const n_poly_t B, nmod_t mod);

FLINT_DLL void n_poly_mod_mulmod(n_poly_t res, const n_poly_t poly1,
                           const n_poly_t poly2, const n_poly_t f, nmod_t mod);

FLINT_DLL int n_poly_mod_invmod(n_poly_t A, const n_poly_t B, const n_poly_t P,
                                                                   nmod_t mod);

FLINT_DLL void n_poly_mod_gcd(n_poly_t G, const n_poly_t A, const n_poly_t B,
                                                                   nmod_t mod);

FLINT_DLL void n_poly_mod_xgcd(n_poly_t G, n_poly_t S, n_poly_t T,
                               const n_poly_t A, const n_poly_t B, nmod_t mod);

FLINT_DLL void n_poly_mod_inv_series(n_poly_t Qinv, const n_poly_t Q, slong n,
                                                                   nmod_t mod);

FLINT_DLL void n_poly_mod_div_series(n_poly_t Q, const n_poly_t A,
                                    const n_poly_t B, slong order, nmod_t ctx);

/*****************************************************************************/

FLINT_DLL void n_fq_print_pretty(const mp_limb_t * a, const fq_nmod_ctx_t ctx);

N_POLY_INLINE
int _n_fq_is_zero(const mp_limb_t * a, slong d)
{
    do {
        if (a[--d] != 0)
            return 0;
    } while (d > 0);
    return 1;
}

N_POLY_INLINE
void _n_fq_zero(mp_limb_t * a, slong d)
{
    slong i;
    for (i = 0; i < d; i++)
        a[i] = 0;
} 

N_POLY_INLINE
int _n_fq_is_one(const mp_limb_t * a, slong d)
{
    slong i;
    if (a[0] != 1)
        return 0;
    for (i = 1; i < d; i++)
        if (a[i] != 0)
            return 0;
    return 1;
}

N_POLY_INLINE
void _n_fq_one(mp_limb_t * a, slong d)
{
    slong i;
    a[0] = 1;
    for (i = 1; i < d; i++)
        a[i] = 0;
} 

N_POLY_INLINE
void _n_fq_set(
    mp_limb_t * a,          /* length d */
    const mp_limb_t * b,    /* length d */
    slong d)
{
    slong i = 0;
    do {
        a[i] = b[i];
        i++;
    } while (i < d);
}

N_POLY_INLINE
void _n_fq_swap(
    mp_limb_t * a,      /* length d */
    mp_limb_t * b,      /* length d */
    slong d)
{
    slong i = 0;
    do {
        MP_LIMB_SWAP(a[i], b[i]);
        i++;
    } while (i < d);
}

N_POLY_INLINE
int _n_fq_equal(
    mp_limb_t * a,          /* length d */
    const mp_limb_t * b,    /* length d */
    slong d)
{
    slong i = 0;
    do {
        if (a[i] != b[i])
            return 0;
    } while (++i < d);
    return 1;
}

FLINT_DLL int n_fq_is_canonical(
    const mp_limb_t * a,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_fq_print_pretty(
    const mp_limb_t * a,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_fq_get_fq_nmod(
    fq_nmod_t a,
    const mp_limb_t * b,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_fq_set_fq_nmod(
    mp_limb_t * a,
    const fq_nmod_t b,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_fq_add_si(
    mp_limb_t * a,
    const mp_limb_t * b,
    slong c,
    const fq_nmod_ctx_t ctx);

N_POLY_INLINE
void n_fq_add(
    mp_limb_t * a,          /* length d */
    const mp_limb_t * b,    /* length d */
    const mp_limb_t * c,    /* length d */
    const fq_nmod_ctx_t ctx)
{
    slong d = ctx->modulus->length - 1;
    FLINT_ASSERT(d > 0);
    _nmod_vec_add(a, b, c, d, ctx->modulus->mod);
}

N_POLY_INLINE
void n_fq_sub(
    mp_limb_t * a,          /* length d */
    const mp_limb_t * b,    /* length d */
    const mp_limb_t * c,    /* length d */
    const fq_nmod_ctx_t ctx)
{
    slong d = ctx->modulus->length - 1;
    FLINT_ASSERT(d > 0);
    _nmod_vec_sub(a, b, c, d, ctx->modulus->mod);
}

N_POLY_INLINE
void _n_fq_add(
    mp_limb_t * a,          /* length d */
    const mp_limb_t * b,    /* length d */
    const mp_limb_t * c,    /* length d */
    const fq_nmod_ctx_t ctx)
{
    slong d = ctx->modulus->length - 1;
    FLINT_ASSERT(d > 0);
    _nmod_vec_add(a, b, c, d, ctx->modulus->mod);
}

N_POLY_INLINE
void _n_fq_sub(
    mp_limb_t * a,          /* length d */
    const mp_limb_t * b,    /* length d */
    const mp_limb_t * c,    /* length d */
    const fq_nmod_ctx_t ctx)
{
    slong d = ctx->modulus->length - 1;
    FLINT_ASSERT(d > 0);
    _nmod_vec_sub(a, b, c, d, ctx->modulus->mod);
}

FLINT_DLL void _n_fq_madd2(
    mp_limb_t * a,          /* length 2d-1 */
    const mp_limb_t * b,    /* length d */
    const mp_limb_t * c,    /* length d */
    const fq_nmod_ctx_t ctx,
    mp_limb_t * t);         /* length 2d */

FLINT_DLL void _n_fq_mul2(
    mp_limb_t * t,          /* length 2d-1 */
    const mp_limb_t * b,    /* length d */
    const mp_limb_t * c,    /* length d */
    const fq_nmod_ctx_t ctx);

#define N_FQ_REDUCE_ITCH 2
FLINT_DLL void _n_fq_reduce(
    mp_limb_t * a,
    mp_limb_t * b, slong blen,
    const fq_nmod_ctx_t ctx,
    mp_limb_t * t);

/* same itch as reduce */
N_POLY_INLINE
void _n_fq_reduce2(
    mp_limb_t * a,          /* length d */
    mp_limb_t * b,          /* length 2d-1 */
    const fq_nmod_ctx_t ctx,
    mp_limb_t * t)          /* length 2d */
{
    slong blen = 2*fq_nmod_ctx_degree(ctx) - 1;

    FLINT_ASSERT(a != b);

    while (blen > 0 && b[blen - 1] == 0)
        blen--;

    _n_fq_reduce(a, b, blen, ctx, t);
}


#define N_FQ_MUL_ITCH 4
N_POLY_INLINE
void _n_fq_mul(
    mp_limb_t * a,          /* length d */
    const mp_limb_t * b,    /* length d */
    const mp_limb_t * c,    /* length d */
    const fq_nmod_ctx_t ctx,
    mp_limb_t * t)          /* length 4d */
{
    slong d = fq_nmod_ctx_degree(ctx);
    _n_fq_mul2(t, b, c, ctx);
    _n_fq_reduce2(a, t, ctx, t + 2*d);
}


#define N_FQ_INV_ITCH 1
FLINT_DLL void _n_fq_inv(
    mp_limb_t * a,
    const mp_limb_t * b,
    const fq_nmod_ctx_t ctx,
    mp_limb_t * t);

FLINT_DLL void _n_fq_pow_ui(
    mp_limb_t * a,
    const mp_limb_t * b,
    ulong e,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_fq_mul(
    mp_limb_t * a,
    const mp_limb_t * b,
    const mp_limb_t * c,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_fq_inv(
    mp_limb_t * a,
    const mp_limb_t * b,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_fq_pow_ui(
    mp_limb_t * a,
    const mp_limb_t * b,
    ulong e,
    const fq_nmod_ctx_t ctx);

/*****************************************************************************/

#define N_POLY_FQ_DIVREM_DIVCONQUER_CUTOFF 20

FLINT_DLL void _n_poly_fq_one(n_poly_t A, slong d);

N_POLY_INLINE
void n_poly_fq_one(n_poly_t A, const fq_nmod_ctx_t ctx)
{
    _n_poly_fq_one(A, fq_nmod_ctx_degree(ctx));
}

FLINT_DLL int n_poly_fq_is_one(n_poly_t A, const fq_nmod_ctx_t ctx);

FLINT_DLL int n_poly_fq_is_canonical(
    const n_poly_t a,
    const fq_nmod_ctx_t ctx);

N_POLY_INLINE
void _n_poly_fq_normalise(
    n_poly_t A,
    slong d)
{
    while (A->length > 0 && _n_fq_is_zero(A->coeffs + d*(A->length - 1), d))
        A->length--;
}

FLINT_DLL void n_poly_fq_print_pretty(
    const n_poly_t A,
    const char * x,
    const fq_nmod_ctx_t ctx);

FLINT_DLL int n_poly_fq_equal(
    const n_poly_t A,
    const n_poly_t B,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_poly_fq_set(
    n_poly_t A,
    const n_poly_t B,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_poly_fq_set_n_fq(
    n_poly_t A,
    const mp_limb_t * c,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_poly_fq_randtest(
    n_poly_t A,
    flint_rand_t state,
    slong len,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_poly_fq_make_monic(
    n_poly_t A,
    const n_poly_t B,
    const fq_nmod_ctx_t ctx);

FLINT_DLL ulong n_poly_fq_remove(
    n_poly_t f,
    const n_poly_t g,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_poly_fq_get_coeff_fq_nmod(
    fq_nmod_t c,
    const n_poly_t A,
    slong e,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_poly_fq_set_coeff_n_fq(
    n_poly_t A,
    slong j,
    const mp_limb_t * c,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_poly_fq_set_coeff_fq_nmod(
    n_poly_t A,
    slong j,
    const fq_nmod_t c,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_poly_fq_scalar_mul_n_fq(
    n_poly_t A,
    const n_poly_t B,
    const mp_limb_t * c,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_poly_fq_scalar_mul_ui(
    n_poly_t A,
    const n_poly_t B,
    ulong c,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_poly_fq_shift_left_scalar_submul(
    n_poly_t A,
    slong k,
    const mp_limb_t * c,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_poly_fq_evaluate_n_fq(
    mp_limb_t * e,
    const n_poly_t A,
    const mp_limb_t * c,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_poly_fq_get_fq_nmod_poly(
    fq_nmod_poly_t A,
    const n_poly_t B,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_poly_fq_set_fq_nmod_poly(
    n_poly_t A,
    const fq_nmod_poly_t B,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_poly_fq_set_n_fq(
    n_poly_t A,
    const mp_limb_t * c,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_poly_fq_set_fq_nmod(
    n_poly_t A,
    const fq_nmod_t c,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_poly_fq_shift_right(
    n_poly_t A,
    const n_poly_t B,
    slong n,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_poly_fq_shift_left(
    n_poly_t A,
    const n_poly_t B,
    slong n,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_poly_fq_truncate(
    n_poly_t A,
    slong len,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_poly_fq_add(
    n_poly_t A,
    const n_poly_t B,
    const n_poly_t C,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_poly_fq_sub(
    n_poly_t A,
    const n_poly_t B,
    const n_poly_t C,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_poly_fq_neg(
    n_poly_t A,
    const n_poly_t B,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_poly_fq_add_si(
    n_poly_t A,
    const n_poly_t B,
    slong c,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void _n_poly_fq_mul_(
    mp_limb_t * A,
    const mp_limb_t * B, slong Blen,
    const mp_limb_t * C, slong Clen,
    const fq_nmod_ctx_t ctx,
    n_poly_stack_t St);

FLINT_DLL void n_poly_fq_mul_(
    n_poly_t A,
    const n_poly_t B,
    const n_poly_t C,
    const fq_nmod_ctx_t ctx,
    n_poly_stack_t St);

FLINT_DLL void n_poly_fq_mul(
    n_poly_t A,
    const n_poly_t B,
    const n_poly_t C,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_poly_fq_pow(
    n_poly_t A,
    const n_poly_t B,
    ulong e,
    const fq_nmod_ctx_t ctx);

FLINT_DLL ulong n_poly_fq_remove(
    n_poly_t f,
    const n_poly_t g,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_poly_fq_divrem_divconquer_(
    n_poly_t Q,
    n_poly_t R,
    const n_poly_t A,
    const n_poly_t B,
    const fq_nmod_ctx_t ctx,
    n_poly_stack_t St);

FLINT_DLL void n_poly_fq_divrem(
    n_poly_t Q,
    n_poly_t R,
    const n_poly_t A,
    const n_poly_t B,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_poly_fq_gcd(
    n_poly_t G,
    const n_poly_t A,
    const n_poly_t B,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_poly_fq_xgcd(
    n_poly_t G,
    n_poly_t S,
    n_poly_t T,
    const n_poly_t B,
    const n_poly_t C,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_poly_fq_mulmod(
    n_poly_t A,
    const n_poly_t B,
    const n_poly_t C,
    const n_poly_t M,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_poly_fq_rem(
    n_poly_t R,
    const n_poly_t A,
    const n_poly_t B,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_poly_fq_mullow(
    n_poly_t A,
    const n_poly_t B,
    const n_poly_t C,
    slong order,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_poly_fq_inv_series(
    n_poly_t A,
    const n_poly_t B,
    slong order,
    const fq_nmod_ctx_t ctx);


/*****************************************************************************/

N_POLY_INLINE
void n_bpoly_init(n_bpoly_t A)
{
    A->coeffs = NULL;
    A->alloc = 0;
    A->length = 0;
}

FLINT_DLL void n_bpoly_clear(n_bpoly_t A);

N_POLY_INLINE
void n_bpoly_swap(n_bpoly_t A, n_bpoly_t B)
{
    n_bpoly_struct t = *A;
    *A = *B;
    *B = t;
}

FLINT_DLL void n_bpoly_print_pretty(const n_bpoly_t A,
                                        const char * xvar, const char * yvar);

N_POLY_INLINE
void n_bpoly_normalise(n_bpoly_t A)
{
    while (A->length > 0 && n_poly_is_zero(A->coeffs + A->length - 1))
        A->length--;
}

FLINT_DLL void n_bpoly_realloc(n_bpoly_t A, slong len);

N_POLY_INLINE
void n_bpoly_fit_length(n_bpoly_t A, slong len)
{
    if (len > A->alloc)
        n_bpoly_realloc(A, len);
}

N_POLY_INLINE
void n_bpoly_zero(n_bpoly_t A)
{
    A->length = 0;
}

N_POLY_INLINE
int n_bpoly_is_zero(const n_bpoly_t A)
{
    return A->length == 0;
}

FLINT_DLL void _n_bpoly_set(n_bpoly_t A, const n_bpoly_t B);

N_POLY_INLINE
void n_bpoly_set(n_bpoly_t A, const n_bpoly_t B)
{
    if (A != B)
        _n_bpoly_set(A, B);
}

FLINT_DLL void n_bpoly_one(n_bpoly_t A);

FLINT_DLL int n_bpoly_equal(const n_bpoly_t A, const n_bpoly_t B);

FLINT_DLL void n_bpoly_set_coeff(n_bpoly_t A, slong e0, slong e1, mp_limb_t c);

FLINT_DLL void n_bpoly_set_coeff_nonzero(n_bpoly_t A, slong e0, slong e1,
                                                                  mp_limb_t c);

FLINT_DLL void n_bpoly_mod_derivative(n_bpoly_t A, const n_bpoly_t B,
                                                                   nmod_t ctx);

N_POLY_INLINE
mp_limb_t n_bpoly_get_coeff(const n_bpoly_t A, slong e0, slong e1)
{
    if (e0 >= A->length)
        return 0;
    else
        return n_poly_get_coeff(A->coeffs + e0, e1);
}

N_POLY_INLINE
slong n_bpoly_degree0(const n_bpoly_t A)
{
    return A->length - 1;
}

FLINT_DLL slong n_bpoly_degree1(const n_bpoly_t A);

FLINT_DLL void n_bpoly_set_poly_var1(n_bpoly_t A, const n_poly_t B);

FLINT_DLL void n_bpoly_set_poly_var0(n_bpoly_t A, const n_poly_t B);

/*****************************************************************************/

FLINT_DLL int n_bpoly_mod_is_canonical(const n_bpoly_t A, nmod_t mod);

N_POLY_INLINE
ulong n_bpoly_bidegree(const n_bpoly_t A)
{
    ulong x, y;
    FLINT_ASSERT(A->length > 0);
    x = A->length - 1;
    y = A->coeffs[x].length - 1;
    return (x << (FLINT_BITS/2)) + y;
}

FLINT_DLL void n_bpoly_scalar_mul_nmod(n_bpoly_t A, mp_limb_t c, nmod_t ctx);

FLINT_DLL void n_bpoly_mod_content_last(n_poly_t g, const n_bpoly_t A, nmod_t ctx);

FLINT_DLL void n_bpoly_mod_divexact_last(n_bpoly_t A, const n_poly_t b, nmod_t ctx);

FLINT_DLL void n_bpoly_mod_mul_last(n_bpoly_t A, const n_poly_t b, nmod_t ctx);


FLINT_DLL void n_bpoly_mod_taylor_shift_var1(n_bpoly_t A, const n_bpoly_t B,
                                                      mp_limb_t c, nmod_t ctx);

FLINT_DLL void n_bpoly_mod_taylor_shift_var0(n_bpoly_t A, mp_limb_t c,
                                                                   nmod_t ctx);

FLINT_DLL void n_bpoly_mod_add(n_bpoly_t A, const n_bpoly_t B,
                                                const n_bpoly_t C, nmod_t ctx);

FLINT_DLL void n_bpoly_mod_sub(n_bpoly_t A, const n_bpoly_t B,
                                                const n_bpoly_t C, nmod_t ctx);

FLINT_DLL void n_bpoly_mod_make_primitive(n_poly_t g, n_bpoly_t A, nmod_t ctx);

FLINT_DLL void n_bpoly_mod_mul(n_bpoly_t A, const n_bpoly_t B,
                                                const n_bpoly_t C, nmod_t ctx);

FLINT_DLL int n_bpoly_mod_divides(n_bpoly_t Q, const n_bpoly_t A,
                                                const n_bpoly_t B, nmod_t ctx);

FLINT_DLL void n_bpoly_mod_mul_series(n_bpoly_t A, const n_bpoly_t B,
                                   const n_bpoly_t C, slong order, nmod_t ctx);

FLINT_DLL void n_bpoly_mod_divrem_series(n_bpoly_t Q, n_bpoly_t R,
                const n_bpoly_t A, const n_bpoly_t B, slong order, nmod_t ctx);

FLINT_DLL void n_bpoly_mod_interp_reduce_2sm_poly(n_poly_t Ap, n_poly_t Am,
                             const n_bpoly_t A, n_poly_t alphapow, nmod_t mod);

FLINT_DLL void n_bpoly_mod_interp_lift_2sm_poly(slong * deg1, n_bpoly_t T,
              const n_poly_t A, const n_poly_t B, mp_limb_t alpha, nmod_t mod);

FLINT_DLL int n_bpoly_mod_interp_crt_2sm_poly(slong * deg1, n_bpoly_t F,
                 n_bpoly_t T, n_poly_t A, n_poly_t B, const n_poly_t modulus,
                                                n_poly_t alphapow, nmod_t mod);

FLINT_DLL int n_bpoly_mod_gcd_brown_smprime(n_bpoly_t G, n_bpoly_t Abar,
                                    n_bpoly_t Bbar, n_bpoly_t A, n_bpoly_t B,
                                          nmod_t ctx, n_poly_bpoly_stack_t Sp);

/*****************************************************************************/

FLINT_DLL void n_bpoly_fq_get_coeff_fq_nmod(
    fq_nmod_t c,
    const n_bpoly_t A,
    slong e0,
    slong e1,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_bpoly_fq_set_fq_nmod_poly_var0(
    n_bpoly_t A,
    const fq_nmod_poly_t B,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_bpoly_fq_set_n_poly_fq_var0(
    n_bpoly_t A,
    const n_poly_t B,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_bpoly_fq_set_fq_nmod_poly_var1(
    n_bpoly_t A,
    const fq_nmod_poly_t B,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_bpoly_fq_set_n_poly_fq_var1(
    n_bpoly_t A,
    const n_poly_t B,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_bpoly_fq_taylor_shift_var1_fq_nmod(
    n_bpoly_t A,
    const n_bpoly_t B,
    const fq_nmod_t c_,
    const fq_nmod_ctx_t ctx);

FLINT_DLL void n_bpoly_fq_taylor_shift_var0_fq_nmod(
    n_bpoly_t A,
    const fq_nmod_t alpha,
    const fq_nmod_ctx_t ctx);

/*****************************************************************************/

N_POLY_INLINE
void n_tpoly_init(n_tpoly_t A)
{
    A->coeffs = NULL;
    A->alloc = 0;
    A->length = 0;
}

N_POLY_INLINE
void n_tpoly_swap(n_tpoly_t A, n_tpoly_t B)
{
    n_tpoly_struct t = *A;
    *A = *B;
    *B = t;
}

FLINT_DLL void n_tpoly_fit_length(n_tpoly_t A, slong len);

FLINT_DLL void n_tpoly_clear(n_tpoly_t A);

/*****************************************************************************/

N_POLY_INLINE
void n_polyu_init(n_polyu_t A)
{
    A->coeffs = NULL;
    A->exps = NULL;
    A->length = 0;
    A->alloc = 0;
}

FLINT_DLL void n_polyu_clear(n_polyu_t A);

FLINT_DLL void n_polyu_realloc(n_polyu_t A, slong len);

N_POLY_INLINE
void n_polyu_fit_length(n_polyu_t A, slong len)
{
    FLINT_ASSERT(A->alloc >= 0);
    if (len > A->alloc)
        n_polyu_realloc(A, len);
}

N_POLY_INLINE
void n_polyu_swap(n_polyu_t A, n_polyu_t B)
{
    n_polyu_struct T = *B;
    *B = *A;
    *A = T;
}

FLINT_DLL void n_polyu3_print_pretty(const n_polyu_t A, const char * var0,
                                         const char * var1, const char * var2);

FLINT_DLL void n_polyu3_degrees(slong * deg0, slong * deg1, slong * deg2,
                                                            const n_polyu_t A);

/*****************************************************************************/

N_POLY_INLINE
void n_polyun_init(n_polyun_t A)
{
    A->terms = NULL;
    A->length = 0;
    A->alloc = 0;
}

FLINT_DLL void n_polyun_clear(n_polyun_t A);

FLINT_DLL void n_polyun_realloc(n_polyun_t A, slong len);

N_POLY_INLINE
void n_polyun_fit_length(n_polyun_t A, slong len)
{
    if (len > A->alloc)
        n_polyun_realloc(A, len);
}

N_POLY_INLINE
void n_polyun_term_swap(n_polyun_term_struct * A, n_polyun_term_struct * B)
{
    n_polyun_term_struct T = *A;
    *A = *B;
    *B = T;
}

N_POLY_INLINE
void n_polyun_swap(n_polyun_t A, n_polyun_t B)
{
    n_polyun_struct t = *B;
    *B = *A;
    *A = t;
}

FLINT_DLL void n_polyu2n_print_pretty(const n_polyun_t A, const char * var0,
                                      const char * var1, const char * varlast);

FLINT_DLL void n_polyu3n_print_pretty(const n_polyun_t A, const char * var0,
                   const char * var1, const char * var2, const char * varlast);

/*****************************************************************************/

FLINT_DLL void n_poly_stack_init(n_poly_stack_t S);

FLINT_DLL void n_poly_stack_clear(n_poly_stack_t S);

FLINT_DLL n_poly_struct ** n_poly_stack_fit_request(n_poly_stack_t S, slong k);

N_POLY_INLINE
mp_limb_t * n_poly_stack_vec_init(n_poly_stack_t S, slong len)
{
    n_poly_struct * poly_top;
    poly_top = n_poly_stack_fit_request(S, 1)[0];
    S->top += 1;
    n_poly_fit_length(poly_top, len);
    return poly_top->coeffs;
}

N_POLY_INLINE
void n_poly_stack_vec_clear(n_poly_stack_t S)
{
    FLINT_ASSERT(S->top >= 1);
    S->top -= 1;
}

N_POLY_INLINE
n_poly_struct ** n_poly_stack_request(n_poly_stack_t S, slong k)
{
    n_poly_struct ** poly_top;
    poly_top = n_poly_stack_fit_request(S, k);
    S->top += k;
    return poly_top;
}

N_POLY_INLINE
n_poly_struct * n_poly_stack_take_top(n_poly_stack_t S)
{
    /* assume the request for 1 has already been fitted */
    n_poly_struct ** poly_top;
    FLINT_ASSERT(S->top + 1 <= S->alloc);
    poly_top = S->array + S->top;
    S->top += 1;
    return poly_top[0];
}

N_POLY_INLINE
void n_poly_stack_give_back(n_poly_stack_t S, slong k)
{
    FLINT_ASSERT(S->top >= k);
    S->top -= k;
}

N_POLY_INLINE
slong n_poly_stack_size(const n_poly_stack_t S)
{
    return S->top;
}

/*****************************************************************************/

FLINT_DLL void n_bpoly_stack_init(n_bpoly_stack_t S);

FLINT_DLL void n_bpoly_stack_clear(n_bpoly_stack_t S);

FLINT_DLL n_bpoly_struct ** n_bpoly_stack_fit_request(n_bpoly_stack_t S, slong k);

N_POLY_INLINE
n_bpoly_struct ** n_bpoly_stack_request(n_bpoly_stack_t S, slong k)
{
    n_bpoly_struct ** bpoly_top;
    bpoly_top = n_bpoly_stack_fit_request(S, k);
    S->top += k;
    return bpoly_top;
}

N_POLY_INLINE
n_bpoly_struct * n_bpoly_stack_take_top(n_bpoly_stack_t S)
{
    /* assume the request for 1 has already been fitted */
    n_bpoly_struct ** bpoly_top;
    FLINT_ASSERT(S->top + 1 <= S->alloc);
    bpoly_top = S->array + S->top;
    S->top += 1;
    return bpoly_top[0];
}

N_POLY_INLINE
void n_bpoly_stack_give_back(n_bpoly_stack_t S, slong k)
{
    FLINT_ASSERT(S->top >= k);
    S->top -= k;
}

N_POLY_INLINE
slong n_bpoly_stack_size(const n_bpoly_stack_t S)
{
    return S->top;
}

/*****************************************************************************/

typedef struct {
    ulong key;
    slong up;
    slong left;
    slong right;
    int color;
} mpoly_rbnode_ui_struct;

typedef struct {
    slong length;
    mpoly_rbnode_ui_struct * nodes;
    slong node_alloc;
    char * data;
    slong data_alloc;
} mpoly_rbtree_ui_struct;

typedef mpoly_rbtree_ui_struct mpoly_rbtree_ui_t[1];

FLINT_DLL void mpoly_rbtree_ui_init(mpoly_rbtree_ui_t T);

FLINT_DLL void mpoly_rbtree_ui_clear(mpoly_rbtree_ui_t T);

FLINT_DLL void * mpoly_rbtree_ui_lookup(
    mpoly_rbtree_ui_t T,
    int * its_new,
    ulong rcx,
    slong dsize);


#ifdef __cplusplus
}
#endif

#endif

