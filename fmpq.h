/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPQ_H
#define FMPQ_H

#ifdef FMPQ_INLINES_C
#define FMPQ_INLINE FLINT_DLL
#else
#define FMPQ_INLINE static __inline__
#endif

#undef ulong
#define ulong ulongxx /* interferes with system includes */
#include <stdio.h>
#undef ulong

#include <gmp.h>
#include <mpfr.h>
#define ulong mp_limb_t
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"


#ifdef __cplusplus
 extern "C" {
#endif

typedef struct
{
    fmpz num;
    fmpz den;
}
fmpq;

typedef fmpq fmpq_t[1];

#define fmpq_numref(__x) (&(__x)->num)
#define fmpq_denref(__y) (&(__y)->den)

FMPQ_INLINE void fmpq_init(fmpq_t x)
{
    x->num = WORD(0);
    x->den = WORD(1);
}

FMPQ_INLINE void fmpq_clear(fmpq_t x)
{
    fmpz_clear(fmpq_numref(x));
    fmpz_clear(fmpq_denref(x));
}

FMPQ_INLINE void fmpq_zero(fmpq_t res)
{
    fmpz_zero(fmpq_numref(res));
    fmpz_one(fmpq_denref(res));
}

FMPQ_INLINE void fmpq_one(fmpq_t res)
{
    fmpz_one(fmpq_numref(res));
    fmpz_one(fmpq_denref(res));
}

FMPQ_INLINE int fmpq_equal(const fmpq_t x, const fmpq_t y)
{
    return fmpz_equal(fmpq_numref(x), fmpq_numref(y)) &&
           fmpz_equal(fmpq_denref(x), fmpq_denref(y));
}

FMPQ_INLINE int fmpq_sgn(const fmpq_t x)
{
    return fmpz_sgn(fmpq_numref(x));
}

FMPQ_INLINE int fmpq_is_zero(const fmpq_t x)
{
    return fmpz_is_zero(fmpq_numref(x));
}

FMPQ_INLINE int fmpq_is_one(const fmpq_t x)
{
    return fmpz_is_one(fmpq_numref(x)) && fmpz_is_one(fmpq_denref(x));
}

FMPQ_INLINE int fmpq_is_pm1(const fmpq_t x)
{
    return fmpz_is_pm1(fmpq_numref(x)) && fmpz_is_one(fmpq_denref(x));
}

FMPQ_INLINE void fmpq_set(fmpq_t dest, const fmpq_t src)
{
    fmpz_set(fmpq_numref(dest), fmpq_numref(src));
    fmpz_set(fmpq_denref(dest), fmpq_denref(src));
}

FMPQ_INLINE void fmpq_swap(fmpq_t op1, fmpq_t op2)
{
    fmpz_swap(fmpq_numref(op1), fmpq_numref(op2));
    fmpz_swap(fmpq_denref(op1), fmpq_denref(op2));
}

FMPQ_INLINE void fmpq_neg(fmpq_t dest, const fmpq_t src)
{
    fmpz_neg(fmpq_numref(dest), fmpq_numref(src));
    fmpz_set(fmpq_denref(dest), fmpq_denref(src));
}

FMPQ_INLINE void fmpq_abs(fmpq_t dest, const fmpq_t src)
{
    fmpz_abs(fmpq_numref(dest), fmpq_numref(src));
    fmpz_set(fmpq_denref(dest), fmpq_denref(src));
}

FLINT_DLL int _fmpq_cmp(const fmpz_t p, const fmpz_t q, const fmpz_t r, const fmpz_t s);

FLINT_DLL int fmpq_cmp(const fmpq_t x, const fmpq_t y);

FLINT_DLL int _fmpq_cmp_fmpz(const fmpz_t p, const fmpz_t q, const fmpz_t r);

FLINT_DLL int fmpq_cmp_fmpz(const fmpq_t x, const fmpz_t y);

FLINT_DLL int _fmpq_cmp_ui(const fmpz_t p, const fmpz_t q, ulong c);

FLINT_DLL int fmpq_cmp_ui(const fmpq_t x, ulong c);

FLINT_DLL int _fmpq_cmp_si(const fmpz_t p, const fmpz_t q, slong c);

FLINT_DLL int fmpq_cmp_si(const fmpq_t x, slong c);

FLINT_DLL void _fmpq_canonicalise(fmpz_t num, fmpz_t den);

FLINT_DLL void fmpq_canonicalise(fmpq_t res);

FLINT_DLL int _fmpq_is_canonical(const fmpz_t num, const fmpz_t den);

FLINT_DLL int fmpq_is_canonical(const fmpq_t x);

FLINT_DLL void _fmpq_set_ui(fmpz_t rnum, fmpz_t rden, ulong p, ulong q);

FLINT_DLL void fmpq_set_ui(fmpq_t res, ulong p, ulong q);

FLINT_DLL void _fmpq_set_si(fmpz_t rnum, fmpz_t rden, slong p, ulong q);

FLINT_DLL void fmpq_set_si(fmpq_t res, slong p, ulong q);

FMPQ_INLINE int fmpq_equal_ui(fmpq_t q, slong n)
{
    return fmpz_equal_ui(fmpq_numref(q), n) && q->den == WORD(1);
}

FMPQ_INLINE int fmpq_equal_si(fmpq_t q, slong n)
{
    return fmpz_equal_si(fmpq_numref(q), n) && q->den == WORD(1);
}

FLINT_DLL void fmpq_set_fmpz_frac(fmpq_t res, const fmpz_t p, const fmpz_t q);

FLINT_DLL int fmpq_set_str(fmpq_t res, const char * str, int base);

FMPQ_INLINE void fmpq_set_mpq(fmpq_t dest, const mpq_t src)
{
    fmpz_set_mpz(fmpq_numref(dest), mpq_numref(src));
    fmpz_set_mpz(fmpq_denref(dest), mpq_denref(src));
}

FMPQ_INLINE void fmpq_get_mpq(mpq_t dest, const fmpq_t src)
{
    fmpz_get_mpz(mpq_numref(dest), fmpq_numref(src));
    fmpz_get_mpz(mpq_denref(dest), fmpq_denref(src));
}

FLINT_DLL double fmpq_get_d(const fmpq_t a);

FLINT_DLL int fmpq_get_mpfr(mpfr_t r, const fmpq_t x, mpfr_rnd_t rnd);

FLINT_DLL void fmpq_get_mpz_frac(mpz_t a, mpz_t b, fmpq_t c);

FLINT_DLL void flint_mpq_init_set_readonly(mpq_t z, const fmpq_t f);

FLINT_DLL void flint_mpq_clear_readonly(mpq_t z);

FLINT_DLL void fmpq_init_set_readonly(fmpq_t f, const mpq_t z);

FLINT_DLL void fmpq_clear_readonly(fmpq_t f);

FLINT_DLL void fmpq_init_set_mpz_frac_readonly(fmpq_t z, const mpz_t num, const mpz_t den);

FLINT_DLL char * _fmpq_get_str(char * str, int b, const fmpz_t num, const fmpz_t den);

FLINT_DLL char * fmpq_get_str(char * str, int b, const fmpq_t x);

FLINT_DLL int _fmpq_fprint(FILE * file, const fmpz_t num, const fmpz_t den);

FLINT_DLL int fmpq_fprint(FILE * file, const fmpq_t x);

FMPQ_INLINE int _fmpq_print(const fmpz_t num, const fmpz_t den)
{
    return _fmpq_fprint(stdout, num, den);
}

FMPQ_INLINE int fmpq_print(const fmpq_t x)
{
    return fmpq_fprint(stdout, x);
}

FLINT_DLL void _fmpq_randtest(fmpz_t num, fmpz_t den, flint_rand_t state, flint_bitcnt_t bits);

FLINT_DLL void fmpq_randtest(fmpq_t res, flint_rand_t state, flint_bitcnt_t bits);

FLINT_DLL void fmpq_randtest_not_zero(fmpq_t res, flint_rand_t state, flint_bitcnt_t bits);

FLINT_DLL void _fmpq_randbits(fmpz_t num, fmpz_t den, flint_rand_t state, flint_bitcnt_t bits);

FLINT_DLL void fmpq_randbits(fmpq_t res, flint_rand_t state, flint_bitcnt_t bits);

FLINT_DLL void _fmpq_add_small(fmpz_t rnum, fmpz_t rden, slong p1, ulong q1, slong p2, ulong q2);

FLINT_DLL void _fmpq_mul_small(fmpz_t rnum, fmpz_t rden, slong p1, ulong q1, slong p2, ulong q2);

FLINT_DLL void _fmpq_add(fmpz_t rnum, fmpz_t rden, const fmpz_t op1num,
    const fmpz_t op1den, const fmpz_t op2num, const fmpz_t op2den);

FLINT_DLL void fmpq_add(fmpq_t res, const fmpq_t op1, const fmpq_t op2);

FLINT_DLL void _fmpq_add_si(fmpz_t rnum, fmpz_t rden, const fmpz_t p, 
                                                      const fmpz_t q, slong r);

FLINT_DLL void fmpq_add_si(fmpq_t res, const fmpq_t op1, slong c);

FLINT_DLL void _fmpq_add_ui(fmpz_t rnum, fmpz_t rden, const fmpz_t p,
		                                      const fmpz_t q, ulong r);

FLINT_DLL void fmpq_add_ui(fmpq_t res, const fmpq_t op1, ulong c);

FLINT_DLL void _fmpq_add_fmpz(fmpz_t rnum, fmpz_t rden, const fmpz_t p, 
                                               const fmpz_t q, const fmpz_t r);

FLINT_DLL void fmpq_add_fmpz(fmpq_t res, const fmpq_t op1, const fmpz_t c);

FLINT_DLL void _fmpq_sub(fmpz_t rnum, fmpz_t rden, const fmpz_t op1num,
    const fmpz_t op1den, const fmpz_t op2num, const fmpz_t op2den);

FLINT_DLL void fmpq_sub(fmpq_t res, const fmpq_t op1, const fmpq_t op2);

FLINT_DLL void _fmpq_sub_si(fmpz_t rnum, fmpz_t rden, const fmpz_t p, 
                                                      const fmpz_t q, slong r);

FLINT_DLL void fmpq_sub_si(fmpq_t res, const fmpq_t op1, slong c);

FLINT_DLL void _fmpq_sub_ui(fmpz_t rnum, fmpz_t rden, const fmpz_t p,
		                                      const fmpz_t q, ulong r);

FLINT_DLL void fmpq_sub_ui(fmpq_t res, const fmpq_t op1, ulong c);

FLINT_DLL void _fmpq_sub_fmpz(fmpz_t rnum, fmpz_t rden, const fmpz_t p, 
                                               const fmpz_t q, const fmpz_t r);

FLINT_DLL void fmpq_sub_fmpz(fmpq_t res, const fmpq_t op1, const fmpz_t c);

FLINT_DLL void _fmpq_mul_si(fmpz_t rnum, fmpz_t rden, const fmpz_t p,
                                                      const fmpz_t q, slong r);

FLINT_DLL void fmpq_mul_si(fmpq_t res, const fmpq_t op1, slong c);

FLINT_DLL void _fmpq_mul_ui(fmpz_t rnum, fmpz_t rden, const fmpz_t p,
                                                      const fmpz_t q, ulong r);

FLINT_DLL void fmpq_mul_ui(fmpq_t res, const fmpq_t op1, ulong c);

FLINT_DLL void _fmpq_mul(fmpz_t rnum, fmpz_t rden, const fmpz_t op1num,
    const fmpz_t op1den, const fmpz_t op2num, const fmpz_t op2den);

FLINT_DLL void fmpq_mul(fmpq_t res, const fmpq_t op1, const fmpq_t op2);


FLINT_DLL void fmpq_mul_fmpz(fmpq_t res, const fmpq_t op, const fmpz_t x);

FLINT_DLL void _fmpq_pow_si(fmpz_t rnum, fmpz_t rden, 
                  const fmpz_t opnum, const fmpz_t opden, slong e);

FLINT_DLL void fmpq_pow_si(fmpq_t rop, const fmpq_t op, slong e);

FLINT_DLL int fmpq_pow_fmpz(fmpq_t a, const fmpq_t b, const fmpz_t e);

FLINT_DLL void _fmpq_addmul(fmpz_t rnum, fmpz_t rden, const fmpz_t op1num,
    const fmpz_t op1den, const fmpz_t op2num, const fmpz_t op2den);

FLINT_DLL void fmpq_addmul(fmpq_t res, const fmpq_t op1, const fmpq_t op2);


FLINT_DLL void _fmpq_submul(fmpz_t rnum, fmpz_t rden, const fmpz_t op1num,
    const fmpz_t op1den, const fmpz_t op2num, const fmpz_t op2den);

FLINT_DLL void fmpq_submul(fmpq_t res, const fmpq_t op1, const fmpq_t op2);


FLINT_DLL void fmpq_inv(fmpq_t dest, const fmpq_t src);

FLINT_DLL void _fmpq_div(fmpz_t rnum, fmpz_t rden, const fmpz_t op1num,
                const fmpz_t op1den, const fmpz_t op2num, const fmpz_t op2den);

FLINT_DLL void fmpq_div(fmpq_t res, const fmpq_t op1, const fmpq_t op2);

FLINT_DLL void fmpq_div_fmpz(fmpq_t res, const fmpq_t op, const fmpz_t x);


FLINT_DLL void fmpq_mul_2exp(fmpq_t res, const fmpq_t x, flint_bitcnt_t exp);

FLINT_DLL void fmpq_div_2exp(fmpq_t res, const fmpq_t x, flint_bitcnt_t exp);


FLINT_DLL int _fmpq_mod_fmpz(fmpz_t res, const fmpz_t num, const fmpz_t den, const fmpz_t mod);

FLINT_DLL int fmpq_mod_fmpz(fmpz_t res, const fmpq_t x, const fmpz_t mod);

FMPQ_INLINE void
_fmpq_gcd(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q,
            const fmpz_t r, const fmpz_t s)
{
   fmpz_t a, b;
   fmpz_init(a); fmpz_init(b);
   fmpz_mul(a, p, s);
   fmpz_mul(b, q, r);
   fmpz_gcd(rnum, a, b);
   fmpz_mul(rden, q, s);
   _fmpq_canonicalise(rnum, rden);
   fmpz_clear(a); fmpz_clear(b);
}

FMPQ_INLINE void 
fmpq_gcd(fmpq_t res, const fmpq_t op1, const fmpq_t op2)
{
    _fmpq_gcd(fmpq_numref(res), fmpq_denref(res), fmpq_numref(op1), 
              fmpq_denref(op1), fmpq_numref(op2), fmpq_denref(op2));
}

FLINT_DLL int _fmpq_reconstruct_fmpz(fmpz_t num, fmpz_t den, const fmpz_t a, const fmpz_t m);

FLINT_DLL int fmpq_reconstruct_fmpz(fmpq_t res, const fmpz_t a, const fmpz_t m);

FLINT_DLL int _fmpq_reconstruct_fmpz_2_naive(fmpz_t n, fmpz_t d,
               const fmpz_t a, const fmpz_t m, const fmpz_t N, const fmpz_t D);

FLINT_DLL int _fmpq_reconstruct_fmpz_2(fmpz_t n, fmpz_t d,
               const fmpz_t a, const fmpz_t m, const fmpz_t N, const fmpz_t D);

FLINT_DLL int fmpq_reconstruct_fmpz_2(fmpq_t res,
               const fmpz_t a, const fmpz_t m, const fmpz_t N, const fmpz_t D);

FLINT_DLL flint_bitcnt_t fmpq_height_bits(const fmpq_t x);

FLINT_DLL void fmpq_height(fmpz_t height, const fmpq_t x);

FLINT_DLL void _fmpq_next_calkin_wilf(fmpz_t rnum, fmpz_t rden,
    const fmpz_t num, const fmpz_t den);

FLINT_DLL void fmpq_next_calkin_wilf(fmpq_t res, const fmpq_t x);

FLINT_DLL void _fmpq_next_signed_calkin_wilf(fmpz_t rnum, fmpz_t rden,
    const fmpz_t num, const fmpz_t den);

FLINT_DLL void fmpq_next_signed_calkin_wilf(fmpq_t res, const fmpq_t x);

FLINT_DLL void _fmpq_next_minimal(fmpz_t rnum, fmpz_t rden,
    const fmpz_t num, const fmpz_t den);

FLINT_DLL void fmpq_next_minimal(fmpq_t res, const fmpq_t x);

FLINT_DLL void _fmpq_next_signed_minimal(fmpz_t rnum, fmpz_t rden,
    const fmpz_t num, const fmpz_t den);

FLINT_DLL void fmpq_next_signed_minimal(fmpq_t res, const fmpq_t x);

FLINT_DLL void fmpq_farey_neighbors(fmpq_t left, fmpq_t right,
                                             const fmpq_t mid, const fmpz_t Q);

FLINT_DLL void fmpq_simplest_between(fmpq_t mid, const fmpq_t l, const fmpq_t r);

FLINT_DLL void _fmpq_simplest_between(fmpz_t mid_num, fmpz_t mid_den,
                                       const fmpz_t l_num, const fmpz_t l_den,
                                       const fmpz_t r_num, const fmpz_t r_den);

FLINT_DLL slong fmpq_get_cfrac_naive(fmpz * c, fmpq_t rem, const fmpq_t x, slong n);

FLINT_DLL slong fmpq_get_cfrac(fmpz * c, fmpq_t rem, const fmpq_t x, slong n);

FLINT_DLL void fmpq_set_cfrac(fmpq_t x, const fmpz * c, slong n);

FLINT_DLL slong fmpq_cfrac_bound(const fmpq_t x);

FLINT_DLL void fmpq_dedekind_sum_naive(fmpq_t s, const fmpz_t h, const fmpz_t k);

FLINT_DLL void fmpq_dedekind_sum(fmpq_t s, const fmpz_t h, const fmpz_t k);

FLINT_DLL void _fmpq_harmonic_ui(fmpz_t num, fmpz_t den, ulong n);

FLINT_DLL void fmpq_harmonic_ui(fmpq_t x, ulong n);

FLINT_DLL fmpq * _fmpq_vec_init(slong len);

FMPQ_INLINE
void _fmpq_vec_clear(fmpq * vec, slong len)
{
    _fmpz_vec_clear((fmpz *) vec, 2 * len);
}

/*********************** 2x2 integer matrix **********************************/

typedef struct {
    fmpz_t _11, _12, _21, _22;
    int det;    /* 0,1,or,-1: 0 -> don't know, 1 -> 1, -1 -> -1 */
} _fmpz_mat22_struct;

typedef _fmpz_mat22_struct _fmpz_mat22_t[1];

typedef struct {
    mp_limb_t _11, _12, _21, _22;
    int det;    /* ditto */
} _ui_mat22_struct;

typedef _ui_mat22_struct _ui_mat22_t[1];

FLINT_DLL void _fmpz_mat22_init(_fmpz_mat22_t M);

FLINT_DLL void _fmpz_mat22_clear(_fmpz_mat22_t M);

FLINT_DLL void _fmpz_mat22_one(_fmpz_mat22_t M);

FLINT_DLL int _fmpz_mat22_is_one(_fmpz_mat22_t M);

FLINT_DLL flint_bitcnt_t _fmpz_mat22_bits(const _fmpz_mat22_t N);

FLINT_DLL void _fmpz_mat22_rmul(_fmpz_mat22_t M, const _fmpz_mat22_t N);

FLINT_DLL void _fmpz_mat22_rmul_ui(_fmpz_mat22_t M, const _ui_mat22_t N);

FLINT_DLL void _fmpz_mat22_rmul_inv_ui(_fmpz_mat22_t M, const _ui_mat22_t N);

FLINT_DLL void _fmpz_mat22_rmul_elem(_fmpz_mat22_t M, const fmpz_t q);

FLINT_DLL void _fmpz_mat22_rmul_inv_elem(_fmpz_mat22_t M, const fmpz_t q);

FLINT_DLL void _fmpz_mat22_lmul_elem(_fmpz_mat22_t M, const fmpz_t q);

/******** resizable integer vector specific to cfrac functionality ***********/

typedef struct
{
    fmpz * array;
    slong length;
    slong alloc;
    slong limit;
    fmpz_t alt_sum;
    int want_alt_sum;
} _fmpq_cfrac_list_struct;

typedef _fmpq_cfrac_list_struct _fmpq_cfrac_list_t[1];

FLINT_DLL void _fmpq_cfrac_list_init(_fmpq_cfrac_list_t v);

FLINT_DLL void _fmpq_cfrac_list_clear(_fmpq_cfrac_list_t v);

FLINT_DLL void _fmpq_cfrac_list_fit_length(_fmpq_cfrac_list_t v, slong len);

FLINT_DLL void _fmpq_cfrac_list_push_back(_fmpq_cfrac_list_t v, const fmpz_t a);

FLINT_DLL void _fmpq_cfrac_list_push_back_zero(_fmpq_cfrac_list_t v);

FLINT_DLL void _fmpq_cfrac_list_append_ui(_fmpq_cfrac_list_t v, const ulong * a, slong n);

FMPQ_INLINE void _fmpq_cfrac_list_swap(_fmpq_cfrac_list_t a, _fmpq_cfrac_list_t b)
{
    _fmpq_cfrac_list_struct t = *a;
    *a = *b;
    *b = t;
}

/*************** ball for closed interval [left, right] **********************/

typedef struct {
    fmpz_t left_num, left_den, right_num, right_den;
    int exact;
} _fmpq_ball_struct;

typedef _fmpq_ball_struct _fmpq_ball_t[1];

FLINT_DLL void _fmpq_ball_init(_fmpq_ball_t x);

FLINT_DLL void _fmpq_ball_clear(_fmpq_ball_t x);

FMPQ_INLINE void _fmpq_ball_swap(_fmpq_ball_t x, _fmpq_ball_t y)
{
   _fmpq_ball_struct t = *x;
   *x = *y;
   *y = t;
}

FLINT_DLL int _fmpq_ball_gt_one(const _fmpq_ball_t x);

FLINT_DLL void _fmpq_hgcd(_fmpq_cfrac_list_t s, _fmpz_mat22_t M,
                                                   fmpz_t x_num, fmpz_t x_den);

FLINT_DLL void _fmpq_ball_get_cfrac(_fmpq_cfrac_list_t s, _fmpz_mat22_t M,
                                                    int needM, _fmpq_ball_t x);

/* Inlines *******************************************************************/

FLINT_DLL void fmpq_numerator(fmpz_t n, const fmpq_t q);
FLINT_DLL void fmpq_denominator(fmpz_t n, const fmpq_t q);
FLINT_DLL fmpz * fmpq_numerator_ptr(fmpq_t q);
FLINT_DLL fmpz * fmpq_denominator_ptr(fmpq_t q);
FLINT_DLL int fmpq_equal_fmpz(fmpq_t q, fmpz_t n);

#ifdef __cplusplus
}
#endif

#endif
