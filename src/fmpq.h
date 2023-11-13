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
#define FMPQ_INLINE
#else
#define FMPQ_INLINE static __inline__
#endif

#include "fmpz.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Memory management *********************************************************/

FMPQ_INLINE void fmpq_init(fmpq_t x) { x->num = WORD(0); x->den = WORD(1); }
FMPQ_INLINE void fmpq_clear(fmpq_t x) { fmpz_clear(fmpq_numref(x)); fmpz_clear(fmpq_denref(x)); }

/* inlines, for speed, use fmpq_numref and fmpq_denref defined in flint.h */
void fmpq_numerator(fmpz_t n, const fmpq_t q);
void fmpq_denominator(fmpz_t n, const fmpq_t q);
fmpz * fmpq_numerator_ptr(fmpq_t q);
fmpz * fmpq_denominator_ptr(fmpq_t q);

void flint_mpq_init_set_readonly(mpq_t z, const fmpq_t f);
void flint_mpq_clear_readonly(mpq_t z);

void fmpq_init_set_readonly(fmpq_t f, const mpq_t z);
void fmpq_clear_readonly(fmpq_t f);

void fmpq_init_set_mpz_frac_readonly(fmpq_t z, const mpz_t num, const mpz_t den);

void _fmpq_canonicalise(fmpz_t num, fmpz_t den);
void fmpq_canonicalise(fmpq_t res);

int _fmpq_is_canonical(const fmpz_t num, const fmpz_t den);
int fmpq_is_canonical(const fmpq_t x);

/* Randomisation *************************************************************/

void _fmpq_randbits(fmpz_t num, fmpz_t den, flint_rand_t state, flint_bitcnt_t bits);
void fmpq_randbits(fmpq_t res, flint_rand_t state, flint_bitcnt_t bits);

void _fmpq_randtest(fmpz_t num, fmpz_t den, flint_rand_t state, flint_bitcnt_t bits);
void fmpq_randtest(fmpq_t res, flint_rand_t state, flint_bitcnt_t bits);
void fmpq_randtest_not_zero(fmpq_t res, flint_rand_t state, flint_bitcnt_t bits);

/* Assignments and conversions ***********************************************/

FMPQ_INLINE void fmpq_zero(fmpq_t res) { fmpz_zero(fmpq_numref(res)); fmpz_one(fmpq_denref(res)); }
FMPQ_INLINE void fmpq_one(fmpq_t res) { fmpz_one(fmpq_numref(res)); fmpz_one(fmpq_denref(res)); }

FMPQ_INLINE void fmpq_swap(fmpq_t op1, fmpq_t op2) { fmpz_swap(fmpq_numref(op1), fmpq_numref(op2)); fmpz_swap(fmpq_denref(op1), fmpq_denref(op2)); }

void _fmpq_set_ui(fmpz_t rnum, fmpz_t rden, ulong p, ulong q);
void fmpq_set_ui(fmpq_t res, ulong p, ulong q);
void _fmpq_set_si(fmpz_t rnum, fmpz_t rden, slong p, ulong q);
void fmpq_set_si(fmpq_t res, slong p, ulong q);
FMPQ_INLINE void fmpq_set_fmpz(fmpq_t q, const fmpz_t n) { fmpz_set(fmpq_numref(q), n); fmpz_one(fmpq_denref(q)); }
void fmpq_set_fmpz_frac(fmpq_t res, const fmpz_t p, const fmpz_t q);
FMPQ_INLINE void fmpq_set(fmpq_t dest, const fmpq_t src) { fmpz_set(fmpq_numref(dest), fmpq_numref(src)); fmpz_set(fmpq_denref(dest), fmpq_denref(src)); }

FMPQ_INLINE void fmpq_set_mpq(fmpq_t dest, const mpq_t src) { fmpz_set_mpz(fmpq_numref(dest), mpq_numref(src)); fmpz_set_mpz(fmpq_denref(dest), mpq_denref(src)); }
FMPQ_INLINE void fmpq_get_mpq(mpq_t dest, const fmpq_t src) { fmpz_get_mpz(mpq_numref(dest), fmpq_numref(src)); fmpz_get_mpz(mpq_denref(dest), fmpq_denref(src)); }

void fmpq_get_mpz_frac(mpz_t a, mpz_t b, fmpq_t c);

double fmpq_get_d(const fmpq_t a);

#ifdef __MPFR_H
int fmpq_get_mpfr(mpfr_t r, const fmpq_t x, mpfr_rnd_t rnd);
#endif

/* Comparisons ***************************************************************/

FMPQ_INLINE int fmpq_is_zero(const fmpq_t x) { return fmpz_is_zero(fmpq_numref(x)); }
FMPQ_INLINE int fmpq_is_one(const fmpq_t x) { return fmpz_is_one(fmpq_numref(x)) && fmpz_is_one(fmpq_denref(x)); }
FMPQ_INLINE int fmpq_is_pm1(const fmpq_t x) { return fmpz_is_pm1(fmpq_numref(x)) && fmpz_is_one(fmpq_denref(x)); }

FMPQ_INLINE int fmpq_equal_ui(fmpq_t q, ulong n) { return fmpz_equal_ui(fmpq_numref(q), n) && fmpz_is_one(fmpq_denref(q)); }
FMPQ_INLINE int fmpq_equal_si(fmpq_t q, slong n) { return fmpz_equal_si(fmpq_numref(q), n) && fmpz_is_one(fmpq_denref(q)); }
int fmpq_equal_fmpz(fmpq_t q, fmpz_t n);
FMPQ_INLINE int fmpq_equal(const fmpq_t x, const fmpq_t y) { return fmpz_equal(fmpq_numref(x), fmpq_numref(y)) && fmpz_equal(fmpq_denref(x), fmpq_denref(y)); }

int _fmpq_cmp_si(const fmpz_t p, const fmpz_t q, slong c);
int fmpq_cmp_si(const fmpq_t x, slong c);
int _fmpq_cmp_ui(const fmpz_t p, const fmpz_t q, ulong c);
int fmpq_cmp_ui(const fmpq_t x, ulong c);
int _fmpq_cmp_fmpz(const fmpz_t p, const fmpz_t q, const fmpz_t r);
int fmpq_cmp_fmpz(const fmpq_t x, const fmpz_t y);
int _fmpq_cmp(const fmpz_t p, const fmpz_t q, const fmpz_t r, const fmpz_t s);
int fmpq_cmp(const fmpq_t x, const fmpq_t y);

FMPQ_INLINE int fmpq_sgn(const fmpq_t x) { return fmpz_sgn(fmpq_numref(x)); }

flint_bitcnt_t fmpq_height_bits(const fmpq_t x);
void fmpq_height(fmpz_t height, const fmpq_t x);

/* I/O ***********************************************************************/

int fmpq_set_str(fmpq_t res, const char * str, int base);

char * _fmpq_get_str(char * str, int b, const fmpz_t num, const fmpz_t den);
char * fmpq_get_str(char * str, int b, const fmpq_t x);

int _fmpq_print(const fmpz_t num, const fmpz_t den);
int fmpq_print(const fmpq_t x);

#ifdef FLINT_HAVE_FILE
int _fmpq_fprint(FILE * file, const fmpz_t num, const fmpz_t den);
int fmpq_fprint(FILE * file, const fmpq_t x);
#endif

/* Basic arithmetic **********************************************************/

void _fmpq_add_si(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q, slong r);
void fmpq_add_si(fmpq_t res, const fmpq_t op1, slong c);
void _fmpq_add_ui(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q, ulong r);
void fmpq_add_ui(fmpq_t res, const fmpq_t op1, ulong c);
void _fmpq_add_fmpz(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q, const fmpz_t r);
void fmpq_add_fmpz(fmpq_t res, const fmpq_t op1, const fmpz_t c);
void _fmpq_add_small(fmpz_t rnum, fmpz_t rden, slong p1, ulong q1, slong p2, ulong q2);
void _fmpq_add(fmpz_t rnum, fmpz_t rden, const fmpz_t op1num, const fmpz_t op1den, const fmpz_t op2num, const fmpz_t op2den);
void fmpq_add(fmpq_t res, const fmpq_t op1, const fmpq_t op2);

void _fmpq_sub_si(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q, slong r);
void fmpq_sub_si(fmpq_t res, const fmpq_t op1, slong c);
void _fmpq_sub_ui(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q, ulong r);
void fmpq_sub_ui(fmpq_t res, const fmpq_t op1, ulong c);
void _fmpq_sub_fmpz(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q, const fmpz_t r);
void fmpq_sub_fmpz(fmpq_t res, const fmpq_t op1, const fmpz_t c);
void _fmpq_sub(fmpz_t rnum, fmpz_t rden, const fmpz_t op1num, const fmpz_t op1den, const fmpz_t op2num, const fmpz_t op2den);
void fmpq_sub(fmpq_t res, const fmpq_t op1, const fmpq_t op2);

FMPQ_INLINE void fmpq_neg(fmpq_t dest, const fmpq_t src) { fmpz_neg(fmpq_numref(dest), fmpq_numref(src)); fmpz_set(fmpq_denref(dest), fmpq_denref(src)); }
FMPQ_INLINE void fmpq_abs(fmpq_t dest, const fmpq_t src) { fmpz_abs(fmpq_numref(dest), fmpq_numref(src)); fmpz_set(fmpq_denref(dest), fmpq_denref(src)); }

void _fmpq_mul_si(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q, slong r);
void fmpq_mul_si(fmpq_t res, const fmpq_t op1, slong c);
void _fmpq_mul_ui(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q, ulong r);
void fmpq_mul_ui(fmpq_t res, const fmpq_t op1, ulong c);
void fmpq_mul_fmpz(fmpq_t res, const fmpq_t op, const fmpz_t x);
void _fmpq_mul_small(fmpz_t rnum, fmpz_t rden, slong p1, ulong q1, slong p2, ulong q2);
void _fmpq_mul(fmpz_t rnum, fmpz_t rden, const fmpz_t op1num, const fmpz_t op1den, const fmpz_t op2num, const fmpz_t op2den);
void fmpq_mul(fmpq_t res, const fmpq_t op1, const fmpq_t op2);
void fmpq_mul_2exp(fmpq_t res, const fmpq_t x, flint_bitcnt_t exp);

void fmpq_div_fmpz(fmpq_t res, const fmpq_t op, const fmpz_t x);
void _fmpq_div(fmpz_t rnum, fmpz_t rden, const fmpz_t op1num, const fmpz_t op1den, const fmpz_t op2num, const fmpz_t op2den);
void fmpq_div(fmpq_t res, const fmpq_t op1, const fmpq_t op2);
void fmpq_div_2exp(fmpq_t res, const fmpq_t x, flint_bitcnt_t exp);

void fmpq_inv(fmpq_t dest, const fmpq_t src);

void _fmpq_addmul(fmpz_t rnum, fmpz_t rden, const fmpz_t op1num, const fmpz_t op1den, const fmpz_t op2num, const fmpz_t op2den);
void fmpq_addmul(fmpq_t res, const fmpq_t op1, const fmpq_t op2);

void _fmpq_submul(fmpz_t rnum, fmpz_t rden, const fmpz_t op1num, const fmpz_t op1den, const fmpz_t op2num, const fmpz_t op2den);
void fmpq_submul(fmpq_t res, const fmpq_t op1, const fmpq_t op2);

void _fmpq_pow_si(fmpz_t rnum, fmpz_t rden, const fmpz_t opnum, const fmpz_t opden, slong e);
void fmpq_pow_si(fmpq_t rop, const fmpq_t op, slong e);
int fmpq_pow_fmpz(fmpq_t a, const fmpq_t b, const fmpz_t e);

/* Modular reduction and rational reconstruction *****************************/

int _fmpq_mod_fmpz(fmpz_t res, const fmpz_t num, const fmpz_t den, const fmpz_t mod);
int fmpq_mod_fmpz(fmpz_t res, const fmpq_t x, const fmpz_t mod);

int _fmpq_reconstruct_fmpz(fmpz_t num, fmpz_t den, const fmpz_t a, const fmpz_t m);
int fmpq_reconstruct_fmpz(fmpq_t res, const fmpz_t a, const fmpz_t m);

int _fmpq_reconstruct_fmpz_2_naive(fmpz_t n, fmpz_t d, const fmpz_t a, const fmpz_t m, const fmpz_t N, const fmpz_t D);
int _fmpq_reconstruct_fmpz_2(fmpz_t n, fmpz_t d, const fmpz_t a, const fmpz_t m, const fmpz_t N, const fmpz_t D);
int fmpq_reconstruct_fmpz_2(fmpq_t res, const fmpz_t a, const fmpz_t m, const fmpz_t N, const fmpz_t D);

/* GCD ***********************************************************************/

void _fmpq_gcd(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q, const fmpz_t r, const fmpz_t s);
void fmpq_gcd(fmpq_t res, const fmpq_t op1, const fmpq_t op2);

void _fmpq_gcd_cofactors(fmpz_t ng, fmpz_t dg, fmpz_t A, fmpz_t B, const fmpz_t na, const fmpz_t da, const fmpz_t nb, const fmpz_t db);
void fmpq_gcd_cofactors(fmpq_t g, fmpz_t A, fmpz_t B, const fmpq_t a, const fmpq_t b);

/* Rational enumeration ******************************************************/

void _fmpq_next_calkin_wilf(fmpz_t rnum, fmpz_t rden, const fmpz_t num, const fmpz_t den);
void fmpq_next_calkin_wilf(fmpq_t res, const fmpq_t x);

void _fmpq_next_signed_calkin_wilf(fmpz_t rnum, fmpz_t rden, const fmpz_t num, const fmpz_t den);
void fmpq_next_signed_calkin_wilf(fmpq_t res, const fmpq_t x);

void _fmpq_next_minimal(fmpz_t rnum, fmpz_t rden, const fmpz_t num, const fmpz_t den);
void fmpq_next_minimal(fmpq_t res, const fmpq_t x);

void _fmpq_next_signed_minimal(fmpz_t rnum, fmpz_t rden, const fmpz_t num, const fmpz_t den);
void fmpq_next_signed_minimal(fmpq_t res, const fmpq_t x);

void fmpq_farey_neighbors(fmpq_t left, fmpq_t right, const fmpq_t mid, const fmpz_t Q);

void _fmpq_simplest_between(fmpz_t mid_num, fmpz_t mid_den, const fmpz_t l_num, const fmpz_t l_den, const fmpz_t r_num, const fmpz_t r_den);
void fmpq_simplest_between(fmpq_t mid, const fmpq_t l, const fmpq_t r);

/* Continued fractions *******************************************************/

slong fmpq_get_cfrac_naive(fmpz * c, fmpq_t rem, const fmpq_t x, slong n);
slong fmpq_get_cfrac(fmpz * c, fmpq_t rem, const fmpq_t x, slong n);

void fmpq_set_cfrac(fmpq_t x, const fmpz * c, slong n);

slong fmpq_cfrac_bound(const fmpq_t x);

/* Special functions *********************************************************/

void _fmpq_harmonic_ui(fmpz_t num, fmpz_t den, ulong n);
void fmpq_harmonic_ui(fmpq_t x, ulong n);

/* Dedekind sums *************************************************************/

void fmpq_dedekind_sum_naive(fmpq_t s, const fmpz_t h, const fmpz_t k);
void fmpq_dedekind_sum(fmpq_t s, const fmpz_t h, const fmpz_t k);

#ifdef __cplusplus
}
#endif

#endif
