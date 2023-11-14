/*
    Copyright (C) 2009, 2010, 2011 Sebastian Pancratz
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_POLY_Q_H
#define FMPZ_POLY_Q_H

#ifdef FMPZ_POLY_Q_INLINES_C
#define FMPZ_POLY_Q_INLINE
#else
#define FMPZ_POLY_Q_INLINE static inline
#endif

#include "fmpz_poly.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Accessing numerator and denominator ***************************************/

#define fmpz_poly_q_numref(op)  ((op)->num)

#define fmpz_poly_q_denref(op)  ((op)->den)

void fmpz_poly_q_canonicalise(fmpz_poly_q_t rop);

int fmpz_poly_q_is_canonical(const fmpz_poly_q_t op);

/* Memory management *********************************************************/

void fmpz_poly_q_init(fmpz_poly_q_t rop);

void fmpz_poly_q_clear(fmpz_poly_q_t rop);

/* Randomisation *************************************************************/

void fmpz_poly_q_randtest(fmpz_poly_q_t poly, flint_rand_t state,
                          slong len1, flint_bitcnt_t bits1,
                          slong len2, flint_bitcnt_t bits2);

void fmpz_poly_q_randtest_not_zero(fmpz_poly_q_t poly, flint_rand_t state,
                                   slong len1, flint_bitcnt_t bits1,
                                   slong len2, flint_bitcnt_t bits2);

/* Assignment ****************************************************************/

void fmpz_poly_q_set(fmpz_poly_q_t rop, const fmpz_poly_q_t op);

void fmpz_poly_q_set_si(fmpz_poly_q_t rop, slong op);

void fmpz_poly_q_swap(fmpz_poly_q_t op1, fmpz_poly_q_t op2);

FMPZ_POLY_Q_INLINE
void fmpz_poly_q_zero(fmpz_poly_q_t rop)
{
    fmpz_poly_zero(rop->num);
    fmpz_poly_set_si(rop->den, 1);
}

FMPZ_POLY_Q_INLINE
void fmpz_poly_q_one(fmpz_poly_q_t rop)
{
    fmpz_poly_set_si(rop->num, 1);
    fmpz_poly_set_si(rop->den, 1);
}

FMPZ_POLY_Q_INLINE
void fmpz_poly_q_neg(fmpz_poly_q_t rop, const fmpz_poly_q_t op)
{
    fmpz_poly_neg(rop->num, op->num);
    fmpz_poly_set(rop->den, op->den);
}

void fmpz_poly_q_inv(fmpz_poly_q_t rop, const fmpz_poly_q_t op);

/* Comparison ****************************************************************/

FMPZ_POLY_Q_INLINE
int fmpz_poly_q_is_zero(const fmpz_poly_q_t op)
{
    return fmpz_poly_is_zero(op->num);
}

FMPZ_POLY_Q_INLINE
int fmpz_poly_q_is_one(const fmpz_poly_q_t op)
{
    return fmpz_poly_is_one(op->num) && fmpz_poly_is_one(op->den);
}

FMPZ_POLY_Q_INLINE
int fmpz_poly_q_equal(const fmpz_poly_q_t op1, const fmpz_poly_q_t op2)
{
    return fmpz_poly_equal(op1->num, op2->num) && fmpz_poly_equal(op1->den, op2->den);
}

/* Addition and subtraction **************************************************/

void fmpz_poly_q_add_in_place(fmpz_poly_q_t rop, const fmpz_poly_q_t op);
void fmpz_poly_q_sub_in_place(fmpz_poly_q_t rop, const fmpz_poly_q_t op);

void fmpz_poly_q_add(fmpz_poly_q_t rop, const fmpz_poly_q_t op1, const fmpz_poly_q_t op2);
void fmpz_poly_q_sub(fmpz_poly_q_t rop, const fmpz_poly_q_t op1, const fmpz_poly_q_t op2);

void fmpz_poly_q_addmul(fmpz_poly_q_t rop, const fmpz_poly_q_t op1, const fmpz_poly_q_t op2);
void fmpz_poly_q_submul(fmpz_poly_q_t rop, const fmpz_poly_q_t op1, const fmpz_poly_q_t op2);

/* Scalar multiplication and division ****************************************/

void fmpz_poly_q_scalar_mul_si(fmpz_poly_q_t rop, const fmpz_poly_q_t op, slong x);
void fmpz_poly_q_scalar_mul_fmpz(fmpz_poly_q_t rop, const fmpz_poly_q_t op, const fmpz_t x);
void fmpz_poly_q_scalar_mul_fmpq(fmpz_poly_q_t rop, const fmpz_poly_q_t op, const fmpq_t x);

void fmpz_poly_q_scalar_div_si(fmpz_poly_q_t rop, const fmpz_poly_q_t op, slong x);
void fmpz_poly_q_scalar_div_fmpz(fmpz_poly_q_t rop, const fmpz_poly_q_t op, const fmpz_t x);
void fmpz_poly_q_scalar_div_fmpq(fmpz_poly_q_t rop, const fmpz_poly_q_t op, const fmpq_t x);

/* Multiplication and division ***********************************************/

void fmpz_poly_q_mul(fmpz_poly_q_t rop,
                     const fmpz_poly_q_t op1, const fmpz_poly_q_t op2);

void fmpz_poly_q_div(fmpz_poly_q_t rop,
                     const fmpz_poly_q_t op1, const fmpz_poly_q_t op2);

/* Powering ******************************************************************/

void fmpz_poly_q_pow(fmpz_poly_q_t rop, const fmpz_poly_q_t op, ulong exp);

/* Derivative ****************************************************************/

void fmpz_poly_q_derivative(fmpz_poly_q_t rop, const fmpz_poly_q_t op);

/* Evaluation ****************************************************************/

int fmpz_poly_q_evaluate_fmpq(fmpq_t rop, const fmpz_poly_q_t f, const fmpq_t a);

/* Input and output **********************************************************/

int fmpz_poly_q_set_str(fmpz_poly_q_t rop, const char *s);

char * fmpz_poly_q_get_str(const fmpz_poly_q_t op);
char * fmpz_poly_q_get_str_pretty(const fmpz_poly_q_t op, const char *x);

int fmpz_poly_q_print(const fmpz_poly_q_t op);
int fmpz_poly_q_print_pretty(const fmpz_poly_q_t op, const char *x);

/* Declare old functions dead *************************************************/

#define fmpz_poly_q_scalar_mul_mpz _Pragma("GCC error \"'fmpz_poly_q_scalar_mul_mpz' is deprecated. Use 'fmpz_poly_q_scalar_mul_fmpz' instead.\"")
#define fmpz_poly_q_scalar_mul_mpq _Pragma("GCC error \"'fmpz_poly_q_scalar_mul_mpq' is deprecated. Use 'fmpz_poly_q_scalar_mul_fmpq' instead.\"")
#define fmpz_poly_q_scalar_div_mpz _Pragma("GCC error \"'fmpz_poly_q_scalar_div_mpz' is deprecated. Use 'fmpz_poly_q_scalar_div_fmpz' instead.\"")
#define fmpz_poly_q_scalar_div_mpq _Pragma("GCC error \"'fmpz_poly_q_scalar_div_mpq' is deprecated. Use 'fmpz_poly_q_scalar_div_fmpq' instead.\"")
#define fmpz_poly_q_evaluate _Pragma("GCC error \"'fmpz_poly_q_evaluate' is deprecated. Use 'fmpz_poly_q_evaluate_fmpq' instead.\"")

#ifdef __cplusplus
}
#endif

#endif
