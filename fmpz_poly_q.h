/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2009, 2010, 2011 Sebastian Pancratz

******************************************************************************/

#ifndef FMPZ_POLY_Q_H
#define FMPZ_POLY_Q_H

#ifdef FMPZ_POLY_Q_INLINES_C
#define FMPZ_POLY_Q_INLINE FLINT_DLL
#else
#define FMPZ_POLY_Q_INLINE static __inline__
#endif

#undef ulong
#define ulong ulongxx /* interferes with system includes */
#include <stdlib.h>
#undef ulong
#include <gmp.h>
#define ulong mp_limb_t

#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

#ifdef __cplusplus
 extern "C" {
#endif

typedef struct
{
    fmpz_poly_struct *num;
    fmpz_poly_struct *den;
} fmpz_poly_q_struct;

typedef fmpz_poly_q_struct fmpz_poly_q_t[1];

/* Accessing numerator and denominator ***************************************/

#define fmpz_poly_q_numref(op)  ((op)->num)

#define fmpz_poly_q_denref(op)  ((op)->den)

FLINT_DLL void fmpz_poly_q_canonicalise(fmpz_poly_q_t rop);

FLINT_DLL int fmpz_poly_q_is_canonical(const fmpz_poly_q_t op);

/* Memory management *********************************************************/

FLINT_DLL void fmpz_poly_q_init(fmpz_poly_q_t rop);

FLINT_DLL void fmpz_poly_q_clear(fmpz_poly_q_t rop);

/* Randomisation *************************************************************/

FLINT_DLL void fmpz_poly_q_randtest(fmpz_poly_q_t poly, flint_rand_t state,
                          slong len1, mp_bitcnt_t bits1, 
                          slong len2, mp_bitcnt_t bits2);

FLINT_DLL void fmpz_poly_q_randtest_not_zero(fmpz_poly_q_t poly, flint_rand_t state, 
                                   slong len1, mp_bitcnt_t bits1, 
                                   slong len2, mp_bitcnt_t bits2);

/* Assignment ****************************************************************/

FLINT_DLL void fmpz_poly_q_set(fmpz_poly_q_t rop, const fmpz_poly_q_t op);

FLINT_DLL void fmpz_poly_q_set_si(fmpz_poly_q_t rop, slong op);

FLINT_DLL void fmpz_poly_q_swap(fmpz_poly_q_t op1, fmpz_poly_q_t op2);

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

FLINT_DLL void fmpz_poly_q_inv(fmpz_poly_q_t rop, const fmpz_poly_q_t op);

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

FLINT_DLL void fmpz_poly_q_add_in_place(fmpz_poly_q_t rop, const fmpz_poly_q_t op);
FLINT_DLL void fmpz_poly_q_sub_in_place(fmpz_poly_q_t rop, const fmpz_poly_q_t op);

FLINT_DLL void fmpz_poly_q_add(fmpz_poly_q_t rop, const fmpz_poly_q_t op1, const fmpz_poly_q_t op2);
FLINT_DLL void fmpz_poly_q_sub(fmpz_poly_q_t rop, const fmpz_poly_q_t op1, const fmpz_poly_q_t op2);

FLINT_DLL void fmpz_poly_q_addmul(fmpz_poly_q_t rop, const fmpz_poly_q_t op1, const fmpz_poly_q_t op2);
FLINT_DLL void fmpz_poly_q_submul(fmpz_poly_q_t rop, const fmpz_poly_q_t op1, const fmpz_poly_q_t op2);

/* Scalar multiplication and division ****************************************/

FLINT_DLL void fmpz_poly_q_scalar_mul_si(fmpz_poly_q_t rop, const fmpz_poly_q_t op, slong x);
FLINT_DLL void fmpz_poly_q_scalar_mul_mpz(fmpz_poly_q_t rop, const fmpz_poly_q_t op, const mpz_t x);
FLINT_DLL void fmpz_poly_q_scalar_mul_mpq(fmpz_poly_q_t rop, const fmpz_poly_q_t op, const mpq_t x);

FLINT_DLL void fmpz_poly_q_scalar_div_si(fmpz_poly_q_t rop, const fmpz_poly_q_t op, slong x);
FLINT_DLL void fmpz_poly_q_scalar_div_mpz(fmpz_poly_q_t rop, const fmpz_poly_q_t op, const mpz_t x);
FLINT_DLL void fmpz_poly_q_scalar_div_mpq(fmpz_poly_q_t rop, const fmpz_poly_q_t op, const mpq_t x);

/* Multiplication and division ***********************************************/

FLINT_DLL void fmpz_poly_q_mul(fmpz_poly_q_t rop, 
                     const fmpz_poly_q_t op1, const fmpz_poly_q_t op2);

FLINT_DLL void fmpz_poly_q_div(fmpz_poly_q_t rop, 
                     const fmpz_poly_q_t op1, const fmpz_poly_q_t op2);

/* Powering ******************************************************************/

FLINT_DLL void fmpz_poly_q_pow(fmpz_poly_q_t rop, const fmpz_poly_q_t op, ulong exp);

/* Derivative ****************************************************************/

FLINT_DLL void fmpz_poly_q_derivative(fmpz_poly_q_t rop, const fmpz_poly_q_t op);

/* Evaluation ****************************************************************/

FLINT_DLL int fmpz_poly_q_evaluate(mpq_t rop, const fmpz_poly_q_t f, const mpq_t a);

/* Input and output **********************************************************/

FLINT_DLL int fmpz_poly_q_set_str(fmpz_poly_q_t rop, const char *s);

char * fmpz_poly_q_get_str(const fmpz_poly_q_t op);
char * fmpz_poly_q_get_str_pretty(const fmpz_poly_q_t op, const char *x);

FLINT_DLL int fmpz_poly_q_print(const fmpz_poly_q_t op);
FLINT_DLL int fmpz_poly_q_print_pretty(const fmpz_poly_q_t op, const char *x);

#ifdef __cplusplus
}
#endif

#endif

