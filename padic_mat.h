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

    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#ifndef PADIC_MAT_H
#define PADIC_MAT_H

#undef ulong /* interferes with system includes */
#include <stdio.h>
#define ulong unsigned long

#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpq_mat.h"
#include "padic.h"

typedef struct
{
    fmpz_mat_struct mat;
    long val;
} padic_mat_struct;

typedef padic_mat_struct padic_mat_t[1];

/* Macros  *******************************************************************/

#define padic_mat(A)             (&((A)->mat))
#define padic_mat_unit(A, i, j)  ((A)->mat.rows[i] + (j))
#define padic_mat_val(A)         ((A)->val)

#define padic_mat_nrows(A)       (((A)->mat).r)
#define padic_mat_ncols(A)       (((A)->mat).c)

/* Memory management  ********************************************************/

void padic_mat_init(padic_mat_t A, long r, long c);
void padic_mat_clear(padic_mat_t A);

void _padic_mat_canonicalise(padic_mat_t A, const padic_ctx_t ctx);
void _padic_mat_reduce(padic_mat_t A, const padic_ctx_t ctx);
void padic_mat_reduce(padic_mat_t A, const padic_ctx_t ctx);

static __inline__ int
padic_mat_is_empty(const padic_mat_t A)
{
    return fmpz_mat_is_empty(padic_mat(A));
}

static __inline__ int
padic_mat_is_square(const padic_mat_t A)
{
    return fmpz_mat_is_square(padic_mat(A));
}

static __inline__ int 
_padic_mat_is_canonical(const padic_mat_t A, const fmpz_t p)
{
    if (fmpz_mat_is_zero(padic_mat(A)))
    {
        return (padic_mat_val(A) == 0);
    }
    else
    {
        long i, j;
        int canonical = 0;

        for (i = 0; i < padic_mat(A)->r; i++)
            for (j = 0; j < padic_mat(A)->c; j++)
                if (!fmpz_divisible(padic_mat_unit(A, i, j), p))
                    canonical = 1;

        return canonical;
    }
}

/* Basic assignment **********************************************************/

void padic_mat_set(padic_mat_t B, const padic_mat_t A);

void padic_mat_swap(padic_mat_t A, padic_mat_t B);

void padic_mat_zero(padic_mat_t A);

void _padic_mat_one(padic_mat_t A);

void padic_mat_one(padic_mat_t A, const padic_ctx_t ctx);

/* Conversions ***************************************************************/

void padic_mat_set_fmpq_mat(padic_mat_t B, 
                            const fmpq_mat_t A, const padic_ctx_t ctx);

void padic_mat_get_fmpq_mat(fmpq_mat_t B, 
                            const padic_mat_t A, const padic_ctx_t ctx);

/* Entries *******************************************************************/

void padic_mat_get_entry_padic(padic_t rop, 
                               const padic_mat_t op, long i, long j, 
                               const padic_ctx_t ctx);

void padic_mat_set_entry_padic(padic_mat_t rop, long i, long j, 
                               const padic_t op, const padic_ctx_t ctx);

/* Comparison ****************************************************************/

int padic_mat_equal(const padic_mat_t A, const padic_mat_t B);

int padic_mat_is_zero(const padic_mat_t A);

/* Input and output  *********************************************************/

int padic_mat_fprint(FILE * file, 
                     const padic_mat_t A, const padic_ctx_t ctx);

int padic_mat_fprint_pretty(FILE * file, const padic_mat_t A, 
                                         const padic_ctx_t ctx);

static __inline__
int padic_mat_print(const padic_mat_t A, const padic_ctx_t ctx)
{
    return padic_mat_fprint(stdout, A, ctx);
}

static __inline__
int padic_mat_print_pretty(const padic_mat_t A, const padic_ctx_t ctx)
{
    return padic_mat_fprint_pretty(stdout, A, ctx);
}

/* Random matrix generation  *************************************************/

void padic_mat_randtest(padic_mat_t mat, flint_rand_t state, 
                        const padic_ctx_t ctx);

/* Transpose *****************************************************************/

void padic_mat_transpose(padic_mat_t B, const padic_mat_t A);

/* Addition and subtraction **************************************************/

void _padic_mat_add(padic_mat_t C, const padic_mat_t A, const padic_mat_t B, 
                                   const padic_ctx_t ctx);
void padic_mat_add(padic_mat_t C, const padic_mat_t A, const padic_mat_t B, 
                                  const padic_ctx_t ctx);

void _padic_mat_sub(padic_mat_t C, const padic_mat_t A, const padic_mat_t B, 
                                   const padic_ctx_t ctx);
void padic_mat_sub(padic_mat_t C, const padic_mat_t A, const padic_mat_t B, 
                                  const padic_ctx_t ctx);

void _padic_mat_neg(padic_mat_t B, const padic_mat_t A);
void padic_mat_neg(padic_mat_t B, const padic_mat_t A, const padic_ctx_t ctx);

/* Scalar operations *********************************************************/

void _padic_mat_scalar_mul_padic(padic_mat_t B, 
                                 const padic_mat_t A, const padic_t c, 
                                 const padic_ctx_t ctx);
void padic_mat_scalar_mul_padic(padic_mat_t B, 
                                const padic_mat_t A, const padic_t c, 
                                const padic_ctx_t ctx);

void _padic_mat_scalar_mul_fmpz(padic_mat_t B, 
                                const padic_mat_t A, const fmpz_t c, 
                                const padic_ctx_t ctx);
void padic_mat_scalar_mul_fmpz(padic_mat_t B, 
                               const padic_mat_t A, const fmpz_t c, 
                               const padic_ctx_t ctx);

void padic_mat_scalar_div_fmpz(padic_mat_t B, 
                               const padic_mat_t A, const fmpz_t c, 
                               const padic_ctx_t ctx);

/* Multiplication ************************************************************/

void _padic_mat_mul(padic_mat_t C, const padic_mat_t A, const padic_mat_t B, 
                                   const padic_ctx_t ctx);
void padic_mat_mul(padic_mat_t C, const padic_mat_t A, const padic_mat_t B, 
                                  const padic_ctx_t ctx);

#endif

