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

    Copyright (C) 2011, 2012, 2013 Sebastian Pancratz

******************************************************************************/

#ifndef PADIC_MAT_H
#define PADIC_MAT_H

#ifdef PADIC_MAT_INLINES_C
#define PADIC_MAT_INLINE FLINT_DLL
#else
#define PADIC_MAT_INLINE static __inline__
#endif

#undef ulong
#define ulong ulongxx /* interferes with system includes */
#include <stdio.h>
#undef ulong
#include <gmp.h>
#define ulong mp_limb_t

#include "flint.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpq_mat.h"
#include "padic.h"

#ifdef __cplusplus
 extern "C" {
#endif


typedef struct
{
    fmpz_mat_struct mat;
    slong val;
    slong N;
} padic_mat_struct;

typedef padic_mat_struct padic_mat_t[1];

/* Macros  *******************************************************************/

PADIC_MAT_INLINE
fmpz_mat_struct * padic_mat(const padic_mat_t A) 
{
   return (fmpz_mat_struct *)(&(A->mat));
}

PADIC_MAT_INLINE
fmpz * padic_mat_entry(const padic_mat_t A, slong i, slong j)
{
   return A->mat.rows[i] + j;
}

#define padic_mat_val(A) ((A)->val)
#define padic_mat_prec(A) ((A)->N)

PADIC_MAT_INLINE
slong padic_mat_get_val(const padic_mat_t A)
{
   return A->val;
}

PADIC_MAT_INLINE
slong padic_mat_get_prec(const padic_mat_t A)
{
   return A->N;
}

PADIC_MAT_INLINE
slong padic_mat_nrows(const padic_mat_t A)
{
   return (A->mat).r;
}

PADIC_MAT_INLINE
slong padic_mat_ncols(const padic_mat_t A)
{
   return (A->mat).c;
}

/* Memory management  ********************************************************/

FLINT_DLL void padic_mat_init(padic_mat_t A, slong r, slong c);

FLINT_DLL void padic_mat_init2(padic_mat_t A, slong r, slong c, slong prec);

FLINT_DLL void padic_mat_clear(padic_mat_t A);

FLINT_DLL void _padic_mat_canonicalise(padic_mat_t A, const padic_ctx_t ctx);

FLINT_DLL void _padic_mat_reduce(padic_mat_t A, const padic_ctx_t ctx);

FLINT_DLL void padic_mat_reduce(padic_mat_t A, const padic_ctx_t ctx);

PADIC_MAT_INLINE int
padic_mat_is_empty(const padic_mat_t A)
{
    return fmpz_mat_is_empty(padic_mat(A));
}

PADIC_MAT_INLINE int
padic_mat_is_square(const padic_mat_t A)
{
    return fmpz_mat_is_square(padic_mat(A));
}

PADIC_MAT_INLINE int 
padic_mat_is_canonical(const padic_mat_t A, const padic_ctx_t ctx)
{
    if (fmpz_mat_is_zero(padic_mat(A)))
    {
        return (padic_mat_val(A) == 0);
    }
    else
    {
        slong i, j;
        int canonical = 0;

        for (i = 0; i < padic_mat(A)->r; i++)
            for (j = 0; j < padic_mat(A)->c; j++)
                if (!fmpz_divisible(padic_mat_entry(A, i, j), ctx->p))
                    canonical = 1;
        return canonical;
    }
}

PADIC_MAT_INLINE int 
padic_mat_is_reduced(const padic_mat_t A, const padic_ctx_t ctx)
{
    if (padic_mat_is_empty(A))
    {
        return 1;
    }
    else if (fmpz_mat_is_zero(padic_mat(A)))
    {
        return (padic_mat_val(A) == 0);
    }
    else if (padic_mat_is_canonical(A, ctx))
    {
        const slong v = padic_mat_val(A);
        const slong N = padic_mat_prec(A);

        if (v >= N)
        {
            return 0;
        }
        else
        {
            slong i, j;
            fmpz_t pN;
            int reduced = 1;
            int alloc = _padic_ctx_pow_ui(pN, N - v, ctx);

            for (i = 0; (i < padic_mat_nrows(A)) && reduced; i++)
                for (j = 0; (j < padic_mat_ncols(A)) && reduced; j++)
                    reduced = (fmpz_cmp(padic_mat_entry(A, i, j), pN) < 0);

            if (alloc)
                fmpz_clear(pN);

            return reduced;
        }
    }
    else
    {
        return 0;
    }
}

/* Basic assignment **********************************************************/

FLINT_DLL void padic_mat_set(padic_mat_t B, const padic_mat_t A, const padic_ctx_t ctx);

FLINT_DLL void padic_mat_swap(padic_mat_t A, padic_mat_t B);

FLINT_DLL void padic_mat_zero(padic_mat_t A);

FLINT_DLL void padic_mat_one(padic_mat_t A);

/* Conversions ***************************************************************/

FLINT_DLL void padic_mat_set_fmpq_mat(padic_mat_t B, 
                            const fmpq_mat_t A, const padic_ctx_t ctx);

FLINT_DLL void padic_mat_get_fmpq_mat(fmpq_mat_t B, 
                            const padic_mat_t A, const padic_ctx_t ctx);

/* Entries *******************************************************************/

FLINT_DLL void padic_mat_get_entry_padic(padic_t rop, 
                               const padic_mat_t op, slong i, slong j, 
                               const padic_ctx_t ctx);

FLINT_DLL void padic_mat_set_entry_padic(padic_mat_t rop, slong i, slong j, 
                               const padic_t op, const padic_ctx_t ctx);

/* Comparison ****************************************************************/

FLINT_DLL int padic_mat_equal(const padic_mat_t A, const padic_mat_t B);

FLINT_DLL int padic_mat_is_zero(const padic_mat_t A);

/* Input and output  *********************************************************/

FLINT_DLL int padic_mat_fprint(FILE * file, 
                     const padic_mat_t A, const padic_ctx_t ctx);

FLINT_DLL int padic_mat_fprint_pretty(FILE * file, const padic_mat_t A, 
                                         const padic_ctx_t ctx);

PADIC_MAT_INLINE
int padic_mat_print(const padic_mat_t A, const padic_ctx_t ctx)
{
    return padic_mat_fprint(stdout, A, ctx);
}

PADIC_MAT_INLINE
int padic_mat_print_pretty(const padic_mat_t A, const padic_ctx_t ctx)
{
    return padic_mat_fprint_pretty(stdout, A, ctx);
}

/* Random matrix generation  *************************************************/

FLINT_DLL void padic_mat_randtest(padic_mat_t mat, flint_rand_t state, 
                        const padic_ctx_t ctx);

/* Transpose *****************************************************************/

FLINT_DLL void padic_mat_transpose(padic_mat_t B, const padic_mat_t A);

/* Addition and subtraction **************************************************/

FLINT_DLL void _padic_mat_add(padic_mat_t C, const padic_mat_t A, const padic_mat_t B, 
                                   const padic_ctx_t ctx);
FLINT_DLL void padic_mat_add(padic_mat_t C, const padic_mat_t A, const padic_mat_t B, 
                                  const padic_ctx_t ctx);

FLINT_DLL void _padic_mat_sub(padic_mat_t C, const padic_mat_t A, const padic_mat_t B, 
                                   const padic_ctx_t ctx);
FLINT_DLL void padic_mat_sub(padic_mat_t C, const padic_mat_t A, const padic_mat_t B, 
                                  const padic_ctx_t ctx);

FLINT_DLL void _padic_mat_neg(padic_mat_t B, const padic_mat_t A);
FLINT_DLL void padic_mat_neg(padic_mat_t B, const padic_mat_t A, const padic_ctx_t ctx);

/* Scalar operations *********************************************************/

FLINT_DLL void _padic_mat_scalar_mul_padic(padic_mat_t B, 
                                 const padic_mat_t A, const padic_t c, 
                                 const padic_ctx_t ctx);
FLINT_DLL void padic_mat_scalar_mul_padic(padic_mat_t B, 
                                const padic_mat_t A, const padic_t c, 
                                const padic_ctx_t ctx);

FLINT_DLL void _padic_mat_scalar_mul_fmpz(padic_mat_t B, 
                                const padic_mat_t A, const fmpz_t c, 
                                const padic_ctx_t ctx);
FLINT_DLL void padic_mat_scalar_mul_fmpz(padic_mat_t B, 
                               const padic_mat_t A, const fmpz_t c, 
                               const padic_ctx_t ctx);

FLINT_DLL void padic_mat_scalar_div_fmpz(padic_mat_t B, 
                               const padic_mat_t A, const fmpz_t c, 
                               const padic_ctx_t ctx);

/* Multiplication ************************************************************/

FLINT_DLL void padic_mat_mul(padic_mat_t C, const padic_mat_t A, const padic_mat_t B, 
                                  const padic_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif

