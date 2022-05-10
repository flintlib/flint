/*
    Copyright (C) 2011, 2012, 2013 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef PADIC_MAT_H
#define PADIC_MAT_H

#ifdef PADIC_MAT_INLINES_C
#define PADIC_MAT_INLINE FLINT_DLL
#else
#define PADIC_MAT_INLINE static __inline__
#endif

#include "fmpz_mini.h"

#ifdef __cplusplus
 extern "C" {
#endif


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

#define FMPZ_MAT_IS_EMPTY(mat) ((mat->r == 0) || (mat->c == 0))
#define FMPZ_MAT_IS_SQUARE(mat) (mat->r == mat->c)

PADIC_MAT_INLINE int
padic_mat_is_empty(const padic_mat_t A)
{
    return FMPZ_MAT_IS_EMPTY(padic_mat(A));
}

PADIC_MAT_INLINE int
padic_mat_is_square(const padic_mat_t A)
{
    return FMPZ_MAT_IS_SQUARE(padic_mat(A));
}

FLINT_DLL int padic_mat_is_canonical(const padic_mat_t A, const padic_ctx_t ctx);

FLINT_DLL int padic_mat_is_reduced(const padic_mat_t A, const padic_ctx_t ctx);

/* Basic assignment **********************************************************/

FLINT_DLL void padic_mat_set(padic_mat_t B, const padic_mat_t A, const padic_ctx_t ctx);

FLINT_DLL void padic_mat_swap(padic_mat_t A, padic_mat_t B);

PADIC_MAT_INLINE void
padic_mat_swap_entrywise(padic_mat_t mat1, padic_mat_t mat2)
{
    slong i, j;

    for (i = 0; i < padic_mat_nrows(mat1); i++)
        for (j = 0; j < padic_mat_ncols(mat1); j++)
            fmpz_swap(padic_mat_entry(mat2, i, j), padic_mat_entry(mat1, i, j));
}

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

#if defined (H_STDIO)               \
  || defined (_H_STDIO)             \
  || defined (_STDIO_H)             \
  || defined (_STDIO_H_)            \
  || defined (__STDIO_H)            \
  || defined (__STDIO_H__)          \
  || defined (_STDIO_INCLUDED)      \
  || defined (__dj_include_stdio_h_)\
  || defined (_FILE_DEFINED)        \
  || defined (__STDIO__)            \
  || defined (_MSL_STDIO_H)         \
  || defined (_STDIO_H_INCLUDED)    \
  || defined (_ISO_STDIO_ISO_H)     \
  || defined (__STDIO_LOADED)       \
  || defined (_STDIO)

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

#endif

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
