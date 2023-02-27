/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef D_MAT_H
#define D_MAT_H

#ifdef D_MAT_INLINES_C
#define D_MAT_INLINE FLINT_DLL
#else
#define D_MAT_INLINE static __inline__
#endif

#undef ulong
#define ulong ulongxx /* interferes with system includes */
#include <stdio.h>
#undef ulong

#include <gmp.h>
#define ulong mp_limb_t
#include "flint.h"
#include "d_vec.h"

#ifdef __cplusplus
 extern "C" {
#endif

typedef struct
{
    double * entries;
    slong r;
    slong c;
    double ** rows;
} d_mat_struct;

typedef d_mat_struct d_mat_t[1];

#define d_mat_entry(mat,i,j) (*((mat)->rows[i] + (j)))

D_MAT_INLINE
double * d_mat_entry_ptr(const d_mat_t mat, slong i, slong j)
{
   return mat->rows[i] + j;
}

D_MAT_INLINE
double d_mat_get_entry(const d_mat_t mat, slong i, slong j)
{
   return mat->rows[i][j];
}

D_MAT_INLINE
slong d_mat_nrows(const d_mat_t mat)
{
    return mat->r;
}

D_MAT_INLINE
slong d_mat_ncols(const d_mat_t mat)
{
    return mat->c;
}

/* Memory management  ********************************************************/

FLINT_DLL void d_mat_init(d_mat_t mat, slong rows, slong cols);

FLINT_DLL void d_mat_swap(d_mat_t mat1, d_mat_t mat2);

D_MAT_INLINE void
d_mat_swap_entrywise(d_mat_t mat1, d_mat_t mat2)
{
    slong i, j;
    for (i = 0; i < d_mat_nrows(mat1); i++)
    {
       double * row1 = mat1->rows[i];
       double * row2 = mat2->rows[i];
       for (j = 0; j < d_mat_ncols(mat1); j++)
          DOUBLE_SWAP(row1[j], row2[j]);
    }
}

FLINT_DLL void d_mat_set(d_mat_t mat1, const d_mat_t mat2);

FLINT_DLL void d_mat_clear(d_mat_t mat);

FLINT_DLL int d_mat_equal(const d_mat_t mat1, const d_mat_t mat2);

FLINT_DLL int d_mat_approx_equal(const d_mat_t mat1, const d_mat_t mat2, double eps);

FLINT_DLL int d_mat_is_zero(const d_mat_t mat);

FLINT_DLL int d_mat_is_approx_zero(const d_mat_t mat, double eps);

D_MAT_INLINE
int d_mat_is_empty(const d_mat_t mat)
{
    return (mat->r == 0) || (mat->c == 0);
}

D_MAT_INLINE
int d_mat_is_square(const d_mat_t mat)
{
    return (mat->r == mat->c);
}

FLINT_DLL void d_mat_zero(d_mat_t mat);

FLINT_DLL void d_mat_one(d_mat_t mat);


/* Input and output  *********************************************************/

FLINT_DLL void d_mat_print(const d_mat_t mat);

/* Random matrix generation  *************************************************/

FLINT_DLL void d_mat_randtest(d_mat_t mat, flint_rand_t state, slong minexp,
                    slong maxexp);

/* Transpose */

FLINT_DLL void d_mat_transpose(d_mat_t B, const d_mat_t A);

/* Multiplication */

FLINT_DLL void d_mat_mul_classical(d_mat_t C, const d_mat_t A, const d_mat_t B);

/* Permutations */

D_MAT_INLINE
void d_mat_swap_rows(d_mat_t mat, slong r, slong s)
{
    if (r != s)
    {
        double * u;

        u = mat->rows[s];
        mat->rows[s] = mat->rows[r];
        mat->rows[r] = u; 
    }
}

/* Gram-Schmidt Orthogonalisation and QR Decomposition  ********************************************************/

FLINT_DLL void d_mat_gso(d_mat_t B, const d_mat_t A);

FLINT_DLL void d_mat_qr(d_mat_t Q, d_mat_t R, const d_mat_t A);

#ifdef __cplusplus
}
#endif

#endif

