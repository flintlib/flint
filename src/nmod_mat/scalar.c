/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2014 Ashish Kedia
    Copyright (C) 2023 Albin Ahlbäck
    Copyright (C) 2024 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint-mparam.h"
#include "ulong_extras.h"
#include "nmod_mat.h"
#include "nmod.h"
#include "fmpz.h"
#include "nmod_vec.h"

/*
    Windows (stride above the number of columns) go row by row through the
    vector functions from NMOD_VEC_SIMD_MIN_LEN columns where a kernel applies
    to n, from NMOD_MAT_SCALAR_ROWS_FULLWORD_MIN_LEN columns at 64 bits, and
    otherwise from NMOD_MAT_SCALAR_ROWS_MIN_LEN columns. Without a kernel,
    below 64 bits, each row computes Shoup's constant again, a division,
    which cost more than the row loop saved at 8-12 columns (63 bits with
    AVX2: 0.64-0.88x the loops on the entries); from 15 columns, 0.98-1.26x.
    At 64 bits the rows need no such constant (Montgomery's reduction from 8
    entries in scalar_mul): 1.03-1.73x the loops at 8-15 columns, except
    0.89-0.95x at 8 columns with clang on Intel. Nov 2026.
*/
#define NMOD_MAT_SCALAR_ROWS_MIN_LEN 16
#define NMOD_MAT_SCALAR_ROWS_FULLWORD_MIN_LEN 8

/* whether windows of len columns go row by row (add: scalar_addmul) */
#define NMOD_MAT_SCALAR_ROWS(len, mod, add)                               \
    ((len) >= NMOD_MAT_SCALAR_ROWS_MIN_LEN                                \
     || ((len) >= NMOD_VEC_SIMD_MIN_LEN                                   \
         && NMOD_BITS(mod) <= NMOD_VEC_SIMD_SCALAR_MAX_BITS(add))         \
     || ((len) >= NMOD_MAT_SCALAR_ROWS_FULLWORD_MIN_LEN                   \
         && NMOD_BITS(mod) == FLINT_BITS))

/*******************/
/*  scalar_addmul  */
/*******************/

void
_nmod_mat_scalar_addmul_ui_generic(nmod_mat_t C, const nmod_mat_t A,
                                          const nmod_mat_t B, const ulong c)
{
    for (slong i = 0; i < A->r; i++)
    {
        for (slong j = 0; j < A->c; j++)
        {
            nmod_mat_entry(C, i, j) =
                n_addmod(nmod_mat_entry(A, i, j),
                         nmod_mul(nmod_mat_entry(B, i, j), c, B->mod),
                         A->mod.n);
        }
    }
}

void
_nmod_mat_scalar_addmul_ui_precomp(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B,
                                    const ulong c, const ulong c_pr)
{
    for (slong i = 0; i < A->r; i++)
    {
        for (slong j = 0; j < A->c; j++)
        {
            nmod_mat_entry(C, i, j) =
                n_addmod(nmod_mat_entry(A, i, j),
                         n_mulmod_shoup(c, nmod_mat_entry(B, i, j),
                                          c_pr, B->mod.n),
                         A->mod.n);
        }
    }
}

void
nmod_mat_scalar_addmul_ui(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B,
                          ulong c)
{
    c = nmod_set_ui(c, C->mod);

    if (c == UWORD(0))
        nmod_mat_set(C, A);
    else if (c == UWORD(1))
        nmod_mat_add(C, A, B);
    else if (c == A->mod.n - UWORD(1))
        nmod_mat_sub(C, A, B);
    else
    {
        /* C = A, then C += c B (C = c B, then C += A, when C is B only): in
           one call on all the entries when the rows follow each other
           (stride = number of columns), else row by row (see above); below,
           the loops on the entries (calls per row of 4 to 6 entries:
           0.41-0.98x, Nov 2026) */
        slong nrows = A->r, len = A->c;
        int add = (C != B || C == A);

        if (A->stride == len && B->stride == len && C->stride == len)
        {
            len *= nrows;
            nrows = (len != 0);
        }
        else if (!NMOD_MAT_SCALAR_ROWS(len, A->mod, add))
        {
            if (A->r * A->c >= FLINT_MULMOD_SHOUP_THRESHOLD && A->mod.norm > 0)
                _nmod_mat_scalar_addmul_ui_precomp(C, A, B, c,
                                        _nmod_mulmod_precomp_shoup(c, A->mod));
            else
                _nmod_mat_scalar_addmul_ui_generic(C, A, B, c);
            return;
        }

        for (slong i = 0; i < nrows; i++)
        {
            nn_ptr Ci = nmod_mat_entry_ptr(C, i, 0);
            nn_srcptr Ai = nmod_mat_entry_ptr(A, i, 0);
            nn_srcptr Bi = nmod_mat_entry_ptr(B, i, 0);

            if (!add)
            {
                _nmod_vec_scalar_mul_nmod(Ci, Bi, len, c, A->mod);
                _nmod_vec_add(Ci, Ci, Ai, len, A->mod);
            }
            else
            {
                if (C != A)
                    _nmod_vec_set(Ci, Ai, len);
                _nmod_vec_scalar_addmul_nmod(Ci, Bi, len, c, A->mod);
            }
        }
    }
}

/****************/
/*  scalar_mul  */
/****************/

void
_nmod_mat_scalar_mul_generic(nmod_mat_t B, const nmod_mat_t A, ulong c)
{
    for (slong i = 0; i < A->r; i++)
        for (slong j = 0; j < A->c; j++)
            nmod_mat_entry(B, i, j) = nmod_mul(nmod_mat_entry(A, i, j), c, A->mod);
}

void
_nmod_mat_scalar_mul_precomp(nmod_mat_t B, const nmod_mat_t A, ulong c, ulong c_pr)
{
    for (slong i = 0; i < A->r; i++)
        for (slong j = 0; j < A->c; j++)
            nmod_mat_entry(B, i, j) = n_mulmod_shoup(
                c, nmod_mat_entry(A, i, j), c_pr, A->mod.n);
}

void
nmod_mat_scalar_mul(nmod_mat_t B, const nmod_mat_t A, ulong c)
{
    if (c == UWORD(0))
        nmod_mat_zero(B);
    else if (c == UWORD(1))
        nmod_mat_set(B, A);
    else if (c == A->mod.n - UWORD(1))
        nmod_mat_neg(B, A);
    else if (A->stride == A->c && B->stride == A->c)
    {
        /* in one call on all the entries when the rows follow each other */
        if (A->r * A->c != 0)
            _nmod_vec_scalar_mul_nmod(B->entries, A->entries, A->r * A->c,
                                      c, A->mod);
    }
    else if (NMOD_MAT_SCALAR_ROWS(A->c, A->mod, 0))
    {
        /* windows: one call per row (see above; with a kernel, from 8 x 8,
           1.07-4.81x on average) */
        for (slong i = 0; i < A->r; i++)
            _nmod_vec_scalar_mul_nmod(nmod_mat_entry_ptr(B, i, 0),
                                      nmod_mat_entry_ptr(A, i, 0), A->c,
                                      c, A->mod);
    }
    /* below, the loops on the entries (calls per row: 0.83-0.98x at 4 x 4
       and 6 x 6, Nov 2026) */
    else if (A->r * A->c >= FLINT_MULMOD_SHOUP_THRESHOLD && A->mod.norm > 0)
    {
        const ulong c_pr = _nmod_mulmod_precomp_shoup(c, A->mod);
        _nmod_mat_scalar_mul_precomp(B, A, c, c_pr);
    }
    else
        _nmod_mat_scalar_mul_generic(B, A, c);
}


void
nmod_mat_scalar_mul_fmpz(nmod_mat_t B, const nmod_mat_t A, const fmpz_t c)
{
    nmod_mat_scalar_mul(B, A, fmpz_get_nmod(c, B->mod));
}
