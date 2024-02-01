/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* NOTE: This module is only used for fmpz_lll and should probably be
 * deprecated. Corresponding source file can be found in fmpz_lll/mpf-impl.c. */

#ifndef FLINT_MPF_H
#define FLINT_MPF_H

#include "flint.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef __mpf_struct mpf;

typedef struct
{
    mpf * entries;
    slong r;
    slong c;
    flint_bitcnt_t prec;
    mpf ** rows;
} mpf_mat_struct;

typedef mpf_mat_struct mpf_mat_t[1];

mpf * _mpf_vec_init(slong len, flint_bitcnt_t prec);
void _mpf_vec_clear(mpf * vec, slong len);
void _mpf_vec_set_fmpz_vec(mpf * appv, const fmpz * vec, slong len);
int _mpf_vec_dot2(mpf_t res, const mpf * vec1, const mpf * vec2, slong len2, flint_bitcnt_t prec);
void _mpf_vec_norm2(mpf_t res, const mpf * vec, slong len, flint_bitcnt_t prec);
void _mpf_vec_norm(mpf_t res, const mpf * vec, slong len);

FLINT_FORCE_INLINE
mpf * mpf_mat_entry(const mpf_mat_t mat, slong i, slong j)
{
   return mat->rows[i] + j;
}
void mpf_mat_init(mpf_mat_t mat, slong rows, slong cols, flint_bitcnt_t prec);
void mpf_mat_clear(mpf_mat_t mat);

#ifdef __cplusplus
}
#endif

#endif
