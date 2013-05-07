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

    Copyright (C) 2010 William Hart

******************************************************************************/

#ifndef MFPR_MAT_H
#define MPFR_MAT_H

#include <gmp.h>
#include <mpfr.h> 

#ifdef __cplusplus
 extern "C" {
#endif

typedef struct
{
    __mpfr_struct * entries;
    long r;
    long c;
    mp_bitcnt_t prec;
    __mpfr_struct ** rows;
} mpfr_mat_struct;

/* fmpz_mat_t allows reference-like semantics for fmpz_mat_struct */
typedef mpfr_mat_struct mpfr_mat_t[1];

void mpfr_mat_init(mpfr_mat_t mat, long rows, long cols, mpfr_prec_t prec);

void mpfr_mat_clear(mpfr_mat_t mat);

#ifdef __cplusplus
}
#endif

#endif

