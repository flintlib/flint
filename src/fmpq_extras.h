/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPQ_EXTRAS_H
#define FMPQ_EXTRAS_H

#include "flint.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    fmpz_t _11, _12, _21, _22;
    int det;    /* 0,1,or,-1: 0 -> don't know, 1 -> 1, -1 -> -1 */
} _fmpz_mat22_struct;

typedef _fmpz_mat22_struct _fmpz_mat22_t[1];

typedef struct
{
    mp_limb_t _11, _12, _21, _22;
    int det;    /* ditto */
} _ui_mat22_struct;

typedef _ui_mat22_struct _ui_mat22_t[1];

/* ball for closed interval [left, right] */

typedef struct
{
    fmpz_t left_num, left_den, right_num, right_den;
    int exact;
} _fmpq_ball_struct;

typedef _fmpq_ball_struct _fmpq_ball_t[1];

typedef struct
{
    fmpz * array;
    slong length;
    slong alloc;
    slong limit;
    fmpz_t alt_sum;
    int want_alt_sum;
} _fmpq_cfrac_list_struct;

typedef _fmpq_cfrac_list_struct _fmpq_cfrac_list_t[1];

/* 2x2 matrices **************************************************************/

void _fmpz_mat22_init(_fmpz_mat22_t M);
void _fmpz_mat22_clear(_fmpz_mat22_t M);

void _fmpz_mat22_one(_fmpz_mat22_t M);
int _fmpz_mat22_is_one(_fmpz_mat22_t M);

flint_bitcnt_t _fmpz_mat22_bits(const _fmpz_mat22_t N);

void _fmpz_mat22_rmul(_fmpz_mat22_t M, const _fmpz_mat22_t N);

void _fmpz_mat22_addmul_inv_vec(fmpz_t ya, fmpz_t yb, _fmpz_mat22_t N, fmpz_t xa, fmpz_t xb);

void _fmpz_mat22_addmul_inv_mat(fmpz_t A11, fmpz_t A12, fmpz_t A21, fmpz_t A22, _fmpz_mat22_t M, fmpz_t B11, fmpz_t B12, fmpz_t B21, fmpz_t B22);

void _fmpz_mat22_rmul_ui(_fmpz_mat22_t M, const _ui_mat22_t N);
void _fmpz_mat22_rmul_inv_ui(_fmpz_mat22_t M, const _ui_mat22_t N);

void _fmpz_mat22_rmul_elem(_fmpz_mat22_t M, const fmpz_t q);
void _fmpz_mat22_rmul_inv_elem(_fmpz_mat22_t M, const fmpz_t q);

void _fmpz_mat22_lmul_elem(_fmpz_mat22_t M, const fmpz_t q);

/* Continued fractions *******************************************************/

void _fmpq_cfrac_list_init(_fmpq_cfrac_list_t v);
void _fmpq_cfrac_list_clear(_fmpq_cfrac_list_t v);

void _fmpq_cfrac_list_fit_length(_fmpq_cfrac_list_t v, slong len);

void _fmpq_cfrac_list_push_back(_fmpq_cfrac_list_t v, const fmpz_t a);
void _fmpq_cfrac_list_push_back_zero(_fmpq_cfrac_list_t v);

void _fmpq_cfrac_list_append_ui(_fmpq_cfrac_list_t v, const ulong * a, slong n);

void _fmpq_cfrac_list_swap(_fmpq_cfrac_list_t a, _fmpq_cfrac_list_t b);

/* Balls *********************************************************************/

void _fmpq_ball_init(_fmpq_ball_t x);
void _fmpq_ball_clear(_fmpq_ball_t x);

void _fmpq_ball_swap(_fmpq_ball_t x, _fmpq_ball_t y);

int _fmpq_ball_gt_one(const _fmpq_ball_t x);

void _fmpq_hgcd(_fmpq_cfrac_list_t s, _fmpz_mat22_t M, fmpz_t x_num, fmpz_t x_den);

void _fmpq_ball_get_cfrac(_fmpq_cfrac_list_t s, _fmpz_mat22_t M, int needM, _fmpq_ball_t x);

#ifdef __cplusplus
}
#endif

#endif
