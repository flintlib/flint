/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

static void
sp2gz_randtest_trig(fmpz_mat_t mat, flint_rand_t state, slong bits)
{
    slong g = sp2gz_dim(mat);
    fmpz_mat_t b, bt;

    fmpz_mat_init(b, g, g);
    fmpz_mat_init(bt, g, g);
    bits = FLINT_MAX(bits, 1);

    fmpz_mat_randbits(b, state, bits);
    fmpz_mat_transpose(bt, b);
    fmpz_mat_add(b, b, bt);
    fmpz_mat_scalar_tdiv_q_2exp(b, b, 1);
    sp2gz_trig(mat, b);

    fmpz_mat_clear(b);
    fmpz_mat_clear(bt);
}

static void
sp2gz_randtest_block_diag(fmpz_mat_t mat, flint_rand_t state, slong bits)
{
    slong g = sp2gz_dim(mat);
    fmpz_mat_t u;

    fmpz_mat_init(u, g, g);
    bits = FLINT_MAX(bits, 1);

    fmpz_mat_one(u);
    fmpz_mat_randops(u, state, 2 * bits * g);
    sp2gz_block_diag(mat, u);

    fmpz_mat_clear(u);
}

void
sp2gz_randtest(fmpz_mat_t mat, flint_rand_t state, slong bits)
{
    slong g = sp2gz_dim(mat);
    slong b = bits/5 + 1;
    fmpz_mat_t aux;

    fmpz_mat_init(aux, 2 * g, 2 * g);

    sp2gz_randtest_block_diag(mat, state, b);
    sp2gz_randtest_trig(aux, state, b);
    fmpz_mat_mul(mat, mat, aux);
    sp2gz_j(aux);
    fmpz_mat_mul(mat, mat, aux);
    sp2gz_randtest_trig(aux, state, b);
    fmpz_mat_mul(mat, mat, aux);
    sp2gz_j(aux);
    fmpz_mat_mul(mat, mat, aux);
    sp2gz_randtest_trig(aux, state, b);
    fmpz_mat_mul(mat, mat, aux);
    sp2gz_j(aux);
    fmpz_mat_mul(mat, mat, aux);
    sp2gz_randtest_trig(aux, state, b);
    fmpz_mat_mul(mat, mat, aux);

    fmpz_mat_clear(aux);
}
