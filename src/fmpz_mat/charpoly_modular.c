/*
    Copyright (C) 2012, 2013 Sebastian Pancratz
    Copyrigth (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "nmod_mat.h"
#include "nmod_poly.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "fmpz_poly.h"
#include "thread_pool.h"
#include "thread_support.h"

typedef struct
{
    const fmpz_mat_struct * op;
    const nmod_mat_struct * op_mod;
    nn_ptr primes;
    nn_ptr coeffs;
} _charpoly_mod_arg_t;

static void _charpoly_mod_worker(slong i, void * arg_ptr)
{
    _charpoly_mod_arg_t *arg = (_charpoly_mod_arg_t *) arg_ptr;

    slong n = (arg->op)->r;

    nmod_poly_t poly;
    poly->coeffs = arg->coeffs + i * (n + 1);
    poly->length = poly->alloc = n + 1;

    if (arg->op_mod == NULL)
    {
        nmod_mat_t mat;
        nmod_mat_init(mat, n, n, arg->primes[i]);
        fmpz_mat_get_nmod_mat(mat, arg->op);
        nmod_mat_charpoly(poly, mat);
        nmod_mat_clear(mat);
    }
    else
    {
        nmod_mat_charpoly(poly, arg->op_mod + i);
    }
}

void _fmpz_mat_charpoly_modular(fmpz * rop, const fmpz_mat_t op)
{
    const slong n = op->r;

    if (n == 0)
    {
        fmpz_one(rop);
        return;
    }

    slong bound;
    slong pbits = FLINT_BITS - 1;
    ulong p = (UWORD(1) << pbits);

    fmpz_t b;
    fmpz_init(b);
    fmpz_mat_charpoly_bound(b, op);
    bound = fmpz_bits(b);
    fmpz_clear(b);

    /* Account for signs */
    bound++;

    slong num_primes = (bound + pbits - 1) / pbits;
    ulong *primes = flint_malloc(num_primes * sizeof(ulong));
    for (slong i = 0; i < num_primes; i++)
    {
        p = n_nextprime(p, 0);
        primes[i] = p;
    }

    nn_ptr coeff_buf = flint_malloc(num_primes * (n + 1) * sizeof(ulong));

    _charpoly_mod_arg_t args;
    args.op = op;
    args.op_mod = NULL;
    args.primes = primes;
    args.coeffs = coeff_buf;

    /* Want asymptotically fast multimodular reduction */
    if (num_primes > FMPZ_MAT_MOD_PRIMES_COMB_CUTOFF)
    {
        nmod_mat_struct * Amod;

        Amod = flint_malloc(sizeof(nmod_mat_struct) * num_primes);
        for (slong i = 0; i < num_primes; i++)
            nmod_mat_init(Amod + i, n, n, primes[i]);

        fmpz_mat_multi_mod_ui((nmod_mat_t *) Amod, num_primes, op);
        args.op_mod = Amod;

        flint_parallel_do(_charpoly_mod_worker, &args, num_primes, 0, FLINT_PARALLEL_UNIFORM);

        for (slong i = 0; i < num_primes; i++)
            nmod_mat_clear(Amod + i);
        flint_free(Amod);
    }
    else
    {
        /* Just do the mod p reductions inline with the mod p charpolys */
        flint_parallel_do(_charpoly_mod_worker, &args, num_primes, 0, FLINT_PARALLEL_UNIFORM);
    }

    nn_srcptr *residues = flint_malloc(num_primes * sizeof(nn_srcptr));
    for (slong i = 0; i < num_primes; i++)
        residues[i] = coeff_buf + i * (n + 1);

    _fmpz_vec_multi_CRT_ui(rop, residues, n + 1, primes, num_primes, 1);

    flint_free(residues);
    flint_free(coeff_buf);
    flint_free(primes);
}

void fmpz_mat_charpoly_modular(fmpz_poly_t cp, const fmpz_mat_t mat)
{
     fmpz_poly_fit_length(cp, mat->r + 1);
    _fmpz_poly_set_length(cp, mat->r + 1);
    _fmpz_mat_charpoly_modular(cp->coeffs, mat);
}

