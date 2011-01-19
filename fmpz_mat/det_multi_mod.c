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

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "nmod_mat.h"
#include "nmod_vec.h"


void
fmpz_mat_det_multi_mod(fmpz_t det, const fmpz_mat_t A, int proved)
{
    fmpz_t prod, det_new;
    mp_limb_t prime;
    nmod_mat_t Amod;
    long i, j;

    fmpz_init(prod);
    fmpz_init(det_new);

    prime = n_nextprime(1UL<<(FLINT_BITS-1), proved);

    nmod_mat_init(Amod, A->r, A->c, prime);
    fmpz_set_ui(prod, prime);
    for (i = 0; i < A->r; i++)
        for (j = 0; j < A->c; j++)
            Amod->rows[i][j] = fmpz_fdiv_ui(A->rows[i]+j, prime);
    fmpz_set_ui(det, nmod_mat_det(Amod));

    /* TODO: support proved = 1, implementing the Hadamard bound */
    while (1)
    {
        /* XXX: add a method for changing nmod_mat modulus */
        nmod_t m;
        prime = n_nextprime(prime, proved);
        nmod_init(&m, prime);
        Amod->mod = m;

        for (i = 0; i < A->r; i++)
            for (j = 0; j < A->c; j++)
                Amod->rows[i][j] = fmpz_fdiv_ui(A->rows[i]+j, prime);

        fmpz_CRT_ui(det_new, det, prod, nmod_mat_det(Amod), prime);

        if (!proved && fmpz_equal(det_new, det))
            break;

        fmpz_mul_ui(prod, prod, prime);
        fmpz_set(det, det_new);
    }

    nmod_mat_clear(Amod);
    fmpz_clear(prod);
    fmpz_clear(det_new);
}
