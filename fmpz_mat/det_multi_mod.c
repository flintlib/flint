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
    mp_limb_t prime, detmod;
    nmod_mat_t Amod;
    long dim = A->r;

    if (dim < 1)
    {
        fmpz_set_ui(det, 1UL);
        return;
    }

    fmpz_init(prod);
    fmpz_init(det_new);

    /* Pick moduli allowing fast single-limb Gaussian elimination */
    prime = n_nextprime(n_sqrt((-1UL) / (dim+1)) * 0.98, proved);

    nmod_mat_init(Amod, A->r, A->c, prime);
    fmpz_set_ui(prod, prime);

    fmpz_mat_get_nmod_mat(Amod, A);
    fmpz_set_ui(det, _nmod_mat_det_rowreduce(Amod));

    /* TODO: support proved = 1, implementing the Hadamard bound */
    while (1)
    {
        prime = n_nextprime(prime, proved);
        _nmod_mat_set_mod(Amod, prime);

        fmpz_mat_get_nmod_mat(Amod, A);
        detmod = _nmod_mat_det_rowreduce(Amod);
        fmpz_CRT_ui(det_new, det, prod, detmod, prime);

        if (!proved && fmpz_equal(det_new, det))
            break;

        fmpz_mul_ui(prod, prod, prime);
        fmpz_set(det, det_new);
    }

    nmod_mat_clear(Amod);
    fmpz_clear(prod);
    fmpz_clear(det_new);
}
