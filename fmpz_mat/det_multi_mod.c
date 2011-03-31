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
    fmpz_t prod, stable_prod, det_new, bound;
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
    fmpz_init(bound);
    fmpz_init(stable_prod);

    fmpz_mat_det_bound(bound, A);
    fmpz_mul_ui(bound, bound, 2UL);  /* signed */

    prime = _nmod_mat_fast_rowreduce_modulus(dim, dim, proved);

    /* 2 as a modulus would not let us distinguish between -1 and 1 */
    if (prime == 2UL)
        prime = 3UL;

    nmod_mat_init(Amod, A->r, A->c, prime);
    fmpz_set_ui(prod, prime);

    fmpz_mat_get_nmod_mat(Amod, A);
    detmod = _nmod_mat_det_rowreduce(Amod);

    fmpz_set_ui(det, detmod);
    /* May be signed */
    fmpz_sub_ui(det_new, det, prime);
    if (fmpz_cmpabs(det_new, det) < 0)
        fmpz_set(det, det_new);

    fmpz_set_ui(stable_prod, detmod == 0UL ? prime : 1UL);

    while (fmpz_cmp(prod, bound) < 0)
    {
        prime = n_nextprime(prime, proved);
        _nmod_mat_set_mod(Amod, prime);

        fmpz_mat_get_nmod_mat(Amod, A);
        detmod = _nmod_mat_det_rowreduce(Amod);
        fmpz_CRT_ui(det_new, det, prod, detmod, prime);

        if (fmpz_equal(det_new, det))
        {
            fmpz_mul_ui(stable_prod, stable_prod, prime);
            if (!proved && fmpz_bits(stable_prod) > 100)
                break;
        }
        else
        {
            fmpz_set_ui(stable_prod, prime);
        }

        /* printf("prime,bits %lu, %ld\n", prime, fmpz_bits(prod)); */

        fmpz_mul_ui(prod, prod, prime);
        fmpz_set(det, det_new);
    }

    nmod_mat_clear(Amod);
    fmpz_clear(prod);
    fmpz_clear(stable_prod);
    fmpz_clear(det_new);
    fmpz_clear(bound);
}
