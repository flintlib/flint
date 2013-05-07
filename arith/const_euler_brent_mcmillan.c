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

#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq.h"
#include "arith.h"

void
mpfr_const_euler_brent_mcmillan(mpfr_t res, mpfr_rnd_t rnd)
{
    fmpq_bsplit_t sum;
    fmpq *ab, *cd, *pq;
    mpfr_t S0, K0, I0, t, u, v;
    long bits, wp, n, nterms1, nterms2, k;

    bits = mpfr_get_prec(res) + 20;
    n = 0.08665 * bits + 1;
    nterms1 = 4.9706258 * n + 1;
    nterms2 = 2 * n + 1;
    wp = bits + FLINT_BIT_COUNT(n);

    fmpq_bsplit_init(sum);

    /* Compute S0, I0 */
    ab = _fmpq_vec_init(nterms1);
    cd = _fmpq_vec_init(nterms1);
    pq = _fmpq_vec_init(nterms1);

    for (k = 0; k < nterms1; k++) fmpq_set_si(ab + k, 1, 1);
    for (k = 0; k < nterms1; k++) fmpq_set_si(cd + k, 1, k + 1);
    for (k = 0; k < nterms1; k++)
    {
        fmpq_set_si(pq + k, n, k + 1);
        fmpq_mul(pq + k, pq + k, pq + k);
    }

    fmpq_bsplit_sum_abcdpq(sum, ab, cd, pq, 0, nterms1);

    mpfr_init2(S0, wp);
    mpfr_init2(K0, wp);
    mpfr_init2(I0, wp);
    mpfr_init2(t, wp);
    mpfr_init2(u, wp);
    mpfr_init2(v, wp);

    fmpq_bsplit_get_mpfr(S0, sum);

    fmpz_zero(sum->D);  /* hack */
    fmpq_bsplit_get_mpfr(I0, sum);
    mpfr_add_ui(I0, I0, 1, MPFR_RNDN);

    fmpq_bsplit_clear(sum);
    fmpq_bsplit_init(sum);

    /* Compute K0 */
    fmpq_set_si(pq + 0, 1, 1);

    for (k = 1; k < nterms2; k++)
    {
        fmpz_set_si(fmpq_numref(pq + k), 1 - 2*k);
        fmpz_pow_ui(fmpq_numref(pq + k), fmpq_numref(pq + k), 3);
        fmpz_neg(fmpq_numref(pq + k), fmpq_numref(pq + k));
        fmpz_set_ui(fmpq_denref(pq + k), 32 * k);
        fmpz_mul_ui(fmpq_denref(pq + k), fmpq_denref(pq + k), n);
        fmpz_mul_ui(fmpq_denref(pq + k), fmpq_denref(pq + k), n);
    }

    fmpq_bsplit_sum_pq(sum, pq, 0, nterms2);
    fmpq_bsplit_get_mpfr(K0, sum);
    mpfr_div_ui(K0, K0, 4 * n, MPFR_RNDN);

    mpfr_mul(S0, S0, I0, MPFR_RNDN);
    mpfr_sub(S0, S0, K0, MPFR_RNDN);
    mpfr_mul(I0, I0, I0, MPFR_RNDN);
    mpfr_div(t, S0, I0, MPFR_RNDN);

    mpfr_set_ui(u, n, MPFR_RNDN);
    mpfr_log(u, u, MPFR_RNDN);
    mpfr_sub(res, t, u, rnd);

    mpfr_clear(S0);
    mpfr_clear(K0);
    mpfr_clear(I0);
    mpfr_clear(t);
    mpfr_clear(u);
    mpfr_clear(v);

    fmpq_bsplit_clear(sum);
    _fmpq_vec_clear(ab, nterms1);
    _fmpq_vec_clear(cd, nterms1);
    _fmpq_vec_clear(pq, nterms1);
}
