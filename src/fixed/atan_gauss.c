/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "mpn_extras.h"
#include "arb.h"
#include "fixed.h"

/* Thread-local cache of the angles theta_0 = pi/2, theta_j =
   2 arg(pi_j) of the first _fixed_atan_gauss_num nonreal Gaussian
   primes, for the diophantine (multi-prime) trigonometric reduction
   (sin_cos_diophantine.c), laid out exactly like the logarithm cache
   of log_primes.c: _fixed_atan_gauss_n limbs per entry, a unit limb
   on top (theta < pi < 4), the lowest fraction limb a guard, each
   entry the exact floor of theta B^(n-1) so that shorter reads are
   exact floors too.  The values come from _fixed_atan_gauss_vec_fball
   (machin_bsplit.c: the Machin-type sets of machin_tab.c by binary
   splitting in fball arithmetic, single arctangents of small quotients
   between neighbouring primes beyond, arb's static table of the first
   13 at small precisions). */

/* real and imaginary parts of the first 64 nonreal Gaussian primes */
const signed char _fixed_gaussian_primes[2 * FIXED_ATAN_GAUSS_MAX] = {
    1, 1, 1, 2, 2, 3, 1, 4, 2, 5, 1, 6, 4, 5, 2, 7, 5, 6, 3, 8, 5, 8, 4, 9,
    1, 10, 3, 10, 7, 8, 4, 11, 7, 10, 6, 11, 2, 13, 9, 10, 7, 12, 1, 14, 2, 15, 8, 13,
    4, 15, 1, 16, 10, 13, 9, 14, 5, 16, 2, 17, 12, 13, 11, 14, 9, 16, 5, 18, 8, 17, 7, 18,
    10, 17, 6, 19, 1, 20, 3, 20, 14, 15, 12, 17, 7, 20, 4, 21, 10, 19, 5, 22, 11, 20, 10, 21,
    14, 19, 13, 20, 1, 24, 8, 23, 5, 24, 17, 18, 16, 19, 4, 25, 13, 22, 6, 25, 12, 23, 1, 26,
    5, 26, 15, 22, 2, 27, 9, 26
};

FLINT_TLS_PREFIX nn_ptr _fixed_atan_gauss = NULL;
FLINT_TLS_PREFIX slong _fixed_atan_gauss_n = 0;
FLINT_TLS_PREFIX slong _fixed_atan_gauss_num = 0;
static FLINT_TLS_PREFIX int _fixed_atan_gauss_cleanup_registered = 0;

static void
_fixed_atan_gauss_cleanup(void)
{
    flint_free(_fixed_atan_gauss);
    _fixed_atan_gauss = NULL;
    _fixed_atan_gauss_n = 0;
    _fixed_atan_gauss_num = 0;
    _fixed_atan_gauss_cleanup_registered = 0;
}

/* res_j = theta_j = 2 arg(pi_j), j < num, as arb balls at prec bits */
void
_fixed_atan_gauss_vec(arb_ptr res, slong num, slong prec)
{
    fball_struct * v;
    slong j;

    v = flint_malloc(num * sizeof(fball_struct));
    for (j = 0; j < num; j++)
        fball_init(v + j);

    _fixed_atan_gauss_vec_fball(v, num, (prec + FLINT_BITS - 1) / FLINT_BITS + 2);

    for (j = 0; j < num; j++)
    {
        fball_get_arb(res + j, v + j);
        arb_set_round(res + j, res + j, prec);
        fball_clear(v + j);
    }
    flint_free(v);
}

void
_fixed_atan_gauss_ensure(slong num, slong nv)
{
    slong nc, j, guard;
    fball_struct * v;

    FLINT_ASSERT(num >= 1 && num <= FIXED_ATAN_GAUSS_MAX);

    if (nv + 2 <= _fixed_atan_gauss_n && num <= _fixed_atan_gauss_num)
        return;

    nc = FLINT_MAX(nv + 2, _fixed_atan_gauss_n);
    num = FLINT_MAX(num, _fixed_atan_gauss_num);

    flint_free(_fixed_atan_gauss);
    _fixed_atan_gauss = flint_malloc(num * nc * sizeof(ulong));
    _fixed_atan_gauss_n = nc;
    _fixed_atan_gauss_num = num;

    v = flint_malloc(num * sizeof(fball_struct));
    for (j = 0; j < num; j++)
        fball_init(v + j);

    /* undetermined floors (astronomically rare) are retried at higher
       precision */
    for (guard = 2; ; guard *= 2)
    {
        _fixed_atan_gauss_vec_fball(v, num, nc + guard);
        if (_fixed_store_floors(_fixed_atan_gauss, nc, v, num))
            break;
    }

    for (j = 0; j < num; j++)
        fball_clear(v + j);
    flint_free(v);

    if (!_fixed_atan_gauss_cleanup_registered)
    {
        flint_register_cleanup_function(_fixed_atan_gauss_cleanup);
        _fixed_atan_gauss_cleanup_registered = 1;
    }
}

nn_srcptr
_fixed_atan_gauss_entry(slong j, slong nv)
{
    FLINT_ASSERT(j >= 0 && j < _fixed_atan_gauss_num);
    FLINT_ASSERT(nv >= 1 && nv + 2 <= _fixed_atan_gauss_n);
    return _fixed_atan_gauss + j * _fixed_atan_gauss_n
        + (_fixed_atan_gauss_n - (nv + 1));
}

void
_fixed_atan_gauss_clear(void)
{
    flint_free(_fixed_atan_gauss);
    _fixed_atan_gauss = NULL;
    _fixed_atan_gauss_n = 0;
    _fixed_atan_gauss_num = 0;
}
