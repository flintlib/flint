/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_theta.h"

/* Bolza to Mumford:
   0 -> {0,5}
   2 -> {2,5}
   4 -> {4,5}
   5 -> empty set
   Then: {i,j} -> S = {i+1,j+1} */

/* Mumford to characteristics:
   emptyset -> 0
   {1,2} -> 2
   {1,4} -> 9
   {1,6} -> 12
   {2,3} -> 8
   {2,5} -> 15
   {3,4} -> 3
   {3,6} -> 6
   {4,5} -> 4
   {5,6} -> 1 */

/* See Bolza, "Darstellung von Invarianten durch \theta-Functionen, p.493 */
static void
bolza_E(acb_t E, acb_srcptr th, slong prec)
{
    acb_ptr R;
    acb_ptr v;
    acb_ptr cmp;
    acb_t P;
    slong k;

    R = _acb_vec_init(15);
    v = _acb_vec_init(16);
    cmp = _acb_vec_init(15);
    acb_init(P);

    for (k = 0; k < 16; k++)
    {
        acb_pow_ui(&v[k], &th[k], 4, prec);
    }

    acb_sub(&R[0], &v[2], &v[6], prec);
    acb_sub(&cmp[0], &v[1], &v[9], prec);
    acb_sub(&R[1], &v[8], &v[12], prec);
    acb_sub(&cmp[1], &v[1], &v[3], prec);
    acb_sub(&R[2], &v[0], &v[4], prec);
    acb_add(&cmp[2], &v[9], &v[3], prec);
    acb_sub(&R[3], &v[4], &v[12], prec);
    acb_sub(&cmp[3], &v[2], &v[3], prec);
    acb_sub(&R[4], &v[0], &v[8], prec);
    acb_add(&cmp[4], &v[6], &v[3], prec);
    acb_sub(&R[5], &v[4], &v[6], prec);
    acb_sub(&cmp[5], &v[8], &v[9], prec);
    acb_sub(&R[6], &v[0], &v[2], prec);
    acb_add(&cmp[6], &v[12], &v[9], prec);
    acb_add(&R[7], &v[12], &v[6], prec);
    acb_sub(&cmp[7], &v[0], &v[1], prec);
    acb_sub(&R[8], &v[4], &v[2], prec);
    acb_sub(&cmp[8], &v[8], &v[1], prec);
    acb_add(&R[9], &v[8], &v[2], prec);
    acb_add(&cmp[9], &v[4], &v[1], prec);
    acb_sub(&R[10], &v[0], &v[6], prec);
    acb_add(&cmp[10], &v[12], &v[1], prec);
    acb_add(&R[11], &v[12], &v[2], prec);
    acb_sub(&cmp[11], &v[0], &v[9], prec);
    acb_sub(&R[12], &v[4], &v[8], prec);
    acb_sub(&cmp[12], &v[2], &v[1], prec);
    acb_add(&R[13], &v[6], &v[8], prec);
    acb_sub(&cmp[13], &v[0], &v[3], prec);
    acb_sub(&R[14], &v[0], &v[12], prec);
    acb_add(&cmp[14], &v[2], &v[9], prec);

    acb_one(P);
    for (k = 0; k < 16; k++)
    {
        if (acb_theta_char_is_even(k, 2))
        {
            acb_mul(P, P, &th[k], prec);
        }
    }
    acb_one(E);
    for (k = 0; k < 15; k++)
    {
        acb_mul(E, E, &R[k], prec);
    }
    acb_mul(E, E, P, prec); /* prod (theta) * prod(Ri) */

    _acb_vec_clear(R, 15);
    _acb_vec_clear(v, 16);
    _acb_vec_clear(cmp, 15);
    acb_clear(P);
}

/* See Igusa, "Modular forms and projective invariants" p. 848 */
void
acb_theta_g2_chi35(acb_t res, acb_srcptr th, slong prec)
{
    acb_t t;
    acb_init(t);

    bolza_E(t, th, prec);
    acb_neg(res, t);
    acb_mul_2exp_si(res, res, -37);

    acb_clear(t);
}
