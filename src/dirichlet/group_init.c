/*
    Copyright (C) 2015 Jonathan Bober
    Copyright (C) 2016 Fredrik Johansson
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "dirichlet.h"

/* see https://arxiv.org/abs/2206.14193 */
static ulong
primitive_root_p_and_p2(ulong p)
{
    if (p == 40487)
        return 10;

#if FLINT_BITS == 64
    if (p == UWORD(6692367337))
        return 7;

    if (p > UWORD(10000000000000000))
        flint_throw(FLINT_ERROR, "p > 10^16 not implemented in %s\n", __func__);
#endif

    return n_primitive_root_prime(p);
}

void
dirichlet_prime_group_init(dirichlet_prime_group_struct * P, ulong p, int e)
{
    P->p = p;
    P->e = e;
    if (p == 2 || p == 4)
    {
        P->p = 2;
        nmod_init(&P->pe, UWORD(1) << e);
        if (p == 2)
        {
            P->e = 2;
            nmod_init(&P->phi, 2);
            P->g = (UWORD(1) << e) - 1;
        }
        else
        {
            nmod_init(&P->phi, UWORD(1) << (e - 2));
            P->g = 5;
        }
    }
    else
    {
        ulong pe1;
        pe1 = n_pow(p, e - 1);
        nmod_init(&P->phi, (p-1) * pe1);
        nmod_init(&P->pe, p * pe1);
        P->g = primitive_root_p_and_p2(p);
    }
    P->dlog = NULL;
}

static void
dirichlet_group_lift_generators(dirichlet_group_t G)
{
    slong k;
    dirichlet_prime_group_struct * P = G->P;

    G->expo = G->phi_q = 1;
    if (G->neven)
    {
        G->phi_q = G->q_even / 2;
        G->expo = P[G->neven - 1].phi.n;
    }
    for (k = G->neven; k < G->num; k++)
    {
        G->phi_q *= P[k].phi.n;
        G->expo *= P[k].phi.n / n_gcd(G->expo, P[k].p - 1);
    }

    for (k = 0; k < G->num; k++)
    {
        nmod_t pe;
        ulong qpe, v;
        G->PHI[k] = G->expo / G->P[k].phi.n;
        /* lift generators mod q */
        /* u * p^e + v * q/p^e = 1 -> g mod q = 1 + (g-1) * v*(q/p^e) */
        pe = G->P[k].pe;
        qpe = G->q / pe.n;
        if (G->q < G->P[k].pe.n)
        {
            flint_printf("lift generator %wu from %wu to %wu e=%wu\n",
                G->P[k].g, G->P[k].pe.n, G->q, G->P[k].e);
        }
        v = nmod_inv(qpe % pe.n, pe);
        /* v * qpe < q */
        G->generators[k] = (1 + nmod_mul(G->P[k].g-1, v * qpe, G->mod)) % G->q;
    }
}

int
dirichlet_group_init(dirichlet_group_t G, ulong q)
{
    slong k;
    ulong e2;
    n_factor_t fac;

    G->q = q;
    nmod_init(&G->mod, q);


    e2 = n_remove(&q, 2);
    G->q_even = UWORD(1) << e2;
    /* number of components at p=2 */
    G->neven = (e2 >= 3) ? 2 : (e2 == 2) ? 1 : 0;

    /* warning: only factor odd part */
    n_factor_init(&fac);
    n_factor(&fac, q, 1);

#if FLINT_BITS == 64
    for (k = 0; k < fac.num; k++)
        if (fac.p[k] > UWORD(10000000000000000))
            return 0;
#endif

    G->num = fac.num + G->neven;
    G->P = flint_malloc(G->num * sizeof(dirichlet_prime_group_struct));
    G->generators = flint_malloc(G->num * sizeof(ulong));
    G->PHI = flint_malloc(G->num * sizeof(ulong));

    /* even part */
    if (G->neven >= 1)
        dirichlet_prime_group_init(&G->P[0], 2, e2);
    if (G->neven == 2)
        dirichlet_prime_group_init(&G->P[1], 4, e2);

    /* odd part */
    G->rad_q = 1;
    for (k = G->neven; k < G->num; k++)
    {
        ulong p, e;
        p = fac.p[k - G->neven];
        e = fac.exp[k - G->neven];
        G->rad_q *= p;
        dirichlet_prime_group_init(&G->P[k], p, e);
    }
    dirichlet_group_lift_generators(G);

    return 1;
}

void
dirichlet_subgroup_init(dirichlet_group_t H, const dirichlet_group_t G, ulong h)
{
    int s[15];    /* selected components */
    slong k, j, e2;

    H->q = h;
    nmod_init(&H->mod, h);

    /* even components */

    e2 = n_remove(&h, 2);
    H->q_even = UWORD(1) << e2;

    s[0] = s[1] = 0;
    j = 0;
    if (e2 >= 2)
        s[j++] = 2;
    if (e2 >= 3)
        s[j++] = e2;

    H->neven = j;

    /* odd components */
    for (k = G->neven; k < G->num; k++)
    {
        ulong p = G->P[k].p;

        s[k] = n_remove(&h, p);

        if (s[k] > 0)
        {
            j++;
            H->rad_q *= p;
        }
    }

    H->num = j;
    H->P = flint_malloc(j * sizeof(dirichlet_prime_group_struct));
    H->generators = flint_malloc(j * sizeof(ulong));
    H->PHI = flint_malloc(j * sizeof(ulong));

    j = 0;
    for (k = 0; k < H->neven; k++)
    {
            H->P[j] = G->P[k];
            if (H->q_even < G->q_even)
            {
                nmod_init(&H->P[j].pe, H->q_even);
                H->P[j].e = s[k];
                if (k == 0)
                    H->P[j].g = H->q_even - 1;
                else
                    nmod_init(&H->P[j].phi, H->q_even / 4);
            }
            j++;
    }

    for (k = G->neven; k < G->num; k++)
    {
        if (s[k])
        {
            H->P[j] = G->P[k];
            if (s[k] < G->P[k].e)
            {
                ulong pf, p = H->P[j].p;
                H->P[j].e = s[k];
                pf = n_pow(p, s[k]);
                nmod_init(&H->P[j].pe, pf);
                nmod_init(&H->P[j].phi, (p-1) * pf / p);
            }
            j++;
        }
    }

    dirichlet_group_lift_generators(H);
}
