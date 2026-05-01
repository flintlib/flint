/*
    Copyright (C) 2025, Vincent Neiger, Éric Schost, Mael Hostettler
    Copyright (C) 2026, Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "impl.h"

/* specialized for moduli that support n_mulmod_shoup */
void _nmod_geometric_progression_evaluate_init_nonfullword(nmod_geometric_progression_t G,
                                                           ulong r, slong len, nmod_t mod,
                                                           ulong q, ulong q_pr_quo, ulong q_pr_rem,
                                                           ulong inv_r,
                                                           ulong inv_q, ulong inv_q_pr_quo, ulong inv_q_pr_rem)
{
    /* G->ev_f = sum_{0 <= i < 2*len - 1} q**(i*i/2) * x**i */
    /* G->ev_s[i] = 1 / q**(i*i/2)                          */
    nmod_poly_init2_preinv(G->ev_f, mod.n, mod.ninv, 2*len - 1);
    G->ev_s = _nmod_vec_init(len);

    G->ev_f->length = 2*len - 1;

    /* precomputations for Shoup multiplication */
    ulong r_pow_i2 = r;
    ulong r_pr_quo = n_mulmod_precomp_shoup(r_pow_i2, mod.n);

    G->ev_f->coeffs[0] = 1;
    for (slong i = 1; i < 2*len - 1; i++)
    {
        G->ev_f->coeffs[i] = n_mulmod_shoup(r_pow_i2, G->ev_f->coeffs[i - 1], r_pr_quo, mod.n);
        n_mulmod_and_precomp_shoup(&r_pow_i2, &r_pr_quo, q, r_pow_i2, q_pr_quo, q_pr_rem, r_pr_quo, mod.n);
    }

    r_pr_quo = n_mulmod_precomp_shoup(inv_r, mod.n);
    G->ev_s[0] = 1;
    for (slong i = 1; i < len; i++)
    {
        G->ev_s[i] = n_mulmod_shoup(inv_r, G->ev_s[i - 1], r_pr_quo, mod.n);
        n_mulmod_and_precomp_shoup(&inv_r, &r_pr_quo, inv_q, inv_r, inv_q_pr_quo, inv_q_pr_rem, r_pr_quo, mod.n);
    }
}

/* general variant */
void _nmod_geometric_progression_evaluate_init(nmod_geometric_progression_t G,
                                               ulong r, slong len, nmod_t mod,
                                               ulong q, ulong inv_r, ulong inv_q)
{
    /* G->ev_f = sum_{0 <= i < 2*len - 1} q**(i*i/2) * x**i */
    /* G->ev_s[i] = 1 / q**(i*i/2)                          */
    nmod_poly_init2_preinv(G->ev_f, mod.n, mod.ninv, 2*len - 1);
    G->ev_s = _nmod_vec_init(len);

    G->ev_f->length = 2*len - 1;
    G->ev_f->coeffs[0] = 1;
    for (slong i = 1; i < 2*len - 1; i++)
    {
        G->ev_f->coeffs[i] = nmod_mul(G->ev_f->coeffs[i - 1], r, mod);
        r = nmod_mul(r, q, mod);  /* r**(2*i+1) */
    }

    G->ev_s[0] = 1;
    for (slong i = 1; i < len; i++)
    {
        G->ev_s[i] = nmod_mul(G->ev_s[i - 1], inv_r, mod);
        inv_r = nmod_mul(inv_r, inv_q, mod);
    }
}

void _nmod_geometric_progression_interpolate_init_nonfullword(nmod_geometric_progression_t G,
                                                              slong len, nmod_t mod,
                                                              ulong q, ulong q_pr_quo, ulong q_pr_rem,
                                                              ulong inv_q, ulong inv_q_pr_quo, ulong inv_q_pr_rem)
{
    /* quantities for Newton interpolation/evaluation/change-of-basis */
    /* see [Bostan - Schost, J.Complexity 2005, Section 5.1]          */
    /* write u_i for prod_{1 <= k <= i} (q**k - 1),                   */
    /*   and q_i for q**(i*(i-1)/2)                                   */

    /* coeff(G->int_f1, i) = (-1)**i * q_i / u_i */
    /* coeff(G->int_f2, i) = q_i / u_i, i.e. G->int_f2 == G->int_f1(-x) */
    nmod_poly_init2_preinv(G->int_f1, mod.n, mod.ninv, len);
    nmod_poly_init2_preinv(G->int_f2, mod.n, mod.ninv, len);
    G->int_f1->length = len;
    G->int_f2->length = len;

    /* G->int_s2[i] = (-1)**i * u_i / q_i = inverse of coeff(G->int_f1, i) */
    /* G->int_s3[i] = (-1)**i / u_i */
    /* G->int_s1[i] = 1 / u_i */
    G->int_s1 = _nmod_vec_init(len);
    G->int_s2 = _nmod_vec_init(len);
    G->int_s3 = _nmod_vec_init(len);

    G->int_f1->coeffs[0] = 1;
    G->int_f2->coeffs[0] = 1;
    G->int_s1[0] = 1;
    G->int_s2[0] = 1;
    G->int_s3[0] = 1;

    ulong q_pow_i = 1;
    ulong inv_q_pow_i = 1;
    ulong inv_q_i = 1;
    ulong prod_diff = 1;

    /* TODO improve with Shoup */
    for (slong i = 1; i < len; i++)
    {
        inv_q_i = nmod_mul(inv_q_i, inv_q_pow_i, mod);      /* 1 / q_i */
        inv_q_pow_i = nmod_mul(inv_q_pow_i, inv_q, mod);    /* 1 / q**i */
        G->int_f2->coeffs[i] = nmod_mul(G->int_f2->coeffs[i-1], q_pow_i, mod);  /* q_i */
        q_pow_i = nmod_mul(q_pow_i, q, mod);                /* q**i */
        G->int_s3[i] = q_pow_i - 1;                         /* temporarily, G->int_s3[i] = q**i - 1 */
        prod_diff = nmod_mul(q_pow_i - 1, prod_diff, mod);  /* u_i */
        if (i % 2)  /* i is odd */
            G->int_s2[i] = nmod_mul(mod.n - prod_diff, inv_q_i, mod);  /* (-1)**i * u_i / q_i */
        else  /* i is even */
            G->int_s2[i] = nmod_mul(prod_diff, inv_q_i, mod);          /* (-1)**i * u_i / q_i */
    }

    G->int_s1[len-1] = nmod_inv(prod_diff, mod);  /* 1 / u_{len-1} */
    for (slong i = len - 1; i > 0; i--)
    {
        ulong w_i = G->int_s1[i];                           /* 1 / u_i */
        ulong tmp = nmod_mul(G->int_f2->coeffs[i], w_i, mod);
        G->int_f2->coeffs[i] = tmp;                         /* q_i / u_i */
        G->int_s1[i-1] = nmod_mul(G->int_s3[i], w_i, mod);  /* 1 / u_{i-1} */
        if (i % 2)  /* i is odd */
        {
            G->int_s3[i] = mod.n - w_i;                  /* (-1)**i * G->int_s1[i] */
            G->int_f1->coeffs[i] = mod.n - tmp;          /* (-1)**i * G->int_f2[i] */
        }
        else  /* i is even */
        {
            G->int_s3[i] = w_i;                          /* (-1)**i * G->int_s1[i] */
            G->int_f1->coeffs[i] = tmp;                  /* (-1)**i * G->int_f2[i] */
        }
    }
}

void _nmod_geometric_progression_interpolate_init(nmod_geometric_progression_t G,
                                                  slong len, nmod_t mod,
                                                  ulong q, ulong inv_q)
{
    /* quantities for Newton interpolation/evaluation/change-of-basis */
    /* see [Bostan - Schost, J.Complexity 2005, Section 5.1]          */
    /* write u_i for prod_{1 <= k <= i} (q**k - 1),                   */
    /*   and q_i for q**(i*(i-1)/2)                                   */

    /* coeff(G->int_f1, i) = (-1)**i * q_i / u_i */
    /* coeff(G->int_f2, i) = q_i / u_i, i.e. G->int_f2 == G->int_f1(-x) */
    nmod_poly_init2_preinv(G->int_f1, mod.n, mod.ninv, len);
    nmod_poly_init2_preinv(G->int_f2, mod.n, mod.ninv, len);
    G->int_f1->length = len;
    G->int_f2->length = len;

    /* G->int_s2[i] = (-1)**i * u_i / q_i = inverse of coeff(G->int_f1, i) */
    /* G->int_s3[i] = (-1)**i / u_i */
    /* G->int_s1[i] = 1 / u_i */
    G->int_s1 = _nmod_vec_init(len);
    G->int_s2 = _nmod_vec_init(len);
    G->int_s3 = _nmod_vec_init(len);

    G->int_f1->coeffs[0] = 1;
    G->int_f2->coeffs[0] = 1;
    G->int_s1[0] = 1;
    G->int_s2[0] = 1;
    G->int_s3[0] = 1;

    ulong q_pow_i = 1;
    ulong inv_q_pow_i = 1;
    ulong inv_q_i = 1;
    ulong prod_diff = 1;

    /* TODO improve with Shoup */
    for (slong i = 1; i < len; i++)
    {
        inv_q_i = nmod_mul(inv_q_i, inv_q_pow_i, mod);      /* 1 / q_i */
        inv_q_pow_i = nmod_mul(inv_q_pow_i, inv_q, mod);    /* 1 / q**i */
        G->int_f2->coeffs[i] = nmod_mul(G->int_f2->coeffs[i-1], q_pow_i, mod);  /* q_i */
        q_pow_i = nmod_mul(q_pow_i, q, mod);                /* q**i */
        G->int_s3[i] = q_pow_i - 1;                         /* temporarily, G->int_s3[i] = q**i - 1 */
        prod_diff = nmod_mul(q_pow_i - 1, prod_diff, mod);  /* u_i */
        if (i % 2)  /* i is odd */
            G->int_s2[i] = nmod_mul(mod.n - prod_diff, inv_q_i, mod);  /* (-1)**i * u_i / q_i */
        else  /* i is even */
            G->int_s2[i] = nmod_mul(prod_diff, inv_q_i, mod);          /* (-1)**i * u_i / q_i */
    }

    G->int_s1[len-1] = nmod_inv(prod_diff, mod);  /* 1 / u_{len-1} */
    for (slong i = len - 1; i > 0; i--)
    {
        ulong w_i = G->int_s1[i];                           /* 1 / u_i */
        ulong tmp = nmod_mul(G->int_f2->coeffs[i], w_i, mod);
        G->int_f2->coeffs[i] = tmp;                         /* q_i / u_i */
        G->int_s1[i-1] = nmod_mul(G->int_s3[i], w_i, mod);  /* 1 / u_{i-1} */
        if (i % 2)  /* i is odd */
        {
            G->int_s3[i] = mod.n - w_i;                  /* (-1)**i * G->int_s1[i] */
            G->int_f1->coeffs[i] = mod.n - tmp;          /* (-1)**i * G->int_f2[i] */
        }
        else  /* i is even */
        {
            G->int_s3[i] = w_i;                          /* (-1)**i * G->int_s1[i] */
            G->int_f1->coeffs[i] = tmp;                  /* (-1)**i * G->int_f2[i] */
        }
    }
}

void nmod_geometric_progression_init(nmod_geometric_progression_t G,
                                     ulong r, slong len, nmod_t mod)
{
    G->len = len;
    G->mod = mod;

    if (NMOD_CAN_USE_SHOUP(mod))
    {
        const ulong q = nmod_mul(r, r, mod);
        const ulong inv_r = nmod_inv(r, mod);
        const ulong inv_q = nmod_mul(inv_r, inv_r, mod);

        ulong q_pr_rem;
        ulong q_pr_quo;
        ulong inv_q_pr_rem;
        ulong inv_q_pr_quo;
        n_mulmod_precomp_shoup_quo_rem(&q_pr_quo, &q_pr_rem, q, mod.n);
        n_mulmod_precomp_shoup_quo_rem(&inv_q_pr_quo, &inv_q_pr_rem, inv_q, mod.n);

        _nmod_geometric_progression_evaluate_init_nonfullword(G, r, len, mod, q, q_pr_quo, q_pr_rem, inv_r, inv_q, inv_q_pr_quo, inv_q_pr_rem);
        _nmod_geometric_progression_interpolate_init_nonfullword(G, len, mod, q, q_pr_quo, q_pr_rem, inv_q, inv_q_pr_quo, inv_q_pr_rem);
    }
    else
    {
        const ulong q = nmod_mul(r, r, mod);
        const ulong inv_r = nmod_inv(r, mod);
        const ulong inv_q = nmod_mul(inv_r, inv_r, mod);

        _nmod_geometric_progression_evaluate_init(G, r, len, mod, q, inv_r, inv_q);
        _nmod_geometric_progression_interpolate_init(G, len, mod, q, inv_q);
    }
}

void nmod_geometric_progression_clear(nmod_geometric_progression_t G)
{
    _nmod_vec_clear(G->ev_s);
    nmod_poly_clear(G->ev_f);
    _nmod_vec_clear(G->int_s1);
    _nmod_vec_clear(G->int_s2);
    _nmod_vec_clear(G->int_s3);
    nmod_poly_clear(G->int_f1);
    nmod_poly_clear(G->int_f2);
}
