/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ACB_MODULAR_H
#define ACB_MODULAR_H

#ifdef ACB_MODULAR_INLINES_C
#define ACB_MODULAR_INLINE
#else
#define ACB_MODULAR_INLINE static inline
#endif

#include "fmpz.h"
#include "acb_types.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    fmpz a;
    fmpz b;
    fmpz c;
    fmpz d;
}
psl2z_struct;

typedef psl2z_struct psl2z_t[1];

ACB_MODULAR_INLINE void
psl2z_init(psl2z_t g)
{
    fmpz_init(&g->a);
    fmpz_init(&g->b);
    fmpz_init(&g->c);
    fmpz_init(&g->d);
    fmpz_one(&g->a);
    fmpz_one(&g->d);
}

ACB_MODULAR_INLINE void
psl2z_clear(psl2z_t g)
{
    fmpz_clear(&g->a);
    fmpz_clear(&g->b);
    fmpz_clear(&g->c);
    fmpz_clear(&g->d);
}

ACB_MODULAR_INLINE void
psl2z_swap(psl2z_t f, psl2z_t g)
{
    psl2z_struct h = *f;
    *f = *g;
    *g = h;
}

ACB_MODULAR_INLINE void
psl2z_set(psl2z_t h, const psl2z_t g)
{
    fmpz_set(&h->a, &g->a);
    fmpz_set(&h->b, &g->b);
    fmpz_set(&h->c, &g->c);
    fmpz_set(&h->d, &g->d);
}

ACB_MODULAR_INLINE void
psl2z_one(psl2z_t g)
{
    fmpz_one(&g->a);
    fmpz_zero(&g->b);
    fmpz_zero(&g->c);
    fmpz_one(&g->d);
}

#ifdef FLINT_HAVE_FILE
void psl2z_fprint(FILE * file, const psl2z_t g);
#endif

void psl2z_print(const psl2z_t g);

ACB_MODULAR_INLINE int
psl2z_equal(const psl2z_t f, const psl2z_t g)
{
    return fmpz_equal(&f->a, &g->a)
        && fmpz_equal(&f->b, &g->b)
        && fmpz_equal(&f->c, &g->c)
        && fmpz_equal(&f->d, &g->d);
}

void psl2z_mul(psl2z_t h, const psl2z_t f, const psl2z_t g);

void psl2z_inv(psl2z_t h, const psl2z_t g);

int psl2z_is_one(const psl2z_t g);

int psl2z_is_correct(const psl2z_t g);

void psl2z_randtest(psl2z_t g, flint_rand_t state, slong bits);

void acb_modular_transform(acb_t w, const psl2z_t g, const acb_t z, slong prec);

void acb_modular_fundamental_domain_approx_d(psl2z_t g,
    double x, double y, double one_minus_eps);

void acb_modular_fundamental_domain_approx_arf(psl2z_t g,
    const arf_t xx, const arf_t yy, const arf_t one_minus_eps, slong prec);

void acb_modular_fundamental_domain_approx(acb_t w, psl2z_t g, const acb_t z,
        const arf_t one_minus_eps, slong prec);

int acb_modular_is_in_fundamental_domain(const acb_t z, const arf_t tol, slong prec);

void acb_modular_addseq_theta(slong * exponents, slong * aindex, slong * bindex, slong num);

void acb_modular_addseq_eta(slong * exponents, slong * aindex, slong * bindex, slong num);

void acb_modular_fill_addseq(slong * tab, slong len);

void acb_modular_theta_transform(int * R, int * S, int * C, const psl2z_t g);

void acb_modular_theta_const_sum(acb_t theta2, acb_t theta3, acb_t theta4,
    const acb_t q, slong prec);

void acb_modular_theta_const_sum_basecase(acb_t theta2, acb_t theta3, acb_t theta4,
    const acb_t q, slong N, slong prec);

void acb_modular_theta_const_sum_rs(acb_t theta2, acb_t theta3, acb_t theta4,
    const acb_t q, slong N, slong prec);

void acb_modular_theta_sum(acb_ptr theta1, acb_ptr theta2,
        acb_ptr theta3, acb_ptr theta4,
    const acb_t w, int w_is_unit, const acb_t q, slong len, slong prec);

void acb_modular_theta_notransform(acb_t theta1, acb_t theta2,
    acb_t theta3, acb_t theta4, const acb_t z, const acb_t tau,
    slong prec);

void acb_modular_theta(acb_t theta1, acb_t theta2,
    acb_t theta3, acb_t theta4, const acb_t z, const acb_t tau,
    slong prec);

void acb_modular_theta_jet_notransform(acb_ptr theta1, acb_ptr theta2,
    acb_ptr theta3, acb_ptr theta4, const acb_t z, const acb_t tau,
    slong len, slong prec);

void acb_modular_theta_jet(acb_ptr theta1, acb_ptr theta2,
    acb_ptr theta3, acb_ptr theta4, const acb_t z, const acb_t tau,
    slong len, slong prec);

void _acb_modular_theta_series(acb_ptr theta1, acb_ptr theta2, acb_ptr theta3, acb_ptr theta4,
    acb_srcptr z, slong zlen, const acb_t tau, slong len, slong prec);

void acb_modular_theta_series(acb_poly_t theta1, acb_poly_t theta2,
    acb_poly_t theta3, acb_poly_t theta4, const acb_poly_t z, const acb_t tau,
        slong len, slong prec);

void acb_modular_j(acb_t z, const acb_t tau, slong prec);

int acb_modular_epsilon_arg(const psl2z_t g);

void acb_modular_eta_sum(acb_t eta, const acb_t q, slong prec);

void acb_modular_eta(acb_t z, const acb_t tau, slong prec);

void acb_modular_lambda(acb_t r, const acb_t tau, slong prec);

void acb_modular_delta(acb_t r, const acb_t tau, slong prec);

void acb_modular_eisenstein(acb_ptr r, const acb_t tau, slong len, slong prec);

void acb_modular_elliptic_p(acb_t r, const acb_t z, const acb_t tau, slong prec);

void acb_modular_elliptic_p_zpx(acb_ptr r, const acb_t z, const acb_t tau, slong len, slong prec);

void acb_modular_elliptic_k(acb_t k, const acb_t m, slong prec);

void acb_modular_elliptic_k_cpx(acb_ptr w, const acb_t m, slong len, slong prec);

void acb_modular_elliptic_e(acb_t res, const acb_t m, slong prec);

void acb_modular_hilbert_class_poly(fmpz_poly_t res, slong D);

/* this is a performance hack until the main arb/acb functions improve */
void _acb_modular_mul(acb_t z, acb_t tmp1, acb_t tmp2, const acb_t x, const acb_t y, slong wprec, slong prec);

#ifdef __cplusplus
}
#endif

#endif

