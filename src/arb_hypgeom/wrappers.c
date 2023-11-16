/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "arb_hypgeom.h"
#include "acb_hypgeom.h"

void
arb_hypgeom_erfi(arb_t res, const arb_t z, slong prec)
{
    if (!arb_is_finite(z))
    {
        arb_indeterminate(res);
    }
    else
    {
        acb_t t;
        acb_init(t);
        arb_set(acb_realref(t), z);
        acb_hypgeom_erfi(t, t, prec);
        arb_swap(res, acb_realref(t));
        acb_clear(t);
    }
}

void
arb_hypgeom_fresnel(arb_t res1, arb_t res2, const arb_t z, int normalized, slong prec)
{
    if (!arb_is_finite(z))
    {
        if (res1 != NULL) arb_indeterminate(res1);
        if (res2 != NULL) arb_indeterminate(res2);
    }
    else
    {
        acb_t t, u;
        acb_init(t);
        acb_init(u);
        arb_set(acb_realref(t), z);
        acb_hypgeom_fresnel(res1 ? t : NULL, res2 ? u : NULL, t, normalized, prec);
        if (res1 != NULL) arb_swap(res1, acb_realref(t));
        if (res2 != NULL) arb_swap(res2, acb_realref(u));
        acb_clear(t);
        acb_clear(u);
    }
}

void
arb_hypgeom_ei(arb_t res, const arb_t z, slong prec)
{
    if (!arb_is_finite(z))
    {
        arb_indeterminate(res);
    }
    else
    {
        acb_t t;
        acb_init(t);
        arb_set(acb_realref(t), z);
        acb_hypgeom_ei(t, t, prec);
        arb_swap(res, acb_realref(t));
        acb_clear(t);
    }
}

void
arb_hypgeom_shi(arb_t res, const arb_t z, slong prec)
{
    if (!arb_is_finite(z))
    {
        arb_indeterminate(res);
    }
    else
    {
        acb_t t;
        acb_init(t);
        arb_set(acb_realref(t), z);
        acb_hypgeom_shi(t, t, prec);
        arb_swap(res, acb_realref(t));
        acb_clear(t);
    }
}

void
arb_hypgeom_chi(arb_t res, const arb_t z, slong prec)
{
    if (!arb_is_finite(z) || !arb_is_positive(z))
    {
        arb_indeterminate(res);
    }
    else
    {
        acb_t t;
        acb_init(t);
        arb_set(acb_realref(t), z);
        acb_hypgeom_chi(t, t, prec);
        arb_swap(res, acb_realref(t));
        acb_clear(t);
    }
}

void
arb_hypgeom_li(arb_t res, const arb_t z, int offset, slong prec)
{
    if (!arb_is_finite(z) || !arb_is_nonnegative(z))
    {
        arb_indeterminate(res);
    }
    else
    {
        acb_t t;
        acb_init(t);
        arb_set(acb_realref(t), z);
        acb_hypgeom_li(t, t, offset, prec);
        arb_swap(res, acb_realref(t));
        acb_clear(t);
    }
}

void
arb_hypgeom_0f1(arb_t res, const arb_t a, const arb_t z, int regularized, slong prec)
{
    acb_t t, u;
    acb_init(t);
    acb_init(u);
    arb_set(acb_realref(t), a);
    arb_set(acb_realref(u), z);
    acb_hypgeom_0f1(t, t, u, regularized, prec);
    if (acb_is_finite(t) && acb_is_real(t))
        arb_swap(res, acb_realref(t));
    else
        arb_indeterminate(res);
    acb_clear(t);
    acb_clear(u);
}

void
arb_hypgeom_m(arb_t res, const arb_t a, const arb_t b, const arb_t z, int regularized, slong prec)
{
    acb_t t, u, v;
    acb_init(t);
    acb_init(u);
    acb_init(v);
    arb_set(acb_realref(t), a);
    arb_set(acb_realref(u), b);
    arb_set(acb_realref(v), z);
    acb_hypgeom_m(t, t, u, v, regularized, prec);
    if (acb_is_finite(t) && acb_is_real(t))
        arb_swap(res, acb_realref(t));
    else
        arb_indeterminate(res);
    acb_clear(t);
    acb_clear(u);
    acb_clear(v);
}

void
arb_hypgeom_1f1(arb_t res, const arb_t a, const arb_t b, const arb_t z, int regularized, slong prec)
{
    arb_hypgeom_m(res, a, b, z, regularized, prec);
}

void
arb_hypgeom_u(arb_t res, const arb_t a, const arb_t b, const arb_t z, slong prec)
{
    acb_t t, u, v;
    acb_init(t);
    acb_init(u);
    acb_init(v);
    arb_set(acb_realref(t), a);
    arb_set(acb_realref(u), b);
    arb_set(acb_realref(v), z);
    acb_hypgeom_u(t, t, u, v, prec);
    if (acb_is_finite(t) && acb_is_real(t))
        arb_swap(res, acb_realref(t));
    else
        arb_indeterminate(res);
    acb_clear(t);
    acb_clear(u);
    acb_clear(v);
}

void
arb_hypgeom_2f1(arb_t res, const arb_t a, const arb_t b, const arb_t c, const arb_t z, int regularized, slong prec)
{
    acb_t t, u, v, w;
    acb_init(t);
    acb_init(u);
    acb_init(v);
    acb_init(w);
    arb_set(acb_realref(t), a);
    arb_set(acb_realref(u), b);
    arb_set(acb_realref(v), c);
    arb_set(acb_realref(w), z);
    acb_hypgeom_2f1(t, t, u, v, w, regularized, prec);
    if (acb_is_finite(t) && acb_is_real(t))
        arb_swap(res, acb_realref(t));
    else
        arb_indeterminate(res);
    acb_clear(t);
    acb_clear(u);
    acb_clear(v);
    acb_clear(w);
}

void
arb_hypgeom_pfq(arb_t res, arb_srcptr a, slong p, arb_srcptr b, slong q, const arb_t z, int regularized, slong prec)
{
    acb_ptr t;
    slong i;
    t = _acb_vec_init(p + q + 1);
    for (i = 0; i < p; i++)
        arb_set(acb_realref(t + i), a + i);
    for (i = 0; i < q; i++)
        arb_set(acb_realref(t + p + i), b + i);
    arb_set(acb_realref(t + p + q), z);
    acb_hypgeom_pfq(t, t, p, t + p, q, t + p + q, regularized, prec);
    if (acb_is_finite(t) && acb_is_real(t))
        arb_swap(res, acb_realref(t));
    else
        arb_indeterminate(res);
    _acb_vec_clear(t, p + q + 1);
}

void
arb_hypgeom_bessel_j(arb_t res, const arb_t nu, const arb_t z, slong prec)
{
    acb_t t, u;
    acb_init(t);
    acb_init(u);
    arb_set(acb_realref(t), nu);
    arb_set(acb_realref(u), z);
    acb_hypgeom_bessel_j(t, t, u, prec);
    if (acb_is_finite(t) && acb_is_real(t))
        arb_swap(res, acb_realref(t));
    else
        arb_indeterminate(res);
    acb_clear(t);
    acb_clear(u);
}

void
arb_hypgeom_bessel_y(arb_t res, const arb_t nu, const arb_t z, slong prec)
{
    acb_t t, u;
    acb_init(t);
    acb_init(u);
    arb_set(acb_realref(t), nu);
    arb_set(acb_realref(u), z);
    acb_hypgeom_bessel_y(t, t, u, prec);
    if (acb_is_finite(t) && acb_is_real(t))
        arb_swap(res, acb_realref(t));
    else
        arb_indeterminate(res);
    acb_clear(t);
    acb_clear(u);
}

void
arb_hypgeom_bessel_jy(arb_t res1, arb_t res2, const arb_t nu, const arb_t z, slong prec)
{
    acb_t t, u;
    acb_init(t);
    acb_init(u);
    arb_set(acb_realref(t), nu);
    arb_set(acb_realref(u), z);
    acb_hypgeom_bessel_jy(t, u, t, u, prec);
    if (acb_is_finite(t) && acb_is_real(t))
        arb_swap(res1, acb_realref(t));
    else
        arb_indeterminate(res1);
    if (acb_is_finite(u) && acb_is_real(u))
        arb_swap(res2, acb_realref(u));
    else
        arb_indeterminate(res2);
    acb_clear(t);
    acb_clear(u);
}

void
arb_hypgeom_bessel_i(arb_t res, const arb_t nu, const arb_t z, slong prec)
{
    acb_t t, u;
    acb_init(t);
    acb_init(u);
    arb_set(acb_realref(t), nu);
    arb_set(acb_realref(u), z);
    acb_hypgeom_bessel_i(t, t, u, prec);
    if (acb_is_finite(t) && acb_is_real(t))
        arb_swap(res, acb_realref(t));
    else
        arb_indeterminate(res);
    acb_clear(t);
    acb_clear(u);
}

void
arb_hypgeom_bessel_i_scaled(arb_t res, const arb_t nu, const arb_t z, slong prec)
{
    acb_t t, u;
    acb_init(t);
    acb_init(u);
    arb_set(acb_realref(t), nu);
    arb_set(acb_realref(u), z);
    acb_hypgeom_bessel_i_scaled(t, t, u, prec);
    if (acb_is_finite(t) && acb_is_real(t))
        arb_swap(res, acb_realref(t));
    else
        arb_indeterminate(res);
    acb_clear(t);
    acb_clear(u);
}

void
arb_hypgeom_bessel_k(arb_t res, const arb_t nu, const arb_t z, slong prec)
{
    acb_t t, u;
    acb_init(t);
    acb_init(u);
    arb_set(acb_realref(t), nu);
    arb_set(acb_realref(u), z);
    acb_hypgeom_bessel_k(t, t, u, prec);
    if (acb_is_finite(t) && acb_is_real(t))
        arb_swap(res, acb_realref(t));
    else
        arb_indeterminate(res);
    acb_clear(t);
    acb_clear(u);
}

void
arb_hypgeom_bessel_k_scaled(arb_t res, const arb_t nu, const arb_t z, slong prec)
{
    acb_t t, u;
    acb_init(t);
    acb_init(u);
    arb_set(acb_realref(t), nu);
    arb_set(acb_realref(u), z);
    acb_hypgeom_bessel_k_scaled(t, t, u, prec);
    if (acb_is_finite(t) && acb_is_real(t))
        arb_swap(res, acb_realref(t));
    else
        arb_indeterminate(res);
    acb_clear(t);
    acb_clear(u);
}

void
arb_hypgeom_expint(arb_t res, const arb_t s, const arb_t z, slong prec)
{
    acb_t t, u;
    acb_init(t);
    acb_init(u);
    arb_set(acb_realref(t), s);
    arb_set(acb_realref(u), z);
    acb_hypgeom_expint(t, t, u, prec);
    if (acb_is_finite(t) && acb_is_real(t))
        arb_swap(res, acb_realref(t));
    else
        arb_indeterminate(res);
    acb_clear(t);
    acb_clear(u);
}

void
arb_hypgeom_gamma_lower(arb_t res, const arb_t s, const arb_t z, int regularized, slong prec)
{
    acb_t t, u;
    acb_init(t);
    acb_init(u);
    arb_set(acb_realref(t), s);
    arb_set(acb_realref(u), z);
    acb_hypgeom_gamma_lower(t, t, u, regularized, prec);
    if (acb_is_finite(t) && acb_is_real(t))
        arb_swap(res, acb_realref(t));
    else
        arb_indeterminate(res);
    acb_clear(t);
    acb_clear(u);
}

void
arb_hypgeom_gamma_upper(arb_t res, const arb_t s, const arb_t z, int regularized, slong prec)
{
    acb_t t, u;
    acb_init(t);
    acb_init(u);
    arb_set(acb_realref(t), s);
    arb_set(acb_realref(u), z);
    acb_hypgeom_gamma_upper(t, t, u, regularized, prec);
    if (acb_is_finite(t) && acb_is_real(t))
        arb_swap(res, acb_realref(t));
    else
        arb_indeterminate(res);
    acb_clear(t);
    acb_clear(u);
}

void
arb_hypgeom_beta_lower(arb_t res, const arb_t a, const arb_t b, const arb_t z, int regularized, slong prec)
{
    acb_t t, u, v;
    acb_init(t);
    acb_init(u);
    acb_init(v);
    arb_set(acb_realref(t), a);
    arb_set(acb_realref(u), b);
    arb_set(acb_realref(v), z);
    acb_hypgeom_beta_lower(t, t, u, v, regularized, prec);
    if (acb_is_finite(t) && acb_is_real(t))
        arb_swap(res, acb_realref(t));
    else
        arb_indeterminate(res);
    acb_clear(t);
    acb_clear(u);
    acb_clear(v);
}

void
arb_hypgeom_chebyshev_t(arb_t res, const arb_t nu, const arb_t z, slong prec)
{
    acb_t t, u;
    acb_init(t);
    acb_init(u);
    arb_set(acb_realref(t), nu);
    arb_set(acb_realref(u), z);
    acb_hypgeom_chebyshev_t(t, t, u, prec);
    if (acb_is_finite(t) && acb_is_real(t))
        arb_swap(res, acb_realref(t));
    else
        arb_indeterminate(res);
    acb_clear(t);
    acb_clear(u);
}

void
arb_hypgeom_chebyshev_u(arb_t res, const arb_t nu, const arb_t z, slong prec)
{
    acb_t t, u;
    acb_init(t);
    acb_init(u);
    arb_set(acb_realref(t), nu);
    arb_set(acb_realref(u), z);
    acb_hypgeom_chebyshev_u(t, t, u, prec);
    if (acb_is_finite(t) && acb_is_real(t))
        arb_swap(res, acb_realref(t));
    else
        arb_indeterminate(res);
    acb_clear(t);
    acb_clear(u);
}

void
arb_hypgeom_jacobi_p(arb_t res, const arb_t n, const arb_t a, const arb_t b, const arb_t z, slong prec)
{
    acb_t t, u, v, w;
    acb_init(t);
    acb_init(u);
    acb_init(v);
    acb_init(w);
    arb_set(acb_realref(t), n);
    arb_set(acb_realref(u), a);
    arb_set(acb_realref(v), b);
    arb_set(acb_realref(w), z);
    acb_hypgeom_jacobi_p(t, t, u, v, w, prec);
    if (acb_is_finite(t) && acb_is_real(t))
        arb_swap(res, acb_realref(t));
    else
        arb_indeterminate(res);
    acb_clear(t);
    acb_clear(u);
    acb_clear(v);
    acb_clear(w);
}

void
arb_hypgeom_gegenbauer_c(arb_t res, const arb_t n, const arb_t m, const arb_t z, slong prec)
{
    acb_t t, u, v;
    acb_init(t);
    acb_init(u);
    acb_init(v);
    arb_set(acb_realref(t), n);
    arb_set(acb_realref(u), m);
    arb_set(acb_realref(v), z);
    acb_hypgeom_gegenbauer_c(t, t, u, v, prec);
    if (acb_is_finite(t) && acb_is_real(t))
        arb_swap(res, acb_realref(t));
    else
        arb_indeterminate(res);
    acb_clear(t);
    acb_clear(u);
    acb_clear(v);
}

void
arb_hypgeom_laguerre_l(arb_t res, const arb_t n, const arb_t m, const arb_t z, slong prec)
{
    acb_t t, u, v;
    acb_init(t);
    acb_init(u);
    acb_init(v);
    arb_set(acb_realref(t), n);
    arb_set(acb_realref(u), m);
    arb_set(acb_realref(v), z);
    acb_hypgeom_laguerre_l(t, t, u, v, prec);
    if (acb_is_finite(t) && acb_is_real(t))
        arb_swap(res, acb_realref(t));
    else
        arb_indeterminate(res);
    acb_clear(t);
    acb_clear(u);
    acb_clear(v);
}

void
arb_hypgeom_hermite_h(arb_t res, const arb_t nu, const arb_t z, slong prec)
{
    acb_t t, u;
    acb_init(t);
    acb_init(u);
    arb_set(acb_realref(t), nu);
    arb_set(acb_realref(u), z);
    acb_hypgeom_hermite_h(t, t, u, prec);
    if (acb_is_finite(t) && acb_is_real(t))
        arb_swap(res, acb_realref(t));
    else
        arb_indeterminate(res);
    acb_clear(t);
    acb_clear(u);
}

void
arb_hypgeom_legendre_q(arb_t res, const arb_t n, const arb_t m, const arb_t z, int type, slong prec)
{
    acb_t t, u, v;
    acb_init(t);
    acb_init(u);
    acb_init(v);
    arb_set(acb_realref(t), n);
    arb_set(acb_realref(u), m);
    arb_set(acb_realref(v), z);
    acb_hypgeom_legendre_q(t, t, u, v, type, prec);
    if (acb_is_finite(t) && acb_is_real(t))
        arb_swap(res, acb_realref(t));
    else
        arb_indeterminate(res);
    acb_clear(t);
    acb_clear(u);
    acb_clear(v);
}

void
arb_hypgeom_dilog(arb_t res, const arb_t z, slong prec)
{
    acb_t t;
    acb_init(t);
    arb_set(acb_realref(t), z);
    acb_hypgeom_dilog(t, t, prec);
    if (acb_is_finite(t) && acb_is_real(t))
        arb_swap(res, acb_realref(t));
    else
        arb_indeterminate(res);
    acb_clear(t);
}
