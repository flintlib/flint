/*
    Copyright (C) 2026 Steve Clapper

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "acb_hypgeom.h"

/*
    Mittag-Leffler function E_{α,β}(z) for complex z, real α > 0, β ∈ ℝ.
    
    Three computational methods:
    1. Series expansion for small |z|
    2. Optimal Parabolic Contour (OPC) integration for mid-range |z|
    3. Asymptotic expansion for large |z|
    
    References:
    - R. Garrappa, "Numerical evaluation of two and three parameter Mittag-Leffler functions",
      SIAM J. Numer. Anal. 53(3), pp. 1350-1369, 2015.
      DOI: 10.1137/140971191
    - R. Gorenflo, A. Kilbas, F. Mainardi, S. Rogosin, "Mittag-Leffler Functions, 
      Related Topics and Applications", Springer, 2014.
*/

typedef struct
{
    acb_t s;        /* Singularity location */
    arb_t phi;      /* φ(s) = Re(s) + |s|/2 */
    int isorigin;   /* Flag for s = 0 */
} ml_sing_t;

typedef struct  
{
    slong regionindex;
    int unboundedright;
    arb_t mu;
    arb_t h;
    slong N;
    arb_t phileft;
    arb_t phiright;
} ml_region_t;

/* Initialize/clear structures */
static void ml_sing_init(ml_sing_t * x)
{
    acb_init(x->s);
    arb_init(x->phi);
    x->isorigin = 0;
}

static void ml_sing_clear(ml_sing_t * x)
{
    acb_clear(x->s);
    arb_clear(x->phi);
}

static void ml_region_init(ml_region_t * r)
{
    arb_init(r->mu);
    arb_init(r->h);
    arb_init(r->phileft);
    arb_init(r->phiright);
    r->regionindex = -1;
    r->unboundedright = 0;
    r->N = 0;
}

static void ml_region_clear(ml_region_t * r)
{
    arb_clear(r->mu);
    arb_clear(r->h);
    arb_clear(r->phileft);
    arb_clear(r->phiright);
}

/* φ(s) = Re(s) + |s|/2 */
static void ml_phi(arb_t out, const acb_t s, slong prec)
{
    arb_t re, ab;
    arb_init(re);
    arb_init(ab);
    
    acb_get_real(re, s);
    acb_abs(ab, s, prec);
    arb_add(out, re, ab, prec);
    arb_mul_2exp_si(out, out, -1);
    
    arb_clear(re);
    arb_clear(ab);
}

/* Comparison function for qsort */
static int ml_cmp_phi_qsort(const void *a, const void *b)
{
    const ml_sing_t *x = (const ml_sing_t *)a;
    const ml_sing_t *y = (const ml_sing_t *)b;
    return arf_cmp(arb_midref(x->phi), arb_midref(y->phi));
}

/* exp(i*theta) as acb */
static void ml_expi_theta(acb_t out, const arb_t theta, slong prec)
{
    arb_t s, c;
    arb_init(s);
    arb_init(c);
    
    arb_sin_cos(s, c, theta, prec);
    arb_set(acb_realref(out), c);
    arb_set(acb_imagref(out), s);
    
    arb_clear(s);
    arb_clear(c);
}

/* Real power x^y for x,y real */
static void ml_arb_pow(arb_t out, const arb_t x, const arb_t y, slong prec)
{
    arb_pow(out, x, y, prec);
}

/* f^(-1/p) for real f>0, p>0: exp(-log(f)/p) */
static void ml_f_to_minus_invp(arb_t out, const arb_t f, const arb_t p, slong prec)
{
    arb_t t;
    arb_init(t);
    
    arb_log(t, f, prec);
    arb_div(t, t, p, prec);
    arb_neg(t, t);
    arb_exp(out, t, prec);
    
    arb_clear(t);
}

/* ceil(upper bound of x) as slong - safe side */
static slong ml_ceil_ub_slong(const arb_t x, slong prec)
{
    arf_t ub;
    double d;
    slong n;
    
    arf_init(ub);
    arb_get_ubound_arf(ub, x, prec);
    d = arf_get_d(ub, ARF_RND_UP);
    arf_clear(ub);
    
    if (!isfinite(d) || d > (double)FLINT_SLONG_MAX / 4)
        return FLINT_SLONG_MAX / 4;
    
    n = (slong)ceil(d);
    if (n < 0)
        n = 0;
    
    return n;
}

/* Find singularities s_j = λ^(1/α) * exp(i*(θ + 2πj)/α) */
static ml_sing_t * ml_find_singularities(slong *count_out, const acb_t lambda, 
                                          const arb_t alpha, slong prec)
{
    arb_t theta, pi, twopi, tmp, jlower, jupper;
    arb_t abslambda, invalpha, radius, arg;
    arf_t lb, ub;
    fmpz_t jminz, jmaxz;
    slong jmin, jmax, npoles, i;
    ml_sing_t *sing;
    
    *count_out = 0;
    
    /* λ = 0: no poles, but include origin singularity */
    if (acb_is_zero(lambda))
    {
        sing = (ml_sing_t *)malloc(sizeof(ml_sing_t));
        if (!sing)
            flint_abort();
        
        ml_sing_init(&sing[0]);
        acb_zero(sing[0].s);
        arb_zero(sing[0].phi);
        sing[0].isorigin = 1;
        *count_out = 1;
        return sing;
    }
    
    arb_init(theta);
    arb_init(pi);
    arb_init(twopi);
    arb_init(tmp);
    arb_init(jlower);
    arb_init(jupper);
    arb_init(abslambda);
    arb_init(invalpha);
    arb_init(radius);
    arb_init(arg);
    arf_init(lb);
    arf_init(ub);
    fmpz_init(jminz);
    fmpz_init(jmaxz);
    
    /* theta = Arg(λ), -π < theta ≤ π */
    acb_arg(theta, lambda, prec);
    arb_const_pi(pi, prec);
    arb_mul_2exp_si(twopi, pi, 1);
    
    /* jlower = -α/2 - θ/(2π) */
    arb_div(tmp, theta, twopi, prec);
    arb_neg(jlower, alpha);
    arb_mul_2exp_si(jlower, jlower, -1);
    arb_sub(jlower, jlower, tmp, prec);
    
    /* jupper = α/2 - θ/(2π) */
    arb_set(jupper, alpha);
    arb_mul_2exp_si(jupper, jupper, -1);
    arb_sub(jupper, jupper, tmp, prec);
    
    /* Integer bounds */
    arb_get_lbound_arf(lb, jlower, prec);
    arb_get_ubound_arf(ub, jupper, prec);
    arf_get_fmpz(jminz, lb, ARF_RND_CEIL);
    arf_get_fmpz(jmaxz, ub, ARF_RND_FLOOR);
    
    if (fmpz_cmp(jmaxz, jminz) < 0)
    {
        /* No poles in main sheet: still include origin singularity */
        sing = (ml_sing_t *)malloc(sizeof(ml_sing_t));
        if (!sing)
            flint_abort();
        
        ml_sing_init(&sing[0]);
        acb_zero(sing[0].s);
        arb_zero(sing[0].phi);
        sing[0].isorigin = 1;
        *count_out = 1;
        goto cleanup_bounds;
    }
    
    jmin = fmpz_get_si(jminz);
    jmax = fmpz_get_si(jmaxz);
    npoles = jmax - jmin + 1;
    
    /* Allocate: origin + poles */
    sing = (ml_sing_t *)malloc((1 + npoles) * sizeof(ml_sing_t));
    if (!sing)
        flint_abort();
    
    for (i = 0; i < 1 + npoles; i++)
        ml_sing_init(&sing[i]);
    
    /* Origin singularity */
    acb_zero(sing[0].s);
    arb_zero(sing[0].phi);
    sing[0].isorigin = 1;
    
    /* radius = |λ|^(1/α) */
    acb_abs(abslambda, lambda, prec);
    arb_inv(invalpha, alpha, prec);
    ml_arb_pow(radius, abslambda, invalpha, prec);
    
    for (i = 0; i < npoles; i++)
    {
        slong j = jmin + i;
        
        /* arg = (θ + 2πj)/α */
        arb_mul_si(tmp, twopi, j, prec);
        arb_add(arg, theta, tmp, prec);
        arb_div(arg, arg, alpha, prec);
        
        /* exp(i*arg) */
        ml_expi_theta(sing[1+i].s, arg, prec);
        
        /* multiply by radius */
        acb_mul_arb(sing[1+i].s, sing[1+i].s, radius, prec);
        
        /* φ */
        ml_phi(sing[1+i].phi, sing[1+i].s, prec);
        sing[1+i].isorigin = 0;
    }
    
    /* Sort by φ increasing */
    qsort(sing, (size_t)(1 + npoles), sizeof(ml_sing_t), ml_cmp_phi_qsort);
    
    *count_out = 1 + npoles;
    
cleanup_bounds:
    arb_clear(theta);
    arb_clear(pi);
    arb_clear(twopi);
    arb_clear(tmp);
    arb_clear(jlower);
    arb_clear(jupper);
    arb_clear(abslambda);
    arb_clear(invalpha);
    arb_clear(radius);
    arb_clear(arg);
    arf_clear(lb);
    arf_clear(ub);
    fmpz_clear(jminz);
    fmpz_clear(jmaxz);
    
    return sing;
}

/* Series expansion: E_{α,β}(z) = Σ z^k / Γ(αk + β) */
static void ml_series(acb_t out, const acb_t z, const arb_t alpha, 
                      const arb_t beta, slong prec)
{
    acb_t sum, term, zpow;
    arb_t akb, gam;
    mag_t tol, mterm;
    slong k;
    
    acb_init(sum);
    acb_init(term);
    acb_init(zpow);
    arb_init(akb);
    arb_init(gam);
    mag_init(tol);
    mag_init(mterm);
    
    /* tol = 2^(-prec) */
    mag_set_ui_2exp_si(tol, 1, -prec);
    
    /* sum = 1/Γ(β) */
    arb_gamma(gam, beta, prec);
    arb_inv(gam, gam, prec);
    acb_set_arb(sum, gam);
    acb_set(zpow, z);
    
    for (k = 1; k < 200000; k++)
    {
        /* akb = αk + β */
        arb_mul_si(akb, alpha, k, prec);
        arb_add(akb, akb, beta, prec);
        
        arb_gamma(gam, akb, prec);
        
        acb_div_arb(term, zpow, gam, prec);
        
        acb_get_mag(mterm, term);
        if (mag_cmp(mterm, tol) < 0)
            break;
        
        acb_add(sum, sum, term, prec);
        acb_mul(zpow, zpow, z, prec);
    }
    
    acb_set(out, sum);
    
    acb_clear(sum);
    acb_clear(term);
    acb_clear(zpow);
    arb_clear(akb);
    arb_clear(gam);
    mag_clear(tol);
    mag_clear(mterm);
}

/* Asymptotic expansion */
static void ml_asymptotic(acb_t out, const acb_t z, const arb_t alpha, 
                          const arb_t beta, slong prec)
{
    acb_t sum, z1a, expz1a, leading, zma, powterm, term;
    arb_t inva, one, oneminusbeta, expo, betam, rg;
    mag_t tol, mterm;
    slong k;
    
    acb_init(sum);
    acb_init(z1a);
    acb_init(expz1a);
    acb_init(leading);
    acb_init(zma);
    acb_init(powterm);
    acb_init(term);
    arb_init(inva);
    arb_init(one);
    arb_init(oneminusbeta);
    arb_init(expo);
    arb_init(betam);
    arb_init(rg);
    mag_init(tol);
    mag_init(mterm);
    
    mag_set_ui_2exp_si(tol, 1, -prec);
    arb_one(one);
    
    /* inva = 1/α */
    arb_inv(inva, alpha, prec);
    
    /* z^(1/α) */
    acb_pow_arb(z1a, z, inva, prec);
    acb_exp(expz1a, z1a, prec);
    
    /* leading = (1/α) z^((1-β)/α) exp(z^(1/α)) */
    arb_sub(oneminusbeta, one, beta, prec);
    arb_mul(expo, oneminusbeta, inva, prec);
    acb_pow_arb(leading, z, expo, prec);
    acb_mul(leading, leading, expz1a, prec);
    acb_mul_arb(leading, leading, inva, prec);
    
    acb_set(sum, leading);
    
    /* z^(-α) */
    arb_neg(expo, alpha);
    acb_pow_arb(zma, z, expo, prec);
    
    acb_set(powterm, zma);
    
    for (k = 1; k < 2000; k++)
    {
        /* term = -z^(-αk) / Γ(β - αk) */
        arb_mul_si(betam, alpha, k, prec);
        arb_sub(betam, beta, betam, prec);
        
        /* Use reciprocal gamma */
        arb_rgamma(rg, betam, prec);
        
        acb_mul_arb(term, powterm, rg, prec);
        acb_neg(term, term);
        
        acb_get_mag(mterm, term);
        if (mag_cmp(mterm, tol) < 0)
            break;
        
        acb_add(sum, sum, term, prec);
        acb_mul(powterm, powterm, zma, prec);
    }
    
    acb_set(out, sum);
    
    acb_clear(sum);
    acb_clear(z1a);
    acb_clear(expz1a);
    acb_clear(leading);
    acb_clear(zma);
    acb_clear(powterm);
    acb_clear(term);
    arb_clear(inva);
    arb_clear(one);
    arb_clear(oneminusbeta);
    arb_clear(expo);
    arb_clear(betam);
    arb_clear(rg);
    mag_clear(tol);
    mag_clear(mterm);
}

/* Parabola: z(u) = μ(i*u - 1)^2 */
static void ml_parabola(acb_t z, const arb_t mu, const arb_t u, slong prec)
{
    acb_t w;
    acb_init(w);
    
    /* w = i*u + 1 */
    acb_set_arb(w, u);
    acb_mul_onei(w, w);
    acb_add_ui(w, w, 1, prec);
    
    acb_sqr(z, w, prec);
    acb_mul_arb(z, z, mu, prec);
    
    acb_clear(w);
}

/* Parabola derivative: z'(u) = 2μ(i - u) */
static void ml_parabola_deriv(acb_t zp, const arb_t mu, const arb_t u, slong prec)
{
    acb_t w;
    acb_init(w);
    
    acb_onei(w);
    acb_sub_arb(w, w, u, prec);
    acb_mul_arb(zp, w, mu, prec);
    acb_mul_2exp_si(zp, zp, 1);  /* × 2 */
    
    acb_clear(w);
}

/* g(u) = exp(z(u)*t) * z(u)^(α-β) * z'(u) / (z(u)^α - λ) */
static void ml_g(acb_t out, const arb_t u, const arb_t mu, const arb_t t,
                 const acb_t lambda, const arb_t alpha, const arb_t beta, slong prec)
{
    acb_t z, zp, expzt, za, zab, num, den;
    arb_t aminusb;
    
    acb_init(z);
    acb_init(zp);
    acb_init(expzt);
    acb_init(za);
    acb_init(zab);
    acb_init(num);
    acb_init(den);
    arb_init(aminusb);
    
    ml_parabola(z, mu, u, prec);
    ml_parabola_deriv(zp, mu, u, prec);
    
    acb_mul_arb(expzt, z, t, prec);
    acb_exp(expzt, expzt, prec);
    
    acb_pow_arb(za, z, alpha, prec);
    
    arb_sub(aminusb, alpha, beta, prec);
    acb_pow_arb(zab, z, aminusb, prec);
    
    acb_sub(den, za, lambda, prec);
    
    acb_mul(num, expzt, zab, prec);
    acb_mul(num, num, zp, prec);
    
    acb_div(out, num, den, prec);
    
    acb_clear(z);
    acb_clear(zp);
    acb_clear(expzt);
    acb_clear(za);
    acb_clear(zab);
    acb_clear(num);
    acb_clear(den);
    arb_clear(aminusb);
}

/* Residue at pole s: (1/α) s^(1-β) exp(s*t) */
static void ml_residue(acb_t out, const acb_t s, const arb_t t,
                       const arb_t alpha, const arb_t beta, slong prec)
{
    acb_t p, e;
    arb_t one, oneminusbeta;
    
    acb_init(p);
    acb_init(e);
    arb_init(one);
    arb_init(oneminusbeta);
    
    arb_one(one);
    arb_sub(oneminusbeta, one, beta, prec);
    
    acb_pow_arb(p, s, oneminusbeta, prec);
    
    acb_mul_arb(e, s, t, prec);
    acb_exp(e, e, prec);
    
    acb_mul(out, p, e, prec);
    acb_div_arb(out, out, alpha, prec);
    
    acb_clear(p);
    acb_clear(e);
    arb_clear(one);
    arb_clear(oneminusbeta);
}

/* Compute bounded region parameters via Garrappa balancing */
static int ml_region_params_bounded(ml_region_t *R, const arb_t phiL, const arb_t phiR,
                                   const arb_t t, const arb_t eps, slong prec)
{
    arb_t sqrtL, sqrtR, d, s, fmin, fmax, ftar, f;
    arb_t epsbar, logeps, w, one, two, twopi;
    arb_t A, denom, tmp1, tmp2, Nh;
    
    arb_init(sqrtL);
    arb_init(sqrtR);
    arb_init(d);
    arb_init(s);
    arb_init(fmin);
    arb_init(fmax);
    arb_init(ftar);
    arb_init(f);
    arb_init(epsbar);
    arb_init(logeps);
    arb_init(w);
    arb_init(one);
    arb_init(two);
    arb_init(twopi);
    arb_init(A);
    arb_init(denom);
    arb_init(tmp1);
    arb_init(tmp2);
    arb_init(Nh);
    
    arb_one(one);
    arb_set_si(two, 2);
    
    if (arb_le(phiR, phiL))
        goto fail;
    
    arb_sqrt(sqrtL, phiL, prec);
    arb_sqrt(sqrtR, phiR, prec);
    arb_sub(d, sqrtR, sqrtL, prec);
    
    if (arb_is_zero(d) || arb_contains_zero(d))
        goto fail;
    
    arb_add(s, sqrtR, sqrtL, prec);
    arb_div(fmin, s, d, prec);
    
    arb_set_d(fmax, 10.0);
    arb_set_d(ftar, 5.0);
    
    if (arb_gt(fmin, fmax))
        goto fail;
    
    arb_set(f, ftar);
    if (arb_lt(f, fmin))
        arb_set(f, fmin);
    
    arb_div(epsbar, eps, f, prec);
    arb_log(logeps, epsbar, prec);
    
    if (arb_is_zero(logeps) || arb_contains_zero(logeps))
        goto fail;
    
    /* w = -(φ_R * t) / log(ε̄) */
    arb_mul(tmp1, phiR, t, prec);
    arb_div(w, tmp1, logeps, prec);
    arb_neg(w, w);
    
    /* A = (1 + w)√φ_L + √φ_R */
    arb_add(tmp1, one, w, prec);
    arb_mul(A, tmp1, sqrtL, prec);
    arb_add(A, A, sqrtR, prec);
    
    /* denom = 2 + w */
    arb_add(denom, two, w, prec);
    
    if (arb_is_zero(denom) || arb_contains_zero(denom))
        goto fail;
    
    /* μ = √A / (2 + w) */
    arb_sqrt(tmp1, A, prec);
    arb_div(R->mu, tmp1, denom, prec);
    
    /* h = (2πi / log(ε̄)) * (√φ_R - √φ_L) / A */
    arb_const_pi(twopi, prec);
    arb_mul_2exp_si(twopi, twopi, 1);
    arb_div(tmp1, twopi, logeps, prec);
    arb_neg(tmp1, tmp1);
    arb_mul(tmp1, tmp1, d, prec);
    arb_div(R->h, tmp1, A, prec);
    
    /* N = ceil(√((t*μ/log(ε̄)) + 1) / h) */
    arb_mul(tmp1, t, R->mu, prec);
    arb_div(tmp2, logeps, tmp1, prec);
    arb_neg(tmp2, tmp2);
    arb_add(tmp2, tmp2, one, prec);
    arb_sqrt(tmp2, tmp2, prec);
    arb_div(Nh, tmp2, R->h, prec);
    
    R->N = ml_ceil_ub_slong(Nh, prec) + 1;
    if (R->N < 8)
        R->N = 8;
    
    arb_set(R->phileft, phiL);
    arb_set(R->phiright, phiR);
    
    arb_clear(sqrtL);
    arb_clear(sqrtR);
    arb_clear(d);
    arb_clear(s);
    arb_clear(fmin);
    arb_clear(fmax);
    arb_clear(ftar);
    arb_clear(f);
    arb_clear(epsbar);
    arb_clear(logeps);
    arb_clear(w);
    arb_clear(one);
    arb_clear(two);
    arb_clear(twopi);
    arb_clear(A);
    arb_clear(denom);
    arb_clear(tmp1);
    arb_clear(tmp2);
    arb_clear(Nh);
    
    return 1;
    
fail:
    arb_clear(sqrtL);
    arb_clear(sqrtR);
    arb_clear(d);
    arb_clear(s);
    arb_clear(fmin);
    arb_clear(fmax);
    arb_clear(ftar);
    arb_clear(f);
    arb_clear(epsbar);
    arb_clear(logeps);
    arb_clear(w);
    arb_clear(one);
    arb_clear(two);
    arb_clear(twopi);
    arb_clear(A);
    arb_clear(denom);
    arb_clear(tmp1);
    arb_clear(tmp2);
    arb_clear(Nh);
    
    return 0;
}

/* OPC mid-range method */
static int ml_opc_mid(acb_t out, const acb_t lambda, const arb_t alpha,
                      const arb_t beta, const arb_t eps, slong prec)
{
    ml_sing_t *sing = NULL;
    slong ns = 0, J, j, bestj, bestN;
    arb_t t;
    acb_t integral, residues, tmpc, gval;
    arb_t u, twopi;
    ml_region_t R;
    
    arb_init(t);
    arb_one(t);
    
    sing = ml_find_singularities(&ns, lambda, alpha, prec);
    if (!sing || ns < 1)
        return 0;
    
    J = ns - 1;
    if (J < 0)
        J = 0;
    
    bestj = -1;
    bestN = FLINT_SLONG_MAX;
    
    for (j = 0; j < J - 1; j++)
    {
        ml_region_init(&R);
        
        if (ml_region_params_bounded(&R, sing[j].phi, sing[j+1].phi, t, eps, prec))
        {
            if (R.N > 0 && R.N < bestN)
            {
                bestN = R.N;
                bestj = j;
            }
        }
        
        ml_region_clear(&R);
    }
    
    if (bestj < 0)
    {
        for (j = 0; j < ns; j++)
            ml_sing_clear(&sing[j]);
        free(sing);
        arb_clear(t);
        return 0;
    }
    
    ml_region_init(&R);
    if (!ml_region_params_bounded(&R, sing[bestj].phi, sing[bestj+1].phi, t, eps, prec))
    {
        ml_region_clear(&R);
        for (j = 0; j < ns; j++)
            ml_sing_clear(&sing[j]);
        free(sing);
        arb_clear(t);
        return 0;
    }
    
    R.regionindex = bestj;
    
    acb_init(residues);
    acb_zero(residues);
    acb_init(tmpc);
    
    for (j = bestj + 1; j < J; j++)
    {
        if (!sing[j].isorigin)
        {
            ml_residue(tmpc, sing[j].s, t, alpha, beta, prec);
            acb_add(residues, residues, tmpc, prec);
        }
    }
    
    acb_init(integral);
    acb_zero(integral);
    acb_init(gval);
    arb_init(u);
    arb_init(twopi);
    
    for (slong k = -R.N; k <= R.N; k++)
    {
        arb_set_si(u, k);
        arb_mul(u, u, R.h, prec);
        
        ml_g(gval, u, R.mu, t, lambda, alpha, beta, prec);
        
        if (k != -R.N && k != R.N)
            acb_mul_2exp_si(gval, gval, -1);
        
        acb_add(integral, integral, gval, prec);
    }
    
    acb_mul_arb(integral, integral, R.h, prec);
    
    /* Divide by 2πi */
    arb_const_pi(twopi, prec);
    arb_mul_2exp_si(twopi, twopi, 1);
    acb_div_arb(integral, integral, twopi, prec);
    acb_div_onei(integral, integral);
    
    acb_add(out, residues, integral, prec);
    
    acb_clear(integral);
    acb_clear(residues);
    acb_clear(tmpc);
    acb_clear(gval);
    arb_clear(u);
    arb_clear(twopi);
    ml_region_clear(&R);
    
    for (j = 0; j < ns; j++)
        ml_sing_clear(&sing[j]);
    free(sing);
    arb_clear(t);
    
    return 1;
}

/* Main public API */
void acb_hypgeom_mittag_leffler_e(acb_t out, const acb_t z, const arb_t alpha,
                                  const arb_t beta, double epsd, slong prec)
{
    arb_t eps, absz, seriesthr, asympthr, gam;
    
    arb_init(eps);
    arb_init(absz);
    arb_init(seriesthr);
    arb_init(asympthr);
    arb_init(gam);
    
    arb_set_d(eps, epsd);
    
    /* Default thresholds */
    arb_set_d(seriesthr, 0.3);
    arb_set_d(asympthr, 50.0);
    
    /* E_{α,β}(0) = 1/Γ(β) */
    if (acb_is_zero(z))
    {
        arb_gamma(gam, beta, prec);
        arb_inv(gam, gam, prec);
        acb_set_arb(out, gam);
        goto done;
    }
    
    /* E_{1,1}(z) = exp(z) */
    if (arb_is_one(alpha) && arb_is_one(beta))
    {
        acb_exp(out, z, prec);
        goto done;
    }
    
    acb_abs(absz, z, prec);
    
    if (arb_lt(absz, seriesthr))
    {
        ml_series(out, z, alpha, beta, prec);
    }
    else if (arb_gt(absz, asympthr))
    {
        ml_asymptotic(out, z, alpha, beta, prec);
    }
    else
    {
        if (!ml_opc_mid(out, z, alpha, beta, eps, prec))
        {
            /* Fallback to series */
            ml_series(out, z, alpha, beta, prec);
        }
    }
    
done:
    arb_clear(eps);
    arb_clear(absz);
    arb_clear(seriesthr);
    arb_clear(asympthr);
    arb_clear(gam);
}
