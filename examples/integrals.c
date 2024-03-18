/*
  Rigorous numerical integration (with fast convergence for piecewise
  holomorphic functions) using Gauss-Legendre quadrature and adaptive
  subdivision.

  Author: Fredrik Johansson.
  This file is in the public domain.
*/

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "profiler.h"
#include "arb_hypgeom.h"
#include "acb_hypgeom.h"
#include "acb_dirichlet.h"
#include "acb_modular.h"
#include "acb_calc.h"

/* ------------------------------------------------------------------------- */
/*  Example integrands                                                       */
/* ------------------------------------------------------------------------- */

/* f(z) = sin(z) */
int
f_sin(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_sin(res, z, prec);

    return 0;
}

/* f(z) = floor(z) */
int
f_floor(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_real_floor(res, z, order != 0, prec);

    return 0;
}

/* f(z) = sqrt(1-z^2) */
int
f_circle(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_one(res);
    acb_submul(res, z, z, prec);
    acb_real_sqrtpos(res, res, order != 0, prec);

    return 0;
}

/* f(z) = 1/(1+z^2) */
int
f_atanderiv(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_mul(res, z, z, prec);
    acb_add_ui(res, res, 1, prec);
    acb_inv(res, res, prec);

    return 0;
}

/* f(z) = sin(z + exp(z)) -- Rump's oscillatory example */
int
f_rump(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_exp(res, z, prec);
    acb_add(res, res, z, prec);
    acb_sin(res, res, prec);

    return 0;
}

/* f(z) = |z^4 + 10z^3 + 19z^2 - 6z - 6| exp(z)   (for real z) --
   Helfgott's integral on MathOverflow */
int
f_helfgott(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_add_si(res, z, 10, prec);
    acb_mul(res, res, z, prec);
    acb_add_si(res, res, 19, prec);
    acb_mul(res, res, z, prec);
    acb_add_si(res, res, -6, prec);
    acb_mul(res, res, z, prec);
    acb_add_si(res, res, -6, prec);

    acb_real_abs(res, res, order != 0, prec);

    if (acb_is_finite(res))
    {
        acb_t t;
        acb_init(t);
        acb_exp(t, z, prec);
        acb_mul(res, res, t, prec);
        acb_clear(t);
    }

    return 0;
}

/* f(z) = zeta(z) */
int
f_zeta(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_zeta(res, z, prec);

    return 0;
}

/* f(z) = z sin(1/z), assume on real interval */
int
f_essing2(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    if ((order == 0) && acb_is_real(z) && arb_contains_zero(acb_realref(z)))
    {
        /* todo: arb_zero_pm_one, arb_unit_interval? */
        acb_zero(res);
        mag_one(arb_radref(acb_realref(res)));
    }
    else
    {
        acb_inv(res, z, prec);
        acb_sin(res, res, prec);
    }

    acb_mul(res, res, z, prec);

    return 0;
}

/* f(z) = sin(1/z), assume on real interval */
int
f_essing(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    if ((order == 0) && acb_is_real(z) && arb_contains_zero(acb_realref(z)))
    {
        /* todo: arb_zero_pm_one, arb_unit_interval? */
        acb_zero(res);
        mag_one(arb_radref(acb_realref(res)));
    }
    else
    {
        acb_inv(res, z, prec);
        acb_sin(res, res, prec);
    }

    return 0;
}

/* f(z) = exp(-z) z^1000 */
int
f_factorial1000(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    acb_t t;

    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_init(t);
    acb_pow_ui(t, z, 1000, prec);
    acb_neg(res, z);
    acb_exp(res, res, prec);
    acb_mul(res, res, t, prec);
    acb_clear(t);

    return 0;
}

/* f(z) = gamma(z) */
int
f_gamma(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_gamma(res, z, prec);

    return 0;
}

/* f(z) = sin(z) + exp(-200-z^2) */
int
f_sin_plus_small(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    acb_t t;

    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_init(t);
    acb_mul(t, z, z, prec);
    acb_add_ui(t, t, 200, prec);
    acb_neg(t, t);
    acb_exp(t, t, prec);
    acb_sin(res, z, prec);
    acb_add(res, res, t, prec);
    acb_clear(t);

    return 0;
}

/* f(z) = exp(z) */
int
f_exp(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_exp(res, z, prec);

    return 0;
}

/* f(z) = exp(-z^2) */
int
f_gaussian(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_mul(res, z, z, prec);
    acb_neg(res, res);
    acb_exp(res, res, prec);

    return 0;
}

int
f_monster(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    acb_t t;

    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_init(t);

    acb_exp(t, z, prec);
    acb_real_floor(res, t, order != 0, prec);

    if (acb_is_finite(res))
    {
        acb_sub(res, t, res, prec);
        acb_add(t, t, z, prec);
        acb_sin(t, t, prec);
        acb_mul(res, res, t, prec);
    }

    acb_clear(t);

    return 0;
}

/* f(z) = sech(10(x-0.2))^2 + sech(100(x-0.4))^4 + sech(1000(x-0.6))^6 */
int
f_spike(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    acb_t a, b, c;

    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_init(a);
    acb_init(b);
    acb_init(c);

    acb_mul_ui(a, z, 10, prec);
    acb_sub_ui(a, a, 2, prec);
    acb_sech(a, a, prec);
    acb_pow_ui(a, a, 2, prec);

    acb_mul_ui(b, z, 100, prec);
    acb_sub_ui(b, b, 40, prec);
    acb_sech(b, b, prec);
    acb_pow_ui(b, b, 4, prec);

    acb_mul_ui(c, z, 1000, prec);
    acb_sub_ui(c, c, 600, prec);
    acb_sech(c, c, prec);
    acb_pow_ui(c, c, 6, prec);

    acb_add(res, a, b, prec);
    acb_add(res, res, c, prec);

    acb_clear(a);
    acb_clear(b);
    acb_clear(c);

    return 0;
}

/* f(z) = sech(z) */
int
f_sech(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_sech(res, z, prec);

    return 0;
}

/* f(z) = sech^3(z) */
int
f_sech3(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_sech(res, z, prec);
    acb_cube(res, res, prec);

    return 0;
}

/* f(z) = -log(z) / (1 + z) */
int
f_log_div1p(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    acb_t t;

    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_init(t);
    acb_add_ui(t, z, 1, prec);
    acb_log(res, z, prec);
    acb_div(res, res, t, prec);
    acb_neg(res, res);

    acb_clear(t);

    return 0;
}

/* f(z) = z exp(-z) / (1 + exp(-z)) */
int
f_log_div1p_transformed(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    acb_t t;

    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_init(t);
    acb_neg(t, z);
    acb_exp(t, t, prec);
    acb_add_ui(res, t, 1, prec);
    acb_div(res, t, res, prec);
    acb_mul(res, res, z, prec);

    acb_clear(t);

    return 0;
}

int
f_elliptic_p_laurent_n(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    slong n;
    acb_t tau;

    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    n = ((slong *)(param))[0];

    acb_init(tau);
    acb_onei(tau);

    acb_modular_elliptic_p(res, z, tau, prec);

    acb_pow_si(tau, z, -n - 1, prec);
    acb_mul(res, res, tau, prec);

    acb_clear(tau);

    return 0;
}

/* f(z) = zeta'(z) / zeta(z) */
int
f_zeta_frac(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    acb_struct t[2];

    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_init(t);
    acb_init(t + 1);

    acb_dirichlet_zeta_jet(t, z, 0, 2, prec);
    acb_div(res, t + 1, t, prec);

    acb_clear(t);
    acb_clear(t + 1);

    return 0;
}

int
f_lambertw(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    acb_t t;

    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_init(t);

    prec = FLINT_MIN(prec, acb_rel_accuracy_bits(z) + 10);

    if (order != 0)
    {
        /* check for branch cut */
        arb_const_e(acb_realref(t), prec);
        acb_inv(t, t, prec);
        acb_add(t, t, z, prec);

        if (arb_contains_zero(acb_imagref(t)) &&
            arb_contains_nonpositive(acb_realref(t)))
        {
            acb_indeterminate(t);
        }
    }

    if (acb_is_finite(t))
    {
        fmpz_t k;
        fmpz_init(k);
        acb_lambertw(res, z, k, 0, prec);
        fmpz_clear(k);
    }
    else
    {
        acb_indeterminate(res);
    }

    acb_clear(t);

    return 0;
}

/* f(z) = max(sin(z), cos(z)) */
int
f_max_sin_cos(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    acb_t s, c;

    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_init(s);
    acb_init(c);

    acb_sin_cos(s, c, z, prec);
    acb_real_max(res, s, c, order != 0, prec);

    acb_clear(s);
    acb_clear(c);

    return 0;
}

/* f(z) = erf(z/sqrt(0.0002)*0.5 +1.5)*exp(-z), example provided by Silviu-Ioan Filip */
int
f_erf_bent(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    acb_t t;

    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_init(t);

    acb_set_ui(t, 1250);
    acb_sqrt(t, t, prec);
    acb_mul(t, t, z, prec);
    acb_set_d(res, 1.5);
    acb_add(res, res, t, prec);
    acb_hypgeom_erf(res, res, prec);

    acb_neg(t, z);
    acb_exp(t, t, prec);
    acb_mul(res, res, t, prec);

    acb_clear(t);

    return 0;
}

/* f(z) = Ai(z) */
int
f_airy_ai(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */


    acb_hypgeom_airy(res, NULL, NULL, NULL, z, prec);

    return 0;
}

int
f_horror(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    acb_t s, t;

    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_init(s);
    acb_init(t);

    acb_real_floor(res, z, order != 0, prec);

    if (acb_is_finite(res))
    {
        acb_sub(res, z, res, prec);
        acb_set_d(t, 0.5);
        acb_sub(res, res, t, prec);
        acb_sin_cos(s, t, z, prec);
        acb_real_max(s, s, t, order != 0, prec);
        acb_mul(res, res, s, prec);
    }

    acb_clear(s);
    acb_clear(t);

    return 0;
}

/* f(z) = sqrt(z) */
int
f_sqrt(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_sqrt_analytic(res, z, order != 0, prec);

    return 0;
}

/* f(z) = 1/sqrt(z) */
int
f_rsqrt(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_rsqrt_analytic(res, z, order != 0, prec);

    return 0;
}

/* f(z) = rgamma(z) */
int
f_rgamma(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_rgamma(res, z, prec);

    return 0;
}

/* f(z) = exp(-z^2+iz) */
int
f_gaussian_twist(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_mul_onei(res, z);
    acb_submul(res, z, z, prec);
    acb_exp(res, res, prec);

    return 0;
}

/* f(z) = exp(-z) Ai(-z) */
int
f_exp_airy(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    acb_t t;

    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_init(t);
    acb_neg(t, z);
    acb_hypgeom_airy(res, NULL, NULL, NULL, t, prec);
    acb_exp(t, t, prec);
    acb_mul(res, res, t, prec);
    acb_clear(t);

    return 0;
}

/* f(z) = z sin(z) / (1 + cos(z)^2) */
int
f_sin_cos_frac(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    acb_t s, c;

    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_init(s);
    acb_init(c);

    acb_sin_cos(s, c, z, prec);
    acb_mul(c, c, c, prec);
    acb_add_ui(c, c, 1, prec);
    acb_mul(s, s, z, prec);
    acb_div(res, s, c, prec);

    acb_clear(s);
    acb_clear(c);

    return 0;
}

/* f(z) = sin((1/1000 + (1-z)^2)^(-3/2)), example from Mioara Joldes' thesis
                                          (suggested by Nicolas Brisebarre) */
int
f_sin_near_essing(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    acb_t t, u;

    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_init(t);
    acb_init(u);

    acb_sub_ui(t, z, 1, prec);
    acb_neg(t, t);
    acb_mul(t, t, t, prec);
    acb_one(u);
    acb_div_ui(u, u, 1000, prec);
    acb_add(t, t, u, prec);
    acb_set_d(u, -1.5);
    acb_pow_analytic(t, t, u, order != 0, prec);
    acb_sin(res, t, prec);

    acb_clear(t);
    acb_clear(u);

    return 0;
}

/* f(z) = exp(-z) (I_0(z/k))^k, from Bruno Salvy */
int
f_scaled_bessel(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    acb_t nu;

    ulong k = ((ulong *) param)[0];

    acb_init(nu);
    acb_div_ui(res, z, k, prec);
    acb_hypgeom_bessel_i_scaled(res, nu, res, prec);
    acb_pow_ui(res, res, k, prec);
    acb_clear(nu);

    return 0;
}

/*
Bound for scaled Bessel function: 2/(2 pi x)^(1/2)
Bound for tail of integral: 2 N (k / (pi N))^(k / 2) / (k - 2).
*/
void
scaled_bessel_tail_bound(arb_t b, ulong k, const arb_t N, slong prec)
{
    arb_const_pi(b, prec);
    arb_mul(b, b, N, prec);
    arb_ui_div(b, k, b, prec);
    arb_sqrt(b, b, prec);
    arb_pow_ui(b, b, k, prec);
    arb_mul(b, b, N, prec);
    arb_mul_ui(b, b, 2, prec);
    arb_div_ui(b, b, k - 2, prec);
}

void
scaled_bessel_select_N(arb_t N, ulong k, slong prec)
{
    slong e;
    double f = log(k/3.14159265358979)/log(2);

    e = 1;
    while ((k / 2.0) * (e - f) - e < prec + 5)
        e++;

    arb_one(N);
    arb_mul_2exp_si(N, N, e);
}

/* ------------------------------------------------------------------------- */
/*  Main test program                                                        */
/* ------------------------------------------------------------------------- */

#define NUM_INTEGRALS 37

const char * descr[NUM_INTEGRALS] =
{
    "int_0^100 sin(x) dx",
    "4 int_0^1 1/(1+x^2) dx",
    "2 int_0^{inf} 1/(1+x^2) dx   (using domain truncation)",
    "4 int_0^1 sqrt(1-x^2) dx",
    "int_0^8 sin(x+exp(x)) dx",
    "int_1^101 floor(x) dx",
    "int_0^1 |x^4+10x^3+19x^2-6x-6| exp(x) dx",
    "1/(2 pi i) int zeta(s) ds  (closed path around s = 1)",
    "int_0^1 sin(1/x) dx  (slow convergence, use -heap and/or -tol)",
    "int_0^1 x sin(1/x) dx  (slow convergence, use -heap and/or -tol)",
    "int_0^10000 x^1000 exp(-x) dx",
    "int_1^{1+1000i} gamma(x) dx",
    "int_{-10}^{10} sin(x) + exp(-200-x^2) dx",
    "int_{-1020}^{-1010} exp(x) dx  (use -tol 0 for relative error)",
    "int_0^{inf} exp(-x^2) dx   (using domain truncation)",
    "int_0^1 sech(10(x-0.2))^2 + sech(100(x-0.4))^4 + sech(1000(x-0.6))^6 dx",
    "int_0^8 (exp(x)-floor(exp(x))) sin(x+exp(x)) dx  (use higher -eval)",
    "int_0^{inf} sech(x) dx   (using domain truncation)",
    "int_0^{inf} sech^3(x) dx   (using domain truncation)",
    "int_0^1 -log(x)/(1+x) dx   (using domain truncation)",
    "int_0^{inf} x exp(-x)/(1+exp(-x)) dx   (using domain truncation)",
    "int_C wp(x)/x^(11) dx   (contour for 10th Laurent coefficient of Weierstrass p-function)",
    "N(1000) = count zeros with 0 < t <= 1000 of zeta(s) using argument principle",
    "int_0^{1000} W_0(x) dx",
    "int_0^pi max(sin(x), cos(x)) dx",
    "int_{-1}^1 erf(x/sqrt(0.0002)*0.5+1.5)*exp(-x) dx",
    "int_{-10}^10 Ai(x) dx",
    "int_0^10 (x-floor(x)-1/2) max(sin(x),cos(x)) dx",
    "int_{-1-i}^{-1+i} sqrt(x) dx",
    "int_0^{inf} exp(-x^2+ix) dx   (using domain truncation)",
    "int_0^{inf} exp(-x) Ai(-x) dx   (using domain truncation)",
    "int_0^pi x sin(x) / (1 + cos(x)^2) dx",
    "int_0^3 sin(0.001 + (1-x)^2)^(-3/2)) dx  (slow convergence, use higher -eval)",
    "int_0^{inf} exp(-x) I_0(x/3)^3 dx   (using domain truncation)",
    "int_0^{inf} exp(-x) I_0(x/15)^{15} dx   (using domain truncation)",
    "int_{-1-i}^{-1+i} 1/sqrt(x) dx",
    "int_0^{inf} 1/gamma(x) dx   (using domain truncation)",
};

int main(int argc, char *argv[])
{
    acb_t s, t, a, b;
    mag_t tol;
    slong num_threads;
    slong prec, goal;
    slong N;
    ulong k;
    int integral, ifrom, ito;
    int i, twice, havegoal, havetol;
    acb_calc_integrate_opt_t options;

    ifrom = ito = -1;

    for (i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "-i"))
        {
            if (!strcmp(argv[i+1], "all"))
            {
                ifrom = 0;
                ito = NUM_INTEGRALS - 1;
            }
            else
            {
                ifrom = ito = atol(argv[i+1]);
                if (ito < 0 || ito >= NUM_INTEGRALS)
                    flint_abort();
            }
        }
    }

    if (ifrom == -1)
    {
        flint_printf("Compute integrals using acb_calc_integrate.\n");
        flint_printf("Usage: integrals -i n [-prec p] [-tol eps] [-twice] [...]\n\n");
        flint_printf("-i n       - compute integral n (0 <= n <= %d), or \"-i all\"\n", NUM_INTEGRALS - 1);
        flint_printf("-prec p    - precision in bits (default p = 64)\n");
        flint_printf("-goal p    - approximate relative accuracy goal (default p)\n");
        flint_printf("-tol eps   - approximate absolute error goal (default 2^-p)\n");
        flint_printf("-twice     - run twice (to see overhead of computing nodes)\n");
        flint_printf("-heap      - use heap for subinterval queue\n");
        flint_printf("-verbose   - show information\n");
        flint_printf("-verbose2  - show more information\n");
        flint_printf("-deg n     - use quadrature degree up to n\n");
        flint_printf("-eval n    - limit number of function evaluations to n\n");
        flint_printf("-depth n   - limit subinterval queue size to n\n");
        flint_printf("-threads n - use parallel computation with n threads\n\n");
        flint_printf("Implemented integrals:\n");
        for (integral = 0; integral < NUM_INTEGRALS; integral++)
            flint_printf("I%d = %s\n", integral, descr[integral]);
        flint_printf("\n");
        return 1;
    }

    acb_calc_integrate_opt_init(options);

    prec = 64;
    twice = 0;
    goal = 0;
    havetol = havegoal = 0;
    num_threads = 1;

    acb_init(a);
    acb_init(b);
    acb_init(s);
    acb_init(t);
    mag_init(tol);

    for (i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "-prec"))
        {
            prec = atol(argv[i+1]);
        }
        else if (!strcmp(argv[i], "-twice"))
        {
            twice = 1;
        }
        else if (!strcmp(argv[i], "-goal"))
        {
            goal = atol(argv[i+1]);
            if (goal < 0)
            {
                flint_printf("expected goal >= 0\n");
                return 1;
            }
            havegoal = 1;
        }
        else if (!strcmp(argv[i], "-tol"))
        {
            arb_t x;
            arb_init(x);
            arb_set_str(x, argv[i+1], 10);
            arb_get_mag(tol, x);
            arb_clear(x);
            havetol = 1;
        }
        else if (!strcmp(argv[i], "-deg"))
        {
            options->deg_limit = atol(argv[i+1]);
        }
        else if (!strcmp(argv[i], "-eval"))
        {
            options->eval_limit = atol(argv[i+1]);
        }
        else if (!strcmp(argv[i], "-depth"))
        {
            options->depth_limit = atol(argv[i+1]);
        }
        else if (!strcmp(argv[i], "-verbose"))
        {
            options->verbose = 1;
        }
        else if (!strcmp(argv[i], "-verbose2"))
        {
            options->verbose = 2;
        }
        else if (!strcmp(argv[i], "-heap"))
        {
            options->use_heap = 1;
        }
        else if (!strcmp(argv[i], "-threads"))
        {
            num_threads = atol(argv[i+1]);
        }
    }

    if (num_threads >= 2)
        flint_set_num_threads(num_threads);

    if (!havegoal)
        goal = prec;

    if (!havetol)
        mag_set_ui_2exp_si(tol, 1, -prec);

    for (integral = ifrom; integral <= ito; integral++)
    {
        flint_printf("I%d = %s ...\n", integral, descr[integral]);

        for (i = 0; i < 1 + twice; i++)
        {
            TIMEIT_ONCE_START
            switch (integral)
            {
            case 0:
                acb_set_d(a, 0);
                acb_set_d(b, 100);
                acb_calc_integrate(s, f_sin, NULL, a, b, goal, tol, options, prec);
                break;

            case 1:
                acb_set_d(a, 0);
                acb_set_d(b, 1);
                acb_calc_integrate(s, f_atanderiv, NULL, a, b, goal, tol, options, prec);
                acb_mul_2exp_si(s, s, 2);
                break;

            case 2:
                acb_set_d(a, 0);
                acb_one(b);
                acb_mul_2exp_si(b, b, goal);
                acb_calc_integrate(s, f_atanderiv, NULL, a, b, goal, tol, options, prec);
                arb_add_error_2exp_si(acb_realref(s), -goal);
                acb_mul_2exp_si(s, s, 1);
                break;

            case 3:
                acb_set_d(a, 0);
                acb_set_d(b, 1);
                acb_calc_integrate(s, f_circle, NULL, a, b, goal, tol, options, prec);
                acb_mul_2exp_si(s, s, 2);
                break;

            case 4:
                acb_set_d(a, 0);
                acb_set_d(b, 8);
                acb_calc_integrate(s, f_rump, NULL, a, b, goal, tol, options, prec);
                break;

            case 5:
                acb_set_d(a, 1);
                acb_set_d(b, 101);
                acb_calc_integrate(s, f_floor, NULL, a, b, goal, tol, options, prec);
                break;

            case 6:
                acb_set_d(a, 0);
                acb_set_d(b, 1);
                acb_calc_integrate(s, f_helfgott, NULL, a, b, goal, tol, options, prec);
                break;

            case 7:
                acb_zero(s);

                acb_set_d_d(a, -1.0, -1.0);
                acb_set_d_d(b, 2.0, -1.0);
                acb_calc_integrate(t, f_zeta, NULL, a, b, goal, tol, options, prec);
                acb_add(s, s, t, prec);

                acb_set_d_d(a, 2.0, -1.0);
                acb_set_d_d(b, 2.0, 1.0);
                acb_calc_integrate(t, f_zeta, NULL, a, b, goal, tol, options, prec);
                acb_add(s, s, t, prec);

                acb_set_d_d(a, 2.0, 1.0);
                acb_set_d_d(b, -1.0, 1.0);
                acb_calc_integrate(t, f_zeta, NULL, a, b, goal, tol, options, prec);
                acb_add(s, s, t, prec);

                acb_set_d_d(a, -1.0, 1.0);
                acb_set_d_d(b, -1.0, -1.0);
                acb_calc_integrate(t, f_zeta, NULL, a, b, goal, tol, options, prec);
                acb_add(s, s, t, prec);

                acb_const_pi(t, prec);
                acb_div(s, s, t, prec);
                acb_mul_2exp_si(s, s, -1);
                acb_div_onei(s, s);
                break;

            case 8:
                acb_set_d(a, 0);
                acb_set_d(b, 1);
                acb_calc_integrate(s, f_essing, NULL, a, b, goal, tol, options, prec);
                break;

            case 9:
                acb_set_d(a, 0);
                acb_set_d(b, 1);
                acb_calc_integrate(s, f_essing2, NULL, a, b, goal, tol, options, prec);
                break;

            case 10:
                acb_set_d(a, 0);
                acb_set_d(b, 10000);
                acb_calc_integrate(s, f_factorial1000, NULL, a, b, goal, tol, options, prec);
                break;

            case 11:
                acb_set_d_d(a, 1.0, 0.0);
                acb_set_d_d(b, 1.0, 1000.0);
                acb_calc_integrate(s, f_gamma, NULL, a, b, goal, tol, options, prec);
                break;

            case 12:
                acb_set_d(a, -10.0);
                acb_set_d(b, 10.0);
                acb_calc_integrate(s, f_sin_plus_small, NULL, a, b, goal, tol, options, prec);
                break;

            case 13:
                acb_set_d(a, -1020.0);
                acb_set_d(b, -1010.0);
                acb_calc_integrate(s, f_exp, NULL, a, b, goal, tol, options, prec);
                break;

            case 14:
                acb_set_d(a, 0);
                acb_set_d(b, ceil(sqrt(goal * 0.693147181) + 1.0));
                acb_calc_integrate(s, f_gaussian, NULL, a, b, goal, tol, options, prec);
                acb_mul(b, b, b, prec);
                acb_neg(b, b);
                acb_exp(b, b, prec);
                arb_add_error(acb_realref(s), acb_realref(b));
                break;

            case 15:
                acb_set_d(a, 0.0);
                acb_set_d(b, 1.0);
                acb_calc_integrate(s, f_spike, NULL, a, b, goal, tol, options, prec);
                break;

            case 16:
                acb_set_d(a, 0.0);
                acb_set_d(b, 8.0);
                acb_calc_integrate(s, f_monster, NULL, a, b, goal, tol, options, prec);
                break;

            case 17:
                acb_set_d(a, 0);
                acb_set_d(b, ceil(goal * 0.693147181 + 1.0));
                acb_calc_integrate(s, f_sech, NULL, a, b, goal, tol, options, prec);
                acb_neg(b, b);
                acb_exp(b, b, prec);
                acb_mul_2exp_si(b, b, 1);
                arb_add_error(acb_realref(s), acb_realref(b));
                break;

            case 18:
                acb_set_d(a, 0);
                acb_set_d(b, ceil(goal * 0.693147181 / 3.0 + 2.0));
                acb_calc_integrate(s, f_sech3, NULL, a, b, goal, tol, options, prec);
                acb_neg(b, b);
                acb_mul_ui(b, b, 3, prec);
                acb_exp(b, b, prec);
                acb_mul_2exp_si(b, b, 3);
                acb_div_ui(b, b, 3, prec);
                arb_add_error(acb_realref(s), acb_realref(b));
                break;

            case 19:
                if (goal < 0)
                    abort();
                /* error bound 2^-N (1+N) when truncated at 2^-N */
                N = goal + FLINT_BIT_COUNT(goal);
                acb_one(a);
                acb_mul_2exp_si(a, a, -N);
                acb_one(b);
                acb_calc_integrate(s, f_log_div1p, NULL, a, b, goal, tol, options, prec);
                acb_set_ui(b, N + 1);
                acb_mul_2exp_si(b, b, -N);
                arb_add_error(acb_realref(s), acb_realref(b));
                break;

           case 20:
                if (goal < 0)
                    abort();
                /* error bound (N+1) exp(-N) when truncated at N */
                N = goal + FLINT_BIT_COUNT(goal);
                acb_zero(a);
                acb_set_ui(b, N);
                acb_calc_integrate(s, f_log_div1p_transformed, NULL, a, b, goal, tol, options, prec);
                acb_neg(b, b);
                acb_exp(b, b, prec);
                acb_mul_ui(b, b, N + 1, prec);
                arb_add_error(acb_realref(s), acb_realref(b));
                break;

            case 21:

                acb_zero(s);

                N = 10;

                acb_set_d_d(a, 0.5, -0.5);
                acb_set_d_d(b, 0.5, 0.5);
                acb_calc_integrate(t, f_elliptic_p_laurent_n, &N, a, b, goal, tol, options, prec);
                acb_add(s, s, t, prec);

                acb_set_d_d(a, 0.5, 0.5);
                acb_set_d_d(b, -0.5, 0.5);
                acb_calc_integrate(t, f_elliptic_p_laurent_n, &N, a, b, goal, tol, options, prec);
                acb_add(s, s, t, prec);

                acb_set_d_d(a, -0.5, 0.5);
                acb_set_d_d(b, -0.5, -0.5);
                acb_calc_integrate(t, f_elliptic_p_laurent_n, &N, a, b, goal, tol, options, prec);
                acb_add(s, s, t, prec);

                acb_set_d_d(a, -0.5, -0.5);
                acb_set_d_d(b, 0.5, -0.5);
                acb_calc_integrate(t, f_elliptic_p_laurent_n, &N, a, b, goal, tol, options, prec);
                acb_add(s, s, t, prec);

                acb_const_pi(t, prec);
                acb_div(s, s, t, prec);
                acb_mul_2exp_si(s, s, -1);
                acb_div_onei(s, s);
                break;

            case 22:

                acb_zero(s);

                N = 1000;

                acb_set_d_d(a, 100.0, 0.0);
                acb_set_d_d(b, 100.0, N);
                acb_calc_integrate(t, f_zeta_frac, NULL, a, b, goal, tol, options, prec);
                acb_add(s, s, t, prec);

                acb_set_d_d(a, 100, N);
                acb_set_d_d(b, 0.5, N);
                acb_calc_integrate(t, f_zeta_frac, NULL, a, b, goal, tol, options, prec);
                acb_add(s, s, t, prec);

                acb_div_onei(s, s);
                arb_zero(acb_imagref(s));

                acb_set_ui(t, N);
                acb_dirichlet_hardy_theta(t, t, NULL, NULL, 1, prec);
                acb_add(s, s, t, prec);

                acb_const_pi(t, prec);
                acb_div(s, s, t, prec);
                acb_add_ui(s, s, 1, prec);
                break;

            case 23:
                acb_set_d(a, 0.0);
                acb_set_d(b, 1000.0);
                acb_calc_integrate(s, f_lambertw, NULL, a, b, goal, tol, options, prec);
                break;

            case 24:
                acb_set_d(a, 0.0);
                acb_const_pi(b, prec);
                acb_calc_integrate(s, f_max_sin_cos, NULL, a, b, goal, tol, options, prec);
                break;

            case 25:
                acb_set_si(a, -1);
                acb_set_si(b, 1);
                acb_calc_integrate(s, f_erf_bent, NULL, a, b, goal, tol, options, prec);
                break;

            case 26:
                acb_set_si(a, -10);
                acb_set_si(b, 10);
                acb_calc_integrate(s, f_airy_ai, NULL, a, b, goal, tol, options, prec);
                break;

            case 27:
                acb_set_si(a, 0);
                acb_set_si(b, 10);
                acb_calc_integrate(s, f_horror, NULL, a, b, goal, tol, options, prec);
                break;

            case 28:
                acb_set_d_d(a, -1, -1);
                acb_set_d_d(b, -1, 1);
                acb_calc_integrate(s, f_sqrt, NULL, a, b, goal, tol, options, prec);
                break;

            case 29:
                acb_set_d(a, 0);
                acb_set_d(b, ceil(sqrt(goal * 0.693147181) + 1.0));
                acb_calc_integrate(s, f_gaussian_twist, NULL, a, b, goal, tol, options, prec);
                acb_mul(b, b, b, prec);
                acb_neg(b, b);
                acb_exp(b, b, prec);
                arb_add_error(acb_realref(s), acb_realref(b));
                arb_add_error(acb_imagref(s), acb_realref(b));
                break;

            case 30:
                acb_set_d(a, 0);
                acb_set_d(b, ceil(goal * 0.693147181 + 1.0));
                acb_calc_integrate(s, f_exp_airy, NULL, a, b, goal, tol, options, prec);
                acb_neg(b, b);
                acb_exp(b, b, prec);
                acb_mul_2exp_si(b, b, 1);
                arb_add_error(acb_realref(s), acb_realref(b));
                break;

            case 31:
                acb_zero(a);
                acb_const_pi(b, prec);
                acb_calc_integrate(s, f_sin_cos_frac, NULL, a, b, goal, tol, options, prec);
                break;

            case 32:
                acb_zero(a);
                acb_set_ui(b, 3);
                acb_calc_integrate(s, f_sin_near_essing, NULL, a, b, goal, tol, options, prec);
                break;

            case 33:
                acb_zero(a);
                acb_zero(b);
                k = 3;
                scaled_bessel_select_N(acb_realref(b), k, prec);
                acb_calc_integrate(s, f_scaled_bessel, &k, a, b, goal, tol, options, prec);
                scaled_bessel_tail_bound(acb_realref(a), k, acb_realref(b), prec);
                arb_add_error(acb_realref(s), acb_realref(a));
                break;

            case 34:
                acb_zero(a);
                acb_zero(b);
                k = 15;
                scaled_bessel_select_N(acb_realref(b), k, prec);
                acb_calc_integrate(s, f_scaled_bessel, &k, a, b, goal, tol, options, prec);
                scaled_bessel_tail_bound(acb_realref(a), k, acb_realref(b), prec);
                arb_add_error(acb_realref(s), acb_realref(a));
                break;

            case 35:
                acb_set_d_d(a, -1, -1);
                acb_set_d_d(b, -1, 1);
                acb_calc_integrate(s, f_rsqrt, NULL, a, b, goal, tol, options, prec);
                break;

            case 36:
                if (goal < 0)
                    abort();
                acb_zero(a);
                acb_set_ui(b, 4 + (goal + 1) / 2);
                acb_calc_integrate(s, f_rgamma, NULL, a, b, goal, tol, options, prec);
                arb_add_error_2exp_si(acb_realref(s), -goal);
                break;

            default:
                abort();
            }

            TIMEIT_ONCE_STOP
        }
        flint_printf("I%d = ", integral);
        acb_printn(s, 3.333 * prec, 0);
        flint_printf("\n\n");
    }

    acb_clear(a);
    acb_clear(b);
    acb_clear(s);
    acb_clear(t);
    mag_clear(tol);

    flint_cleanup_master();
    return 0;
}

