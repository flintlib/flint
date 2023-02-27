/* This file is public domain. Author: Fredrik Johansson. */

#include "acb_calc.h"
#include "flint/profiler.h"

int
sinx(acb_ptr out, const acb_t inp, void * params, slong order, slong prec)
{
    int xlen = FLINT_MIN(2, order);
    acb_set(out, inp);
    if (xlen > 1)
        acb_one(out + 1);
    _acb_poly_sin_series(out, out, xlen, order, prec);
    return 0;
}

int
elliptic(acb_ptr out, const acb_t inp, void * params, slong order, slong prec)
{
    acb_ptr t;
    t = _acb_vec_init(order);
    acb_set(t, inp);
    if (order > 1)
        acb_one(t + 1);
    _acb_poly_sin_series(t, t, FLINT_MIN(2, order), order, prec);
    _acb_poly_mullow(out, t, order, t, order, order, prec);
    _acb_vec_scalar_mul_2exp_si(t, out, order, -1);
    acb_sub_ui(t, t, 1, prec);
    _acb_vec_neg(t, t, order);
    _acb_poly_rsqrt_series(out, t, order, order, prec);
    _acb_vec_clear(t, order);
    return 0;
}

int
bessel(acb_ptr out, const acb_t inp, void * params, slong order, slong prec)
{
    acb_ptr t;
    acb_t z;
    ulong n;

    t = _acb_vec_init(order);
    acb_init(z);

    acb_set(t, inp);
    if (order > 1)
        acb_one(t + 1);

    n = 10;
    arb_set_si(acb_realref(z), 20);
    arb_set_si(acb_imagref(z), 10);

    /* z sin(t) */
    _acb_poly_sin_series(out, t, FLINT_MIN(2, order), order, prec);
    _acb_vec_scalar_mul(out, out, order, z, prec);

    /* t n */
    _acb_vec_scalar_mul_ui(t, t, FLINT_MIN(2, order), n, prec);

    _acb_poly_sub(out, t, FLINT_MIN(2, order), out, order, prec);

    _acb_poly_cos_series(out, out, order, order, prec);

    _acb_vec_clear(t, order);
    acb_clear(z);
    return 0;
}

int main(int argc, char *argv[])
{
    acb_t r, s, a, b;
    arf_t inr, outr;
    slong digits, prec;

    if (argc < 2)
    {
        flint_printf("integrals d\n");
        flint_printf("compute integrals using d decimal digits of precision\n");
        return 1;
    }

    acb_init(r);
    acb_init(s);
    acb_init(a);
    acb_init(b);
    arf_init(inr);
    arf_init(outr);

    arb_calc_verbose = 0;

    digits = atol(argv[1]);
    prec = digits * 3.32193;
    flint_printf("Digits: %wd\n", digits);

    flint_printf("----------------------------------------------------------------\n");
    flint_printf("Integral of sin(t) from 0 to 100.\n");
    arf_set_d(inr, 0.125);
    arf_set_d(outr, 1.0);
    TIMEIT_ONCE_START
    acb_set_si(a, 0);
    acb_set_si(b, 100);
    acb_calc_integrate_taylor(r, sinx, NULL, a, b, inr, outr, prec, 1.1 * prec);
    flint_printf("RESULT:\n");
    acb_printn(r, digits, 0); flint_printf("\n");
    TIMEIT_ONCE_STOP

    flint_printf("----------------------------------------------------------------\n");
    flint_printf("Elliptic integral F(phi, m) = integral of 1/sqrt(1 - m*sin(t)^2)\n");
    flint_printf("from 0 to phi, with phi = 6+6i, m = 1/2. Integration path\n");
    flint_printf("0 -> 6 -> 6+6i.\n");
    arf_set_d(inr, 0.2);
    arf_set_d(outr, 0.5);
    TIMEIT_ONCE_START
    acb_set_si(a, 0);
    acb_set_si(b, 6);
    acb_calc_integrate_taylor(r, elliptic, NULL, a, b, inr, outr, prec, 1.1 * prec);
    acb_set_si(a, 6);
    arb_set_si(acb_realref(b), 6);
    arb_set_si(acb_imagref(b), 6);
    acb_calc_integrate_taylor(s, elliptic, NULL, a, b, inr, outr, prec, 1.1 * prec);
    acb_add(r, r, s, prec);
    flint_printf("RESULT:\n");
    acb_printn(r, digits, 0); flint_printf("\n");
    TIMEIT_ONCE_STOP

    flint_printf("----------------------------------------------------------------\n");
    flint_printf("Bessel function J_n(z) = (1/pi) * integral of cos(t*n - z*sin(t))\n");
    flint_printf("from 0 to pi. With n = 10, z = 20 + 10i.\n");
    arf_set_d(inr, 0.1);
    arf_set_d(outr, 0.5);
    TIMEIT_ONCE_START
    acb_set_si(a, 0);
    acb_const_pi(b, 3 * prec);
    acb_calc_integrate_taylor(r, bessel, NULL, a, b, inr, outr, prec, 3 * prec);
    acb_div(r, r, b, prec);
    flint_printf("RESULT:\n");
    acb_printn(r, digits, 0); flint_printf("\n");
    TIMEIT_ONCE_STOP

    acb_clear(r);
    acb_clear(s);
    acb_clear(a);
    acb_clear(b);
    arf_clear(inr);
    arf_clear(outr);

    flint_cleanup();
    return 0;
}

