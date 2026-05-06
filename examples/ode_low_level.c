#include "acb_types.h"
#include "acb_poly.h"
#include "acb_mat.h"
#include "acb_ode.h"


// XXX get rid of this file?

void
_acb_ode_sum_swap_mat_ordinary(
        acb_mat_t mat,
        const acb_ode_sum_struct * sum, slong s)
{
    for (slong j = 0; j < sum->nsols; j++)
    {
        FLINT_ASSERT(sum->sol[j].nlogs <= 1);
        if (sum->sol[j].nlogs < 1)
            continue;
        for (slong i = 0; i < sum->nder; i++)
            acb_swap(acb_mat_entry(mat, i, j),
                     acb_ode_sol_sum_ptr(sum->sol + j, s, 0, i));
    }
}


void
ordinary(void)
{
    acb_ode_sum_t sum;

    /* (x^2 + 1)*Dx^3 + 7*x
     * = (x^2 + 1)*Tx^3 + (-3*x^2 - 3)*Tx^2 + (2*x^2 + 2)*Tx + 7*x^4 */
    acb_ode_sum_init(sum, 4, 2, 3, 3);
    acb_poly_set_coeff_si(sum->dop + 3, 2, 1);
    acb_poly_set_coeff_si(sum->dop + 3, 0, 1);
    acb_poly_set_coeff_si(sum->dop + 2, 2, -3);
    acb_poly_set_coeff_si(sum->dop + 2, 0, -3);
    acb_poly_set_coeff_si(sum->dop + 1, 2, 2);
    acb_poly_set_coeff_si(sum->dop + 1, 0, 2);
    acb_poly_set_coeff_si(sum->dop + 0, 4, 7);

    /* sum->flags |= ACB_ODE_WANT_SERIES; */
    sum->flags |= ACB_ODE_APPROX;

    acb_ode_sum_set_ordinary(sum);
    acb_ode_sum_set_ini_echelon(sum);

    acb_set_d_d(sum->pts, 0.25, 0.25);
    acb_set_si(sum->pts + 1, 0);
    mag_set_d(&((sum->pts + 1)->real.rad), 1e-8);

    acb_ode_sum_divconquer(sum, 20, 64);

    acb_mat_t mat;
    acb_mat_init(mat, sum->nder, sum->nsols);

    for (slong i = 0; i < sum->npts; i++)
    {
        _acb_ode_sum_swap_mat_ordinary(mat, sum, i);
        acb_mat_printd(mat, 8);
        /* flint_printf("%{acb_mat}\n\n", mat); */
    }

    acb_ode_sum_clear(sum);
    acb_mat_clear(mat);
}


void
series(void)
{
    acb_ode_sum_t sum;

    acb_ode_sum_init(sum, 2, 0, 1, 1);
    acb_poly_set_coeff_si(sum->dop + 1, 0, 1);
    acb_poly_set_coeff_si(sum->dop + 0, 1, -1);

    acb_ode_sum_set_ordinary(sum);
    acb_ode_sum_set_ini_echelon(sum);
    sum->flags |= ACB_ODE_WANT_SERIES;

    slong len = 5;
    acb_ode_sum_divconquer(sum, len, 64);
    /* Clear the high part reserved for the residual. (In this special case, the
     * high part is zero because block_length = 1.) */
    acb_poly_truncate(sum->sol[0].series, len);

    flint_printf("%{acb_poly}\n", sum->sol[0].series);

    acb_ode_sum_clear(sum);
}


void
bessel_j0(void)
{
    acb_ode_sum_t sum;

    slong dop_order = 2;
    acb_ode_sum_init(sum, dop_order + 1, 1, dop_order, dop_order);
    acb_poly_set_coeff_si(sum->dop + 2, 0, 1);
    acb_poly_set_coeff_si(sum->dop + 0, 2, 1);

    sum->group->shifts[0].n = 0;
    sum->group->shifts[0].mult = 2;

    acb_ode_sum_set_ini_echelon(sum);

    acb_set_d(sum->pts, 0.25);

    acb_ode_sum_divconquer(sum, 10, 64);

    for (slong m = 0; m < sum->nsols; m++)
    {
        acb_ode_sol_struct * sol = sum->sol + m;
        /* series in x */
        flint_printf("f%wd =", m);
        for (slong k = 0; k < sol->nlogs; k++)
        {
            flint_printf(" + (%{acb*})*log(%{acb} + x)^%wd/%wd!",
                         acb_ode_sol_sum_ptr(sol, 0, k, 0), sum->nder,
                         sum->pts, k, k);
        }
        flint_printf("\n");
    }

    acb_ode_sum_clear(sum);
}


void
whittaker_m(void)  /* non-integer exponent, no logs */
{
    acb_t kappa, mu2, half;
    acb_poly_t val;

    acb_init(kappa);
    acb_init(mu2);
    acb_init(half);
    acb_poly_init(val);

    acb_ode_sum_t sum;

    acb_set_si(kappa, 2);
    acb_set_si(mu2, 3);
    acb_set_d(half, .5);

    slong prec = 64;
    slong dop_order = 2;
    slong len = dop_order;

    acb_ode_sum_init(sum, dop_order + 1, 1, 1, len);

    acb_poly_set_coeff_si(sum->dop + 2, 0, 4);
    acb_poly_set_coeff_si(sum->dop + 1, 0, -4);
    acb_poly_set_coeff_si(sum->dop + 0, 2, -1);
    acb_mul_si((sum->dop + 0)->coeffs + 1, kappa, 4, prec);
    acb_poly_set_coeff_si(sum->dop + 0, 0, 1);
    acb_addmul_si((sum->dop + 0)->coeffs + 0, mu2, -4, prec);

    acb_sqrt(sum->group->leader, mu2, prec);
    acb_add(sum->group->leader, sum->group->leader, half, prec);  /* other leader = 1/2 - mu */

    sum->group->shifts[0].n = 0;
    sum->group->shifts[0].mult = 1;

    acb_ode_sum_set_ini_echelon(sum);

    acb_set_d(sum->pts, 1.4242);

    acb_ode_sum_divconquer(sum, -1, prec);

    flint_printf("(%{acb} + x)^(%{acb}) * (%{acb*})\n", sum->pts, sum->group->leader,
                 sum->sol->sums, sum->nder);

    acb_ode_sol_jet(val, sum->group->leader, sum->sol, 0, sum->pts, sum->nder, 1, prec);

    flint_printf("M(%{acb} + x) = %{acb_poly} + O(x^%wd)\n", sum->pts, val, len);

    acb_ode_sum_clear(sum);
    acb_clear(kappa);
    acb_clear(mu2);
    acb_clear(half);
    acb_poly_clear(val);
}


int
main(void)
{
    ordinary();
    series();
    bessel_j0();
    whittaker_m();

    flint_cleanup_master();
}
