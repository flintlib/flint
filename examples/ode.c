#include "acb_types.h"
#include "acb_poly.h"
#include "acb_mat.h"
#include "acb_ode.h"
#include "gr.h"
#include "gr_poly.h"


void
_acb_ode_sum_swap_mat_ordinary(
        acb_mat_t mat,
        const acb_ode_sum_context_struct * ctx, slong s)
{
    for (slong j = 0; j < ctx->nsols; j++)
    {
        FLINT_ASSERT(ctx->sol[j].nlogs <= 1);
        if (ctx->sol[j].nlogs < 1)
            continue;
        for (slong i = 0; i < ctx->nder; i++)
            acb_swap(acb_mat_entry(mat, i, j),
                     (ctx->sol[j].sums + s)->coeffs + i);
    }
}

/* TODO
 * - series
 * - singular
 */

void
ordinary(void)
{
    acb_ode_sum_context_t ctx;

    /* (x^2 + 1)*Dx^3 + 7*x
     * = (x^2 + 1)*Tx^3 + (-3*x^2 - 3)*Tx^2 + (2*x^2 + 2)*Tx + 7*x^4 */
    acb_ode_sum_context_init(ctx, 4, 2, 3, 3);
    acb_poly_set_coeff_si(ctx->dop + 3, 2, 1);
    acb_poly_set_coeff_si(ctx->dop + 3, 0, 1);
    acb_poly_set_coeff_si(ctx->dop + 2, 2, -3);
    acb_poly_set_coeff_si(ctx->dop + 2, 0, -3);
    acb_poly_set_coeff_si(ctx->dop + 1, 2, 2);
    acb_poly_set_coeff_si(ctx->dop + 1, 0, 2);
    acb_poly_set_coeff_si(ctx->dop + 0, 4, 7);

    /* ctx->flags |= acb_ode_WANT_SERIES; */

    acb_ode_sum_ordinary(ctx);
    acb_ode_sum_canonical_basis(ctx);

    acb_set_d_d(ctx->pts, 0.25, 0.25);
    acb_set_si(ctx->pts + 1, 0);
    mag_set_d(&((ctx->pts + 1)->real.rad), 1e-8);

    ctx->prec = 64;
    ctx->sums_prec = 64;

    acb_ode_sum_divconquer(ctx, 20);

    acb_mat_t mat;
    acb_mat_init(mat, ctx->nder, ctx->nsols);

    for (slong i = 0; i < ctx->npts; i++)
    {
        _acb_ode_sum_swap_mat_ordinary(mat, ctx, i);
        acb_mat_printd(mat, 8);
        /* flint_printf("%{acb_mat}\n\n", mat); */
    }

    acb_ode_sum_context_clear(ctx);
    acb_mat_clear(mat);
}


void
series(void)
{
    acb_ode_sum_context_t ctx;

    acb_ode_sum_context_init(ctx, 2, 0, 1, 1);
    acb_poly_set_coeff_si(ctx->dop + 1, 0, 1);
    acb_poly_set_coeff_si(ctx->dop + 0, 1, -1);

    acb_ode_sum_ordinary(ctx);
    acb_ode_sum_canonical_basis(ctx);
    ctx->prec = 64;
    ctx->flags |= acb_ode_WANT_SERIES;

    slong len = 5;
    acb_ode_sum_divconquer(ctx, len);
    /* Clear the high part reserved for the residual. (In this special case, the
     * high part is zero because block_length = 1.) */
    acb_poly_truncate(ctx->sol[0].series, len);

    flint_printf("%{acb_poly}\n", ctx->sol[0].series);

    acb_ode_sum_context_clear(ctx);
}


void
bessel_j0(void)
{
    acb_ode_sum_context_t ctx;

    slong dop_order = 2;
    acb_ode_sum_context_init(ctx, dop_order + 1, 1, dop_order, dop_order);
    acb_poly_set_coeff_si(ctx->dop + 2, 0, 1);
    acb_poly_set_coeff_si(ctx->dop + 0, 2, 1);

    ctx->sing_shifts[0].n = 0;
    ctx->sing_shifts[0].mult = 2;

    acb_ode_sum_canonical_basis(ctx);

    acb_set_d(ctx->pts, 0.25);

    ctx->prec = 64;
    ctx->sums_prec = 64;

    acb_ode_sum_divconquer(ctx, 10);

    for (slong m = 0; m < ctx->nsols; m++)
    {
        acb_ode_sol_struct * sol = ctx->sol + m;
        /* series in x */
        flint_printf("f%wd =", m);
        for (slong k = 0; k < sol->nlogs; k++)
        {
            flint_printf(" + (%{acb_poly})*log(%{acb} + x)^%wd/%wd",
                         acb_ode_sol_sum_ptr(sol, 0, k),
                         ctx->pts, k, k);
        }
        flint_printf("\n");
    }

    acb_ode_sum_context_clear(ctx);
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

    acb_ode_sum_context_t ctx;

    acb_set_si(kappa, 2);
    acb_set_si(mu2, 3);
    acb_set_d(half, .5);

    slong prec = 64;
    slong dop_order = 2;
    slong len = dop_order;

    acb_ode_sum_context_init(ctx, dop_order + 1, 1, 1, len);

    acb_poly_set_coeff_si(ctx->dop + 2, 0, 4);
    acb_poly_set_coeff_si(ctx->dop + 1, 0, -4);
    acb_poly_set_coeff_si(ctx->dop + 0, 2, -1);
    acb_mul_si((ctx->dop + 0)->coeffs + 1, kappa, 4, prec);
    acb_poly_set_coeff_si(ctx->dop + 0, 0, 1);
    acb_addmul_si((ctx->dop + 0)->coeffs + 0, mu2, -4, prec);

    acb_sqrt(ctx->expo, mu2, prec);
    acb_add(ctx->expo, ctx->expo, half, prec);  /* other expo = mu2 - 1/2 */

    ctx->sing_shifts[0].n = 0;
    ctx->sing_shifts[0].mult = 1;

    acb_ode_sum_canonical_basis(ctx);

    acb_set_d(ctx->pts, 1.4242);

    ctx->prec = prec;
    ctx->sums_prec = prec;

    acb_ode_sum_divconquer(ctx, prec);

    acb_poly_struct * f = acb_ode_sol_sum_ptr(ctx->sol, 0, 0);

    flint_printf("(%{acb} + x)^(%{acb}) * (%{acb_poly})\n",
                 ctx->pts, ctx->expo, f);

    _acb_ode_sol_value(val, ctx->expo, f, ctx->sol[0].nlogs, ctx->pts,
                       ctx->nder, 1, ctx->prec);

    flint_printf("M(%{acb} + x) = %{acb_poly} + O(x^%wd)\n", ctx->pts, val, len);

    acb_ode_sum_context_clear(ctx);
    acb_clear(kappa);
    acb_clear(mu2);
    acb_clear(half);
    acb_poly_clear(val);
}


void
fundamental_matrix(const char * dop_str,
                   const acb_ode_exponents_struct * expos,
                   double pt_d)
{
    gr_ctx_t CC, Pol, Dop;
    gr_ptr dop;

    slong prec = 30;

    int status = GR_SUCCESS;

    gr_ctx_init_complex_acb(CC, prec);
    gr_ctx_init_gr_poly(Pol, CC);
    gr_ctx_init_gr_poly(Dop, Pol);  /* should be Ore poly */

    GR_TMP_INIT(dop, Dop);

    status |= gr_ctx_set_gen_name(Pol, "z");
    status |= gr_ctx_set_gen_name(Dop, "Tz");
    status |= gr_set_str(dop, dop_str, Dop);
    status |= gr_println(dop, Dop);
    GR_MUST_SUCCEED(status);

    slong dop_order = gr_poly_length(dop, Dop) - 1;

    acb_mat_t mat;
    acb_mat_init(mat, dop_order, dop_order);

    acb_t pt;
    acb_init(pt);
    acb_set_d(pt, pt_d);

    acb_ode_fundamental_matrix(mat, dop, Dop, expos, pt, 1, 0, 8, prec);

    flint_printf("%{acb_mat}\n", mat);

    acb_mat_clear(mat);
    GR_TMP_CLEAR(dop, Dop);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(CC);
    acb_clear(pt);
}


void
apery(void)
{
    acb_ode_shift_struct shift[1] = {{ .n = 0, .mult = 3 }};
    acb_ode_group_struct grp[1] = {{ .nshifts = 1, .shifts = shift }};
    acb_init(grp->expo);
    acb_zero(grp->expo);
    acb_ode_exponents_struct expos[1] = {{ .len = 1, .grps = grp }};

    fundamental_matrix(
            "(z^2 - 34*z + 1)*Tz^3 + (3*z^2 - 51*z)*Tz^2 + (3*z^2 - 27*z)*Tz + z^2 - 5*z",
            expos,
            0.015625);

    acb_clear(grp->expo);
}


void
multiple_shifts(void)
{
    /* const char * dop = "Tz^4 - 4*Tz^3 + 3*Tz^2 - z"; */
    const char * dop = "Tz^6 - 6*Tz^5 + 12*Tz^4 - 10*Tz^3 + 3*Tz^2 + z^2";

    acb_ode_shift_struct shift[3] = {
        { .n = 0, .mult = 2 },
        { .n = 1, .mult = 3 },
        { .n = 3, .mult = 1 },
    };
    acb_ode_group_struct grp[1] = {{ .nshifts = 3, .shifts = shift }};
    acb_init(grp->expo);
    acb_zero(grp->expo);
    acb_ode_exponents_struct expos[1] = {{ .len = 1, .grps = grp }};

    fundamental_matrix(dop, expos, 2.);

    acb_clear(grp->expo);
}


int
main(void)
{
    /* ordinary(); */
    /* series(); */
    /* bessel_j0(); */
    /* whittaker_m(); */

    /* apery(); */
    multiple_shifts();

    flint_cleanup_master();
}
