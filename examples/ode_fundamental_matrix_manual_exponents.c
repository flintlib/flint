#include "acb_types.h"
#include "acb.h"
#include "acb_mat.h"
#include "acb_ode.h"
#include "gr.h"
#include "gr_ore_poly.h"

// XXX also manual sing?

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
    gr_ctx_init_gr_ore_poly(Dop, Pol, 0, ORE_ALGEBRA_EULER_DERIVATIVE);

    GR_TMP_INIT(dop, Dop);

    status |= gr_ctx_set_gen_name(Pol, "z");
    status |= gr_ctx_set_gen_name(Dop, "Tz");
    status |= gr_ore_poly_set_str(dop, dop_str, Dop);

    status |= gr_println(dop, Dop);
    flint_printf("\n");

    GR_MUST_SUCCEED(status);

    slong dop_order = gr_ore_poly_length(dop, Dop) - 1;

    acb_mat_t mat;
    acb_mat_init(mat, dop_order, dop_order);

    acb_t pt;
    acb_init(pt);
    acb_set_d(pt, pt_d);

    GR_MUST_SUCCEED(acb_ode_fundamental_matrix(mat, dop, Dop, expos, NULL, pt, 0, prec));

    flint_printf("%{acb_mat}\n\n", mat);
    flint_printf("--------\n\n");

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
    acb_init(grp->leader);
    acb_zero(grp->leader);
    acb_ode_exponents_struct expos[1] = {{ .ngroups = 1, .groups = grp }};

    fundamental_matrix(
            "(z^2 - 34*z + 1)*Tz^3 + (3*z^2 - 51*z)*Tz^2 + (3*z^2 - 27*z)*Tz + z^2 - 5*z",
            expos,
            0.015625);

    acb_clear(grp->leader);
}


void
multiple_shifts(void)
{
    const char * dop = "Tz^6 - 6*Tz^5 + 12*Tz^4 - 10*Tz^3 + 3*Tz^2 + z^2";

    acb_ode_shift_struct shift[3] = {
        { .n = 0, .mult = 2 },
        { .n = 1, .mult = 3 },
        { .n = 3, .mult = 1 },
    };
    acb_ode_group_struct grp[1] = {{ .nshifts = 3, .shifts = shift }};
    acb_init(grp->leader);
    acb_zero(grp->leader);
    acb_ode_exponents_struct expos[1] = {{ .ngroups = 1, .groups = grp }};

    fundamental_matrix(dop, expos, 2.);

    acb_clear(grp->leader);
}


void
whittaker(void)
{
    slong prec = 30;

    acb_t kappa, mu, half;
    acb_init(kappa);
    acb_init(mu);
    acb_init(half);

    const char * dop = "4*Tz^2 - 4*Tz - z^2 + 8*z - 11";
    acb_set_si(kappa, 2);
    acb_set_si(mu, 3);
    acb_sqrt(mu, mu, prec);
    acb_set_d(half, .5);

    acb_ode_shift_struct shift[1] = { { .n = 0, .mult = 1 } };
    acb_ode_group_struct grp[2] = {
        { .nshifts = 1, .shifts = shift },
        { .nshifts = 1, .shifts = shift },
    };
    acb_init(grp[0].leader);
    acb_sub(grp[0].leader, half, mu, prec);
    acb_init(grp[1].leader);
    acb_add(grp[1].leader, half, mu, prec);
    acb_ode_exponents_struct expos[1] = {{ .ngroups = 2, .groups = grp }};

    fundamental_matrix(dop, expos, 2.);

    acb_clear(grp[0].leader);
    acb_clear(grp[1].leader);
    acb_clear(kappa);
    acb_clear(mu);
    acb_clear(half);
}


int
main(void)
{
    apery();
    multiple_shifts();
    whittaker();

    flint_cleanup_master();
}
