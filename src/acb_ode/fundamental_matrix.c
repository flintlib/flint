#include "acb_types.h"
#include "acb_poly.h"
#include "acb_mat.h"
#include "acb_ode.h"

#include "gr.h"
#include "gr_poly.h"


static void
init_context(acb_ode_sum_context_t ctx,
             const acb_poly_struct * dop, slong dop_len,
             const acb_ode_group_struct * grp,
             acb_srcptr pts, slong npts, slong nder,
             slong prec)
{
    acb_ode_sum_context_init(ctx, dop_len, npts, grp->nshifts, nder);
    _acb_poly_vec_set(ctx->dop, dop, dop_len);
    acb_ode_sum_group(ctx, grp);
    _acb_vec_set(ctx->pts, pts, npts);
    ctx->prec = prec;        /* XXX */
    ctx->sums_prec = prec;
}


static slong
col(const acb_ode_group_struct * grp, slong s, slong k)
{
    slong j = 0;
    for (slong s1 = 0; s1 <= s; s1++)
        j += grp->shifts[s1].mult;
    return j - 1 - k;
}


static void
fill_column(acb_mat_struct * mat,
            const acb_ode_group_struct * grp, slong s, slong k,
            const acb_poly_struct * val)
{
    slong j = col(grp, s, k);
    slong len = FLINT_MIN(acb_mat_nrows(mat), val->length);
    for (slong i = 0; i < len; i++)
    {
        acb_ptr a = acb_mat_entry(mat, i, j);
        acb_swap(a, val->coeffs + i);
    }
}


static void
fix_column_echelon(acb_mat_struct * mat,
                   const acb_ode_group_struct * grp, slong s, slong k,
                   const acb_mat_t extini, slong prec)
{
    slong mult = grp->shifts[s].mult;
    slong delta = mult - 1 - k;
    slong j = col(grp, s, k);

    for (slong s1 = s + 1, j1 = col(grp, s1, 0); s1 < grp->nshifts; s1++)
    {
        slong mult1 =  grp->shifts[s1].mult;
        for (slong k1 = FLINT_MAX(0, mult1 - delta); k1 < mult1; k1++, j1++)
        {
            acb_ptr cc = acb_mat_entry(extini, s1, k1 + delta);
            /* flint_printf("fix: s=%wd k=%wd s1=%wd k1=%wd j=%wd j1=%d cc=%{acb}\n", s, k, s1, k1, j, j1, cc); */
            for (slong i = 0; i < acb_mat_nrows(mat); i++)
            {
                acb_ptr a  = acb_mat_entry(mat, i, j);
                acb_ptr a1 = acb_mat_entry(mat, i, j1);
                /* flint_printf("i=%wd %{acb} - %{acb}", i, a, a1); */
                acb_submul(a, cc, a1, prec);
                /* flint_printf(" = %{acb}\n", a); */
            }
        }
    }

}


static void
fill_group(acb_mat_struct * mat, const acb_ode_sum_context_t ctx,
           const acb_ode_group_struct * grp, slong p,
           acb_ode_basis_t basis, slong prec)
{
    for (slong s = grp->nshifts - 1; s >= 0; s--)
    {
        acb_poly_struct * val = _acb_poly_vec_init(ctx->sol[s].nlogs);

        slong mult = grp->shifts[s].mult;
        acb_ode_sum_value(val, mult, ctx, s, p, prec);
        /* flint_printf("s=%wd mult=%wd val=%{acb_poly}\n\n", s, mult, val); */

        for (slong k = 0; k < mult; k++)
        {
            slong delta = mult - 1 - k;
            fill_column(mat, grp, s, k, val + delta);

            /* flint_printf("s=%wd k=%wd after fill \n%{acb_mat}\n\n", s, k, mat); */

            switch (basis)
            {
                case acb_ode_BASIS_CASCADE:
                    break;
                case acb_ode_BASIS_ECHELON:
                    fix_column_echelon(mat, grp, s, k, ctx->sol[s].extini, prec);
                    break;
                default:
                    FLINT_ASSERT(0);
            }

            /* flint_printf("s=%wd k=%wd after fix \n%{acb_mat}\n\n", s, k, mat); */
        }

        _acb_poly_vec_clear(val, mult);
    }
}


void
_acb_ode_fundamental_matrix(
        acb_mat_struct * mat,
        const acb_poly_struct * dop, slong dop_len,
        const acb_ode_group_struct * groups, slong ngroups,
        acb_srcptr pts, slong npts,
        acb_ode_basis_t basis,
        acb_ode_sum_worker_t sum_worker,
        slong nterms, slong prec)
{
    slong nder = acb_mat_nrows(mat);

    for (slong g = 0, j = 0; g < ngroups; g++)
    {
        slong glen = acb_ode_group_length(groups + g);

        acb_ode_sum_context_t ctx;
        init_context(ctx, dop, dop_len, groups + g, pts, npts, nder, prec);
        /* XXX tmp */
        /* ctx->flags |= acb_ode_WANT_SERIES; */

        sum_worker(ctx, nterms);

        for (slong p = 0; p < npts; p++)
        {
            acb_mat_t win;
            acb_mat_window_init(win, mat + p, j, 0, j + glen, acb_mat_nrows(mat));

            fill_group(win, ctx, groups + g, p, basis, prec);

            acb_mat_window_clear(win);
        }

        acb_ode_sum_context_clear(ctx);

        j += glen;
    }
}


int
acb_ode_fundamental_matrix(
        acb_mat_struct * mat,
        gr_srcptr dop, gr_ctx_t dop_ctx,
        const acb_ode_exponents_struct * expos,
        acb_srcptr pts, slong npts,                   /* XXX gr_vec? */
        acb_ode_basis_t basis,
        slong nterms, slong prec)
{
    gr_ctx_t CC, Pol, Dop;
    gr_ptr _dop;

    if (dop_ctx->which_ring != GR_CTX_GR_POLY)
        return GR_UNABLE;

    ulong which_base = POLYNOMIAL_ELEM_CTX(dop_ctx)->which_ring;
    if (!(which_base == GR_CTX_GR_POLY || which_base == GR_CTX_FMPZ_POLY ||
          which_base == GR_CTX_FMPQ_POLY))
        return GR_UNABLE;

    int status = GR_SUCCESS;

    gr_ctx_init_complex_acb(CC, prec);
    gr_ctx_init_gr_poly(Pol, CC);
    gr_ctx_init_gr_poly(Dop, Pol);  /* XXX should be Ore poly */

    status |= gr_ctx_set_gen_name(Dop, POLYNOMIAL_CTX(dop_ctx)->var);
    status |= gr_ctx_set_gen_name(
            Pol, POLYNOMIAL_CTX(POLYNOMIAL_ELEM_CTX(dop_ctx))->var);

    GR_TMP_INIT(_dop, Dop);
    status |= gr_set_other(_dop, dop, dop_ctx, Dop);

    GR_MUST_SUCCEED(status);

    gr_ptr dop_coeff = gr_poly_entry_ptr(_dop, 0, Dop);
    slong dop_len = gr_poly_length(_dop, Dop);

    _acb_ode_fundamental_matrix(mat, dop_coeff, dop_len,
                                expos->grps, expos->len,
                                pts, npts,
                                basis, acb_ode_sum_divconquer,
                                nterms, prec);

    GR_TMP_CLEAR(_dop, Dop);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(CC);

    return status;
}
