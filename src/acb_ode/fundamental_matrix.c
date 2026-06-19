#include "acb_types.h"
#include "acb_poly.h"
#include "acb_mat.h"
#include "acb_ode.h"
#include "fmpz_vec.h"
#include "gr.h"
#include "gr_poly.h"
#include "gr_ore_poly.h"
#include "gr_vec.h"


static slong
col(const acb_ode_group_struct * grp, slong s, slong k)
{
    slong j = 0;
    for (slong s1 = 0; s1 <= s; s1++)
        j += grp->shifts[s1].mult;
    return j - 1 - k;
}


static void
fill_column(acb_mat_t mat,
            const acb_ode_group_struct * grp, slong s, slong k,
            const acb_poly_struct * val)
{
    slong j = col(grp, s, k);
    slong len = FLINT_MIN(acb_mat_nrows(mat), val->length);

    slong i = 0;
    for (; i < len; i++)
    {
        acb_ptr a = acb_mat_entry(mat, i, j);
        acb_swap(a, val->coeffs + i);
    }
    for (; i < acb_mat_nrows(mat); i++)
        acb_zero(acb_mat_entry(mat, i, j));
}


static void
fix_column_echelon(acb_mat_struct * mat,
                   const acb_ode_group_struct * grp, slong s, slong k,
                   const acb_mat_t extini, slong prec)
{
    slong mult = grp->shifts[s].mult;
    slong delta = mult - 1 - k;
    slong j = col(grp, s, k);

    for (slong s1 = s + 1; s1 < grp->nshifts; s1++)
    {
        slong j1 = col(grp, s1, 0);
        slong mult1 =  grp->shifts[s1].mult;
        for (slong k1 = FLINT_MAX(0, mult1 - delta); k1 < mult1; k1++)
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
            j1--;
        }
    }

}


static void
fill_group(acb_mat_t mat, const acb_ode_sum_t sum,
           const acb_ode_group_struct * grp, slong p,
           acb_ode_basis_t basis, slong prec)
{
    for (slong s = grp->nshifts - 1; s >= 0; s--)
    {
        acb_ode_sol_struct * sol = sum->sol + s;
        slong mult = grp->shifts[s].mult;

        acb_poly_struct * val = _acb_poly_vec_init(mult);

        acb_ode_sol_jet(val, sum->group->leader, sol, p, sum->pts + p,
                        acb_mat_nrows(mat), mult, prec);
        // flint_printf("s=%wd mult=%wd val=%{acb_poly}\n\n", s, mult, val);

        for (slong k = 0; k < mult; k++)
        {
            slong delta = mult - 1 - k;
            fill_column(mat, grp, s, k, val + delta);

            // flint_printf("s=%wd k=%wd after fill \n%{acb_mat}\n\n", s, k, mat);

            switch (basis)
            {
                case ACB_ODE_BASIS_CASCADE:
                    break;
                case ACB_ODE_BASIS_ECHELON:
                    fix_column_echelon(mat, grp, s, k, sum->sol[s].extini, prec);
                    break;
                default:
                    FLINT_ASSERT(0);
            }

            // flint_printf("s=%wd k=%wd after fix \n%{acb_mat}\n\n", s, k, mat);
        }

        _acb_poly_vec_clear(val, mult);
    }
}


void
_acb_ode_fundamental_matrix_vec(
        acb_mat_struct * mat,
        const acb_poly_struct * dop, slong dop_len,
        const acb_ode_exponents_t expos,
        acb_srcptr lcrt,
        acb_ode_bound_t bound,  /* mutable */
        acb_srcptr pts, slong npts,
        acb_ode_basis_t basis,
        acb_ode_sum_worker_t sum_worker,
        void * worker_data,
        slong prec)
{
    mag_t cvrad, mag;

    mag_init(cvrad);
    mag_init(mag);

    slong nder = 0;
    for (slong p = 0; p < npts; p++)
        nder = FLINT_MAX(nder, acb_mat_nrows(mat + p));

    if (nder == 0)
        return;

    _acb_vec_get_mag_lower(cvrad, lcrt, dop[dop_len - 1].length - 1);

    /* XXX maybe not the best place to test this;
       redundant with sum_set_points */
    _acb_vec_get_mag(mag, pts, npts);
    flint_printf("mag=%{mag} cvrad=%{mag}\n", mag, cvrad);
    if (mag_cmp(mag, cvrad) >= 0)
    {
        for (slong p = 0; p < npts; p++)
            acb_mat_indeterminate(mat + p);
        goto cleanup;
    }

    for (slong g = 0, j = 0; g < expos->ngroups; g++)
    {
        acb_ode_sum_t sum;
        acb_ode_group_bound_t gbound;

        acb_ode_group_struct * group = expos->groups + g;
        slong glen = acb_ode_group_length(expos->groups + g);

        flint_printf("GROUP #%ld leader=%{acb} shifts+muls=%{slong*} glen=%ld\n",
                     g, group->leader, group->shifts, 2*group->nshifts, glen);

        acb_ode_group_bound_init(gbound);
        acb_ode_group_bound_precompute(gbound, dop, dop_len, expos, g, bound, bound->prec);

        acb_ode_sum_init(sum, dop_len, npts, group->nshifts, nder);
        acb_ode_sum_set_diffop(sum, dop, dop_len, cvrad);
        acb_ode_sum_set_group(sum, group);
        /* todo: only avoids part of the redundant computations between
           solutions of the same group in the presence of logs */
        acb_ode_sum_set_ini_highest(sum);
        acb_ode_sum_set_points(sum, pts, npts);
        sum->data = worker_data;

        sum_worker(sum, -1, bound, gbound, prec);

        for (slong p = 0; p < npts; p++)
        {
            acb_mat_t win;
            acb_mat_window_init(win, mat + p, 0, j, acb_mat_nrows(mat + p), j + glen);

            fill_group(win, sum, expos->groups + g, p, basis, prec);

            flint_printf("g=%ld p=%ld win=%{acb_mat}\n", g, p, win);

            acb_mat_window_clear(win);
        }

        acb_ode_sum_clear(sum);
        acb_ode_group_bound_clear(gbound);

        j += glen;
    }

cleanup:
    mag_clear(cvrad);
    mag_clear(mag);
}


static slong
x_valuation(const gr_ore_poly_t dop, gr_ctx_t Cst)
{
    slong val = WORD_MAX;

    for (slong i = 0; val > 0 && i < dop->length; i++)
    {
        gr_poly_struct * pol = ((gr_poly_struct *) (dop->coeffs)) + i;
        slong v = 0;
        while (v < val && v < pol->length
               && gr_is_zero(gr_poly_coeff_ptr(pol, v, Cst), Cst) == T_TRUE)
            v++;
        if (v < pol->length)
            val = v;
    }

    return val;
}


int
acb_ode_fundamental_matrix_vec(
        acb_mat_struct * mat,
        const gr_ore_poly_t dop, gr_ore_poly_ctx_t Dop,
        const acb_ode_exponents_t expos,
        acb_srcptr lcrt,
        acb_srcptr pts, slong npts,
        acb_ode_basis_t basis,
        slong prec)
{
    gr_ctx_struct * Scalars, * Pol;
    gr_ctx_t CC, CCx, CCxT, bCC;
    acb_ode_exponents_t _expos;
    gr_ore_poly_struct * _dop;
    gr_vec_t sing;
    acb_struct * _lcrt = NULL;
    acb_ode_bound_t bound;

    if (gr_ore_poly_ctx_over_gr_poly_base_ptrs(&Scalars, &Pol, Dop) != GR_SUCCESS)
    {
        flint_printf("%s:%d UNABLE\n", __FILE__, __LINE__);
        return GR_UNABLE;
    }

    /* todo: support d/dx (handle lcrt) */
    if (GR_ORE_POLY_CTX(Dop)->which_algebra != ORE_ALGEBRA_EULER_DERIVATIVE)
    {
        flint_printf("%s:%d UNABLE\n", __FILE__, __LINE__);
        return GR_UNABLE;
    }

    if (dop->length == 0)
        return GR_DOMAIN;
    if (dop->length == 1)
        return GR_SUCCESS;

    int is_zero = gr_ore_poly_is_zero(dop, Dop);
    if (is_zero != T_FALSE)
        return gr_check(truth_not(is_zero));

    gr_poly_struct * lc = (gr_poly_struct *) (dop->coeffs) + dop->length - 1;
    slong xval = x_valuation(dop, Scalars);

    int is_regular = xval < lc->length
        && truth_not(gr_is_zero(gr_poly_coeff_ptr(lc, xval, Scalars), Scalars));
    if (is_regular != T_TRUE)
        return gr_check(is_regular);
    if (gr_is_zero(lc->coeffs, Scalars) != T_FALSE)
        /* we could handle some cases if we were more careful with lcroots */
        return GR_UNABLE;

    gr_ctx_init_complex_acb(bCC, MAG_BITS);  // XXX prec choice?
    gr_vec_init(sing, 0, bCC);

    int status = GR_SUCCESS;

    gr_ctx_init_complex_acb(CC, prec);
    gr_ctx_init_gr_poly(CCx, CC);
    gr_ore_poly_ctx_init(CCxT, CCx, 0, ORE_ALGEBRA_EULER_DERIVATIVE);

    GR_TMP_INIT(_dop, CCxT);

    char * gen = NULL;
    status |= gr_ctx_gen_name(&gen, 0, Pol);
    status |= gr_ctx_set_gen_name(CCx, gen);
    if(status) flint_printf("%s:%d status=%d\n", __FILE__, __LINE__, status);
    flint_free(gen);

    status |= gr_ctx_set_gen_name(CCxT, GR_ORE_POLY_CTX(Dop)->var);
    if(status) flint_printf("%s:%d status=%d\n", __FILE__, __LINE__, status);

    /* XXX reorganize to do this at the *working* precision */
    /* todo: use gr_set_other once its gr_ore_poly version is powerful enough */
    status |= gr_poly_set_gr_poly_other((gr_ptr) _dop, (gr_ptr) dop, Pol, CCx);
    if(status) flint_printf("%s:%d status=%d\n", __FILE__, __LINE__, status);

    if (expos == NULL)
    {
        acb_ode_exponents_init(_expos);
        status |= acb_ode_exponents(_expos, dop, Dop, CC);
        if(status) flint_printf("%s:%d status=%d\n", __FILE__, __LINE__, status);
        expos = _expos;
    }

    if (lcrt == NULL)
    {
        fmpz_vec_t sing_mult;
        fmpz_vec_init(sing_mult, 0);

        /* todo: should use (multiple) root enclosures supporting polynomials
           with interval coefficients */
        status |= gr_poly_roots_other(sing, sing_mult, lc,
                                      Scalars, 0, bCC);
        if(status) flint_printf("%s:%d status=%d\n", __FILE__, __LINE__, status);

        _lcrt = flint_malloc(sizeof(acb_struct) * (lc->length - 1));
        for (slong j0 = 0; j0 < sing->length; j0++)
        {
            slong m = fmpz_get_si((fmpz *) (sing_mult->entries) + j0);
            for (slong j1 = 0; j1 < m; j1++)
                _lcrt[j0 + j1] = ((acb_ptr) sing->entries)[j0];
        }
        lcrt = _lcrt;

        fmpz_vec_clear(sing_mult);
    }

    if (status != GR_SUCCESS)
        goto cleanup;

    acb_ode_bound_init(bound);
    slong dop_clen = _acb_poly_vec_length(_dop->coeffs, _dop->length);
    slong pol_part_len = dop_clen / 2 + 3;
    acb_ode_bound_precompute(bound, _dop->coeffs, _dop->length, lcrt,
                             pol_part_len, MAG_BITS);

    _acb_ode_fundamental_matrix_vec(mat, _dop->coeffs, _dop->length,
                                    expos, lcrt, bound,
                                    pts, npts,
                                    basis,
                                    acb_ode_sum_divconquer, NULL,
                                    prec);

    acb_ode_bound_clear(bound);

cleanup:
    flint_free(_lcrt);
    gr_vec_clear(sing, bCC);
    if (expos == _expos)
        acb_ode_exponents_clear(_expos);
    GR_TMP_CLEAR(_dop, CCxT);
    gr_ctx_clear(bCC);
    gr_ctx_clear(CCxT);
    gr_ctx_clear(CCx);
    gr_ctx_clear(CC);

    return status;
}


int
acb_ode_fundamental_matrix(
        acb_mat_t mat,
        const gr_ore_poly_t dop, gr_ore_poly_ctx_t Dop,
        const acb_ode_exponents_t expos,
        acb_srcptr lcrt,
        const acb_t pt,
        acb_ode_basis_t basis,
        slong prec)
{
    return acb_ode_fundamental_matrix_vec(mat, dop, Dop, expos, lcrt, pt, 1,
                                          basis, prec);
}
