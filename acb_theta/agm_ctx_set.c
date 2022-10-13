
#include "acb_theta.h"

/* Set projective theta values at tau/2 */

static void
agm_ctx_set_th(acb_theta_agm_ctx_t ctx, slong prec)
{
    slong g = acb_theta_agm_ctx_g(ctx);
    acb_mat_t half;

    acb_mat_init(half, g, g);    
    acb_mat_scalar_mul_2exp_si(half, acb_theta_agm_ctx_tau(ctx), -1);
    
    if (acb_theta_agm_ctx_is_ext(ctx))
    {
        acb_theta_naive_proj(acb_theta_agm_ctx_th(ctx),
                acb_theta_agm_ctx_z(ctx), 2, half, prec);
    }
    else
    {
        acb_theta_naive_const_proj(acb_theta_agm_ctx_th(ctx), half, prec);
    }

    acb_mat_clear(half);
}

/* Candidates for symplectic matrices */

static void
fmpz_mat_Mi(fmpz_mat_t N, slong i)
{
    slong g = fmpz_mat_nrows(N)/2;

    fmpz_mat_one(N);
    fmpz_one(fmpz_mat_entry(N, i, i+g));
    fmpz_set_si(fmpz_mat_entry(N, i+g, i), -1);
    fmpz_zero(fmpz_mat_entry(N, i+g, i+g));
}

static void
fmpz_mat_Nij(fmpz_mat_t N, slong i, slong j)
{  
    slong g = fmpz_mat_nrows(N)/2;

    fmpz_mat_one(N);
    fmpz_one(fmpz_mat_entry(N, i, j+g));
    fmpz_one(fmpz_mat_entry(N, j, i+g));
    fmpz_set_si(fmpz_mat_entry(N, i+g, j), -1);
    fmpz_set_si(fmpz_mat_entry(N, j+g, i), -1);
    fmpz_zero(fmpz_mat_entry(N, i+g, i+g));
    fmpz_zero(fmpz_mat_entry(N, j+g, j+g));
}

static void
agm_ctx_candidates(fmpz_mat_struct* Ni, slong try, slong g)
{
    slong j, u, v, c;
    fmpz_mat_t J;
    flint_rand_t state;

    flint_randinit(state);
    fmpz_mat_init(J, 2*g, 2*g);

    fmpz_mat_J(J);
    
    /* Change state according to try */
    for (j = 0; j < try; j++) n_randint(state, 2);
  
    fmpz_mat_one(&Ni[0]);
    if (g == 1)
    {
        fmpz_mat_J(&Ni[1]);
    }
    else if (g == 2 && try == 0)
    {
        fmpz_mat_Mi(&Ni[1], 0);
        fmpz_mat_Mi(&Ni[2], 1);
        fmpz_mat_Nij(&Ni[3], 0, 1);
    }
    else
    {
        for (j = 1; j < (1<<g); j++)
        {
            /* (JM)^2 for random upper triangular M */
            fmpz_mat_one(&Ni[j]);
            for (u = 0; u < g; u++)
            {
                for (v = u; v < g; v++)
                {
                    c = n_randint(state, 3) - 1;
                    fmpz_set_si(fmpz_mat_entry(&Ni[j], u, v+g), c);
                    fmpz_set_si(fmpz_mat_entry(&Ni[j], v, u+g), c);
                }
            }
            fmpz_mat_mul(&Ni[j], J, &Ni[j]);
            fmpz_mat_mul(&Ni[j], &Ni[j], &Ni[j]);
        }
    }
    flint_randclear(state);
    fmpz_mat_clear(J);
}

/* Collect data for a single matrix. We only keep bad steps for which relative
   distance between roots is > 1/8 */

static void
agm_ctx_set_k2_ab(acb_theta_agm_ctx_t ctx, slong k)
{    
    acb_theta_agm_ctx_k2(ctx, k)
        = acb_theta_k2(acb_theta_agm_ctx_mat(ctx, k));
    acb_theta_agm_ctx_ab(ctx, k)
        = acb_theta_transform_image_char(acb_theta_agm_ctx_eps(ctx, k), 0,
                acb_theta_agm_ctx_mat(ctx, k));
}

static int
agm_ctx_is_good_step(acb_srcptr roots, slong n, int is_ext, slong lowprec,
        slong prec)
{
    arb_t eps1, eps2;
    int res;

    arb_init(eps1);
    arb_init(eps2);
    
    acb_theta_agm_rel_dist(eps1, roots, n, lowprec, prec);
    if (is_ext)
    {
        acb_theta_agm_rel_dist(eps2, roots+n, n, lowprec, prec);
        arb_max(eps1, eps1, eps2, lowprec);
    }
    arb_mul_2exp_si(eps1, eps1, 3);
    arb_sub_si(eps1, eps1, 1, prec);
    res = arb_is_negative(eps1);

    arb_clear(eps1);
    arb_clear(eps2);
    return res;
}

static void
agm_ctx_rescale_roots(acb_theta_agm_ctx_t ctx, slong k, slong prec)
{
    slong g = acb_theta_agm_ctx_g(ctx);
    slong n = 1<<g;
    slong nb_bad = acb_theta_agm_ctx_nb_bad(ctx, k);
    int is_ext = acb_theta_agm_ctx_is_ext(ctx);
    acb_t scal;
    slong i;

    acb_init(scal);

    if (is_ext && (nb_bad > 0))
    {
        acb_inv(scal, &acb_theta_agm_ctx_roots(ctx, k)[n], prec);
        _acb_vec_scalar_mul(acb_theta_agm_ctx_roots(ctx, k),
                acb_theta_agm_ctx_roots(ctx, k), 2*n*nb_bad, scal, prec);
        acb_inv(scal, &acb_theta_agm_ctx_roots(ctx, k)[0], prec);
        for (i = 0; i < nb_bad; i++)
        {
            _acb_vec_scalar_mul(acb_theta_agm_ctx_roots(ctx, k) + 2*n*i,
                    acb_theta_agm_ctx_roots(ctx, k) + 2*n*i, n, scal, prec);
            acb_sqrt(scal, scal, prec);
        }
    }
    else if (nb_bad > 0)
    {
        acb_inv(scal, &acb_theta_agm_ctx_roots(ctx, k)[0], prec);
        _acb_vec_scalar_mul(acb_theta_agm_ctx_roots(ctx, k),
                acb_theta_agm_ctx_roots(ctx, k), n*nb_bad, scal, prec);
    }

    acb_clear(scal);
}

static void
agm_ctx_set_roots(acb_theta_agm_ctx_t ctx, slong k, slong prec)
{
    slong g = acb_theta_agm_ctx_g(ctx);
    slong n = 1<<g;
    int is_ext = acb_theta_agm_ctx_is_ext(ctx);
    slong lowprec = ACB_THETA_AGM_LOWPREC;
    slong nb_th;
    acb_mat_t Ntau;
    acb_ptr Nz;
    acb_ptr roots;
    slong nb1, nb2;
    
    nb_th = 1<<g;
    if (acb_theta_agm_ctx_is_ext(ctx)) nb_th *= 2;

    acb_mat_init(Ntau, g, g);
    Nz = _acb_vec_init(2*g);
    
    /* Transform z, tau to medium precision */
    if (is_ext)
    {
        acb_siegel_transform_ext(Nz, Ntau, acb_theta_agm_ctx_mat(ctx, k),
                acb_theta_agm_ctx_z(ctx), acb_theta_agm_ctx_tau(ctx), prec);
        nb1 = acb_theta_agm_ext_nb_bad_steps(Nz, Ntau, prec);
    }
    else
    {
        acb_siegel_transform(Ntau, acb_theta_agm_ctx_mat(ctx, k),
                acb_theta_agm_ctx_tau(ctx), prec);
        nb1 = acb_theta_agm_nb_bad_steps(Ntau, prec);
    }

    /* Compute all roots to low precision */
    roots = _acb_vec_init(nb1 * nb_th);
    if (is_ext)
    {
        acb_theta_agm_ext_roots(roots, Nz, Ntau, nb1, lowprec);
    }
    else
    {
        acb_theta_agm_roots(roots, Ntau, nb1, lowprec);
    }

    /* Find out how many bad steps we exactly have */
    nb2 = nb1;
    while (nb2 > 0 &&
            agm_ctx_is_good_step(&roots[(nb2-1) * nb_th],
                    n, is_ext, lowprec, prec))
    {
        nb2 = nb2 - 1;
    }

    /* Set bad steps and roots */
    acb_theta_agm_ctx_reset_steps(ctx, k, nb2);
    _acb_vec_set(acb_theta_agm_ctx_roots(ctx, k), roots, nb2 * nb_th);
    agm_ctx_rescale_roots(ctx, k, lowprec);

    acb_mat_clear(Ntau);
    _acb_vec_clear(Nz, 2*g);
    _acb_vec_clear(roots, nb1 * nb_th);
}

static void
agm_ctx_get_mi_Mi(arf_struct* mi, arf_struct* Mi,
        const acb_theta_agm_ctx_t ctx, slong k, slong prec)
{
    slong g = acb_theta_agm_ctx_g(ctx);
    slong nb_bad = acb_theta_agm_ctx_nb_bad(ctx, k);
    slong nb_th;
    slong i;
    arb_t abs;

    nb_th = 1<<g;
    if (acb_theta_agm_ctx_is_ext(ctx)) nb_th *= 2;

    arb_init(abs);

    for (i = 0; i < nb_bad; i++)
    {
        acb_theta_agm_max_abs(abs, acb_theta_agm_ctx_roots(ctx, k) + i * nb_th,
                nb_th, prec);
        arb_sqr(abs, abs, prec);
        arb_get_ubound_arf(&Mi[i], abs, prec);
        acb_theta_agm_min_abs(abs, acb_theta_agm_ctx_roots(ctx, k) + i * nb_th,
                nb_th, prec);
        arb_sqr(abs, abs, prec);
        arb_get_lbound_arf(&mi[i], abs, prec);
    }

    arb_clear(abs);
}

static void
agm_ctx_deform_good(arf_t rad, arf_t min, arf_t max,
        const acb_theta_agm_ctx_t ctx, slong k, slong prec)
{
    slong g = acb_theta_agm_ctx_g(ctx);
    slong n = 1<<g;
    slong nb_bad = acb_theta_agm_ctx_nb_bad(ctx, k);
    acb_ptr a;
    arb_t eps, abs, m;

    a = _acb_vec_init(n);
    arb_init(eps);
    arb_init(abs);
    arb_init(m);

    /* Get start of good steps */
    acb_theta_agm_step_sqrt(a, acb_theta_agm_ctx_roots(ctx, k)
            + (nb_bad - 1) * n, g, prec);

    /* Get absolute dist; accept that as deformation; compute new min, max */
    acb_theta_agm_abs_dist(eps, a, n, prec, prec);
    arb_get_lbound_arf(rad, eps, prec);
    
    acb_abs(abs, &a[0], prec);
    arb_mul_si(eps, eps, 2, prec);
    arb_add(m, abs, eps, prec);
    arb_get_ubound_arf(max, m, prec);

    arb_sub(m, abs, eps, prec);
    arb_get_lbound_arf(min, m, prec);

    _acb_vec_clear(a, n);
    arb_clear(eps);
    arb_clear(abs);
    arb_clear(m);
}

static void
agm_ctx_deform_good_ext(arf_t rad, arf_t min, arf_t max,
        const acb_theta_agm_ctx_t ctx, slong k, slong prec)
{
    slong g = acb_theta_agm_ctx_g(ctx);
    slong n = 1<<g;
    slong nb_bad = acb_theta_agm_ctx_nb_bad(ctx, k);
    slong lowprec = ACB_THETA_AGM_LOWPREC;
    acb_ptr a;
    arb_t abs, x, y;
    arf_t mu, Mu, ms, Ms, eps;
    arf_t c1, c2, r;

    a = _acb_vec_init(2*n);
    arb_init(abs);
    arb_init(x);
    arb_init(y);
    arf_init(mu);
    arf_init(Mu);
    arf_init(ms);
    arf_init(Ms);
    arf_init(eps);
    arf_init(c1);
    arf_init(c2);
    arf_init(r);
    
    /* Get start of good steps */
    acb_theta_agm_ext_step_sqrt(a, acb_theta_agm_ctx_roots(ctx, k)
            + (nb_bad - 1) * 2*n, g, prec);

    /* Accepted deformation is minimum of:
       - 1/2 * min(abs of extended part)
       - absolute distance of pure Borchardt part */

    acb_theta_agm_min_abs(abs, a, n, prec);
    arb_mul_2exp_si(abs, abs, -1);
    acb_theta_agm_abs_dist(x, a+n, n, lowprec, prec);
    arb_min(abs, abs, x, prec);
    arb_get_lbound_arf(rad, abs, prec);

    /* Compute new min, max for extended and Borchardt parts, and relative
       distance for Borchardt part */

    acb_theta_agm_min_abs(x, a, n, prec);
    arb_sub_arf(x, x, rad, prec);
    arb_get_lbound_arf(mu, x, prec);

    acb_theta_agm_max_abs(x, a, n, prec);
    arb_add_arf(x, x, rad, prec);
    arb_get_ubound_arf(Mu, x, prec);

    acb_theta_agm_abs_dist(abs, a+n, n, lowprec, prec);    
    acb_abs(x, &a[0], prec);
    arb_sub(y, x, abs, prec);
    arb_sub_arf(y, y, rad, prec);
    arb_get_lbound_arf(ms, y, prec);

    arb_add(y, x, abs, prec);
    arb_add_arf(y, y, rad, prec);
    arb_get_ubound_arf(Ms, y, prec);

    arb_add_si(x, abs, 1, prec);
    arb_set_arf(abs, rad);
    arb_addmul_si(x, abs, 2, prec);
    arb_one(y);
    arb_sub(y, y, abs, prec);
    arb_div(x, x, y, prec);
    arb_sub_si(x, x, 1, prec);
    arb_get_ubound_arf(eps, x, prec);

    /* Compute minimal convergence rates on whole disk */
    acb_theta_agm_ext_conv_rate(c1, c2, r, eps, mu, Mu, prec);

    /* Deduce maximal, minimal values for extended Borchardt */
    acb_theta_agm_ext_rel_err(eps, c2, r, 1, prec);
    arf_mul(eps, eps, eps, prec, ARF_RND_CEIL);

    arb_set_arf(x, Mu);
    arb_mul_arf(x, x, Ms, prec);
    arb_set_arf(y, ms);
    arb_sqr(y, y, prec);
    arb_div(x, x, y, prec);
    arb_one(y);
    arb_add_arf(y, y, eps, prec);
    arb_mul(x, x, y, prec);
    arb_pow_ui(x, x, nb_bad, prec);
    arb_get_ubound_arf(max, x, prec);
    
    arb_set_arf(x, mu);
    arb_mul_arf(x, x, ms, prec);
    arb_set_arf(y, Ms);
    arb_sqr(y, y, prec);
    arb_div(x, x, y, prec);
    arb_one(y);
    arb_sub_arf(y, y, eps, prec);
    arb_mul(x, x, y, prec);
    arb_pow_ui(x, x, nb_bad, prec);
    arb_get_lbound_arf(min, x, prec);

    /* Compare with maximal, minimal values for regular Borchardt */
    arf_max(max, max, Ms);
    arf_min(min, min, ms);
    
    _acb_vec_clear(a, 2*n);
    arb_clear(abs);
    arb_clear(x);
    arb_clear(y);
    arf_clear(mu);
    arf_clear(Mu);
    arf_clear(ms);
    arf_clear(Ms);
    arf_clear(eps);
    arf_clear(c1);
    arf_clear(c2);
    arf_clear(r);
}

static void
agm_ctx_get_bounds(arf_t rad, arf_t min, arf_t max,
        const acb_theta_agm_ctx_t ctx, slong k, slong prec)
{
    int is_ext = acb_theta_agm_ctx_is_ext(ctx);
    slong lowprec = ACB_THETA_AGM_LOWPREC;
    slong nb_bad = acb_theta_agm_ctx_nb_bad(ctx, k);
    arf_struct* mi;
    arf_struct* Mi;
    slong i;
    
    mi = flint_malloc(nb_bad * sizeof(arf_struct));
    Mi = flint_malloc(nb_bad * sizeof(arf_struct));
    for (i = 0; i < nb_bad; i++)
    {
        arf_init(&mi[k]);
        arf_init(&Mi[k]);
    }
    
    /* Get mi, Mi */
    agm_ctx_get_mi_Mi(mi, Mi, ctx, k, lowprec);
    
    /* Pick an accepted deformation of first good step; deduce min, max */
    if (is_ext) agm_ctx_deform_good_ext(rad, min, max, ctx, k, lowprec);
    else agm_ctx_deform_good(rad, min, max, ctx, k, lowprec);
    
    /* Propagate radius back to projectivized theta values */
    acb_theta_agm_radius(rad, mi, Mi, rad, nb_bad, lowprec);
    
    /* Propagate radius back to projective theta(tau/2) */
    if (is_ext)
    {
        acb_theta_dupl_transform_radius(rad, rad,
            acb_theta_agm_ctx_th(ctx), acb_theta_agm_ctx_mat(ctx, k), lowprec);
    }
    else
    {        
        acb_theta_dupl_transform_radius_const(rad, rad,
            acb_theta_agm_ctx_th(ctx), acb_theta_agm_ctx_mat(ctx, k), lowprec);
    }

    for (i = 0; i < nb_bad; i++)
    {
        arf_clear(&mi[k]);
        arf_clear(&Mi[k]);
    }
    flint_free(mi);
    flint_free(Mi);
}

/* Collect global bounds */

static void
agm_ctx_get_rho_M(arf_t rho, arf_t M, const acb_theta_agm_ctx_t ctx, slong prec)
{
    slong g = acb_theta_agm_ctx_g(ctx);
    slong is_ext = acb_theta_agm_ctx_is_ext(ctx);
    acb_ptr dupl;
    arb_t abs0, abs1, abs;
    arf_t m0, rad, min, max;
    slong k;
    slong nb = acb_theta_agm_ctx_nb(ctx);

    if (is_ext) dupl = _acb_vec_init(1<<(2*g+1));
    else dupl = _acb_vec_init(1<<(2*g));
    arb_init(abs0);
    arb_init(abs1);
    arb_init(abs);
    arf_init(m0);
    arf_init(rad);
    arf_init(min);
    arf_init(max);
    
    if (is_ext) acb_theta_dupl_all(dupl, acb_theta_agm_ctx_th(ctx), g, prec);
    else acb_theta_dupl_all_const(dupl, acb_theta_agm_ctx_th(ctx), g, prec);
    
    agm_ctx_get_bounds(rho, m0, max, ctx, 0, prec);
    arf_zero(M);
    acb_abs(abs0, &dupl[0], prec);
    if (is_ext) acb_abs(abs1, &dupl[1<<(2*g)], prec);
    
    for (k = 1; k < nb; k++)
    {
        agm_ctx_get_bounds(rad, min, max, ctx, k, prec);        
        arf_min(rho, rho, rad);

        acb_abs(abs, &dupl[acb_theta_agm_ctx_ab(ctx, k)], prec);
        arb_mul_arf(abs, abs, max, prec);
        arb_div_arf(abs, abs, min, prec);
        arb_div(abs, abs, abs0, prec);
        arb_get_ubound_arf(max, abs, prec);
        arf_max(M, M, max);

        if (is_ext)
        {
            acb_abs(abs, &dupl[acb_theta_agm_ctx_ab(ctx, k) + (1<<(2*g))], prec);
            arb_mul_arf(abs, abs, max, prec);
            arb_div_arf(abs, abs, min, prec);
            arb_div(abs, abs, abs1, prec);
            arb_get_ubound_arf(max, abs, prec);
            arf_max(M, M, max);
        }
    }
    arf_div(M, M, m0, prec, ARF_RND_CEIL);    

    if (is_ext) _acb_vec_clear(dupl, 1<<(2*g+1));
    else _acb_vec_clear(dupl, 1<<(2*g));
    arb_clear(abs0);
    arb_clear(abs1);
    arb_clear(abs);
    arf_clear(m0);
    arf_clear(rad);
    arf_clear(min);
    arf_clear(max);
}

/* Get B3 */

static void
agm_ctx_get_B3(arf_t B3, const arf_t rho, const arf_t M,
        const acb_theta_agm_ctx_t ctx, slong prec)
{
    slong dim = acb_theta_agm_ctx_dim(ctx);
    arf_t B2;
    fmpz_t e;
    slong exp;
    arb_t eta;
    acb_ptr r;
    acb_mat_t fd;
    arb_t norm, bound;
    slong lowprec = ACB_THETA_AGM_LOWPREC;
    int res;
        
    arf_init(B2);
    fmpz_init(e);
    arb_init(eta);
    r = _acb_vec_init(dim);
    acb_mat_init(fd, dim, dim);
    arb_init(norm);
    arb_init(bound);
        
    /* Evaluate finite difference */
    acb_theta_cauchy(B2, rho, M, 2, dim, lowprec);
    arf_frexp(B2, e, B2);    
    exp = fmpz_get_si(e);
    arf_mul_2exp_si(B2, B2, exp);
    
    arb_one(eta);
    arb_mul_2exp_si(eta, eta, FLINT_MIN(- exp - n_clog(dim, 2), - prec/2));
    acb_theta_newton_fd(r, fd, acb_theta_agm_ctx_th(ctx), eta, ctx, prec);

    flint_printf("Finite diff:\n");
    acb_mat_printd(fd, 10); flint_printf("\n");    
    res = acb_mat_inv(fd, fd, prec);

    if (!res) arb_pos_inf(norm);
    else acb_mat_ninf(norm, fd, lowprec);
    
    flint_printf("Inv, norm:\n");
    acb_mat_printd(fd, 10); flint_printf("\n");
    arb_printd(norm, 10); flint_printf("\n");
      
    /* Is ||FD^-1||*n*B2*eta less than 1? If yes, deduce bound on dF^(-1) */
    arb_mul_arf(bound, norm, B2, lowprec);
    arb_mul_si(bound, bound, dim, lowprec);
    arb_mul(bound, bound, eta, lowprec);
    arb_sub_si(eta, bound, 1, lowprec);
    if (!arb_is_negative(eta))
    {
        arf_pos_inf(B3);
    }
    else
    {
        arb_mul(bound, bound, norm, lowprec);
        arb_add(bound, bound, norm, lowprec);
        arb_get_ubound_arf(B3, bound, lowprec);
    }

    arf_clear(B2);
    fmpz_clear(e);
    arb_clear(eta);
    _acb_vec_clear(r, dim);
    acb_mat_clear(fd);
    arb_clear(norm);
    arb_clear(bound);
}

/* Get logs */

static slong
fmpz_get_si_with_warning(const fmpz_t e)
{
    if (fmpz_cmp_si(e, WORD_MAX) > 0
            || fmpz_cmp_si(e, WORD_MIN) < 0)
    {
        flint_printf("agm_ctx_set: Error (cannot convert to slong)\n");
        fmpz_print(e); flint_printf("\n");
        fflush(stdout);
        flint_abort();
    }
    return fmpz_get_si(e);
}

static void
agm_ctx_set_logs(acb_theta_agm_ctx_t ctx, const arf_t rho, const arf_t M,
        const arf_t B3, slong prec)
{
    slong g = acb_theta_agm_ctx_g(ctx);
    slong dim = acb_theta_agm_ctx_dim(ctx);
    slong lowprec = ACB_THETA_AGM_LOWPREC;
    slong nb_th;
    arf_t c;
    fmpz_t e;
    arb_t abs;
    
    arf_init(c);
    fmpz_init(e);
    arb_init(abs);

    nb_th = 1<<g;
    if (acb_theta_agm_ctx_is_ext(ctx)) nb_th *= 2;
  
    arf_frexp(c, e, rho);
    acb_theta_agm_ctx_log_rho(ctx) = fmpz_get_si_with_warning(e) - 1;
    
    arf_frexp(c, e, M);
    acb_theta_agm_ctx_log_M(ctx) = fmpz_get_si_with_warning(e);
    
    arf_mul_si(c, M, 2*(dim+1), prec, ARF_RND_CEIL);
    arf_div(c, c, rho, prec, ARF_RND_CEIL);
    arf_frexp(c, e, c);
    acb_theta_agm_ctx_log_B1(ctx) = fmpz_get_si_with_warning(e);
    
    arf_mul_si(c, M, 2*(dim+1)*(dim+2), prec, ARF_RND_CEIL);
    arf_div(c, c, rho, prec, ARF_RND_CEIL);
    arf_div(c, c, rho, prec, ARF_RND_CEIL);
    arf_frexp(c, e, c);
    acb_theta_agm_ctx_log_B2(ctx) = fmpz_get_si_with_warning(e);
    
    arf_frexp(c, e, B3);
    acb_theta_agm_ctx_log_B3(ctx) = fmpz_get_si_with_warning(e);

    acb_theta_agm_max_abs(abs, acb_theta_agm_ctx_th(ctx), nb_th, prec);
    arb_get_ubound_arf(c, abs, lowprec);
    arf_frexp(c, e, c);
    acb_theta_agm_ctx_log_th(ctx) = fmpz_get_si_with_warning(e);

    arf_clear(c);
    fmpz_clear(e);
    arb_clear(abs);
}

/* User function */

int
acb_theta_agm_ctx_set(acb_theta_agm_ctx_t ctx, slong prec)
{
    slong g = acb_theta_agm_ctx_g(ctx);
    slong nb = acb_theta_agm_ctx_nb(ctx);
    arf_t rho, M, B3;
    int try = -1;
    int res = 0;
    slong k;

    arf_init(rho);
    arf_init(M);
    arf_init(B3);
    
    agm_ctx_set_th(ctx, prec);

    /* Try different matrix setups, starting with try=0 */
    while (try < ACB_THETA_AGM_NB_MATRIX_SETUPS)
    {
        try++;
        agm_ctx_candidates(ctx->mat, try, g);
        
        /* Set data each matrix */
        for (k = 0; k < nb; k++)
        {
            agm_ctx_set_k2_ab(ctx, k);
            agm_ctx_set_roots(ctx, k, prec);
        }

        /* Deduce global rho, M, B3 */
        agm_ctx_get_rho_M(rho, M, ctx, prec);
        agm_ctx_get_B3(B3, rho, M, ctx, prec);

        /* Set logs if valid, otherwise continue */
        if (arf_is_finite(M) && arf_is_finite(B3) && arf_cmp_si(rho,0) > 0)
        {
            res = 1;
            agm_ctx_set_logs(ctx, rho, M, B3, prec);
            break;
        }
    }
        
    arf_clear(rho);
    arf_clear(M);
    arf_clear(B3);
    return res;
}
