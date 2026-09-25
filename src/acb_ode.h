#ifndef ACB_ODE_H
#define ACB_ODE_H

#ifdef ACB_ODE_INLINES_C
#define ACB_ODE_INLINE
#else
#define ACB_ODE_INLINE static inline
#endif

// /* temporary (for sage) */
// #include "gmp.h"
// typedef mp_limb_signed_t slong;
// /* */

#include "acb_types.h"
#include "gr_types.h"
#include "gr_ore_poly.h"

// #ifdef __cplusplus
// extern "C" {
// #endif

/****************************************************************************/

// probably want:
// - a data structure for storing “local” information (shifted singularities,
// local exponent structure...),
// - a function (acb_ode_focus?) for “moving it” to another point
typedef struct
{
    // XXX sing_gr?
    acb_ptr sing_acb;
    slong sing_len;

    // XXX maybe something to track the precision at which transition matrices
    // are known??
}
acb_ode_context_struct;

typedef acb_ode_context_struct acb_ode_context_t[1];

/****************************************************************************/

typedef struct
{
    slong n;
    slong mult;
}
acb_ode_shift_struct;

typedef acb_ode_shift_struct acb_ode_shift_t[1];


typedef struct
{
    acb_t leader;
    acb_ode_shift_struct * shifts;
    slong nshifts;
}
acb_ode_group_struct;

typedef acb_ode_group_struct acb_ode_group_t[1];


typedef struct
{
    /* Memory managed by the exponent object, partial sharing across groups */
    acb_ode_group_struct * groups;
    slong ngroups;
}
acb_ode_exponents_struct;

typedef acb_ode_exponents_struct acb_ode_exponents_t[1];


typedef enum
{
    ACB_ODE_BASIS_ECHELON = 0,
    ACB_ODE_BASIS_CASCADE,

    ACB_ODE_NUM_BASES
}
acb_ode_basis_t;

void acb_ode_group_init(acb_ode_group_t group, slong len);
void acb_ode_group_clear(acb_ode_group_t group);

void acb_ode_group_set(acb_ode_group_t dest, const acb_ode_group_t src);
slong acb_ode_group_length(const acb_ode_group_t grp);
slong acb_ode_group_multiplicity (const acb_ode_group_t group, slong n);
slong acb_ode_group_nlogs (const acb_ode_group_t group, slong n);

void acb_ode_exponents_init(acb_ode_exponents_t expos);
void acb_ode_exponents_clear(acb_ode_exponents_t expos);

slong acb_ode_exponents_length(const acb_ode_exponents_t expos);
void acb_ode_exponents_println(const acb_ode_exponents_t expos);

void acb_ode_exponents_ordinary(acb_ode_exponents_t expos, slong dop_order);
void acb_ode_exponents_randtest(acb_ode_exponents_t expos, flint_rand_t state, slong len, slong disp, slong prec, slong mag_bits);
int acb_ode_exponents(acb_ode_exponents_t expos, const gr_ore_poly_t dop, gr_ctx_t Dop, gr_ctx_t CC);

void acb_ode_indicial_polynomial_from_exponents(acb_poly_t ind, const acb_ode_exponents_t expos, slong prec);

/****************************************************************************
 * Error bounds
 ****************************************************************************/

typedef struct
{
    slong prec;

    arb_t cst;       /* 1/abs(lc_x(lc_θ(dop))) */
    arb_poly_t den_lbound;
    arb_ptr den_rt;  /* roots(lc_θ(dop)), with multiple roots repeated */
    slong den_rt_len;

    acb_poly_struct * all_nums;
    slong all_nums_len;
    slong pol_part_len;
}
acb_ode_bound_struct;

typedef acb_ode_bound_struct acb_ode_bound_t[1];

typedef struct
{
    slong length;
    struct
    {
        acb_t root;
        mag_t global_lbound;
        slong mult;
        slong n_min;
    }
    * r;
}
acb_ode_ind_lbound_struct;

typedef acb_ode_ind_lbound_struct acb_ode_ind_lbound_t[1];

typedef struct
{
    slong length;
    mag_ptr h;
}
acb_ode_stairs_struct;

typedef acb_ode_stairs_struct acb_ode_stairs_t[1];

typedef struct
{
    acb_poly_t ind;  /* MONIC indicial polynomial shifted so that 0 corresponds
                        to the group leader */
    acb_ode_ind_lbound_t ind_lbound;
    acb_ode_stairs_t stairs;
}
acb_ode_group_bound_struct;

typedef acb_ode_group_bound_struct acb_ode_group_bound_t[1];

void acb_ode_bound_init(acb_ode_bound_t bound);
void acb_ode_bound_clear(acb_ode_bound_t bound);
void acb_ode_bound_precompute(acb_ode_bound_t bound,
        const acb_poly_struct *dop, slong dop_len, acb_srcptr lcrt,
        slong pol_part_len, slong prec);

void acb_ode_group_bound_init(acb_ode_group_bound_t gbound);
void acb_ode_group_bound_clear(acb_ode_group_bound_t gbound);
void acb_ode_group_bound_precompute(acb_ode_group_bound_t gbound,
        const acb_poly_struct * dop, slong dop_len,
        const acb_ode_exponents_t expos, slong grp,
        const acb_ode_bound_t bound,
        slong prec);

void acb_ode_bound_precompute_integrand(arb_poly_t itg_pol, arb_poly_t itg_num,
        const acb_ode_group_t group, const acb_ode_bound_t bound,
        const acb_ode_group_bound_t gbound, slong n0, slong nlogs, slong prec);
void acb_ode_tail_bound_jet_precomp(arb_poly_t res,
        const acb_ode_bound_t bound, slong n,
        const arb_poly_t itg_pol, const arb_poly_t itg_num,
        const arb_poly_t nres_maj,
        const arb_t rad, slong ord, slong prec);
void acb_ode_tail_bound_jet(arb_poly_t res,
        const acb_ode_group_t group, const acb_ode_bound_t bound,
        const acb_ode_group_bound_t gbound,
        slong n, slong nlogs, const arb_poly_t nres_maj,
        const arb_t rad, slong ord, slong prec);

void acb_ode_ind_lbound_init(acb_ode_ind_lbound_t ind_lbound);
void acb_ode_ind_lbound_clear(acb_ode_ind_lbound_t ind_lbound);
void acb_ode_ind_lbound_precompute(acb_ode_ind_lbound_t ind_lbound,
        const acb_ode_exponents_t expos, slong grp, slong prec);
void acb_ode_ind_lbound_eval(mag_t res, const acb_ode_ind_lbound_t ind_lbound,
        slong n, slong prec);

void acb_ode_stairs_init(acb_ode_stairs_t stairs);
void acb_ode_stairs_clear(acb_ode_stairs_t stairs);
void acb_ode_stairs_precompute(acb_ode_stairs_t stairs,
        const acb_poly_struct * num, slong len,
        const acb_poly_t ind, const acb_ode_group_t group,
        const acb_ode_ind_lbound_t ind_lbound, slong prec);

void acb_ode_bound_rat_vec(mag_ptr res,
        const acb_poly_struct * num, slong len,
        const acb_poly_t ind, const acb_ode_group_t group,
        const acb_ode_ind_lbound_t ind_lbound, const acb_ode_stairs_t stairs,
        slong n0, slong ord, slong prec);
void acb_ode_bound_rat_ordinary_vec(mag_ptr res,
        const acb_poly_struct * num, slong len,
        const acb_poly_t ind, const acb_ode_ind_lbound_t ind_lbound,
        slong n0, slong ord, slong prec);
void acb_ode_bound_rat_ref_vec(mag_ptr res,
        const acb_poly_struct *num, slong len,
        const acb_poly_t ind, slong n, slong mult, slong ord, slong prec);

/****************************************************************************/

typedef struct
{
    double neglogterm;
    double cvg_rate;
    slong accuracy;
    double loss_rate;
    double terms_wanted;
    double terms_full_prec;
    slong prec_wanted;
}
acb_ode_cvest_struct;

typedef acb_ode_cvest_struct acb_ode_cvest_t[1];

void acb_ode_cvest_init(acb_ode_cvest_t cvest);
void acb_ode_cvest_clear(acb_ode_cvest_t cvest);
void acb_ode_cvest_update(acb_ode_cvest_t cvest, const acb_ode_cvest_t old, mag_t est, slong accuracy, slong stride, slong prec, slong work_prec);

/****************************************************************************/

/* todo: better vectors of polynomials / polynomial with vector coefficients
   (single buffer, row stride, essentially self-resizing matrices) */

typedef struct
{
    /* Vector of polynomials in x, not necessarily normalized (used as buffers).
       The polynomial of index k corresponds to the series in front of
       log(x)^k/k! in the solution. Its coefficients 0..n-n0-1 are the
       coefficients of x^n0 to x^{n-1} in the series. FIXME: Generally speaking,
       the coefficients of index >= n0 are used to store the coefficients of
       x^n, x^{n+1}, ... in the *residual* of the sum truncated at order n, but
       this is currently not a properly maintained invariant. The values of n0
       and n are provided by the associated sum object.  */
    acb_poly_struct * series;

    /* coefficients of Taylor expansions (jets) of the partial sums
     * dim: sums[npts][alloc_logs][nder] */
    acb_ptr sums;
    /* "extended" initial value matrix; nshifts × logs; *rows* correspond to
       shifts in increasing order; column k contains the coeff of
       log(x)^k/k!*x^(leader+shift); the first mult of the entries of a given
       row are initial values, the remaining entries are filled by series
       expansion functions and used for recombining solutions */
    acb_mat_t extini;

    /* bounds on the tails of successive derivatives, valid uniformly for all
     * coefficients of log(x)^k/k! and all evaluation points */
    // XXX move to cvest?
    mag_ptr tb;

    /* convergence estimates */
    acb_ode_cvest_t cvest;

    slong alloc_logs;   /* vector size, >= #logs in full series */
    slong nlogs;        /* length in log of current partial sum */
    slong future_logs;  /* >= #additional logs in future partial sums */
    slong npts;
    slong nder;

    /* flag indicating that no further terms should be added */
    int done;
}
acb_ode_sol_struct;

typedef acb_ode_sol_struct acb_ode_sol_t[1];

ACB_ODE_INLINE acb_ptr
acb_ode_sol_sum_ptr(const acb_ode_sol_t sol,
                    slong p, slong k, slong i)
{
    return sol->sums + (p * sol->alloc_logs + k) * sol->nder + i;
}

void acb_ode_sol_init(acb_ode_sol_t sol, slong nshifts,
                      slong nlogs, slong npts, slong nder);

void acb_ode_sol_clear(acb_ode_sol_t sol);

void acb_ode_sol_zero(acb_ode_sol_t sol);

void acb_ode_sol_fit_length(acb_ode_sol_t sol, slong len);


void
acb_ode_sol_set_ini(acb_ode_sol_t sol, acb_srcptr ini,
                    const acb_ode_shift_t shifts);  // unused

void _acb_ode_sol_jet(acb_poly_struct * val, const acb_t expo, acb_srcptr sums, slong stride, mag_srcptr tb, const acb_t pt, slong nlogs, slong ord, slong nfrobshifts, slong prec);
void acb_ode_sol_jet(acb_poly_struct * val, const acb_t expo, const acb_ode_sol_t sol, slong p, const acb_t pt, slong ord, slong nfrobshifts, slong prec);

slong
acb_ode_sol_estimate_terms(mag_t est,
                           const acb_ode_sol_t sol, slong off, slong len,
                           const mag_t radpow);
slong acb_ode_sol_estimate_sums(mag_t mag, mag_t rad, const acb_ode_sol_t sol);


/****************************************************************************/

#define ACB_ODE_WANT_SERIES 1
#define ACB_ODE_APPROX      2

/* to be adapted for compatibility with other summation algorithms */

typedef struct
{
    /* operator */
    acb_poly_struct * dop;  /* XXX FJ: move outside the struct? */
    slong dop_len;
    slong dop_clen;
    mag_t cvrad;

    /* solution group */  // XXX full exponent structure???
    acb_ode_group_t group;
    acb_poly_t ind;

    /* solutions */
    /* XXX possible alternative design: move out of sum struct, use sum as a
       context when manipulating sols, replace some sum_* functions by sol_*_vec
       functions, some of which would also update the context object */
    acb_ode_sol_struct * sol;
    slong nsols;
    slong n0;  /* index of the first stored coefficient */
    slong n;   /* current truncation order */

    /* evaluation points */
    acb_ptr pts;   /* evaluation points x[j] for 0 <= j < npts */
    acb_ptr pows;  /* x[j]^{n-i} for i < nder (i-major, cyclic on i);
                      only the block i=0 is really used for nonzero points */
    char * shifted_sums;  /* nonzero if the i-th derivatives of the sums at x[j]
                             stored in sol have extra x[j]^i factors */
    slong npts;
    slong nder;

    /* working precisions */
    slong wp;
    slong sum_wp;

    /* error bounds */
    mag_t mag;
    mag_t magpow;
    arb_poly_t itg_pol, itg_num;  // itg_n ?

    /* options */
    ulong flags;
    void * data;     /* user data for sum workers */

    /* miscellaneous state */
    fmpz * binom_n;  /* binom(n, j) for j < nder */
    int have_precomputed;
}
acb_ode_sum_struct;

typedef acb_ode_sum_struct acb_ode_sum_t[1];

void acb_ode_sum_init(acb_ode_sum_t sum, slong dop_len, slong npts, slong nsols, slong nder);
void acb_ode_sum_clear(acb_ode_sum_t sum);

void acb_ode_sum_fit_length(acb_ode_sum_t sum, slong len);
void acb_ode_sum_discard_head(acb_ode_sum_struct * sum, slong n0);

void acb_ode_sum_set_diffop(acb_ode_sum_t sum, const acb_poly_struct * dop, slong dop_len, const mag_t cvrad);

void acb_ode_sum_set_group(acb_ode_sum_t sum, const acb_ode_group_t grp);
void acb_ode_sum_set_ordinary(acb_ode_sum_t sum);

void acb_ode_sum_set_ini_echelon(acb_ode_sum_t sum);
void acb_ode_sum_set_ini_highest(acb_ode_sum_t sum);

void acb_ode_sum_set_points(acb_ode_sum_t sum, acb_srcptr pts, slong npts);

slong acb_ode_sum_max_nlogs(acb_ode_sum_struct * sum);
slong acb_ode_sum_max_nlogs_xn(acb_ode_sum_struct * sum, slong off);

void acb_ode_sum_precompute(acb_ode_sum_t sum);

void _acb_ode_sum_forward_1(acb_ode_sum_struct * sum);
int _acb_ode_sum_forward_divconquer(acb_ode_sum_struct * sum, slong high, slong block_len, slong stride, acb_ode_bound_t bound, acb_ode_group_bound_t gbound, slong prec);
void _acb_ode_sum_update_residuals(acb_ode_sum_struct * sum, slong low, slong mid, slong high);
void _acb_ode_sum_fix(acb_ode_sum_struct * sum);

int _acb_ode_sum_done(acb_ode_sum_struct * sum, slong stride,
                      acb_ode_bound_t bound, acb_ode_group_bound_t gbound,
                      slong prec);

void acb_ode_sum_divconquer(acb_ode_sum_t sum, slong nterms,
                            acb_ode_bound_t bound, acb_ode_group_bound_t gbound,
                            slong prec);

void _acb_ode_poly_negdivrevhigh(acb_ptr res, acb_srcptr a, acb_srcptr cst,
                             acb_srcptr b, slong len, slong prec);
void acb_ode_poly_taylor_shift_aps_trunc(acb_poly_t g, const acb_poly_t f, acb_srcptr a, slong n, slong len, slong prec);


/****************************************************************************/

typedef void (* acb_ode_sum_worker_t)(acb_ode_sum_t sum, slong nterms, acb_ode_bound_t bound, acb_ode_group_bound_t gbound, slong prec);

void _acb_ode_fundamental_matrix_vec(
        acb_mat_struct * mat,
        const acb_poly_struct * dop, slong dop_len,
        const acb_ode_exponents_t expos,
        acb_srcptr lcrt,
        acb_ode_bound_t bound,  /* mutable */
        acb_srcptr pts, slong npts,
        acb_ode_basis_t basis,
        acb_ode_sum_worker_t sum_worker,
        void * worker_data,
        slong prec);

WARN_UNUSED_RESULT int acb_ode_fundamental_matrix_vec(
        acb_mat_struct * mat,
        const gr_ore_poly_t dop, gr_ore_poly_ctx_t Dop,
        const acb_ode_exponents_t expos,
        acb_srcptr lcrt,
        acb_srcptr pts, slong npts,
        acb_ode_basis_t basis,
        slong prec);

WARN_UNUSED_RESULT int acb_ode_fundamental_matrix(
        acb_mat_t mat,
        const gr_ore_poly_t dop, gr_ore_poly_ctx_t Dop,
        const acb_ode_exponents_t expos,
        acb_srcptr lcrt,
        const acb_t pt,
        acb_ode_basis_t basis,
        slong prec);

/****************************************************************************/

void acb_ode_randtest_acb(acb_poly_struct * dop, acb_struct * lcroots, acb_ode_exponents_t expos, flint_rand_t state, slong disp, slong lcdeg, slong clen, slong len, slong prec);

/****************************************************************************/

void _acb_ode_apply_diffop_basecase_weights(
        acb_ptr weights,
        const acb_poly_struct * dop, slong dop_len,
        acb_srcptr expo, slong offset,
        slong flen, slong nlogs, slong start, slong len,
        slong prec);

void _acb_ode_apply_diffop_basecase_precomp(
        acb_poly_struct * g, slong goff,
        acb_srcptr weights, slong weights_nlogs,
        const acb_poly_struct * f, slong foff, slong flen,
        slong nlogs,
        slong start, slong len,
        slong prec);

void acb_ode_apply_diffop_basecase(
        acb_poly_struct * g,
        const acb_poly_struct * dop, slong dop_len, // gr_ore_poly?
        acb_srcptr expo, slong offset,
        const acb_poly_struct * f,
        slong nlogs,
        slong start, slong len,
        slong prec);

void _acb_ode_apply_diffop_polmul(
        acb_poly_struct * g, slong goff,
        const acb_poly_struct * dop, slong dop_len,
        acb_srcptr expo, slong offset,
        const acb_poly_struct * f, slong foff, slong flen,
        slong nlogs,
        slong start, slong len,
        slong prec);

void acb_ode_apply_diffop_polmul(
        acb_poly_struct * g,
        const acb_poly_struct * dop, slong dop_len, // gr_ore_poly?
        acb_srcptr expo, slong offset,
        const acb_poly_struct * f,
        slong nlogs,
        slong start, slong len,
        slong prec);

// todo: acb_ode_apply_diffop

/****************************************************************************/

void _acb_ode_solution_growth(mag_t order, mag_t base, const acb_poly_struct * dop, slong dop_len);
slong acb_ode_choose_prec(slong * rec_prec, const acb_poly_struct * dop, slong dop_len, mag_srcptr rad, mag_srcptr cvrad, slong tgt_prec);


// #ifdef __cplusplus
// }
// #endif

#endif
