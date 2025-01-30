#ifndef ACB_HOLONOMIC_H
#define ACB_HOLONOMIC_H

/* temporary (for sage) */
#include "gmp.h"
typedef mp_limb_signed_t slong;
/* */

#include "acb_types.h"
#include "gr_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/****************************************************************************/

typedef struct
{
    slong n;
    slong mult;
}
acb_holonomic_shift_struct;


typedef struct
{
    acb_t expo;
    slong nshifts;
    acb_holonomic_shift_struct * shifts;
}
acb_holonomic_group_struct;


typedef struct
{
    slong len;
    acb_holonomic_group_struct * grps;
}
acb_holonomic_exponents_struct;


typedef enum
{
    ACB_HOLONOMIC_BASIS_ECHELON = 0,
    ACB_HOLONOMIC_BASIS_CASCADE  /* XXX: same as Frobenius? */
}
acb_holonomic_basis_t;


slong acb_holonomic_group_length(const acb_holonomic_group_struct * grp);

/****************************************************************************/

typedef struct
{
    acb_mat_t extini;

    slong nlogs;
    /* vector of polynomials in x, entries = coefficients of log(x)^k/k! of
     * a chunk of solution and/or rhs */
    acb_poly_struct * series;
    /* vector of polynomials holding the jets of coefficients wrt log(ξ)^k/k!
     * of the partial sums:
     * sums + j*alloc_logs + k is the jet of the coeff of log^k/k! at the
     * evaluation point of index j */
    /* XXX switch to acb_struct*? we are not really using the polynomial
     * structure */
    acb_poly_struct * sums;

    slong alloc_logs;  /* XXX support automatic growth? */
    slong alloc_pts;
}
acb_holonomic_sol_struct;

typedef acb_holonomic_sol_struct acb_holonomic_sol_t[1];

void acb_holonomic_sol_init(acb_holonomic_sol_struct * sol, slong nshifts,
                            slong nlogs, slong npts);

void acb_holonomic_sol_clear(acb_holonomic_sol_struct * sol);

void acb_holonomic_sol_reset(acb_holonomic_sol_struct * sol);

void acb_holonomic_sol_fit_length(acb_holonomic_sol_struct * sol, slong len,
                                  slong nder);

void acb_holonomic_sol_unit_ini(acb_holonomic_sol_struct * sol, slong i0,
                                const acb_holonomic_shift_struct * shifts);

acb_poly_struct * acb_holonomic_sol_sum_ptr(
        const acb_holonomic_sol_struct * sol, slong j, slong k);

void
_acb_holonomic_sol_add_term(
        acb_holonomic_sol_struct * sol, slong base, slong n,
        acb_poly_struct * rhs, slong rhs_nlogs,
        /* XXX define a separate struct with the next three? */
        slong mult, slong ini_i, const acb_poly_struct * ind_n,
        acb_srcptr pows, const char * shifted, slong npts,
        const fmpz * binom_n,
        slong nder, slong prec, slong sums_prec);

void _acb_holonomic_sol_value(acb_poly_struct * val, acb_srcptr expo,
                              const acb_poly_struct * sums, slong nlogs,
                              acb_srcptr pt, slong nder,
                              slong nshifts, slong prec);

/****************************************************************************/

/* XXX see if this can be reused by other summation algorithms (probably not as
 * is), rename */

#define ACB_HOLONOMIC_WANT_SERIES 1

typedef struct
{
    /* operator */
    acb_poly_struct * dop;  /* XXX FJ: move outside the struct? */
    slong dop_len;
    slong dop_clen;

    /* solution group */
    acb_t expo;
    acb_holonomic_shift_struct * sing_shifts;
    acb_poly_t ind;

    /* solutions */
    acb_holonomic_sol_struct * sol;
    slong nsols;

    /* evaluation points */
    acb_ptr pts;   /* evaluation points x[j] for 0 < j < npts */
    acb_ptr pows;  /* x[j]^{n-i} for i < nder (i-major, cyclic on i);
                      only the block i=0 is really used for nonzero points */
    char * shifted_sums;
    slong npts;
    slong nder;

    /* working precisions */
    slong prec;
    slong sums_prec;
    slong bounds_prec;

    /* TODO: error bounds etc. */

    /* ancillary data */
    fmpz * binom_n;  /* binom(n, j) for j < nder */

    /* options */
    slong block_size;
    ulong flags;

    int have_precomputed;
}
acb_holonomic_sum_context_struct;

typedef acb_holonomic_sum_context_struct acb_holonomic_sum_context_t[1];

void acb_holonomic_sum_context_init(
        acb_holonomic_sum_context_struct * ctx,
        slong dop_len, slong npts, slong nsols, slong nder);
void acb_holonomic_sum_context_clear(acb_holonomic_sum_context_struct * ctx);

void acb_holonomic_sum_ordinary(acb_holonomic_sum_context_struct * ctx);
void acb_holonomic_sum_canonical_basis(acb_holonomic_sum_context_struct * ctx);
void acb_holonomic_sum_highest(acb_holonomic_sum_context_struct * ctx);
void acb_holonomic_sum_group(acb_holonomic_sum_context_struct * ctx, const acb_holonomic_group_struct * grp);


void _acb_holonomic_sum_precompute(acb_holonomic_sum_context_struct * ctx);
void _acb_holonomic_sum_reset(acb_holonomic_sum_context_struct * ctx);

void acb_holonomic_sum_divconquer(
        acb_holonomic_sum_context_struct * ctx, slong nterms);

void acb_holonomic_sum_value(acb_poly_struct * val, slong nfrobshifts,
                             const acb_holonomic_sum_context_struct * ctx,
                             slong i, slong j, slong prec);

/****************************************************************************/

/* XXX void * data or fix args */
typedef void (* acb_holonomic_sum_worker_t)(
        acb_holonomic_sum_context_struct * ctx, slong nterms);

void _acb_holonomic_fundamental_matrix(
        acb_mat_struct * mat,
        const acb_poly_struct * dop, slong dop_len,
        const acb_holonomic_group_struct * groups, slong ngroups,
        acb_srcptr pts, slong npts,
        acb_holonomic_basis_t basis,
        acb_holonomic_sum_worker_t sum_worker,
        slong nterms, slong prec);

int acb_holonomic_fundamental_matrix(
        acb_mat_struct * mat,
        gr_srcptr dop, gr_ctx_t dop_ctx,
        const acb_holonomic_exponents_struct * expos,
        acb_srcptr pts, slong npts,
        acb_holonomic_basis_t basis,
        slong nterms, slong prec);

/****************************************************************************/

void _acb_holonomic_apply_diffop_basecase_weights(
        acb_ptr weights,
        const acb_poly_struct * dop, slong dop_len,
        acb_srcptr expo, slong offset,
        slong flen, slong nlogs, slong start, slong len,
        slong prec);

void _acb_holonomic_apply_diffop_basecase_precomp(
        acb_poly_struct * g, slong goff,
        acb_srcptr weights, slong weights_nlogs,
        const acb_poly_struct * f, slong foff, slong flen,
        slong nlogs,
        slong start, slong len,
        slong prec);

void acb_holonomic_apply_diffop_basecase(
        acb_poly_struct * g,
        const acb_poly_struct * dop, slong dop_len,
        acb_srcptr expo, slong offset,
        const acb_poly_struct * f,
        slong nlogs,
        slong start, slong len,
        slong prec);

void _acb_holonomic_apply_diffop_polmul(
        acb_poly_struct * g, slong goff,
        const acb_poly_struct * dop, slong dop_len,
        acb_srcptr expo, slong offset,
        const acb_poly_struct * f, slong foff, slong flen,
        slong nlogs,
        slong start, slong len,
        slong prec);

void acb_holonomic_apply_diffop_polmul(
        acb_poly_struct * g,
        const acb_poly_struct * dop, slong dop_len,
        acb_srcptr expo, slong offset,
        const acb_poly_struct * f,
        slong nlogs,
        slong start, slong len,
        slong prec);

#ifdef __cplusplus
        }
#endif

#endif
