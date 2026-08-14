#include <stdalign.h>
#include <string.h>
#include "acb.h"
#include "acb_ode.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "gr_poly.h"
#include "gr_vec.h"

void
acb_ode_group_init(acb_ode_group_t group, slong len)
{
    acb_init(group->leader);
    group->shifts = flint_malloc(sizeof(acb_ode_shift_struct) * len);
    group->nshifts = 0;
}

void
acb_ode_group_clear(acb_ode_group_t group)
{
    acb_clear(group->leader);
    flint_free(group->shifts);
}

void
acb_ode_group_set(acb_ode_group_t dest, const acb_ode_group_t src)
{
    acb_set(dest->leader, src->leader);
    memcpy(dest->shifts, src->shifts,
           src->nshifts * sizeof(acb_ode_shift_struct));
    dest->nshifts = src->nshifts;
}

slong
acb_ode_group_length(const acb_ode_group_t grp)
{
    slong len = 0;

    for (slong s = 0; s < grp->nshifts; s++)
        len += grp->shifts[s].mult;

    return len;
}

slong
acb_ode_group_multiplicity (const acb_ode_group_t group, slong n)
{
    for (slong s = 0; s < group->nshifts; s++)
        if (group->shifts[s].n == n)
            return group->shifts[s].mult;
    return 0;
}

slong
acb_ode_group_nlogs (const acb_ode_group_t group, slong n)
{
    slong nlogs = 0;
    for (slong s = 0; s < group->nshifts; s++) {
        if (group->shifts[s].n > n)
            break;
        nlogs += group->shifts[s].mult;
    }
    return nlogs;
}

void
acb_ode_exponents_init(acb_ode_exponents_t expos)
{
    expos->ngroups = 0;
    expos->groups = NULL;
}

void
acb_ode_exponents_clear(acb_ode_exponents_t expos)
{
    for (slong g = 0; g < expos->ngroups; g++)
        /* group objects referenced by exponents objects share data, so we
           cannot call group_clear */
        acb_clear(expos->groups[g].leader);
    flint_free(expos->groups);
}

slong
acb_ode_exponents_length(const acb_ode_exponents_t expos)
{
    slong len = 0;

    for (slong g = 0; g < expos->ngroups; g++)
        len += acb_ode_group_length(expos->groups + g);

    return len;
}

void
acb_ode_exponents_ordinary(acb_ode_exponents_t expos, slong order)
{
    acb_ode_exponents_clear(expos);

    slong off = sizeof(acb_ode_group_struct);
    slong bufsz = off + sizeof(acb_ode_shift_struct) * order;
    unsigned char * buf = flint_malloc(bufsz);
    expos->groups = (acb_ode_group_struct *) buf;
    expos->groups->shifts = (acb_ode_shift_struct *) (buf + off);

    acb_init(expos->groups->leader);
    for (slong n = 0; n < order; n++)
    {
        expos->groups->shifts[n].n = n;
        expos->groups->shifts[n].mult = 1;
    }
    expos->groups->nshifts = order;

    expos->ngroups = 1;
}

void
acb_ode_exponents_println(const acb_ode_exponents_t expos)
{
    flint_printf("[");
    for (slong g = 0; g < expos->ngroups; g++)
    {
        acb_ode_group_struct * grp = expos->groups + g;
        flint_printf("%{acb} ", grp->leader);
        for (slong s = 0; s < grp->nshifts; s++)
        {
            flint_printf("+ %wd", grp->shifts[s].n);
            if (grp->shifts[s].mult != 1)
                flint_printf(" (mult=%wd)", grp->shifts[s].mult);
            if (g < expos->ngroups - 1 || s < grp->nshifts - 1)
                flint_printf(", ");
        }
    }
    flint_printf("]\n");
}



/* Random exponents, exactly len in total, nontrivial shift structure with high
   probability. */
void
acb_ode_exponents_randtest(acb_ode_exponents_t expos, flint_rand_t state,
                           slong len, slong disp, slong prec, slong mag_bits)
{
    arb_t tmp;
    arb_init(tmp);

    acb_ode_exponents_clear(expos);
    void * buf = flint_malloc(len * (sizeof(acb_ode_group_struct) +
                                     sizeof(acb_ode_shift_struct)));
    expos->groups = (acb_ode_group_struct *) buf;
    acb_ode_shift_struct * shift =
        (acb_ode_shift_struct *) (expos->groups + len);

    expos->ngroups = 0;
    int precise = n_randint(state, 64);
    for (acb_ode_group_struct * grp = expos->groups; len > 0; grp++)
    {
        expos->ngroups++;

        acb_init(grp->leader);
        if (precise)
            acb_randtest_precise(grp->leader, state, prec, mag_bits);
        else
            acb_randtest(grp->leader, state, prec, mag_bits);

        /* leaders differing by integers would be unsound, so we perturb one
           and pretend they were close but different */
        for (slong g = 0; g < expos->ngroups - 1; g++)
        {
            if (acb_is_exact(grp->leader) && acb_is_exact(expos->groups[g].leader))
            {
                arb_sub(tmp,
                        acb_realref(grp->leader),
                        acb_realref(expos->groups[g].leader),
                        ARF_PREC_EXACT);
                if (arb_is_int(tmp))
                    mag_set_d(&grp->leader->real.rad, 1e-50);
            }
        }

        grp->shifts = shift;
        grp->nshifts = 0;

        slong grplen = 1 + n_randint(state, len);
        len -= grplen;

        for (slong n = 0;;)
        {
            // flint_printf("len = %wd grplen = %wd n = %wd\n", len, grplen, n);

            grp->nshifts++;

            shift->n = n;
            if (n == disp)
                shift->mult = grplen;
            else
                shift->mult = 1 + n_randint(state, grplen);
            grplen -= shift->mult;

            shift++;

            if (grplen == 0)
                break;

            n += 1 + n_randint(state, disp - n);
        }
    }

    arb_clear(tmp);
}

/* todo: also provide a function for computing the exponents starting from the
   indicial polynomial */
int
acb_ode_exponents(acb_ode_exponents_t expos, const gr_ore_poly_t dop,
                  gr_ctx_t Dop, gr_ctx_t CC)
{
    int status = GR_SUCCESS;

    gr_ctx_struct * Scalars, * Pol;
    gr_ctx_t ZZ, ZZvec;
    gr_poly_t ind, polc;
    gr_vec_t slshifts, slmult, roots;
    gr_poly_vec_t fac, slfac;
    fmpz_vec_t mult, rtmult;

    if (gr_ore_poly_ctx_over_gr_poly_base_ptrs(&Scalars, &Pol, Dop) != GR_SUCCESS)
    {
        // flint_printf("%s:%d UNABLE\n", __FILE__, __LINE__);
        return GR_UNABLE;
    }

    /* For now exponents are of type acb_t */
    if (CC->which_ring != GR_CTX_CC_ACB)
    {
        // flint_printf("%s:%d UNABLE\n", __FILE__, __LINE__);
        return GR_UNABLE;
    }

    gr_poly_init(ind, Scalars);

    /* Compute the local exponents (indicial roots) and organize them into shift
       equivalence classes */

    status |= gr_ore_poly_indicial_polynomial(ind, dop, Dop);
    if(status)
    {
        // flint_printf("%s:%d status=%d\n", __FILE__, __LINE__, status);
        goto cleanup_ind;
    }

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_vector_gr_vec(ZZvec, ZZ);  /* XXX fmpz_vec */

    gr_poly_init(polc, Scalars);
    gr_poly_vec_init(fac, 0, Scalars);
    fmpz_vec_init(mult, 0);
    gr_poly_vec_init(slfac, 0, Scalars);
    gr_vec_init(slshifts, 0, ZZvec);
    gr_vec_init(slmult, 0, ZZvec);
    gr_vec_init(roots, 0, CC);
    fmpz_vec_init(rtmult, 0);

    /* XXX In some cases this should use the irreducible factorization for
       representing the exponents exactly; in others, this should compute the
       roots of the shiftless factors without further decomposing them.

       Quite often it would be possible to continue even if gr_factor fails,
       however the current implementation of shiftless decomposition requires
       factoring. */

    status |= gr_factor(polc /* unused */, (gr_vec_struct *) fac,
                        mult, ind, 0, Pol);
    // if(status) flint_printf("%s:%d status=%d\n", __FILE__, __LINE__, status);
    status |= gr_poly_shiftless_decomposition_from_factors(
            slfac, slshifts, slmult, fac, mult, Scalars);
    // if(status) flint_printf("%s:%d status=%d\n", __FILE__, __LINE__, status);

    /*
    gr_ctx_t Polvec, ZZvecvec;
    gr_ctx_init_vector_gr_vec(Polvec, Pol);
    gr_ctx_init_vector_gr_vec(ZZvecvec, ZZvec);
    flint_printf("sldec: %{gr} %{gr} %{gr}\n", slfac, Polvec, slshifts,
                 ZZvecvec, slmult, ZZvecvec);
    gr_ctx_clear(ZZvecvec);
    gr_ctx_clear(Polvec);
    */

    acb_ode_exponents_clear(expos);

    /* Allocate space for sum(deg(slfac)) exponents and sum(len(shiftset))
       shifts; exponents belonging to the same irred factor will share a shift
       set */

    expos->ngroups = 0;
    for (slong i = 0; i < slfac->length; i++)
    {
        gr_poly_struct * indfactor = gr_poly_vec_entry_ptr(slfac, i, Scalars);
        expos->ngroups += indfactor->length - 1;
    }

    size_t shifts_offset = expos->ngroups * sizeof(acb_ode_group_struct);
    size_t defect = shifts_offset % alignof(acb_ode_shift_struct);
    if (defect != 0)
        shifts_offset += alignof(acb_ode_shift_struct) - defect;

    size_t bufsz = shifts_offset;
    for (slong i = 0; i < slshifts->length; i++)
    {
        gr_vec_struct * shi = GR_ENTRY(slshifts->entries, i, sizeof(gr_vec_struct));
        bufsz += shi->length * sizeof(acb_ode_shift_struct);
    }

    unsigned char * buf = flint_malloc(bufsz);

    expos->groups = (acb_ode_group_struct *) buf;
    acb_ode_shift_struct * shifts = (acb_ode_shift_struct *) (buf + shifts_offset);

    for (slong g = 0; g < expos->ngroups; g++)
        gr_init(expos->groups[g].leader, CC);

    /* Populate the data structure */

    for (slong i = 0, g = 0; i < slfac->length; i++)
    {
        gr_poly_struct * indfactor = gr_poly_vec_entry_ptr(slfac, i, Scalars);
        gr_vec_struct * shi = gr_vec_entry_ptr(slshifts, i, ZZvec);
        gr_vec_struct * mi = gr_vec_entry_ptr(slmult, i, ZZvec);

        status |= gr_poly_roots_other(roots, rtmult /* unused */,
                                      indfactor, Scalars, 0, CC);
        // if(status) flint_printf("%s:%d status=%d\n", __FILE__, __LINE__, status);
        FLINT_ASSERT(roots->length + 1 == indfactor->length || status != GR_SUCCESS);

        slong len = shi->length;
        slong maxshift = fmpz_get_si((fmpz *) shi->entries + len - 1);
        for (slong s = 0; s < len; s++)
        {
            /* The shifts in the shiftless decomposition apply to the variable
               of the indicial polynomial; we want shifts on the roots */
            shifts[s].n = maxshift - fmpz_get_si((fmpz *) shi->entries + len - 1 - s);
            shifts[s].mult = fmpz_get_si((fmpz *) mi->entries + len - 1 - s);
        }

        for (slong j = 0; j < roots->length; j++, g++)
        {
            status |= gr_sub_si(expos->groups[g].leader,
                                GR_ENTRY(roots->entries, j, CC->sizeof_elem),
                                maxshift, CC);
            // if(status) flint_printf("%s:%d status=%d\n", __FILE__, __LINE__, status);
            /* XXX exponents share shift pointers but not shift counts... */
            expos->groups[g].nshifts = shi->length;
            expos->groups[g].shifts = shifts;
        }

        shifts += shi->length;
    }

    FLINT_ASSERT((unsigned char *) shifts == buf + bufsz);

    fmpz_vec_clear(rtmult);
    gr_vec_clear(roots, CC);
    gr_vec_clear(slmult, ZZvec);
    gr_vec_clear(slshifts, ZZvec);
    gr_poly_vec_clear(slfac, Scalars);
    fmpz_vec_clear(mult);
    gr_poly_vec_clear(fac, Scalars);
    gr_poly_clear(polc, Scalars);

    gr_ctx_clear(ZZvec);
    gr_ctx_clear(ZZ);

cleanup_ind:
    gr_poly_clear(ind, Scalars);

    return status;
}
