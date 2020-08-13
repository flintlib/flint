/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/


#include "mpoly.h"

/*
    a and b are arrays of packed monomials
    define
        score(e) = (# cross products in a X b <= e)
    finds a monomial e such that
        lower <= score(e) <= upper
        or is as close as possible and store its score in e_score
    e is written in the same format as a and b
    three arrays find, gind, and hind, each of length a_len
        are need for working space
    the return pointer is one of find, gind, or hind
        the elements of this array are indices into b where the first monomial
        <= e was found
*/
void mpoly_search_monomials(
    slong ** e_ind, ulong * e, slong * e_score,
    slong * find, slong * gind, slong * hind,
    slong lower, slong upper,
    const ulong * a, slong a_len, const ulong * b, slong b_len,
    slong N, const ulong * cmpmask)
{
    slong i, j, x;
    slong maxdiff, maxind;
    /*
        for each i, there is an integer 0 <= find[i] <= blen such that
            a[i] + b[find[i]-1] < fexp <= a[i] + b[find[i]]
        ( If fexp < a[0] + b[blen-1] then find[i] is blen.
          Similaryly if fexp >= a[i] +  b[0], then find[i] is 0 )
        fscore is score(fexp)
        ditto for g and h

        We always maintain paths f, h, g with corresponding exponents
        fexp > hexp > gexp. These paths are non-increasing. Example:

            b_len |
                  |
             g=>  |___________
                  |           _
                  |            _
             h => |_            _______
                  |  ______________
                  |               _____
             f => |_______
                  |       _______
                0 +--------------______
                  0                   a_len
    */
    slong fscore, gscore, hscore, tscore;
    ulong * fexp, * gexp, * hexp, * texp;
    slong * tind;
    ulong * temp_exp;

    FLINT_ASSERT(a_len > 0);
    FLINT_ASSERT(b_len > 0);
    FLINT_ASSERT(lower <= upper);

    /* set f to correspond to an upperbound on all products */
    fscore = a_len * b_len;
    fexp = (ulong *) flint_malloc(N*sizeof(ulong));
    mpoly_monomial_add_mp(fexp, a + 0*N, b + 0*N, N);
    for (i = 0; i < a_len; i++)
        find[i] = 0;

    /* set g to correspond to a lowerbound on all products */
    gscore = 1;
    gexp = (ulong *) flint_malloc(N*sizeof(ulong));
    mpoly_monomial_add_mp(gexp, a + (a_len - 1)*N, b + (b_len - 1)*N, N);
    for (i = 0; i < a_len; i++)
        gind[i] = b_len;
    gind[a_len - 1] = b_len - 1;

    /* just allocate h */
    hexp = (ulong *) flint_malloc(N*sizeof(ulong));

    temp_exp = (ulong *) flint_malloc(N*sizeof(ulong));

    /* early exit */
    if (fscore == gscore)
        goto return_f;

    /* main loop */
    while (gscore < lower && upper < fscore)
    {
        /* find the index 'maxind' where gind[i] - find[i] is largest */
        maxdiff = -1;
        maxind = -1;
        for (i = 0; i < a_len; i++)
        {
            if (maxdiff < gind[i] - find[i])
            {
                maxdiff = gind[i] - find[i];
                maxind = i;
            }
        }

        if (maxdiff == 0)
        {
            /* f and g are the same path */
            break;

        } else if (maxdiff == 1)
        {
            /* there may or may not be another path between */
            maxind = -1;
            for (i = 0; i < a_len; i++)
            {
                if (gind[i] > find[i])
                {
                    mpoly_monomial_add_mp(temp_exp, a + i*N, b + find[i]*N, N);
                    if (mpoly_monomial_equal(temp_exp, fexp, N) == 0)
                    {
                        maxind = i;
                        hind[maxind] = find[i];
                        mpoly_monomial_add_mp(hexp, a + maxind*N,
                                                        b + hind[maxind]*N, N);
                    }
                }
            }
            if (maxind == -1)
                /* there is no path between */
                break;

        } else
        {
            /* there is definitely a path between */
            hind[maxind] = (gind[maxind] + find[maxind])/2;    
        }


        /*
            the point (maxind, hind[maxind)) is now set to a bisector
            get the corresponding monomial into hexp
        */
        mpoly_monomial_add_mp(hexp, a + maxind*N, b + hind[maxind]*N, N);

        FLINT_ASSERT(mpoly_monomial_lt(hexp, fexp, N, cmpmask));
        FLINT_ASSERT(mpoly_monomial_lt(gexp, hexp, N, cmpmask));

        /*
            find new path for h through the point
        */
        hscore = gscore + gind[maxind] - hind[maxind];

        /*
            find new path for h to the right of the point
        */
        for (i = maxind + 1; i < a_len; i++)
        {
            x = find[i];
            for (j = FLINT_MIN(hind[i-1], gind[i]) - 1; j >= find[i]; j--)
            {
                mpoly_monomial_add_mp(temp_exp, a + i*N, b + j*N, N);

                if (mpoly_monomial_lt(hexp, temp_exp, N, cmpmask))
                {
                    x = j + 1;
                    break;
                }
            }
            hind[i] = x;
            hscore += gind[i] - hind[i];            
        }

        /*
            find new path for h to the left of the point
        */
        for (i = maxind - 1; i >= 0; i--)
        {
            x = FLINT_MAX(hind[i+1], find[i]);
            for (j = FLINT_MAX(hind[i+1], find[i]); j < gind[i]; j++)
            {
                mpoly_monomial_add_mp(temp_exp, a + i*N, b + j*N, N);
                if (mpoly_monomial_lt(hexp, temp_exp, N, cmpmask))
                    x = j + 1;
                else
                    break;
            }
            hind[i] = x;
            hscore += gind[i] - hind[i];            
        }

        if (hscore <= upper) 
        {
            tind = gind; tscore = gscore; texp = gexp;
            gind = hind; gscore = hscore; gexp = hexp;
            hind = tind; hscore = tscore; hexp = texp;
        } else {
            tind = find; tscore = fscore; texp = fexp;
            find = hind; fscore = hscore; fexp = hexp;
            hind = tind; hscore = tscore; hexp = texp;
        }
    }


    /* upper and lower bounds are out of range */
    if (fscore <= lower)
        goto return_f;
    else if (gscore >= upper)
        goto return_g;
    /* found something in range */
    else if (fscore <= upper)
        goto return_f;
    else if (gscore >= lower)
        goto return_g;
    /* could not get in range - choose closest one */
    else if (fscore - upper < lower - gscore)
        goto return_f;
    else
        goto return_g;

return_g:
    mpoly_monomial_set(e, gexp, N);
    *e_score = gscore;
    tind = gind;
    goto cleanup;

return_f:
    mpoly_monomial_set(e, fexp, N);
    *e_score = fscore;
    tind = find;

cleanup:
    flint_free(temp_exp);
    flint_free(hexp);
    flint_free(gexp);
    flint_free(fexp);
    * e_ind = tind;
}
