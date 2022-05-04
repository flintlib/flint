/*
    Copyright (C) 2008, 2009, William Hart 
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "nmod_vec.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

/* The better mpn_addmul_1 is, the larger FMPZ_CRT_UI_CUTOFF can be. */
#define FMPZ_CRT_UI_CUTOFF 50

/*
    The better fmpz_fdiv_ui (mpn_mod_1) is, the larger FMPZ_MOD_UI_CUTOFF can
    be. Very poor implementations of mpn_mod_1 can have 10, and anything up to
    100 is reasonable.
*/
#define FMPZ_MOD_UI_CUTOFF 75

/* minimal number of basecase chunks in which to partition */
#define FMPZ_CRT_UI_MULTIPLE_CUTOFF 4
#define FMPZ_MOD_UI_MULTIPLE_CUTOFF 4


void fmpz_comb_temp_init(fmpz_comb_temp_t CT, const fmpz_comb_t C)
{
    CT->Alen = FLINT_MAX(C->crt_klen, C->mod_klen);
    CT->Tlen = FLINT_MAX(C->crt_P->localsize, C->mod_P->localsize);

    CT->A = _fmpz_vec_init(CT->Alen);
    CT->T = _fmpz_vec_init(CT->Tlen);
}


void fmpz_comb_init(fmpz_comb_t C, mp_srcptr m, slong len)
{
    int success;
    slong l, i, j, k, s;
    ulong tt, mm, mt;
    fmpz_poly_t M, Mm; /* only used for resizable fmpz array convenience */

    if (len < 1)
        flint_throw(FLINT_ERROR, "fmpz_comb_init: len should be positive");

    fmpz_poly_init(Mm);
    fmpz_poly_init(M);

    fmpz_multi_CRT_init(C->crt_P);
    C->packed_multipliers = NULL;
    C->step = NULL;
    C->crt_lu = NULL;
    C->crt_lu_alloc = 0;
    C->crt_offsets = NULL;
    C->crt_offsets_alloc = 0;

    fmpz_multi_mod_init(C->mod_P);
    C->mod_lu = NULL;
    C->mod_lu_alloc = 0;
    C->mod_offsets = NULL;
    C->mod_offsets_alloc = 0;

    C->num_primes = len;

    /* set up crt */

    for (l = 0, i = 0, k = 0; l < len; k++)
    {
        if (k + 1 >= C->crt_offsets_alloc)
        {
            C->crt_offsets_alloc = FLINT_MAX(k + 1, C->crt_offsets_alloc*3/2);
            C->crt_offsets = FLINT_ARRAY_REALLOC(C->crt_offsets,
                                                  C->crt_offsets_alloc, slong);
        }

        fmpz_poly_fit_length(M, k + 1);
        fmpz_one(M->coeffs + k);

        for (j = i;
             l < len && fmpz_size(M->coeffs + k) <= FMPZ_CRT_UI_CUTOFF;
             j++)
        {
            if (j + 1 >= C->crt_lu_alloc)
            {
                C->crt_lu_alloc = FLINT_MAX(j + 1, C->crt_lu_alloc*3/2);
                C->crt_lu = FLINT_ARRAY_REALLOC(C->crt_lu, C->crt_lu_alloc,
                                                                crt_lut_entry);
            }

            C->crt_lu[j].i0 = 0;
            C->crt_lu[j].i1 = 0;
            C->crt_lu[j].i2 = 0;

            mm = m[l];
            if (l + 1 < len && !n_mul_checked(&mt, mm, m[l+1]))
            {
                mm = mt;
                if (l + 2 < len && !n_mul_checked(&mt, mm, m[l+2]))
                {
                    mm = mt;
                    success = (1 == n_gcdinv(&C->crt_lu[j].i0,
                                          m[l+1]*m[l+2] % m[l+0], m[l+0])) &&
                              (1 == n_gcdinv(&C->crt_lu[j].i1,
                                          m[l+0]*m[l+2] % m[l+1], m[l+1])) &&
                              (1 == n_gcdinv(&C->crt_lu[j].i2,
                                          m[l+0]*m[l+1] % m[l+2], m[l+2]));
                    if (!success)
                        goto bad_moduli;

                    C->crt_lu[j].i0 *= m[l+1]*m[l+2];
                    C->crt_lu[j].i1 *= m[l+0]*m[l+2];
                    C->crt_lu[j].i2 *= m[l+0]*m[l+1];
                    l += 3;
                }
                else
                {
                    success = (1 == n_gcdinv(&C->crt_lu[j].i0,
                                             m[l+1] % m[l+0], m[l+0])) &&
                              (1 == n_gcdinv(&C->crt_lu[j].i1,
                                             m[l+0] % m[l+1], m[l+1]));
                    if (!success)
                        goto bad_moduli;

                    C->crt_lu[j].i0 *= m[l+1];
                    C->crt_lu[j].i1 *= m[l+0];
                    l += 2;
                }
            }
            else
            {
                C->crt_lu[j].i0 = 1;
                l += 1;
            }

            nmod_init(&C->crt_lu[j].mod, mm);
            fmpz_mul_ui(M->coeffs + k, M->coeffs + k, C->crt_lu[j].mod.n);
        }

        C->crt_offsets[k] = j;

        fmpz_poly_fit_length(Mm, j);

        for ( ; i < j; i++)
        {
            fmpz_divexact_ui(Mm->coeffs + i, M->coeffs + k, C->crt_lu[i].mod.n);
            tt = fmpz_fdiv_ui(Mm->coeffs + i, C->crt_lu[i].mod.n);

            success = (1 == n_gcdinv(&tt, tt, C->crt_lu[i].mod.n));
            if (!success)
                goto bad_moduli;

            C->crt_lu[i].i0 = nmod_mul(tt, C->crt_lu[i].i0, C->crt_lu[i].mod);
            C->crt_lu[i].i1 = nmod_mul(tt, C->crt_lu[i].i1, C->crt_lu[i].mod);
            C->crt_lu[i].i2 = nmod_mul(tt, C->crt_lu[i].i2, C->crt_lu[i].mod);
        }
    }

    /*
        avoid small last chunk
        and have at least FMPZ_CRT_UI_MULTIPLE_CUTOFF chunks or one big chunk
    */
    while (k > 1 && (k < FMPZ_CRT_UI_MULTIPLE_CUTOFF ||
                     fmpz_size(M->coeffs + k - 1) <= FMPZ_CRT_UI_CUTOFF*3/4))
    {
        k--;
        for (i = (k >= 2) ? C->crt_offsets[k - 2] : 0;
             i < C->crt_offsets[k];
             i++)
        {
            int last = (i >= C->crt_offsets[k - 1]);

            fmpz_mul(Mm->coeffs + i, Mm->coeffs + i, M->coeffs + k - last);
            tt = fmpz_fdiv_ui(M->coeffs + k - last, C->crt_lu[i].mod.n);

            success = (1 == n_gcdinv(&tt, tt, C->crt_lu[i].mod.n));
            if (!success)
                goto bad_moduli;

            C->crt_lu[i].i0 = nmod_mul(tt, C->crt_lu[i].i0, C->crt_lu[i].mod);
            C->crt_lu[i].i1 = nmod_mul(tt, C->crt_lu[i].i1, C->crt_lu[i].mod);
            C->crt_lu[i].i2 = nmod_mul(tt, C->crt_lu[i].i2, C->crt_lu[i].mod);            
        }

        C->crt_offsets[k - 1] = C->crt_offsets[k];
        fmpz_mul(M->coeffs + k - 1, M->coeffs + k - 1, M->coeffs + k);
    }

    C->crt_klen = k;

    success = fmpz_multi_CRT_precompute(C->crt_P, M->coeffs, C->crt_klen);
    if (!success)
        goto bad_moduli;

    C->step = FLINT_ARRAY_ALLOC(C->crt_klen, slong);

    l = 0;
    for (k = 0, i = 0; k < C->crt_klen; k++)
    {
        int all_large = 1;

        s = 1;
        for (j = i; j < C->crt_offsets[k]; j++)
        {
            if (C->crt_lu[j].i1 != 0 || C->crt_lu[j].i2 != 0)
                all_large = 0;

            FLINT_ASSERT(fmpz_cmp(Mm->coeffs + j, M->coeffs + k) <= 0);
            s = FLINT_MAX(s, fmpz_size(Mm->coeffs + j));
        }

        if (all_large)
        {
            s = 1;
            for (j = i; j < C->crt_offsets[k]; j++)
            {
                FLINT_ASSERT(C->crt_lu[i].i1 == 0 && C->crt_lu[i].i2 == 0);
                fmpz_mul_ui(Mm->coeffs + j, Mm->coeffs + j, C->crt_lu[j].i0);
                s = FLINT_MAX(s, fmpz_size(Mm->coeffs + j));
            }

            C->step[k] = -s - 1;
        }
        else
        {
            C->step[k] = s;
        }

        for ( ; i < C->crt_offsets[k]; i++)
        {
            l += s;
        }
    }

    C->packed_multipliers = FLINT_ARRAY_ALLOC(l, mp_limb_t);

    l = 0;
    for (k = 0, i = 0; k < C->crt_klen; k++)
    {
        s = C->step[k];
        if (s < 0)
            s = -s - 1;

        for ( ; i < C->crt_offsets[k]; i++)
        {
            fmpz_get_ui_array(C->packed_multipliers + l, s, Mm->coeffs + i);
            l += s;
        }
    }

    /* set up mod */

    for (l = 0, i = 0, k = 0; l < len; k++)
    {
        if (k + 1 >= C->mod_offsets_alloc)
        {
            C->mod_offsets_alloc = FLINT_MAX(k + 1, C->mod_offsets_alloc*3/2);
            C->mod_offsets = FLINT_ARRAY_REALLOC(C->mod_offsets,
                                                  C->mod_offsets_alloc, slong);
        }

        fmpz_poly_fit_length(M, k + 1);
        fmpz_one(M->coeffs + k);

        for (j = i; l < len && fmpz_size(M->coeffs + k) <= FMPZ_MOD_UI_CUTOFF; j++)
        {
            if (j + 1 >= C->mod_lu_alloc)
            {
                C->mod_lu_alloc = FLINT_MAX(j + 1, C->mod_lu_alloc*3/2);
                C->mod_lu = FLINT_ARRAY_REALLOC(C->mod_lu, C->mod_lu_alloc,
                                                                mod_lut_entry);
            }

            C->mod_lu[j].mod0.n = 0;
            C->mod_lu[j].mod1.n = 0;
            C->mod_lu[j].mod2.n = 0;

            mm = m[l];
            if (l + 1 < len && !n_mul_checked(&mt, mm, m[l+1]))
            {
                mm = mt;
                if (l + 2 < len && !n_mul_checked(&mt, mm, m[l+2]))
                {
                    mm = mt;

                    success = (m[l+0] != 0) && (m[l+1] != 0) && (m[l+2] != 0);
                    if (!success)
                            goto zero_moduli;

                    nmod_init(&C->mod_lu[j].mod0, m[l+0]);
                    nmod_init(&C->mod_lu[j].mod1, m[l+1]);
                    nmod_init(&C->mod_lu[j].mod2, m[l+2]);
                    l += 3;
                }
                else
                {
                    success = (m[l+0] != 0) && (m[l+1] != 0);
                    if (!success)
                        goto zero_moduli;

                    nmod_init(&C->mod_lu[j].mod0, m[l+0]);
                    nmod_init(&C->mod_lu[j].mod1, m[l+1]);
                    l += 2;
                }
            }
            else
            {
                success = (m[l+0] != 0);
                if (!success)
                    goto zero_moduli;

                nmod_init(&C->mod_lu[j].mod0, m[l+0]);
                l += 1;
            }

            nmod_init(&C->mod_lu[j].mod, mm);
            fmpz_mul_ui(M->coeffs + k, M->coeffs + k, C->mod_lu[j].mod.n);
        }

        C->mod_offsets[k] = j;

        i = j;
    }

    /*
       avoid small last chunk
       and have at least FMPZ_MOD_UI_MULTIPLE_CUTOFF chunks or one big chunk
    */
    while (k > 1 && (k < FMPZ_MOD_UI_MULTIPLE_CUTOFF ||
                     fmpz_size(M->coeffs + k - 1) < FMPZ_MOD_UI_CUTOFF*3/4))
    {
        k--;
        C->mod_offsets[k - 1] = C->mod_offsets[k];
        fmpz_mul(M->coeffs + k - 1, M->coeffs + k - 1, M->coeffs + k);
    }

    C->mod_klen = k;

    success = fmpz_multi_mod_precompute(C->mod_P, M->coeffs, C->mod_klen);
    if (!success)
        goto zero_moduli;

    fmpz_poly_clear(M);
    fmpz_poly_clear(Mm);

    return;

bad_moduli:

    flint_throw(FLINT_ERROR, "fmpz_comb_init: moduli are not pairwise prime");

zero_moduli:

    flint_throw(FLINT_ERROR, "fmpz_comb_init: moduli are not nonzero");
}

