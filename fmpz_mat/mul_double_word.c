/*
    Copyright (C) 2010,2011,2018 Fredrik Johansson
    Copyright (C) 2016 Aaditya Thakkar
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

/*
    There is a branchy version because the k loop takes time to warmup.
    Once k is big enough, the loop runs fast and there is no point in checking
    for small inputs.
*/
#define FMPZ_MAT_MUL_4_BRANCHLESS_CUTOFF 16

/* 2x2 -> 4 signed addmul */
static void _do_row_22_4_signed_branchy(
    fmpz * CR,
    const mp_limb_t * AR,
    const mp_limb_t * B,
    slong br,
    slong bc)
{
    slong j, k, l;
    mp_limb_t s[4], t3, t2, t1, t0, w3, w2, w1, w0;
    mp_limb_t A0, A1, B0, B1;
    mp_limb_t u2, u1, u0;

    for (j = 0, l = 0; j < bc; j++)
    {
        t3 = t2 = t1 = t0 = 0;
        u2 = u1 = u0 = 0;

        for (k = 0; k < br; k++, l++)
        {
            A0 = AR[2*k + 0];
            A1 = AR[2*k + 1];
            B0 = B[2*l + 0];
            B1 = B[2*l + 1];

            if (FLINT_SIGN_EXT(A0) == A1 && FLINT_SIGN_EXT(B0) == B1)
            {
                smul_ppmm(w1, w0, B0, A0);
                add_sssaaaaaa(u2, u1, u0, u2, u1, u0,
                              FLINT_SIGN_EXT(w1), w1, w0);
            }
            else
            {
                sub_ddmmss(t3, t2, t3, t2, 0, FLINT_SIGN_EXT(A1)&B0);
                sub_ddmmss(t3, t2, t3, t2, 0, FLINT_SIGN_EXT(B1)&A0);

                smul_ppmm(w3, w2, B1, A1);
                add_ssaaaa(t3, t2, t3, t2, w3, w2);

                umul_ppmm(w1, w0, B0, A0);
                add_sssaaaaaa(u2, u1, u0, u2, u1, u0, UWORD(0), w1, w0);

                umul_ppmm(w2, w1, A1, B0);
                add_sssaaaaaa(t3, t2, t1, t3, t2, t1, UWORD(0), w2, w1);

                umul_ppmm(w2, w1, B1, A0);
                add_sssaaaaaa(t3, t2, t1, t3, t2, t1, UWORD(0), w2, w1);
            }
        }

        add_ssssaaaaaaaa(s[3], s[2], s[1], s[0], t3, t2, t1, t0,
                         FLINT_SIGN_EXT(u2), u2, u1, u0);

        fmpz_set_signed_ui_array(CR + j, s, 4);
    }
}

/* 2x2 -> 4 signed addmul */
static void _do_row_22_4_signed(
    fmpz * CR,
    const mp_limb_t * AR,
    const mp_limb_t * B,
    slong br,
    slong bc)
{
    slong j, k, l;
    mp_limb_t s[4], t3, t2, t1, t0, w3, w2, w1, w0;
    mp_limb_t A0, A1, B0, B1;
    mp_limb_t v3, v2, u2, u1, u0;

    for (j = 0, l = 0; j < bc; j++)
    {
        t3 = t2 = t1 = t0 = 0;
        u2 = u1 = u0 = 0;
        v3 = v2 = 0;

        for (k = 0; k < br; k++, l++)
        {
            A0 = AR[2*k + 0];
            A1 = AR[2*k + 1];
            B0 = B[2*l + 0];
            B1 = B[2*l + 1];

            sub_ddmmss(t3, t2, t3, t2, 0, FLINT_SIGN_EXT(A1)&B0);
            sub_ddmmss(t3, t2, t3, t2, 0, FLINT_SIGN_EXT(B1)&A0);

            smul_ppmm(w3, w2, B1, A1);
            add_ssaaaa(t3, t2, t3, t2, w3, w2);

            umul_ppmm(w1, w0, B0, A0);
            add_sssaaaaaa(u2, u1, u0, u2, u1, u0, UWORD(0), w1, w0);

            umul_ppmm(w2, w1, A1, B0);
            add_sssaaaaaa(t3, t2, t1, t3, t2, t1, UWORD(0), w2, w1);

            umul_ppmm(w2, w1, B1, A0);
            add_sssaaaaaa(v3, v2, t1, v3, v2, t1, UWORD(0), w2, w1);
        }

        add_ssssaaaaaaaa(t3, t2, t1, t0, t3, t2, t1, t0,
                         FLINT_SIGN_EXT(u2), u2, u1, u0);

        add_ssaaaa(t3, t2, t3, t2, v3, v2);

        s[3] = t3;
        s[2] = t2;
        s[1] = t1;
        s[0] = t0;

        fmpz_set_signed_ui_array(CR + j, s, 4);
    }
}


/* 2x2 -> 5 signed addmul */
static void _do_row_22_5_signed(
    fmpz * CR,
    const mp_limb_t * AR,
    const mp_limb_t * B,
    slong br,
    slong bc)
{
    slong j, k, l;
    mp_limb_t s[5];
    mp_limb_t A0, A1, B0, B1;
    mp_limb_t s4, s3, s2, s1, s0;
    mp_limb_t p3, p2, p1, p0;
    mp_limb_t u3, u2;
    mp_limb_t v2;

    for (j = 0, l = 0; j < bc; j++)
    {
        s4 = s3 = s2 = s1 = s0 = 0;
        u3 = u2 = 0;
        v2 = 0;

        for (k = 0; k < br; k++, l++)
        {
            A0 = AR[2*k + 0];
            A1 = AR[2*k + 1];
            B0 = B[2*l + 0];
            B1 = B[2*l + 1];

            sub_ddmmss(u3, u2, u3, u2, UWORD(0), FLINT_SIGN_EXT(A1)&B0);
            sub_ddmmss(u3, u2, u3, u2, UWORD(0), FLINT_SIGN_EXT(B1)&A0);

            smul_ppmm(p3, p2, B1, A1);
            add_sssaaaaaa(s4, s3, s2, s4, s3, s2, FLINT_SIGN_EXT(p3), p3, p2);

            umul_ppmm(p2, p1, A1, B0);
            add_sssaaaaaa(u3, u2, s1, u3, u2, s1, UWORD(0), p2, p1);

            umul_ppmm(p1, p0, B0, A0);
            add_sssaaaaaa(v2, s1, s0, v2, s1, s0, UWORD(0), p1, p0);

            umul_ppmm(p2, p1, B1, A0);
            add_sssaaaaaa(u3, u2, s1, u3, u2, s1, UWORD(0), p2, p1);
        }

        add_sssaaaaaa(s4, s3, s2, s4, s3, s2, FLINT_SIGN_EXT(u3), u3, u2);

        add_sssaaaaaa(s4, s3, s2, s4, s3, s2, UWORD(0), UWORD(0), v2);

        s[0] = s0;
        s[1] = s1;
        s[2] = s2;
        s[3] = s3;
        s[4] = s4;

        fmpz_set_signed_ui_array(CR + j, s, 5);
    }
}

/* 2x2 -> 4 unsigned addmul */
static void _do_row_22_4_unsigned(
    fmpz * CR,
    const mp_limb_t * AR,
    const mp_limb_t * B,
    slong br,
    slong bc)
{
    slong j, k, l;
    mp_limb_t s[4];
    mp_limb_t A0, A1, B0, B1;
    mp_limb_t p3, p2, p1, p0;
    mp_limb_t s3, s2, s1, s0;
    mp_limb_t u3, u2, u1;
    mp_limb_t v3, v2;

    for (j = 0, l = 0; j < bc; j++)
    {
        s3 = s2 = s1 = s0 = 0;
        u3 = u2 = u1 = 0;
        v3 = v2 = 0;

        for (k = 0; k < br; k++, l++)
        {
            A0 = AR[2*k + 0];
            A1 = AR[2*k + 1];
            B0 = B[2*l + 0];
            B1 = B[2*l + 1];

            umul_ppmm(p2, p1, A1, B0);
            add_sssaaaaaa(u3, u2, u1, u3, u2, u1, UWORD(0), p2, p1);

            umul_ppmm(p3, p2, A1, B1);
            umul_ppmm(p1, p0, A0, B0);
            add_ssssaaaaaaaa(s3, s2, s1, s0, s3, s2, s1, s0, p3, p2, p1, p0);

            umul_ppmm(p2, p1, A0, B1);
            add_sssaaaaaa(v3, v2, u1, v3, v2, u1, UWORD(0), p2, p1);
        }

        add_ssaaaa(u3, u2, u3, u2, v3, v2);
        add_sssaaaaaa(s3, s2, s1, s3, s2, s1, u3, u2, u1);

        s[0] = s0;
        s[1] = s1;
        s[2] = s2;
        s[3] = s3;

        fmpz_set_ui_array(CR + j, s, 4);
    }
}

/* 2x2 -> 5 unsigned addmul */
static void _do_row_22_5_unsigned(
    fmpz * CR,
    const mp_limb_t * AR,
    const mp_limb_t * B,
    slong br,
    slong bc)
{
    slong j, k, l;
    mp_limb_t s[5];
    mp_limb_t A0, A1, B0, B1;
    mp_limb_t p3, p2, p1, p0;
    mp_limb_t s4, s3, s2, s1, s0;
    mp_limb_t u2, u1;
    mp_limb_t v3, v2;

    for (j = 0, l = 0; j < bc; j++)
    {
        s4 = s3 = s2 = s1 = s0 = 0;
        u2 = u1 = 0;
        v3 = v2 = 0;

        for (k = 0; k < br; k++, l++)
        {
            A0 = AR[2*k + 0];
            A1 = AR[2*k + 1];
            B0 = B[2*l + 0];
            B1 = B[2*l + 1];

            umul_ppmm(p2, p1, A1, B0);
            add_sssaaaaaa(s3, s2, s1, s3, s2, s1, UWORD(0), p2, p1);

            umul_ppmm(p1, p0, B0, A0);
            add_sssaaaaaa(u2, u1, s0, u2, u1, s0, UWORD(0), p1, p0);

            umul_ppmm(p2, p1, B1, A0);
            add_sssaaaaaa(s3, s2, s1, s3, s2, s1, UWORD(0), p2, p1);

            umul_ppmm(p3, p2, B1, A1);
            add_sssaaaaaa(s4, v3, v2, s4, v3, v2, UWORD(0), p3, p2);
        }

        /* s3 is small, so this doesn't overflow */
        add_sssaaaaaa(s3, s2, s1, s3, s2, s1, UWORD(0), u2, u1);

        add_sssaaaaaa(s4, s3, s2, s4, s3, s2, UWORD(0), v3, v2);

        s[0] = s0;
        s[1] = s1;
        s[2] = s2;
        s[3] = s3;
        s[4] = s4;

        fmpz_set_ui_array(CR + j, s, 5);
    }
}


typedef struct {
    slong Astartrow;
    slong Astoprow;
    slong Bstartcol;
    slong Bstopcol;
    slong br;
    slong bc;
    fmpz ** Crows;
    fmpz ** Arows;
    fmpz ** Brows;
    mp_limb_t * BL;
    int sign;
    int words;
} _worker_arg;

static void _red_worker(void * varg)
{
    _worker_arg * arg = (_worker_arg *) varg;
    slong Bstartcol = arg->Bstartcol;
    slong Bstopcol = arg->Bstopcol;
    slong br = arg->br;
    fmpz ** Brows = arg->Brows;
    mp_limb_t * BL = arg->BL;
    int sign = arg->sign;
    slong i, j;

    if (sign)
    {
        for (j = Bstartcol; j < Bstopcol; j++)
            for (i = 0; i < br; i++)
                fmpz_get_signed_uiui(BL + 2*(j*br + i) + 1,
                                     BL + 2*(j*br + i) + 0, &Brows[i][j]);
    }
    else
    {
        for (j = Bstartcol; j < Bstopcol; j++)
            for (i = 0; i < br; i++)
                fmpz_get_uiui(BL + 2*(j*br + i) + 1,
                              BL + 2*(j*br + i) + 0, &Brows[i][j]);
    }
}

static void _mul_worker(void * varg)
{
    _worker_arg * arg = (_worker_arg *) varg;
    slong Astartrow = arg->Astartrow;
    slong Astoprow = arg->Astoprow;
    slong ac = arg->br;
    slong br = arg->br;
    slong bc = arg->bc;
    fmpz ** Crows = arg->Crows;
    fmpz ** Arows = arg->Arows;
    mp_limb_t * BL = arg->BL;
    int sign = arg->sign;
    int words = arg->words;
    mp_limb_t * AL;
    slong i, j;
    TMP_INIT;

    FLINT_ASSERT(sign == 0 || sign == 1);
    FLINT_ASSERT(words == 4 || words == 5);

    TMP_START;

    AL = TMP_ARRAY_ALLOC(2*ac, mp_limb_t);

    if (sign)
    {
        for (i = Astartrow; i < Astoprow; i++)
        {
            for (j = 0; j < ac; j++)
                fmpz_get_signed_uiui(AL + 2*j + 1, AL + 2*j, &Arows[i][j]);

            if (words == 4 && br < FMPZ_MAT_MUL_4_BRANCHLESS_CUTOFF)
                _do_row_22_4_signed_branchy(Crows[i], AL, BL, br, bc);
            else if (words == 4)
                _do_row_22_4_signed(Crows[i], AL, BL, br, bc);
            else
                _do_row_22_5_signed(Crows[i], AL, BL, br, bc);
        }
    }
    else
    {
        for (i = Astartrow; i < Astoprow; i++)
        {
            for (j = 0; j < ac; j++)
                fmpz_get_uiui(AL + 2*j + 1, AL + 2*j, &Arows[i][j]);

            if (words == 4)
                _do_row_22_4_unsigned(Crows[i], AL, BL, br, bc);
            else
                _do_row_22_5_unsigned(Crows[i], AL, BL, br, bc);
        }
    }

    TMP_END;
}


/*
    sign = 1:   max|A|, max|B| < 2^(2*FLINT_BITS - 1)
                max|C| < 2^bits

    sign = 0:   all entries are >= 0 and
                max|A|, max|B| < 2^(2*FLINT_BITS)
                max|C| < 2^bits

    bits + sign <= 5*FLINT_BITS
*/
FLINT_DLL void _fmpz_mat_mul_double_word_internal(
    fmpz_mat_t C,
    const fmpz_mat_t A,
    const fmpz_mat_t B,
    int sign,
    flint_bitcnt_t bits)
{
    slong i;
    slong ar = fmpz_mat_nrows(A);
    slong br = fmpz_mat_nrows(B);
    slong bc = fmpz_mat_ncols(B);
    _worker_arg mainarg;
    thread_pool_handle * handles;
    slong num_workers;
    _worker_arg * args;
    slong limit;
    TMP_INIT;

    FLINT_ASSERT(ar > 0);
    FLINT_ASSERT(br > 0);
    FLINT_ASSERT(bc > 0);

    TMP_START;

    /* limit on number of threads */
    limit = FLINT_MAX(br, bc);
    limit = FLINT_MIN(limit, ar);
    limit = limit < 16 ? 0 : (limit - 16)/16;

    mainarg.Astartrow = 0;
    mainarg.Astoprow = ar;
    mainarg.Bstartcol = 0;
    mainarg.Bstopcol = bc;
    mainarg.br = br;
    mainarg.bc = bc;
    mainarg.Crows = C->rows;
    mainarg.Arows = A->rows;
    mainarg.Brows = B->rows;
    mainarg.BL = TMP_ARRAY_ALLOC(br*bc*2, mp_limb_t);
    mainarg.sign = sign;

    if (bits + sign <= 4*FLINT_BITS)
        mainarg.words = 4;
    else
        mainarg.words = 5;

    if (limit < 2)
    {
use_one_thread:
        _red_worker(&mainarg);
        _mul_worker(&mainarg);
        TMP_END;
        return;
    }

    num_workers = flint_request_threads(&handles, limit);
    if (num_workers < 1)
    {
        flint_give_back_threads(handles, num_workers);
        goto use_one_thread;
    }

    args = FLINT_ARRAY_ALLOC(num_workers, _worker_arg);

    for (i = 0; i < num_workers; i++)
    {
        args[i].Astartrow = (i + 0)*ar/(num_workers + 1);
        args[i].Astoprow = (i + 1)*ar/(num_workers + 1);
        args[i].Bstartcol = (i + 0)*bc/(num_workers + 1);
        args[i].Bstopcol = (i + 1)*bc/(num_workers + 1);
        args[i].br = mainarg.br;
        args[i].bc = mainarg.bc;
        args[i].Crows = mainarg.Crows;
        args[i].Arows = mainarg.Arows;
        args[i].Brows = mainarg.Brows;
        args[i].BL = mainarg.BL;
        args[i].sign = mainarg.sign;
        args[i].words = mainarg.words;
    }

    i = num_workers;
    mainarg.Astartrow = (i + 0)*ar/(num_workers + 1);
    mainarg.Astoprow = (i + 1)*ar/(num_workers + 1);
    mainarg.Bstartcol = (i + 0)*bc/(num_workers + 1);
    mainarg.Bstopcol = (i + 1)*bc/(num_workers + 1);

    for (i = 0; i < num_workers; i++)
        thread_pool_wake(global_thread_pool, handles[i], 0, _red_worker, &args[i]);
    _red_worker(&mainarg);
    for (i = 0; i < num_workers; i++)
        thread_pool_wait(global_thread_pool, handles[i]);

    for (i = 0; i < num_workers; i++)
        thread_pool_wake(global_thread_pool, handles[i], 0, _mul_worker, &args[i]);
    _mul_worker(&mainarg);
    for (i = 0; i < num_workers; i++)
        thread_pool_wait(global_thread_pool, handles[i]);

    flint_give_back_threads(handles, num_workers);
    flint_free(args);

    TMP_END;
    return;
}

void _fmpz_mat_mul_double_word(
    fmpz_mat_t C,
    const fmpz_mat_t A,
    const fmpz_mat_t B)
{
    int sign = 0;
    slong Abits, Bbits;
    flint_bitcnt_t Cbits;

    if (fmpz_mat_is_empty(A) || fmpz_mat_is_empty(B))
    {
        fmpz_mat_zero(C);
        return;
    }

    Abits = fmpz_mat_max_bits(A);
    Bbits = fmpz_mat_max_bits(B);

    if (Abits < 0)
    {
        sign = 1;
        Abits = -Abits;
    }

    if (Bbits < 0)
    {
        sign = 1;
        Bbits = -Bbits;
    }

    FLINT_ASSERT(Abits + sign <= 2*FLINT_BITS);
    FLINT_ASSERT(Bbits + sign <= 2*FLINT_BITS);

    Cbits = FLINT_BIT_COUNT(A->c) + Abits + Bbits;

    _fmpz_mat_mul_double_word_internal(C, A, B, sign, Cbits);
}

