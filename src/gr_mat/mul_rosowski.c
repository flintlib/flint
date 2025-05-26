/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_mat.h"
#include "gr_vec.h"

/*
Reference implementation in Python:

def mul(A, B):
    l = A.rows; n = A.cols; m = B.cols
    if n <= 1 or l <= 1 or m <= 1:
        return A * B
    if n % 2 == 0:
        C = matrix(l, m)
        h = n // 2
        muls = 0
        T = [0] * (h * (m - 1))
        for j in range(1, m):
            for k in range(h):
                T[(j - 1) * h + k] = B[2*k, 0] + B[2*k, j]
        for i in range(l):
            for j in range(1, m):
                C[i, j] = sum((A[i, 2*k] + B[2*k+1, j]) * (A[i, 2*k+1] + T[(j - 1) * h + k]) for k in range(h))
                muls += h
        for j in range(1, m):
            v = sum(B[2*k+1, j] * (B[2*k, 1-1] + B[2*k, j]) for k in range(h))
            muls += h
            for i in range(l):
                C[i, j] -= v
        for i in range(l):
            v = sum(A[i, 2*k] * (B[2*k, 0] + A[i, 2*k+1]) for k in range(h))
            muls += h
            for j in range(1, m):
                C[i, j] -= v
            v += sum(A[i, 2*k+1] * (B[2*k+1, 0] - A[i, 2*k]) for k in range(h))
            muls += h
            C[i, 0] = v
        assert muls == (n*(l*m+l+m-1))//2
    elif m == 2:
        # l x n x m  -> (l x (n - 1) x 2) + (l x 1 x m)
        A1 = A[:,:1]; A2 = A[:,1:]; B1 = B[:1,:]; B2 = B[1:,:]
        C = mul(A2, B2)
        # C += mul(A1, B1)
        for i in range(l):
            for j in range(2):
                C[i, j] += A[i, 0] * B[0, j]
        return C
    else:
        # l x n x m  -> (l x (n - 3) x m) + (l x 3 x m)
        A1 = A[:,:3]; A2 = A[:,3:]; B1 = B[:3,:]; B2 = B[3:,:]
        C = mul(A2, B2)
        # C += mul(A1, B1)
        def a(i, j):
            return A[i-1, j-1]
        def b(i, j):
            return B[i-1, j-1]
        muls = 0
        m0 = 4 if m % 2 else 5
        t0 = b(1,2) * b(2,1); t1 = b(1,3) * b(3,1); t2 = b(2,3) * b(3,2)
        t3 = t0 + t1; t4 = t0 + t2; t5 = t1 + t2
        t6 = b(1,1) - b(1,2) - b(1,3)
        t7 = b(2,2) - b(2,1) - b(2,3)
        t8 = b(3,3) - b(3,1) - b(3,2)
        muls += 3
        if m % 2 == 0:
            u0 = t0 + (b(2,1) - b(2,4)) * (b(1,4) - b(1,2))
            u1 = b(1,4) - b(1,2); u2 = b(2,1) - b(2,4)
            muls += 1
        v = [b(1,j) - b(1,j+1) for j in range(m0, m, 2)]
        w = [b(3,j) - b(3,j+1) for j in range(m0, m, 2)]
        assert len(v) == (m - 3) // 2
        assert len(w) == (m - 3) // 2
        for i in range(1, l+1):
            q0 = (a(i,1) + b(2,1)) * (a(i,2) + b(1,2))
            q1 = (a(i,1) + b(3,1)) * (a(i,3) + b(1,3))
            q2 = (a(i,2) + b(3,2)) * (a(i,3) + b(2,3))
            C[i-1,0] += q0 + q1 + a(i,1) * (t6 - a(i,2) - a(i,3)) - t3
            C[i-1,1] += q0 + q2 + a(i,2) * (t7 - a(i,1) - a(i,3)) - t4
            C[i-1,2] += q1 + q2 + a(i,3) * (t8 - a(i,1) - a(i,2)) - t5
            muls += 6
            if m % 2 == 0:
                C[i-1,3] += q0 + (a(i,1) + u2) * (u1 - a(i,2)) + a(i,3) * b(3,4) - u0
                muls += 2
            r0 = q0 + q1 - t3; r1 = q1 + q2 - t5
            r2 = a(i,1) + b(2,1); r3 = a(i,3) + b(1,3); r4 = a(i,1) + b(3,1)
            r5 = a(i,2) + b(1,2); r6 = a(i,2) + b(3,2); r7 = a(i,3) + b(2,3)
            for j in range(m0, m, 2):
                r8 = (r4 - b(3,j)) * (b(1,j+1) - r3)
                C[i-1,j-1] += r0 + r8 + (r2 - b(2,j)) * (v[(j - 4) // 2] - r5)
                C[i-1,j] += r1 + r8 + (r6 + w[(j - 4) // 2]) * (b(2,j+1) - r7)
                muls += 3
        for j in range(m0, m, 2):
            q0 = (b(3,1) - b(3,j)) * (b(1,3) - b(1,j+1))
            q1 = (b(2,j) - b(2,1)) * (v[(j - 4) // 2] - b(1,2)) + q0
            q2 = (b(3,2) + w[(j - 4) // 2]) * (b(2,3) - b(2,j+1)) + q0
            muls += 3
            for i in range(1, l+1):
                C[i-1,j-1] += q1
                C[i-1,j] += q2
        if m % 2:
            wanted = 3 * (l*m + l + m - 1) // 2
        else:
            wanted = 2*(l-1) + 3*(l*m+m)//2
        assert muls == wanted
    return C
*/

int
gr_mat_mul_rosowski(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    slong Cstride = C->stride;
    slong Astride = A->stride;
    slong Bstride = B->stride;
    gr_ptr cc = C->entries;
    gr_srcptr aa = A->entries;
    gr_srcptr bb = B->entries;
    slong alloc;
    gr_ptr tmp, t, u, v, w;
    slong l, n, m;
    slong h, i, j, k;

    l = gr_mat_nrows(A, ctx);
    n = gr_mat_ncols(A, ctx);
    m = gr_mat_ncols(B, ctx);

    if (l <= 1 || n <= 1 || m <= 1)
        return gr_mat_mul_classical(C, A, B, ctx);

    if (C == A || C == B)
    {
        gr_mat_t T;
        gr_mat_init(T, l, m, ctx);
        status |= gr_mat_mul_rosowski(T, A, B, ctx);
        status |= gr_mat_swap_entrywise(T, C, ctx);
        gr_mat_clear(T, ctx);
        return status;
    }

#define AA(ii, jj) GR_ENTRY(aa, (ii) * (Astride) + (jj), sz)
#define BB(ii, jj) GR_ENTRY(bb, (ii) * (Bstride) + (jj), sz)
#define CC(ii, jj) GR_ENTRY(cc, (ii) * (Cstride) + (jj), sz)

    if (n % 2 == 0)
    {
        h = n / 2;

        alloc = h * (m - 1) + 3;
        GR_TMP_INIT_VEC(tmp, alloc, ctx);
        t = tmp;
        u = GR_ENTRY(t, 1, sz);
        v = GR_ENTRY(u, 1, sz);
        w = GR_ENTRY(v, 1, sz);

#define WW(ii, jj) GR_ENTRY(w, (ii) * (h) + (jj), sz)

        /* Use O(h*m) temp space instead of O(1) to save l*h*(m-1) additions. */
        for (j = 1; j < m; j++)
        {
            for (k = 0; k < h; k++)
            {
                /* W[(j - 1) * h + k] = B[2*k, 0] + B[2*k, j] */
                status |= gr_add(WW(j - 1, k), BB(2 * k, 0), BB(2 * k, j), ctx);
            }
        }

        for (i = 0; i < l; i++)
        {
            for (j = 1; j < m; j++)
            {
                for (k = 0; k < h; k++)
                {
                    status |= gr_add(t, AA(i, 2 * k), BB(2 * k + 1, j), ctx);
                    status |= gr_add(u, AA(i, 2 * k + 1), WW(j - 1, k), ctx);

                    if (k == 0)
                        status |= gr_mul(CC(i, j), t, u, ctx);
                    else
                        status |= gr_addmul(CC(i, j), t, u, ctx);
                }
            }
        }

        for (j = 1; j < m; j++)
        {
            for (k = 0; k < h; k++)
            {
                status |= gr_add(u, BB(2 * k, 0), BB(2 * k, j), ctx);

                if (k == 0)
                    status |= gr_mul(v, BB(2 * k + 1, j), u, ctx);
                else
                    status |= gr_addmul(v, BB(2 * k + 1, j), u, ctx);
            }

            for (i = 0; i < l; i++)
                status |= gr_sub(CC(i, j), CC(i, j), v, ctx);
        }

        for (i = 0; i < l; i++)
        {
            for (k = 0; k < h; k++)
            {
                status |= gr_add(u, BB(2 * k, 0), AA(i, 2 * k + 1), ctx);

                if (k == 0)
                    status |= gr_mul(v, AA(i, 2 * k), u, ctx);
                else
                    status |= gr_addmul(v, AA(i, 2 * k), u, ctx);
            }

            for (j = 1; j < m; j++)
                status |= gr_sub(CC(i, j), CC(i, j), v, ctx);

            for (k = 0; k < h; k++)
            {
                status |= gr_sub(u, BB(2 * k + 1, 0), AA(i, 2 * k), ctx);

                if (k == 0)
                    status |= gr_mul(t, AA(i, 2 * k + 1), u, ctx);
                else
                    status |= gr_addmul(t, AA(i, 2 * k + 1), u, ctx);
            }

            status |= gr_add(CC(i, 0), v, t, ctx);
        }

        GR_TMP_CLEAR_VEC(tmp, alloc, ctx);
    }
    else if (m == 2)
    {
        /* l x n x m  -> (l x (n - 1) x 2) + (l x 1 x m) */

        gr_mat_t A2, B2;
        gr_mat_window_init(A2, A, 0, 1, l, n, ctx);
        gr_mat_window_init(B2, B, 1, 0, n, m, ctx);

        status = gr_mat_mul_rosowski(C, A2, B2, ctx);

        for (i = 0; i < l; i++)
            for (j = 0; j < m; j++)
                status |= gr_addmul(CC(i, j), AA(i, 0), BB(0, j), ctx);
    }
    else
    {
        int even = (m % 2) == 0;
        slong m0 = (m % 2) ? 4 : 5;

        slong alloc1, alloc2, alloc3, alloc4;

        if (n != 3)
        {
            gr_mat_t A2, B2;
            gr_mat_window_init(A2, A, 0, 3, l, n, ctx);
            gr_mat_window_init(B2, B, 3, 0, n, m, ctx);
            status |= gr_mat_mul_rosowski(C, A2, B2, ctx);
        }

        /* We use 1-based indexing in the following section to mirror the
           notation of the paper. To do: consider rewriting the indices. */
#define a(iii, jjj) AA(iii - 1, jjj - 1)
#define b(iii, jjj) BB(iii - 1, jjj - 1)

        gr_ptr t, u;
        gr_ptr t0, t1, t2, t3, t4, t5, t6, t7, t8;
        gr_ptr q0, q1, q2;
        gr_ptr r0, r1, r2, r3, r4, r5, r6, r7, r8;
        gr_ptr u0, u1, u2;   /* even only */
        gr_ptr v, w;  /* (m-3)/2 elements each */

        /* todo: r0...,r8 when m0 < m only */

        alloc1 = 17;
        alloc2 = m0 < m ? 9 : 0;
        alloc3 = even ? 3 : 0;
        alloc4 = 2 * ((m - 3) / 2);
        alloc = alloc1 + alloc2 + alloc3 + alloc4;

        GR_TMP_INIT_VEC(tmp, alloc, ctx);

#define NEXT(ttt) GR_ENTRY(ttt, 1, sz)

        t = tmp;
        u = NEXT(t);
        t0 = NEXT(u);
        t1 = NEXT(t0);
        t2 = NEXT(t1);
        t3 = NEXT(t2);
        t4 = NEXT(t3);
        t5 = NEXT(t4);
        t6 = NEXT(t5);
        t7 = NEXT(t6);
        t8 = NEXT(t7);
        q0 = NEXT(t8);
        q1 = NEXT(q0);
        q2 = NEXT(q1);

        if (m0 < m)
        {
            r0 = GR_ENTRY(tmp, alloc1, sz);
            r1 = NEXT(r0);
            r2 = NEXT(r1);
            r3 = NEXT(r2);
            r4 = NEXT(r3);
            r5 = NEXT(r4);
            r6 = NEXT(r5);
            r7 = NEXT(r6);
            r8 = NEXT(r7);
        }

        if (even)
        {
            u0 = GR_ENTRY(tmp, alloc1 + alloc2, sz);
            u1 = NEXT(u0);
            u2 = NEXT(u1);
        }

        v = GR_ENTRY(tmp, alloc1 + alloc2 + alloc3, sz);
        w = GR_ENTRY(v, (m - 3) / 2, sz);

        /* TODO: spurious temporaries when inner loops not entered? */

        status |= gr_mul(t0, b(1, 2), b(2, 1), ctx);
        status |= gr_mul(t1, b(1, 3), b(3, 1), ctx);
        status |= gr_mul(t2, b(2, 3), b(3, 2), ctx);
        status |= gr_add(t3, t0, t1, ctx);
        status |= gr_add(t4, t0, t2, ctx);
        status |= gr_add(t5, t1, t2, ctx);
        status |= gr_sub(t6, b(1, 1), b(1, 2), ctx);
        status |= gr_sub(t6, t6, b(1, 3), ctx);
        status |= gr_sub(t7, b(2, 2), b(2, 1), ctx);
        status |= gr_sub(t7, t7, b(2, 3), ctx);
        status |= gr_sub(t8, b(3, 3), b(3, 1), ctx);
        status |= gr_sub(t8, t8, b(3, 2), ctx);

        if (even)
        {
            status |= gr_sub(t, b(2, 1), b(2, 4), ctx);
            status |= gr_sub(u, b(1, 4), b(1, 2), ctx);
            status |= gr_mul(u0, t, u, ctx);
            status |= gr_add(u0, u0, t0, ctx);
            status |= gr_sub(u1, b(1, 4), b(1, 2), ctx);
            status |= gr_sub(u2, b(2, 1), b(2, 4), ctx);
        }

        for (j = m0; j < m; j += 2)
        {
            status |= gr_sub(GR_ENTRY(v, (j - m0) / 2, sz), b(1, j), b(1, j + 1), ctx);
            status |= gr_sub(GR_ENTRY(w, (j - m0) / 2, sz), b(3, j), b(3, j + 1), ctx);
        }

        for (i = 1; i < l + 1; i++)
        {
            status |= gr_add(t, a(i, 1), b(2, 1), ctx);
            status |= gr_add(u, a(i, 2), b(1, 2), ctx);
            status |= gr_mul(q0, t, u, ctx);
            status |= gr_add(t, a(i, 1), b(3, 1), ctx);
            status |= gr_add(u, a(i, 3), b(1, 3), ctx);
            status |= gr_mul(q1, t, u, ctx);
            status |= gr_add(t, a(i, 2), b(3, 2), ctx);
            status |= gr_add(u, a(i, 3), b(2, 3), ctx);
            status |= gr_mul(q2, t, u, ctx);

            status |= gr_sub(t, t6, a(i, 2), ctx);
            status |= gr_sub(t, t, a(i, 3), ctx);
            status |= gr_mul(u, a(i, 1), t, ctx);
            status |= gr_add(u, u, q0, ctx);
            status |= gr_add(u, u, q1, ctx);
            if (n == 3)
            {
                status |= gr_sub(CC(i - 1, 0), u, t3, ctx);
            }
            else
            {
                status |= gr_sub(u, u, t3, ctx);
                status |= gr_add(CC(i - 1, 0), CC(i - 1, 0), u, ctx);
            }

            status |= gr_sub(t, t7, a(i, 1), ctx);
            status |= gr_sub(t, t, a(i, 3), ctx);
            status |= gr_mul(u, a(i, 2), t, ctx);
            status |= gr_add(u, u, q0, ctx);
            status |= gr_add(u, u, q2, ctx);
            if (n == 3)
            {
                status |= gr_sub(CC(i - 1, 1), u, t4, ctx);
            }
            else
            {
                status |= gr_sub(u, u, t4, ctx);
                status |= gr_add(CC(i - 1, 1), CC(i - 1, 1), u, ctx);
            }

            status |= gr_sub(t, t8, a(i, 1), ctx);
            status |= gr_sub(t, t, a(i, 2), ctx);
            status |= gr_mul(u, a(i, 3), t, ctx);
            status |= gr_add(u, u, q1, ctx);
            status |= gr_add(u, u, q2, ctx);
            if (n == 3)
            {
                status |= gr_sub(CC(i - 1, 2), u, t5, ctx);
            }
            else
            {
                status |= gr_sub(u, u, t5, ctx);
                status |= gr_add(CC(i - 1, 2), CC(i - 1, 2), u, ctx);
            }

            if (even)
            {
                status |= gr_add(t, a(i, 1), u2, ctx);
                status |= gr_sub(u, u1, a(i, 2), ctx);

                if (n == 3)
                    status |= gr_mul(CC(i - 1, 3), t, u, ctx);
                else
                    status |= gr_addmul(CC(i - 1, 3), t, u, ctx);

                status |= gr_sub(t, q0, u0, ctx);
                status |= gr_addmul(t, a(i, 3), b(3, 4), ctx);
                status |= gr_add(CC(i - 1, 3), CC(i - 1, 3), t, ctx);
            }

            if (m0 < m)
            {
                status |= gr_add(r0, q0, q1, ctx);
                status |= gr_sub(r0, r0, t3, ctx);
                status |= gr_add(r1, q1, q2, ctx);
                status |= gr_sub(r1, r1, t5, ctx);
                status |= gr_add(r2, a(i, 1), b(2, 1), ctx);
                status |= gr_add(r3, a(i, 3), b(1, 3), ctx);
                status |= gr_add(r4, a(i, 1), b(3, 1), ctx);
                status |= gr_add(r5, a(i, 2), b(1, 2), ctx);
                status |= gr_add(r6, a(i, 2), b(3, 2), ctx);
                status |= gr_add(r7, a(i, 3), b(2, 3), ctx);

                for (j = m0; j < m; j += 2)
                {
                    status |= gr_sub(t, r4, b(3, j), ctx);
                    status |= gr_sub(u, b(1, j + 1), r3, ctx);
                    status |= gr_mul(r8, t, u, ctx);

                    if (n == 3)
                    {
                        status |= gr_add(CC(i - 1, j - 1), r0, r8, ctx);
                    }
                    else
                    {
                        status |= gr_add(CC(i - 1, j - 1), CC(i - 1, j - 1), r0, ctx);
                        status |= gr_add(CC(i - 1, j - 1), CC(i - 1, j - 1), r8, ctx);
                    }

                    status |= gr_sub(t, r2, b(2, j), ctx);
                    status |= gr_sub(u, GR_ENTRY(v, (j - 4) / 2, sz), r5, ctx);
                    status |= gr_addmul(CC(i - 1, j - 1), t, u, ctx);

                    if (n == 3)
                    {
                        status |= gr_add(CC(i - 1, j), r1, r8, ctx);
                    }
                    else
                    {
                        status |= gr_add(CC(i - 1, j), CC(i - 1, j), r1, ctx);
                        status |= gr_add(CC(i - 1, j), CC(i - 1, j), r8, ctx);
                    }

                    status |= gr_add(t, r6, GR_ENTRY(w, (j - 4) / 2, sz), ctx);
                    status |= gr_sub(u, b(2, j + 1), r7, ctx);
                    status |= gr_addmul(CC(i - 1, j), t, u, ctx);
                }
            }
        }

        for (j = m0; j < m; j += 2)
        {
            status |= gr_sub(t, b(3, 1), b(3, j), ctx);
            status |= gr_sub(u, b(1, 3), b(1, j + 1), ctx);
            status |= gr_mul(q0, t, u, ctx);

            status |= gr_sub(t, b(2, j), b(2, 1), ctx);
            status |= gr_sub(u, GR_ENTRY(v, (j - 4) / 2, sz), b(1, 2), ctx);
            status |= gr_mul(q1, t, u, ctx);
            status |= gr_add(q1, q1, q0, ctx);

            status |= gr_add(t, b(3, 2), GR_ENTRY(w, (j - 4) / 2, sz), ctx);
            status |= gr_sub(u, b(2, 3), b(2, j + 1), ctx);
            status |= gr_mul(q2, t, u, ctx);
            status |= gr_add(q2, q2, q0, ctx);

            for (i = 1; i < l + 1; i++)
            {
                status |= gr_add(CC(i - 1, j - 1), CC(i - 1, j - 1), q1, ctx);
                status |= gr_add(CC(i - 1, j), CC(i - 1, j), q2, ctx);
            }
        }

        GR_TMP_CLEAR_VEC(tmp, alloc, ctx);
    }

    return status;
}

