#include "flint/profiler.h"
#include "flint/long_extras.h"

#include "gr.h"
#include "gr_vec.h"
#include "gr_mat.h"
#include "gr_poly.h"

typedef int ((*gr_test_function)(gr_ctx_t, flint_rand_t, int));

int
gr_test_binary_op_aliasing(gr_ctx_t R, int (*gr_op)(gr_ptr, gr_srcptr, gr_srcptr, gr_ctx_t), flint_rand_t state, int test_flags)
{
    int status, alias;
    gr_ptr x, y, xy1, xy2;

    GR_TMP_INIT4(x, y, xy1, xy2, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(y, state, R));

    status = GR_SUCCESS;
    status |= gr_op(xy1, x, y, R);

    alias = n_randint(state, 4);
    switch (alias)
    {
        case 0:
            status |= gr_set(xy2, x, R);
            status |= gr_op(xy1, x, y, R);
            status |= gr_op(xy2, xy2, y, R);
            break;
        case 1:
            status |= gr_set(xy2, y, R);
            status |= gr_op(xy1, x, y, R);
            status |= gr_op(xy2, x, xy2, R);
            break;
        case 2:
            status |= gr_set(y, x, R);
            status |= gr_op(xy1, x, y, R);
            status |= gr_op(xy2, x, x, R);
            break;
        default:
            status |= gr_set(y, x, R);
            status |= gr_set(xy2, x, R);
            status |= gr_op(xy1, x, y, R);
            status |= gr_op(xy2, xy2, xy2, R);
    }

    if (status == GR_SUCCESS && gr_equal(xy1, xy2, R) == T_FALSE)
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        printf("\n");
        printf("alias: %d\n", alias);
        printf("x = "); gr_println(x, R);
        printf("y = "); gr_println(y, R);
        printf("y (op) y (1) = "); gr_println(xy1, R);
        printf("x (op) y (2) = "); gr_println(xy2, R);
        printf("\n");
    }

    GR_TMP_CLEAR4(x, y, xy1, xy2, R);

    return status;
}

int
gr_test_set_ui(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr xa, xb, xc, xa_xb;
    ulong a, b, c;

    do {
        a = n_randtest(state);
        b = n_randtest(state);
        c = a + b;
    } while (c < a);

    GR_TMP_INIT4(xa, xb, xc, xa_xb, R);

    GR_MUST_SUCCEED(gr_randtest(xa, state, R));

    status = GR_SUCCESS;
    status |= gr_set_ui(xa, a, R);
    status |= gr_set_ui(xb, b, R);
    status |= gr_set_ui(xc, c, R);
    status |= gr_add(xa_xb, xa, xb, R);

    if (status == GR_SUCCESS && gr_equal(xc, xa_xb, R) == T_FALSE)
        status = GR_TEST_FAIL;
    if (status == GR_SUCCESS && a == 1 && gr_is_one(xa, R) == T_FALSE)
        status = GR_TEST_FAIL;
    if (status == GR_SUCCESS && a == 0 && gr_is_zero(xa, R) == T_FALSE)
        status = GR_TEST_FAIL;
    if (status == GR_SUCCESS && b == 1 && gr_is_one(xb, R) == T_FALSE)
        status = GR_TEST_FAIL;
    if (status == GR_SUCCESS && b == 0 && gr_is_zero(xb, R) == T_FALSE)
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        printf("\n");
        flint_printf("a = %wu\n", a);
        flint_printf("b = %wu\n", b);
        flint_printf("c = %wu\n", c);
        printf("xa = "); gr_println(xa, R);
        printf("xb = "); gr_println(xb, R);
        printf("xc = "); gr_println(xc, R);
        printf("xa + xb = "); gr_println(xa_xb, R);
        printf("\n");
    }

    GR_TMP_CLEAR4(xa, xb, xc, xa_xb, R);

    return status;
}

int
gr_test_set_si(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr xa, xb, xc, xa_xb;
    slong a, b, c;

    do {
        a = z_randtest(state);
        b = z_randtest(state);
    }
    while (z_add_checked(&c, a, b));

    GR_TMP_INIT4(xa, xb, xc, xa_xb, R);

    GR_MUST_SUCCEED(gr_randtest(xa, state, R));

    status = GR_SUCCESS;
    status |= gr_set_si(xa, a, R);
    status |= gr_set_si(xb, b, R);
    status |= gr_set_si(xc, c, R);
    status |= gr_add(xa_xb, xa, xb, R);

    if (status == GR_SUCCESS && gr_equal(xc, xa_xb, R) == T_FALSE)
        status = GR_TEST_FAIL;

    if (status == GR_SUCCESS && a == 1 && gr_is_one(xa, R) == T_FALSE)
        status = GR_TEST_FAIL;
    if (status == GR_SUCCESS && a == 0 && gr_is_zero(xa, R) == T_FALSE)
        status = GR_TEST_FAIL;
    if (status == GR_SUCCESS && b == 1 && gr_is_one(xb, R) == T_FALSE)
        status = GR_TEST_FAIL;
    if (status == GR_SUCCESS && b == 0 && gr_is_zero(xb, R) == T_FALSE)
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        printf("\n");
        flint_printf("a = %wd\n", a);
        flint_printf("b = %wd\n", b);
        flint_printf("c = %wd\n", c);
        printf("xa = "); gr_println(xa, R);
        printf("xb = "); gr_println(xb, R);
        printf("xc = "); gr_println(xc, R);
        printf("xa + xb = "); gr_println(xa_xb, R);
        printf("\n");
    }

    GR_TMP_CLEAR4(xa, xb, xc, xa_xb, R);

    return status;
}

int
gr_test_set_fmpz(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr xa, xb, xc, xa_xb;
    fmpz_t a, b, c;

    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(c);

    fmpz_randtest(a, state, 100);
    fmpz_randtest(b, state, 100);
    fmpz_add(c, a, b);

    GR_TMP_INIT4(xa, xb, xc, xa_xb, R);

    GR_MUST_SUCCEED(gr_randtest(xa, state, R));

    status = GR_SUCCESS;
    status |= gr_set_fmpz(xa, a, R);
    status |= gr_set_fmpz(xb, b, R);
    status |= gr_set_fmpz(xc, c, R);
    status |= gr_add(xa_xb, xa, xb, R);

    if (status == GR_SUCCESS && gr_equal(xc, xa_xb, R) == T_FALSE)
        status = GR_TEST_FAIL;

    if (status == GR_SUCCESS && fmpz_is_one(a) && gr_is_one(xa, R) == T_FALSE)
        status = GR_TEST_FAIL;
    if (status == GR_SUCCESS && fmpz_is_zero(a) && gr_is_zero(xa, R) == T_FALSE)
        status = GR_TEST_FAIL;
    if (status == GR_SUCCESS && fmpz_is_one(b) && gr_is_one(xb, R) == T_FALSE)
        status = GR_TEST_FAIL;
    if (status == GR_SUCCESS && fmpz_is_zero(b) && gr_is_zero(xb, R) == T_FALSE)
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        printf("\n");
        printf("a = "); fmpz_print(a); printf("\n");
        printf("b = "); fmpz_print(b); printf("\n");
        printf("c = "); fmpz_print(c); printf("\n");
        printf("xa = "); gr_println(xa, R);
        printf("xb = "); gr_println(xb, R);
        printf("xc = "); gr_println(xc, R);
        printf("xa + xb = "); gr_println(xa_xb, R);
        printf("\n");
    }

    GR_TMP_CLEAR4(xa, xb, xc, xa_xb, R);

    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(c);

    return status;
}

int
gr_test_set_fmpq(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr xa, xb, xc, xa_xb;
    fmpq_t a, b, c;

    fmpq_init(a);
    fmpq_init(b);
    fmpq_init(c);

    fmpq_randtest(a, state, 100);
    fmpq_randtest(b, state, 100);
    fmpq_add(c, a, b);

    GR_TMP_INIT4(xa, xb, xc, xa_xb, R);

    GR_MUST_SUCCEED(gr_randtest(xa, state, R));

    status = GR_SUCCESS;
    status |= gr_set_fmpq(xa, a, R);
    status |= gr_set_fmpq(xb, b, R);
    status |= gr_set_fmpq(xc, c, R);
    status |= gr_add(xa_xb, xa, xb, R);

    if (status == GR_SUCCESS && gr_equal(xc, xa_xb, R) == T_FALSE)
        status = GR_TEST_FAIL;

    if (status == GR_SUCCESS && fmpq_is_one(a) && gr_is_one(xa, R) == T_FALSE)
        status = GR_TEST_FAIL;
    if (status == GR_SUCCESS && fmpq_is_zero(a) && gr_is_zero(xa, R) == T_FALSE)
        status = GR_TEST_FAIL;
    if (status == GR_SUCCESS && fmpq_is_one(b) && gr_is_one(xb, R) == T_FALSE)
        status = GR_TEST_FAIL;
    if (status == GR_SUCCESS && fmpq_is_zero(b) && gr_is_zero(xb, R) == T_FALSE)
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        printf("\n");
        printf("a = "); fmpq_print(a); printf("\n");
        printf("b = "); fmpq_print(b); printf("\n");
        printf("c = "); fmpq_print(c); printf("\n");
        printf("xa = "); gr_println(xa, R);
        printf("xb = "); gr_println(xb, R);
        printf("xc = "); gr_println(xc, R);
        printf("xa + xb = "); gr_println(xa_xb, R);
        printf("\n");
    }

    GR_TMP_CLEAR4(xa, xb, xc, xa_xb, R);

    fmpq_clear(a);
    fmpq_clear(b);
    fmpq_clear(c);

    return status;
}

int
gr_test_binary_op_type_variants(gr_ctx_t R,
    int (*gr_op)(gr_ptr, gr_srcptr, gr_srcptr, gr_ctx_t),
    int (*gr_op_ui)(gr_ptr, gr_srcptr, ulong, gr_ctx_t),
    int (*gr_op_si)(gr_ptr, gr_srcptr, slong, gr_ctx_t),
    int (*gr_op_fmpz)(gr_ptr, gr_srcptr, const fmpz_t, gr_ctx_t),
    int (*gr_op_fmpq)(gr_ptr, gr_srcptr, const fmpq_t, gr_ctx_t),
    flint_rand_t state, int test_flags)
{
    int status, alias, which;
    gr_ptr x, y, xy1, xy2;
    ulong uy;
    slong sy;
    fmpz_t zy;
    fmpq_t qy;

    GR_TMP_INIT4(x, y, xy1, xy2, R);
    fmpz_init(zy);
    fmpq_init(qy);

    which = 0;

    uy = n_randtest(state);
    sy = (slong) n_randtest(state);
    fmpz_randtest(zy, state, 100);
    fmpq_randtest(qy, state, 100);

    for (which = 0; which < 4; which++)
    {
        GR_MUST_SUCCEED(gr_randtest(x, state, R));
        GR_MUST_SUCCEED(gr_randtest(y, state, R));
        GR_MUST_SUCCEED(gr_randtest(xy1, state, R));
        GR_MUST_SUCCEED(gr_randtest(xy2, state, R));

        status = GR_SUCCESS;
        alias = n_randint(state, 2);

        if (which == 0)
        {
            if (alias)
            {
                status |= gr_set(xy1, x, R);
                status |= gr_op_ui(xy1, xy1, uy, R);
            }
            else
            {
                status |= gr_op_ui(xy1, x, uy, R);
            }
            status |= gr_set_ui(y, uy, R);
            status |= gr_op(xy2, x, y, R);
        }
        else if (which == 1)
        {
            if (alias)
            {
                status |= gr_set(xy1, x, R);
                status |= gr_op_si(xy1, xy1, sy, R);
            }
            else
            {
                status |= gr_op_si(xy1, x, sy, R);
            }
            status |= gr_set_si(y, sy, R);
            status |= gr_op(xy2, x, y, R);
        }
        else if (which == 2)
        {
            if (alias)
            {
                status |= gr_set(xy1, x, R);
                status |= gr_op_fmpz(xy1, xy1, zy, R);
            }
            else
            {
                status |= gr_op_fmpz(xy1, x, zy, R);
            }
            status |= gr_set_fmpz(y, zy, R);
            status |= gr_op(xy2, x, y, R);
        }
        else
        {
            if (alias)
            {
                status |= gr_set(xy1, x, R);
                status |= gr_op_fmpq(xy1, xy1, qy, R);
            }
            else
            {
                status |= gr_op_fmpq(xy1, x, qy, R);
            }
            status |= gr_set_fmpq(y, qy, R);
            status |= gr_op(xy2, x, y, R);
        }

        if (status == GR_SUCCESS && gr_equal(xy1, xy2, R) == T_FALSE)
        {
            status = GR_TEST_FAIL;
            break;
        }
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        printf("\n");
        printf("which: %d\n", which);
        printf("alias: %d\n", alias);
        printf("x = "); gr_println(x, R);
        printf("uy = %lu\n", uy);
        printf("y = "); gr_println(y, R);
        printf("y (op) y (1) = "); gr_println(xy1, R);
        printf("x (op) y (2) = "); gr_println(xy2, R);
        printf("\n");
    }

    GR_TMP_CLEAR4(x, y, xy1, xy2, R);

    fmpz_clear(zy);
    fmpq_clear(qy);

    return status;
}

int
gr_test_binary_op_associative(gr_ctx_t R, int (*gr_op)(gr_ptr, gr_srcptr, gr_srcptr, gr_ctx_t), flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr x, y, z;
    gr_ptr xy, yz, xy_z, x_yz;

    GR_TMP_INIT3(x, y, z, R);
    GR_TMP_INIT4(xy, yz, xy_z, x_yz, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(y, state, R));
    GR_MUST_SUCCEED(gr_randtest(z, state, R));
    GR_MUST_SUCCEED(gr_randtest(xy, state, R));
    GR_MUST_SUCCEED(gr_randtest(yz, state, R));

    status = GR_SUCCESS;
    status |= gr_op(xy, x, y, R);
    status |= gr_op(yz, y, z, R);
    status |= gr_op(xy_z, xy, z, R);
    status |= gr_op(x_yz, x, yz, R);

    if (status == GR_SUCCESS && gr_equal(xy_z, x_yz, R) == T_FALSE)
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        printf("\n");
        printf("x = \n"); gr_println(x, R);
        printf("y = \n"); gr_println(y, R);
        printf("z = \n"); gr_println(z, R);
        printf("x (op) y = \n"); gr_println(xy, R);
        printf("y (op) z = \n"); gr_println(yz, R);
        printf("(x (op) y) (op) z = \n"); gr_println(xy_z, R);
        printf("x (op) (y (op) z) = \n"); gr_println(x_yz, R);
        printf("\n");
    }

    GR_TMP_CLEAR3(x, y, z, R);
    GR_TMP_CLEAR4(xy, yz, xy_z, x_yz, R);

    return status;
}

int
gr_test_binary_op_commutative(gr_ctx_t R, int (*gr_op)(gr_ptr, gr_srcptr, gr_srcptr, gr_ctx_t), flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr x, y, xy, yx;

    GR_TMP_INIT4(x, y, xy, yx, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(y, state, R));

    status = GR_SUCCESS;
    status |= gr_op(xy, x, y, R);
    status |= gr_op(yx, y, x, R);

    if (status == GR_SUCCESS && gr_equal(xy, yx, R) == T_FALSE)
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        printf("\n");
        printf("x = \n"); gr_println(x, R);
        printf("y = \n"); gr_println(y, R);
        printf("y (op) y = \n"); gr_println(xy, R);
        printf("y (op) x = \n"); gr_println(yx, R);
        printf("\n");
    }

    GR_TMP_CLEAR4(x, y, xy, yx, R);

    return status;
}

/*
test x op (y op2 z) = (x op y) op2 (x op z)
*/
int
gr_test_binary_op_left_distributive(gr_ctx_t R,
    int (*gr_op)(gr_ptr, gr_srcptr, gr_srcptr, gr_ctx_t),
    int (*gr_op2)(gr_ptr, gr_srcptr, gr_srcptr, gr_ctx_t),
    flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr x, y, z, yz, x_yz, xy, xz, xy_xz;

    GR_TMP_INIT4(x, y, z, yz, R);
    GR_TMP_INIT4(x_yz, xy, xz, xy_xz, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(y, state, R));
    GR_MUST_SUCCEED(gr_randtest(z, state, R));

    status = GR_SUCCESS;
    status |= gr_op2(yz, y, z, R);
    status |= gr_op(x_yz, x, yz, R);
    status |= gr_op(xy, x, y, R);
    status |= gr_op(xz, x, z, R);
    status |= gr_op2(xy_xz, xy, xz, R);

    if (status == GR_SUCCESS && gr_equal(x_yz, xy_xz, R) == T_FALSE)
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        printf("\n");
        printf("x = \n"); gr_println(x, R);
        printf("y = \n"); gr_println(y, R);
        printf("y (op2) z = \n"); gr_println(yz, R);
        printf("x (op) (y (op2) z) = \n"); gr_println(x_yz, R);
        printf("x (op) y = \n"); gr_println(xy, R);
        printf("x (op) z = \n"); gr_println(xz, R);
        printf("(x op y) (op2) (x op z) = \n"); gr_println(xy_xz, R);
        printf("\n");
    }

    GR_TMP_CLEAR4(x, y, z, yz, R);
    GR_TMP_CLEAR4(x_yz, xy, xz, xy_xz, R);

    return status;
}

/*
test (y op2 z) op x = (y op x) op2 (z op x)
*/
int
gr_test_binary_op_right_distributive(gr_ctx_t R,
    int (*gr_op)(gr_ptr, gr_srcptr, gr_srcptr, gr_ctx_t),
    int (*gr_op2)(gr_ptr, gr_srcptr, gr_srcptr, gr_ctx_t),
    flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr x, y, z, yz, yz_x, yx, zx, yx_zx;

    GR_TMP_INIT4(x, y, z, yz, R);
    GR_TMP_INIT4(yz_x, yx, zx, yx_zx, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(y, state, R));
    GR_MUST_SUCCEED(gr_randtest(z, state, R));

    status = GR_SUCCESS;
    status |= gr_op2(yz, y, z, R);
    status |= gr_op(yz_x, yz, x, R);
    status |= gr_op(yx, y, x, R);
    status |= gr_op(zx, z, x, R);
    status |= gr_op2(yx_zx, yx, zx, R);

    if (status == GR_SUCCESS && gr_equal(yz_x, yx_zx, R) == T_FALSE)
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        printf("\n");
        printf("x = \n"); gr_println(x, R);
        printf("y = \n"); gr_println(y, R);
        printf("y (op2) z = \n"); gr_println(yz, R);
        printf("(y (op2) z) op x = \n"); gr_println(yz_x, R);
        printf("y (op) x = \n"); gr_println(yz, R);
        printf("z (op) x = \n"); gr_println(zx, R);
        printf("(y op x) (op2) (z op x) = \n"); gr_println(yx_zx, R);
        printf("\n");
    }

    GR_TMP_CLEAR4(x, y, z, yz, R);
    GR_TMP_CLEAR4(yz_x, yx, zx, yx_zx, R);

    return status;
}

int
gr_test_init_clear(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr a, b, c, d, e;

    status = GR_SUCCESS;

    GR_TMP_INIT(a, R);
    status |= gr_randtest(a, state, R);
    GR_TMP_CLEAR(a, R);

    GR_TMP_INIT2(a, b, R);
    status |= gr_randtest(a, state, R);
    status |= gr_randtest(b, state, R);
    GR_TMP_CLEAR2(a, b, R);

    GR_TMP_INIT3(a, b, c, R);
    status |= gr_randtest(a, state, R);
    status |= gr_randtest(b, state, R);
    status |= gr_randtest(c, state, R);
    GR_TMP_CLEAR3(a, b, c, R);

    GR_TMP_INIT4(a, b, c, d, R);
    status |= gr_randtest(a, state, R);
    status |= gr_randtest(b, state, R);
    status |= gr_randtest(c, state, R);
    status |= gr_randtest(d, state, R);
    GR_TMP_CLEAR4(a, b, c, d, R);

    GR_TMP_INIT5(a, b, c, d, e, R);
    status |= gr_randtest(a, state, R);
    status |= gr_randtest(b, state, R);
    status |= gr_randtest(c, state, R);
    status |= gr_randtest(d, state, R);
    status |= gr_randtest(e, state, R);
    GR_TMP_CLEAR5(a, b, c, d, e, R);

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    return status;
}

int
gr_test_swap(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr a, b, c, d;
    truth_t equal0, equal1, equal2, equal3, equal4;

    status = GR_SUCCESS;

    GR_TMP_INIT4(a, b, c, d, R);

    status |= gr_randtest(a, state, R);
    status |= gr_randtest(b, state, R);
    status |= gr_set(c, a, R);
    status |= gr_set(d, b, R);
    gr_swap(a, a, R);

    equal0 = gr_equal(a, c, R);

    gr_swap(a, b, R);
    equal1 = gr_equal(b, c, R);
    equal2 = gr_equal(a, d, R);

    gr_swap(a, b, R);
    equal3 = gr_equal(a, c, R);
    equal4 = gr_equal(b, d, R);

    if (status == GR_SUCCESS &&
        (equal0 == T_FALSE || equal1 == T_FALSE ||
         equal2 == T_FALSE || equal3 == T_FALSE ||
         equal4 == T_FALSE))
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    GR_TMP_CLEAR4(a, b, c, d, R);

    return status;
}

int
gr_test_zero_one(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr a;
    truth_t equal;

    status = GR_SUCCESS;

    GR_TMP_INIT(a, R);

    status |= gr_randtest(a, state, R);
    status |= gr_zero(a, R);
    equal = gr_is_zero(a, R);
    if (status == GR_SUCCESS && equal == T_FALSE)
        status = GR_TEST_FAIL;

    status |= gr_randtest(a, state, R);
    status |= gr_one(a, R);
    equal = gr_is_one(a, R);
    if (status == GR_SUCCESS && equal == T_FALSE)
        status = GR_TEST_FAIL;

    status |= gr_randtest(a, state, R);
    status |= gr_one(a, R);
    status |= gr_neg(a, a, R);
    equal = gr_is_neg_one(a, R);
    if (status == GR_SUCCESS && equal == T_FALSE)
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    GR_TMP_CLEAR(a, R);

    return status;
}

int
gr_test_add_associative(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    return gr_test_binary_op_associative(R, gr_add, state, test_flags);
}

int
gr_test_neg(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr x, y, xy;

    GR_TMP_INIT3(x, y, xy, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(y, state, R));
    GR_MUST_SUCCEED(gr_randtest(xy, state, R));

    status = GR_SUCCESS;

    /* check x + (-x) = 0 */
    status |= gr_neg(y, x, R);
    status |= gr_add(xy, x, y, R);

    if (status == GR_SUCCESS && gr_is_zero(xy, R) == T_FALSE)
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        printf("\n");
        printf("x = \n"); gr_println(x, R);
        printf("y = \n"); gr_println(y, R);
        printf("x + y = \n"); gr_println(xy, R);
        printf("\n");
    }

    /* check -(-x) = x, plus aliasing */
    status |= gr_neg(y, y, R);

    if (status == GR_SUCCESS && gr_equal(x, y, R) == T_FALSE)
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        printf("\n");
        printf("x = \n"); gr_println(x, R);
        printf("y = \n"); gr_println(y, R);
        printf("\n");
    }

    GR_TMP_CLEAR3(x, y, xy, R);

    return status;
}

int
gr_test_add_commutative(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    return gr_test_binary_op_commutative(R, gr_add, state, test_flags);
}

int
gr_test_add_aliasing(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    return gr_test_binary_op_aliasing(R, gr_add, state, test_flags);
}

int
gr_test_add_type_variants(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    return gr_test_binary_op_type_variants(R,
        gr_add, gr_add_ui, gr_add_si, gr_add_fmpz, gr_add_fmpq,
            state, test_flags);
}

int
gr_test_sub_equal_neg_add(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr x, y, neg_y, x_sub_y, x_neg_y;

    GR_TMP_INIT5(x, y, neg_y, x_sub_y, x_neg_y, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(y, state, R));
    GR_MUST_SUCCEED(gr_randtest(neg_y, state, R));
    GR_MUST_SUCCEED(gr_randtest(x_sub_y, state, R));
    GR_MUST_SUCCEED(gr_randtest(x_neg_y, state, R));

    status = GR_SUCCESS;
    status |= gr_sub(x_sub_y, x, y, R);
    status |= gr_neg(neg_y, y, R);
    status |= gr_add(x_neg_y, x, neg_y, R);

    if (status == GR_SUCCESS && gr_equal(x_sub_y, x_neg_y, R) == T_FALSE)
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        printf("\n");
        printf("x = \n"); gr_println(x, R);
        printf("y = \n"); gr_println(y, R);
        printf("-y = \n"); gr_println(neg_y, R);
        printf("x - y = \n"); gr_println(x_sub_y, R);
        printf("x + (-y) = \n"); gr_println(x_neg_y, R);
        printf("\n");
    }

    GR_TMP_CLEAR5(x, y, neg_y, x_sub_y, x_neg_y, R);

    return status;
}

int
gr_test_sub_aliasing(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    return gr_test_binary_op_aliasing(R, gr_sub, state, test_flags);
}

int
gr_test_sub_type_variants(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    return gr_test_binary_op_type_variants(R,
        gr_sub, gr_sub_ui, gr_sub_si, gr_sub_fmpz, gr_sub_fmpq,
            state, test_flags);
}

int
gr_test_mul_associative(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    return gr_test_binary_op_associative(R, gr_mul, state, test_flags);
}

int
gr_test_mul_commutative(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    return gr_test_binary_op_commutative(R, gr_mul, state, test_flags);
}

int
gr_test_mul_left_distributive(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    return gr_test_binary_op_left_distributive(R, gr_mul, gr_add, state, test_flags);
}

int
gr_test_mul_right_distributive(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    return gr_test_binary_op_right_distributive(R, gr_mul, gr_add, state, test_flags);
}


int
gr_test_mul_aliasing(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    return gr_test_binary_op_aliasing(R, gr_mul, state, test_flags);
}

int
gr_test_mul_type_variants(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    return gr_test_binary_op_type_variants(R,
        gr_mul, gr_mul_ui, gr_mul_si, gr_mul_fmpz, gr_mul_fmpq,
            state, test_flags);
}

int
gr_test_div_type_variants(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    return gr_test_binary_op_type_variants(R,
        gr_div, gr_div_ui, gr_div_si, gr_div_fmpz, gr_div_fmpq,
            state, test_flags);
}

int
gr_test_inv_involution(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr x, x_inv, x_inv_inv;

    GR_TMP_INIT3(x, x_inv, x_inv_inv, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(x_inv, state, R));
    GR_MUST_SUCCEED(gr_randtest(x_inv_inv, state, R));

    status = GR_SUCCESS;
    status |= gr_inv(x_inv, x, R);
    status |= gr_inv(x_inv_inv, x_inv, R);

    if (status == GR_SUCCESS && gr_equal(x, x_inv_inv, R) == T_FALSE)
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        printf("\n");
        printf("x = \n"); gr_println(x, R);
        printf("x ^ -1 = \n"); gr_println(x_inv, R);
        printf("(x ^ -1) ^ -1 = \n"); gr_println(x_inv_inv, R);
        printf("\n");
    }

    GR_TMP_CLEAR3(x, x_inv, x_inv_inv, R);

    return status;
}

int
gr_test_inv_multiplication(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    truth_t equal1, equal2;
    gr_ptr x, x_inv, x_inv_x, x_x_inv;

    GR_TMP_INIT4(x, x_inv, x_inv_x, x_x_inv, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(x_inv, state, R));
    GR_MUST_SUCCEED(gr_randtest(x_inv_x, state, R));
    GR_MUST_SUCCEED(gr_randtest(x_x_inv, state, R));

    /* todo: split status */
    status = GR_SUCCESS;
    status |= gr_inv(x_inv, x, R);
    status |= gr_mul(x_inv_x, x_inv, x, R);
    status |= gr_mul(x_x_inv, x, x_inv, R);
    equal1 = gr_is_one(x_inv_x, R);
    equal2 = gr_is_one(x_x_inv, R);

    if (status == GR_SUCCESS && (equal1 == T_FALSE || equal2 == T_FALSE))
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        printf("\n");
        printf("x = \n"); gr_println(x, R);
        printf("x ^ -1 = \n"); gr_println(x_inv, R);
        printf("(x ^ -1) * x = \n"); gr_println(x_inv_x, R);
        printf("x * (x ^ -1) = \n"); gr_println(x_x_inv, R);
        printf("\n");
    }

    GR_TMP_CLEAR4(x, x_inv, x_inv_x, x_x_inv, R);

    return status;
}

int
gr_test_div_right_distributive(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    return gr_test_binary_op_right_distributive(R, gr_div, gr_add, state, test_flags);
}

int
gr_test_div_then_mul(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr x, y, xy, xyy;

    GR_TMP_INIT4(x, y, xy, xyy, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(y, state, R));
    GR_MUST_SUCCEED(gr_randtest(xy, state, R));
    GR_MUST_SUCCEED(gr_randtest(xyy, state, R));

    status = GR_SUCCESS;
    status |= gr_div(xy, x, y, R);
    status |= gr_mul(xyy, xy, y, R);

    if (status == GR_SUCCESS && gr_equal(x, xyy, R) == T_FALSE)
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        printf("\n");
        printf("x = \n"); gr_println(x, R);
        printf("y = \n"); gr_println(y, R);
        printf("x / y = \n"); gr_println(xy, R);
        printf("(x / y) * y = \n"); gr_println(xyy, R);
        printf("\n");
    }

    GR_TMP_CLEAR4(x, y, xy, xyy, R);

    return status;
}

int
gr_test_mul_then_div(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    gr_ptr x, y, xy, xyy;

    GR_TMP_INIT4(x, y, xy, xyy, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(y, state, R));
    GR_MUST_SUCCEED(gr_randtest(xy, state, R));
    GR_MUST_SUCCEED(gr_randtest(xyy, state, R));

    status = GR_SUCCESS;
    status |= gr_mul(xy, x, y, R);
    status |= gr_div(xyy, xy, y, R);

    if (status == GR_SUCCESS && gr_equal(x, xyy, R) == T_FALSE)
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        printf("\n");
        printf("x = \n"); gr_println(x, R);
        printf("y = \n"); gr_println(y, R);
        printf("x * y = \n"); gr_println(xy, R);
        printf("(x * y) / y = \n"); gr_println(xyy, R);
        printf("\n");
    }

    GR_TMP_CLEAR4(x, y, xy, xyy, R);

    return status;
}

int
gr_test_pow_ui_exponent_addition(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    ulong a, b;
    gr_ptr x, xa, xb, xab, xaxb;

    GR_TMP_INIT5(x, xa, xb, xab, xaxb, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(xa, state, R));
    GR_MUST_SUCCEED(gr_randtest(xb, state, R));
    GR_MUST_SUCCEED(gr_randtest(xab, state, R));
    GR_MUST_SUCCEED(gr_randtest(xaxb, state, R));

    if (gr_ctx_is_finite(R) == T_TRUE)
    {
        do {
            a = n_randtest(state);
            b = n_randtest(state);
        } while (a + b < a);
    }
    else
    {
        a = n_randtest(state) % 20;
        b = n_randtest(state) % 20;
    }

    status = GR_SUCCESS;

    status |= gr_pow_ui(xa, x, a, R);
    status |= gr_pow_ui(xb, x, b, R);
    status |= gr_pow_ui(xab, x, a + b, R);
    status |= gr_mul(xaxb, xa, xb, R);

    if (status == GR_SUCCESS && gr_equal(xab, xaxb, R) == T_FALSE)
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        printf("\n");
        printf("x = \n"); gr_println(x, R);
        flint_printf("a = %wu\n", a);
        flint_printf("b = %wu\n", b);
        printf("x ^ a = \n"); gr_println(xa, R);
        printf("x ^ b = \n"); gr_println(xb, R);
        printf("x ^ (a + b) = \n"); gr_println(xab, R);
        printf("(x ^ a) * (x ^ b) = \n"); gr_println(xaxb, R);
        printf("\n");
    }

    GR_TMP_CLEAR5(x, xa, xb, xab, xaxb, R);

    return status;
}

int
gr_test_pow_ui_base_scalar_multiplication(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    ulong a;
    slong y;
    gr_ptr x, xa, ya, xya, xaya;

    GR_TMP_INIT3(x, xa, ya, R);
    GR_TMP_INIT2(xya, xaya, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(xa, state, R));
    GR_MUST_SUCCEED(gr_randtest(ya, state, R));

    y = n_randtest(state);

    if (gr_ctx_is_finite(R) == T_TRUE)
        a = n_randtest(state);
    else
        a = n_randtest(state) % 20;

    status = GR_SUCCESS;
    status |= gr_pow_ui(xa, x, a, R);
    status |= gr_set_si(ya, y, R);
    status |= gr_pow_ui(ya, ya, a, R);
    status |= gr_set_si(xya, y, R);
    status |= gr_mul(xya, x, xya, R);   /* todo mul_si */
    status |= gr_pow_ui(xya, xya, a, R);
    status |= gr_mul(xaya, xa, ya, R);

    if (status == GR_SUCCESS && gr_equal(xya, xaya, R) == T_FALSE)
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        printf("\n");
        printf("x = \n"); gr_println(x, R);
        flint_printf("y = %wd\n", y);
        flint_printf("a = %wu\n", a);
        printf("x ^ a = \n"); gr_println(xa, R);
        printf("y ^ a = \n"); gr_println(ya, R);
        printf("(x * y) ^ a = \n"); gr_println(xya, R);
        printf("(x ^ a) * (y ^ a) = \n"); gr_println(xaya, R);
        printf("\n");
    }

    GR_TMP_CLEAR3(x, xa, ya, R);
    GR_TMP_CLEAR2(xya, xaya, R);

    return status;
}

int
gr_test_pow_ui_base_multiplication(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    ulong a;
    gr_ptr x, y, xa, ya, xya, xaya;

    GR_TMP_INIT4(x, y, xa, ya, R);
    GR_TMP_INIT2(xya, xaya, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(y, state, R));
    GR_MUST_SUCCEED(gr_randtest(xa, state, R));
    GR_MUST_SUCCEED(gr_randtest(ya, state, R));

    if (gr_ctx_is_finite(R) == T_TRUE)
        a = n_randtest(state);
    else
        a = n_randtest(state) % 20;

    status = GR_SUCCESS;
    status |= gr_pow_ui(xa, x, a, R);
    status |= gr_pow_ui(ya, y, a, R);
    status |= gr_mul(xya, x, y, R);
    status |= gr_pow_ui(xya, xya, a, R);
    status |= gr_mul(xaya, xa, ya, R);

    if (status == GR_SUCCESS && gr_equal(xya, xaya, R) == T_FALSE)
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        printf("\n");
        printf("x = \n"); gr_println(x, R);
        printf("y = \n"); gr_println(y, R);
        flint_printf("a = %wu\n", a);
        printf("x ^ a = \n"); gr_println(xa, R);
        printf("y ^ a = \n"); gr_println(ya, R);
        printf("(x * y) ^ a = \n"); gr_println(xya, R);
        printf("(x ^ a) * (y ^ a) = \n"); gr_println(xaya, R);
        printf("\n");
    }

    GR_TMP_CLEAR4(x, y, xa, ya, R);
    GR_TMP_CLEAR2(xya, xaya, R);

    return status;
}

int
gr_test_pow_ui_aliasing(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    ulong a;
    gr_ptr x, xa1, xa2;

    GR_TMP_INIT3(x, xa1, xa2, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(xa1, state, R));

    if (gr_ctx_is_finite(R) == T_TRUE)
        a = n_randtest(state);
    else
        a = n_randtest(state) % 20;

    status = GR_SUCCESS;
    status |= gr_pow_ui(xa1, x, a, R);
    status |= gr_set(xa2, x, R);
    status |= gr_pow_ui(xa2, xa2, a, R);

    if (status == GR_SUCCESS && gr_equal(xa1, xa2, R) == T_FALSE)
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        printf("\n");
        printf("x = \n"); gr_println(x, R);
        flint_printf("a = %wu\n", a);
        printf("x ^ a (1) = \n"); gr_println(xa1, R);
        printf("x ^ a (2) = \n"); gr_println(xa2, R);
        printf("\n");
    }

    GR_TMP_CLEAR3(x, xa1, xa2, R);

    return status;
}

int
gr_test_pow_fmpz_exponent_addition(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    fmpz_t a, b, ab;
    gr_ptr x, xa, xb, xab, xaxb;

    GR_TMP_INIT5(x, xa, xb, xab, xaxb, R);

    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(ab);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(xa, state, R));
    GR_MUST_SUCCEED(gr_randtest(xb, state, R));
    GR_MUST_SUCCEED(gr_randtest(xab, state, R));
    GR_MUST_SUCCEED(gr_randtest(xaxb, state, R));

    if (gr_ctx_is_finite(R) == T_TRUE)
    {
        fmpz_randtest(a, state, 100);
        fmpz_randtest(b, state, 100);
    }
    else
    {
        fmpz_randtest(a, state, 4);
        fmpz_randtest(b, state, 4);
    }

    fmpz_add(ab, a, b);

    status = GR_SUCCESS;

    status |= gr_pow_fmpz(xa, x, a, R);
    status |= gr_pow_fmpz(xb, x, b, R);
    status |= gr_pow_fmpz(xab, x, ab, R);
    status |= gr_mul(xaxb, xa, xb, R);

    if (status == GR_SUCCESS && gr_equal(xab, xaxb, R) == T_FALSE)
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        printf("\n");
        printf("x = \n"); gr_println(x, R);
        printf("a = "); fmpz_print(a); printf("\n");
        printf("b = "); fmpz_print(b); printf("\n");
        printf("x ^ a = \n"); gr_println(xa, R);
        printf("x ^ b = \n"); gr_println(xb, R);
        printf("x ^ (a + b) = \n"); gr_println(xab, R);
        printf("(x ^ a) * (x ^ b) = \n"); gr_println(xaxb, R);
        printf("\n");
    }

    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(ab);

    GR_TMP_CLEAR5(x, xa, xb, xab, xaxb, R);

    return status;
}

int
gr_test_sqrt(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status = GR_SUCCESS;
    gr_ptr x, y, y2;
    int perfect;

    GR_TMP_INIT3(x, y, y2, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(y, state, R));

    perfect = n_randint(state, 2);

    if (perfect)
        status |= gr_sqr(x, x, R);

    if (n_randint(state, 2))
    {
        status |= gr_set(y, x, R);
        status |= gr_sqrt(y, y, R);
    }
    else
    {
        status |= gr_sqrt(y, x, R);
    }

    status |= gr_sqr(y2, y, R);

    if (status == GR_SUCCESS && gr_equal(y2, x, R) == T_FALSE)
    {
        status = GR_TEST_FAIL;
    }

    if (status == GR_DOMAIN && perfect)
    {
        status = GR_TEST_FAIL;
    }

    if (status == GR_DOMAIN && perfect && gr_is_square(x, R) == T_FALSE)
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        printf("\n");
        gr_ctx_println(R);
        printf("x = \n"); gr_println(x, R);
        printf("y = \n"); gr_println(y, R);
        printf("y ^ 2 = \n"); gr_println(y2, R);
        printf("\n");
    }

    GR_TMP_CLEAR3(x, y, y2, R);

    return status;
}

int
gr_test_ordered_ring_cmp(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status = GR_SUCCESS;
    gr_ptr x, y, z, xz, yz, zero, xy;
    int cmp1, cmp2, cmp3;

    GR_TMP_INIT5(x, y, z, xz, yz, R);
    GR_TMP_INIT2(zero, xy, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(y, state, R));
    GR_MUST_SUCCEED(gr_randtest(z, state, R));

    /* cmp(x, y) = -cmp(y, x) */
    status |= gr_cmp(&cmp1, x, y, R);
    status |= gr_cmp(&cmp2, y, x, R);

    if (status == GR_SUCCESS && cmp1 != -cmp2)
    {
        status = GR_TEST_FAIL;
    }

    /* x <= y  -->  x + z <= y + z */
    status |= gr_add(xz, x, z, R);
    status |= gr_add(yz, y, z, R);
    status |= gr_cmp(&cmp1, x, y, R);
    status |= gr_cmp(&cmp2, xz, yz, R);

    if (status == GR_SUCCESS && cmp1 != cmp2)
    {
        status = GR_TEST_FAIL;
    }

    /* 0 <= x and 0 <= y --> 0 <= xy */
    status |= gr_cmp(&cmp1, zero, x, R);
    status |= gr_cmp(&cmp2, zero, y, R);
    status |= gr_mul(xy, x, y, R);
    status |= gr_cmp(&cmp3, zero, xy, R);

    if (status == GR_SUCCESS && cmp1 <= 0 && cmp2 <= 0 && cmp3 > 0)
    {
        status = GR_TEST_FAIL;
    }

    if (status & GR_DOMAIN && !(status & GR_UNABLE))
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        printf("\n");
        printf("R = "); gr_ctx_println(R);
        printf("x = \n"); gr_println(x, R);
        printf("y = \n"); gr_println(y, R);
        printf("z = \n"); gr_println(z, R);
        printf("x + z = \n"); gr_println(xz, R);
        printf("y + z = \n"); gr_println(yz, R);
        printf("xy = \n"); gr_println(xy, R);
        printf("cmp = %d, %d, %d\n", cmp1, cmp2, cmp3);
        printf("\n");
    }

    GR_TMP_CLEAR5(x, y, z, xz, yz, R);
    GR_TMP_CLEAR2(zero, xy, R);

    return status;
}

int
gr_test_ordered_ring_cmpabs(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status = GR_SUCCESS;
    gr_ptr x, y, ax, ay;
    int cmp1, cmp2;

    GR_TMP_INIT4(x, y, ax, ay, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(y, state, R));

    status |= gr_abs(ax, x, R);
    status |= gr_abs(ay, y, R);
    status |= gr_cmpabs(&cmp1, x, y, R);
    status |= gr_cmp(&cmp2, ax, ay, R);

    if (status == GR_SUCCESS && cmp1 != cmp2)
    {
        status = GR_TEST_FAIL;
    }

    if (status & GR_DOMAIN && !(status & GR_UNABLE))
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        printf("\n");
        printf("R = "); gr_ctx_println(R);
        printf("x = \n"); gr_println(x, R);
        printf("y = \n"); gr_println(y, R);
        printf("ax = \n"); gr_println(ax, R);
        printf("ay = \n"); gr_println(ay, R);
        printf("cmp = %d, %d\n", cmp1, cmp2);
        printf("\n");
    }

    GR_TMP_CLEAR4(x, y, ax, ay, R);

    return status;
}

int
gr_test_vec_add(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status, aliasing;
    slong i, len;
    gr_ptr x, y, xy1, xy2;

    len = n_randint(state, 10);

    GR_TMP_INIT_VEC(x, len, R);
    GR_TMP_INIT_VEC(y, len, R);
    GR_TMP_INIT_VEC(xy1, len, R);
    GR_TMP_INIT_VEC(xy2, len, R);

    GR_MUST_SUCCEED(_gr_vec_randtest(x, state, len, R));

    GR_MUST_SUCCEED(_gr_vec_randtest(y, state, len, R));
    GR_MUST_SUCCEED(_gr_vec_randtest(xy1, state, len, R));
    GR_MUST_SUCCEED(_gr_vec_randtest(xy2, state, len, R));

    status = GR_SUCCESS;

    aliasing = n_randint(state, 4);
    switch (aliasing)
    {
        case 0:
            status |= _gr_vec_set(xy1, x, len, R);
            status |= _gr_vec_add(xy1, xy1, y, len, R);
            break;
        case 1:
            status |= _gr_vec_set(xy1, y, len, R);
            status |= _gr_vec_add(xy1, x, xy1, len, R);
            break;
        case 2:
            status |= _gr_vec_set(y, x, len, R);
            status |= _gr_vec_set(xy1, x, len, R);
            status |= _gr_vec_add(xy1, xy1, xy1, len, R);
            break;
        default:
            status |= _gr_vec_add(xy1, x, y, len, R);
    }

    for (i = 0; i < len; i++)
        status |= gr_add(GR_ENTRY(xy2, i, R->sizeof_elem),
                         GR_ENTRY(x, i, R->sizeof_elem),
                         GR_ENTRY(y, i, R->sizeof_elem), R);

    if (status == GR_SUCCESS && _gr_vec_equal(xy1, xy2, len, R) == T_FALSE)
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        /* todo: vec print */
        printf("\n");
        printf("aliasing: %d\n", aliasing);
        printf("\n");
    }

    GR_TMP_CLEAR_VEC(x, len, R);
    GR_TMP_CLEAR_VEC(y, len, R);
    GR_TMP_CLEAR_VEC(xy1, len, R);
    GR_TMP_CLEAR_VEC(xy2, len, R);

    return status;
}

/* (AB)C = A(BC) */
int
gr_test_mat_mul_classical_associative(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status;
    gr_mat_t A, B, C, AB, BC, AB_C, A_BC;
    slong m, n, p, q;

    if (gr_ctx_is_finite(R) == T_TRUE)
    {
        m = n_randint(state, 5);
        n = n_randint(state, 5);
        p = n_randint(state, 5);
        q = n_randint(state, 5);
    }
    else
    {
        m = n_randint(state, 3);
        n = n_randint(state, 3);
        p = n_randint(state, 3);
        q = n_randint(state, 3);
    }

    gr_mat_init(A, m, n, R);
    gr_mat_init(B, n, p, R);
    gr_mat_init(C, p, q, R);
    gr_mat_init(AB, m, p, R);
    gr_mat_init(BC, n, q, R);
    gr_mat_init(AB_C, m, q, R);
    gr_mat_init(A_BC, m, q, R);

    GR_MUST_SUCCEED(gr_mat_randtest(A, state, R));
    GR_MUST_SUCCEED(gr_mat_randtest(B, state, R));
    GR_MUST_SUCCEED(gr_mat_randtest(C, state, R));
    GR_MUST_SUCCEED(gr_mat_randtest(AB, state, R));
    GR_MUST_SUCCEED(gr_mat_randtest(BC, state, R));
    GR_MUST_SUCCEED(gr_mat_randtest(AB_C, state, R));
    GR_MUST_SUCCEED(gr_mat_randtest(A_BC, state, R));

    status = GR_SUCCESS;
    status |= gr_mat_mul_classical(AB, A, B, R);
    status |= gr_mat_mul_classical(BC, B, C, R);
    status |= gr_mat_mul_classical(AB_C, AB, C, R);
    status |= gr_mat_mul_classical(A_BC, A, BC, R);

    if (status == GR_SUCCESS && gr_mat_equal(AB_C, A_BC, R) == T_FALSE)
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        /* todo: vec print */
        printf("\n");
        printf("A = \n"); gr_mat_print(A, R); printf("\n");
        printf("B = \n"); gr_mat_print(B, R); printf("\n");
        printf("C = \n"); gr_mat_print(C, R); printf("\n");
        printf("AB = \n"); gr_mat_print(AB, R); printf("\n");
        printf("BC = \n"); gr_mat_print(BC, R); printf("\n");
        printf("AB * C = \n"); gr_mat_print(AB_C, R); printf("\n");
        printf("A * BC = \n"); gr_mat_print(A_BC, R); printf("\n");
        printf("\n");
    }

    gr_mat_clear(A, R);
    gr_mat_clear(B, R);
    gr_mat_clear(C, R);
    gr_mat_clear(AB, R);
    gr_mat_clear(BC, R);
    gr_mat_clear(A_BC, R);
    gr_mat_clear(AB_C, R);

    return status;
}

int
gr_test_integral_domain(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status = GR_SUCCESS;
    gr_ptr x, y, z;

    GR_TMP_INIT3(x, y, z, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));
    GR_MUST_SUCCEED(gr_randtest(y, state, R));
    status |= gr_mul(z, x, y, R);

    if (status == GR_SUCCESS && gr_is_zero(x, R) == T_FALSE && gr_is_zero(y, R) == T_FALSE && gr_is_zero(z, R) == T_TRUE)
    {
        status = GR_TEST_FAIL;
    }

    if ((test_flags & GR_TEST_ALWAYS_ABLE) && (status & GR_UNABLE))
        status = GR_TEST_FAIL;

    if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
    {
        printf("\n");
        printf("x = \n"); gr_println(x, R);
        printf("y = \n"); gr_println(y, R);
        printf("z = \n"); gr_println(z, R);
        printf("\n");
    }

    if (gr_ctx_is_commutative_ring(R) == T_FALSE)
    {
        printf("integral domain is not a commutative ring\n");
        printf("\n");
        status = GR_TEST_FAIL;
    }

    GR_TMP_CLEAR3(x, y, z, R);

    return status;
}

int
gr_test_field(gr_ctx_t R, flint_rand_t state, int test_flags)
{
    int status = GR_SUCCESS;
    gr_ptr x, y, z;

    GR_TMP_INIT3(x, y, z, R);

    GR_MUST_SUCCEED(gr_randtest(x, state, R));

    if (gr_is_zero(x, R) == T_FALSE)
    {
        if (gr_is_invertible(x, R) == T_FALSE)
        {
            status = GR_TEST_FAIL;
        }

        if (gr_inv(y, x, R) == GR_DOMAIN)
        {
            status = GR_TEST_FAIL;
        }

        if (gr_div(z, y, x, R) == GR_DOMAIN)
        {
            status = GR_TEST_FAIL;
        }

        if ((test_flags & GR_TEST_VERBOSE) || status == GR_TEST_FAIL)
        {
            printf("\n");
            printf("x = \n"); gr_println(x, R);
            printf("y = \n"); gr_println(y, R);
            printf("z = \n"); gr_println(z, R);
            printf("\n");
        }
    }

    if (gr_ctx_is_commutative_ring(R) == T_FALSE)
    {
        printf("field is not a commutative ring\n");
        printf("\n");
        status = GR_TEST_FAIL;
    }

    if (gr_ctx_is_integral_domain(R) == T_FALSE)
    {
        printf("field is not an integral domain\n");
        printf("\n");
        status = GR_TEST_FAIL;
    }

    GR_TMP_CLEAR3(x, y, z, R);

    return status;
}

void
gr_test_iter(gr_ctx_t R, flint_rand_t state, const char * descr, gr_test_function func, slong iters, int test_flags)
{
    slong iter, count_success, count_unable, count_domain;
    int status;
    timeit_t timer;

    count_success = 0;
    count_unable = 0;
    count_domain = 0;

    if (test_flags & GR_TEST_VERBOSE)
    {
        printf("%s ... ", descr);
        fflush(stdout);
    }

    timeit_start(timer);

    for (iter = 0; iter < iters; iter++)
    {
        /* flint_printf("iter %ld\n", iter); */
        status = func(R, state, test_flags & ~GR_TEST_VERBOSE);

        if (status == GR_SUCCESS)
            count_success++;

        if (status & GR_UNABLE)
            count_unable++;

        if (status & GR_DOMAIN)
            count_domain++;

        if (status & GR_TEST_FAIL)
        {
            flint_printf("\nFAIL\n");
            abort();
        }
    }

    timeit_stop(timer);

    if (test_flags & GR_TEST_VERBOSE)
    {
        flint_printf("PASS   (%wd successful, %wd domain, %wd unable, 0 wrong, %.3g cpu, %.3g wall)\n",
            count_success, count_domain, count_unable, timer->cpu*0.001, timer->wall*0.001);
    }
}

void
gr_test_ring(gr_ctx_t R, slong iters, int test_flags)
{
    timeit_t timer;
    flint_rand_t state;

    /* test_flags |= GR_TEST_VERBOSE; */

    if (test_flags & GR_TEST_VERBOSE)
    {
        timeit_start(timer);

        flint_printf("===============================================================================\n");
        flint_printf("Testing "); gr_ctx_println(R);
        flint_printf("-------------------------------------------------------------------------------\n");
    }

    flint_randinit(state);

    gr_test_iter(R, state, "init/clear", gr_test_init_clear, iters, test_flags);
    gr_test_iter(R, state, "swap", gr_test_swap, iters, test_flags);
    gr_test_iter(R, state, "zero_one", gr_test_zero_one, iters, test_flags);
    gr_test_iter(R, state, "neg", gr_test_neg, iters, test_flags);

    gr_test_iter(R, state, "set_ui", gr_test_set_ui, iters, test_flags);
    gr_test_iter(R, state, "set_si", gr_test_set_si, iters, test_flags);
    gr_test_iter(R, state, "set_fmpz", gr_test_set_fmpz, iters, test_flags);
    gr_test_iter(R, state, "set_fmpq", gr_test_set_fmpq, iters, test_flags);

    gr_test_iter(R, state, "add: associative", gr_test_add_associative, iters, test_flags);
    gr_test_iter(R, state, "add: commutative", gr_test_add_commutative, iters, test_flags);
    gr_test_iter(R, state, "add: aliasing", gr_test_add_aliasing, iters, test_flags);
    gr_test_iter(R, state, "sub: equal neg add", gr_test_sub_equal_neg_add, iters, test_flags);
    gr_test_iter(R, state, "sub: aliasing", gr_test_sub_aliasing, iters, test_flags);

    gr_test_iter(R, state, "add: ui/si/fmpz/fmpq", gr_test_add_type_variants, iters, test_flags);
    gr_test_iter(R, state, "sub: ui/si/fmpz/fmpq", gr_test_sub_type_variants, iters, test_flags);
    gr_test_iter(R, state, "mul: ui/si/fmpz/fmpq", gr_test_mul_type_variants, iters, test_flags);
    gr_test_iter(R, state, "div: ui/si/fmpz/fmpq", gr_test_div_type_variants, iters, test_flags);

    gr_test_iter(R, state, "mul: associative", gr_test_mul_associative, iters, test_flags);
    if (gr_ctx_is_commutative_ring(R) == T_TRUE)
        gr_test_iter(R, state, "mul: commutative", gr_test_mul_commutative, iters, test_flags);
    gr_test_iter(R, state, "mul: aliasing", gr_test_mul_aliasing, iters, test_flags);
    gr_test_iter(R, state, "mul: left distributive", gr_test_mul_left_distributive, iters, test_flags);
    gr_test_iter(R, state, "mul: right distributive", gr_test_mul_right_distributive, iters, test_flags);

    if (gr_ctx_is_integral_domain(R) == T_TRUE)
        gr_test_iter(R, state, "integral_domain", gr_test_integral_domain, iters, test_flags);

    if (gr_ctx_is_field(R) == T_TRUE)
        gr_test_iter(R, state, "field", gr_test_integral_domain, iters, test_flags);

    gr_test_iter(R, state, "div: distributive", gr_test_div_right_distributive, iters, test_flags);
    gr_test_iter(R, state, "div: div then mul", gr_test_div_then_mul, iters, test_flags);
    gr_test_iter(R, state, "div: mul then div", gr_test_mul_then_div, iters, test_flags);

    gr_test_iter(R, state, "inv: multiplication", gr_test_inv_multiplication, iters, test_flags);
    gr_test_iter(R, state, "inv: involution", gr_test_inv_involution, iters, test_flags);

    gr_test_iter(R, state, "pow_ui: exponent addition", gr_test_pow_ui_exponent_addition, iters, test_flags);
    gr_test_iter(R, state, "pow_ui: base scalar multiplication", gr_test_pow_ui_base_scalar_multiplication, iters, test_flags);

    if (gr_ctx_is_commutative_ring(R) == T_TRUE)
        gr_test_iter(R, state, "pow_ui: base multiplication", gr_test_pow_ui_base_multiplication, iters, test_flags);

    gr_test_iter(R, state, "pow_ui: aliasing", gr_test_pow_ui_exponent_addition, iters, test_flags);
    gr_test_iter(R, state, "pow_fmpz: exponent addition", gr_test_pow_fmpz_exponent_addition, iters, test_flags);

    gr_test_iter(R, state, "sqrt", gr_test_sqrt, iters, test_flags & (~GR_TEST_ALWAYS_ABLE));

    if (gr_ctx_is_ordered_ring(R) == T_TRUE)
    {
        gr_test_iter(R, state, "ordered_ring_cmp", gr_test_ordered_ring_cmp, iters, test_flags);
        gr_test_iter(R, state, "ordered_ring_cmpabs", gr_test_ordered_ring_cmpabs, iters, test_flags);
    }

    gr_test_iter(R, state, "vec_add", gr_test_vec_add, iters, test_flags);

    gr_test_iter(R, state, "mat_mul_classical: associative", gr_test_mat_mul_classical_associative, iters, test_flags);

    flint_randclear(state);

    if (test_flags & GR_TEST_VERBOSE)
    {
        timeit_stop(timer);

        flint_printf("-------------------------------------------------------------------------------\n");
        flint_printf("Tests finished in %.3g cpu, %.3g wall\n", timer->cpu*0.001, timer->wall*0.001);
        flint_printf("===============================================================================\n\n");
    }
}
