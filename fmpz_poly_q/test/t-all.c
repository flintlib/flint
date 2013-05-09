/******************************************************************************

    Copyright (C) 2009, 2010 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "fmpz_poly_q.h"

void test_set(char * in, char * out)
{
    int ans;
    fmpz_poly_q_t rop, op;
    char * res;
    
    fmpz_poly_q_init(op);
    fmpz_poly_q_set_str(op, in);
    
    fmpz_poly_q_init(rop);
    fmpz_poly_q_set(rop, op);
    
    res = fmpz_poly_q_get_str(rop);
    ans = !strcmp(out, res);
    
    if (!ans)
    {
        printf("test_set: failed\n");
        printf("    Expected \"%s\", got \"%s\"\n", out, res);
        abort();
    }
    
    fmpz_poly_q_clear(op);
    fmpz_poly_q_clear(rop);
    flint_free(res);
}

void test_set_si(len_t x, char * out)
{
    int ans;
    fmpz_poly_q_t rop;
    char * res;
    
    fmpz_poly_q_init(rop);
    fmpz_poly_q_set_si(rop, x);
    
    res = fmpz_poly_q_get_str(rop);
    ans = !strcmp(out, res);
    
    if (!ans)
    {
        printf("test_set_si: failed\n");
        printf("    Expected \"%s\", got \"%s\"\n", out, res);
        abort();
    }
    
    fmpz_poly_q_clear(rop);
    flint_free(res);
}

void test_swap(char * in1, char * in2, char * out1, char * out2)
{
    int ans;
    fmpz_poly_q_t op1, op2;
    char * res1;
    char * res2;
    
    fmpz_poly_q_init(op1);
    fmpz_poly_q_set_str(op1, in1);
    fmpz_poly_q_init(op2);
    fmpz_poly_q_set_str(op2, in2);
    
    fmpz_poly_q_swap(op1, op2);
    
    res1 = fmpz_poly_q_get_str(op1);
    res2 = fmpz_poly_q_get_str(op2);
    
    ans = !strcmp(out1, res1) && !strcmp(out2, res2);
    
    if (!ans)
    {
        printf("test_swap: failed\n");
        printf("    Expected \"%s\" \"%s\", got \"%s\" \"%s\"\n", out1, out2, res1, res2);
        abort();
    }
    
    fmpz_poly_q_clear(op1);
    fmpz_poly_q_clear(op2);
    flint_free(res1);
    flint_free(res2);
}

void test_zero(char * in, char * out)
{
    int ans;
    fmpz_poly_q_t op;
    char * res;
    
    fmpz_poly_q_init(op);
    fmpz_poly_q_set_str(op, in);
    
    fmpz_poly_q_zero(op);
    
    res = fmpz_poly_q_get_str(op);
    ans = !strcmp(res, out);
    
    if (!ans)
    {
        printf("test_zero: failed\n");
        printf("    Expected \"%s\", got \"%s\"\n", out, res);
        abort();
    }
    
    fmpz_poly_q_clear(op);
    flint_free(res);
}

void test_neg(char * in, char * out)
{
    int ans;
    fmpz_poly_q_t rop, op;
    char * res;
    
    fmpz_poly_q_init(op);
    fmpz_poly_q_set_str(op, in);
    
    fmpz_poly_q_init(rop);
    fmpz_poly_q_neg(rop, op);
    
    res = fmpz_poly_q_get_str(rop);
    ans = !strcmp(res, out);
    
    if (!ans)
    {
        printf("test_neg: failed\n");
        printf("    Expected \"%s\", got \"%s\"\n", out, res);
        abort();
    }
    
    fmpz_poly_q_clear(op);
    fmpz_poly_q_clear(rop);
    flint_free(res);
}

void test_inv(char * in, char * out)
{
    int ans;
    fmpz_poly_q_t rop, op;
    char * res;
     
    fmpz_poly_q_init(op);
    fmpz_poly_q_set_str(op, in);
    fmpz_poly_q_init(rop);
    fmpz_poly_q_inv(rop, op);
    
    res = fmpz_poly_q_get_str(rop);
    ans = !strcmp(res, out);
    
    if (!ans)
    {
        printf("test_inv: failed\n");
        printf("    Expected \"%s\", got \"%s\"\n", out, res);
        abort();
    }
    
    fmpz_poly_q_clear(op);
    fmpz_poly_q_clear(rop);
    flint_free(res);
}

void test_inv_inplace(char * in, char * out)
{
    int ans;
    fmpz_poly_q_t rop;
    char * res;
     
    fmpz_poly_q_init(rop);
    fmpz_poly_q_set_str(rop, in);
    fmpz_poly_q_inv(rop, rop);
    
    res = fmpz_poly_q_get_str(rop);
    ans = !strcmp(res, out);
    
    if (!ans)
    {
        printf("test_inv_inplace: failed\n");
        printf("    Expected \"%s\", got \"%s\"\n", out, res);
        abort();
    }
    
    fmpz_poly_q_clear(rop);
    flint_free(res);
}

void test_is_zero(char * in, int out)
{
    int ans;
    fmpz_poly_q_t op;
    int res;
    
    fmpz_poly_q_init(op);
    fmpz_poly_q_set_str(op, in);
    
    res = fmpz_poly_q_is_zero(op);
    ans = (res == out);
    
    if (!ans)
    {
        printf("test_equal: failed\n");
        printf("    Expected \"%d\", got \"%d\"\n", out, res);
        abort();
    }
    
    fmpz_poly_q_clear(op);
}

void test_is_one(char * in, int out)
{
    int ans;
    fmpz_poly_q_t op;
    int res;
    
    fmpz_poly_q_init(op);
    fmpz_poly_q_set_str(op, in);
    
    res = fmpz_poly_q_is_one(op);
    ans = (res == out);
    
    if (!ans)
    {
        printf("test_equal: failed\n");
        printf("    Expected \"%d\", got \"%d\"\n", out, res);
        abort();
    }
    
    fmpz_poly_q_clear(op);
}

void test_equal(char * in1, char * in2, int out)
{
    int ans;
    fmpz_poly_q_t op1, op2;
    int res;
    
    fmpz_poly_q_init(op1);
    fmpz_poly_q_set_str(op1, in1);
    
    fmpz_poly_q_init(op2);
    fmpz_poly_q_set_str(op2, in2);
    
    res = fmpz_poly_q_equal(op1, op2);
    ans = (res == out);
    
    if (!ans)
    {
        printf("test_equal: failed\n");
        printf("    Expected \"%d\", got \"%d\"\n", out, res);
        abort();
    }
    
    fmpz_poly_q_clear(op1);
    fmpz_poly_q_clear(op2);
}

void test_add(char * in1, char * in2, char * out)
{
    int ans;
    fmpz_poly_q_t rop, op1, op2;
    char * res;
    
    fmpz_poly_q_init(op1);
    fmpz_poly_q_set_str(op1, in1);
    
    fmpz_poly_q_init(op2);
    fmpz_poly_q_set_str(op2, in2);
    
    fmpz_poly_q_init(rop);
    fmpz_poly_q_add(rop, op1, op2);
    
    res = fmpz_poly_q_get_str(rop);
    ans = !strcmp(out, res);

    if (!ans)
    {
        printf("test_add: failed\n");
        printf("    Expected \"%s\", got \"%s\"\n", out, res);
        abort();
    }
    
    fmpz_poly_q_clear(op1);
    fmpz_poly_q_clear(op2);
    fmpz_poly_q_clear(rop);
    flint_free(res);
}

/* Runs in1 = in1 + in2 */
void test_add_in_place1(char * in1, char * in2, char * out)
{
    int ans;
    fmpz_poly_q_t op1, op2;
    char * res;
    
    fmpz_poly_q_init(op1);
    fmpz_poly_q_set_str(op1, in1);
    
    fmpz_poly_q_init(op2);
    fmpz_poly_q_set_str(op2, in2);
    
    fmpz_poly_q_add(op1, op1, op2);
    
    res = fmpz_poly_q_get_str(op1);
    ans = !strcmp(out, res);

    if (!ans)
    {
        printf("test_add_in_place1: failed\n");
        printf("    Expected \"%s\", got \"%s\"\n", out, res);
        abort();
    }
    
    fmpz_poly_q_clear(op1);
    fmpz_poly_q_clear(op2);
    flint_free(res);
}

/* Runs in2 = in1 + in2 */
void test_add_in_place2(char * in1, char * in2, char * out)
{
    int ans;
    fmpz_poly_q_t op1, op2;
    char * res;
    
    fmpz_poly_q_init(op1);
    fmpz_poly_q_set_str(op1, in1);
    
    fmpz_poly_q_init(op2);
    fmpz_poly_q_set_str(op2, in2);
    
    fmpz_poly_q_add(op2, op1, op2);
    
    res = fmpz_poly_q_get_str(op2);
    ans = !strcmp(out, res);

    if (!ans)
    {
        printf("test_add_in_place2: failed\n");
        printf("    Expected \"%s\", got \"%s\"\n", out, res);
        abort();
    }
    
    fmpz_poly_q_clear(op1);
    fmpz_poly_q_clear(op2);
    flint_free(res);
}

/* Runs out = in + in */
void test_add_in_place3(char * in, char * out)
{
    int ans;
    fmpz_poly_q_t rop, op;
    char * res;
    
    fmpz_poly_q_init(op);
    fmpz_poly_q_set_str(op, in);
    
    fmpz_poly_q_init(rop);
    fmpz_poly_q_add(rop, op, op);
    
    res = fmpz_poly_q_get_str(rop);
    ans = !strcmp(out, res);

    if (!ans)
    {
        printf("test_add_in_place3: failed\n");
        printf("    Expected \"%s\", got \"%s\"\n", out, res);
        abort();
    }
    
    fmpz_poly_q_clear(rop);
    fmpz_poly_q_clear(op);
    flint_free(res);
}

void test_sub(char * in1, char * in2, char * out)
{
    int ans;
    fmpz_poly_q_t rop, op1, op2;
    char * res;
    
    fmpz_poly_q_init(op1);
    fmpz_poly_q_set_str(op1, in1);
    
    fmpz_poly_q_init(op2);
    fmpz_poly_q_set_str(op2, in2);
    
    fmpz_poly_q_init(rop);
    fmpz_poly_q_sub(rop, op1, op2);
    
    res = fmpz_poly_q_get_str(rop);
    ans = !strcmp(out, res);
    
    if (!ans)
    {
        printf("test_sub: failed\n");
        printf("    Expected \"%s\", got \"%s\"\n", out, res);
        abort();
    }
    
    fmpz_poly_q_clear(op1);
    fmpz_poly_q_clear(op2);
    fmpz_poly_q_clear(rop);
    flint_free(res);
}

/* in1 = in1 + in2 */
void test_sub_in_place1(char * in1, char * in2, char * out)
{
    int ans;
    fmpz_poly_q_t op1, op2;
    char * res;
    
    fmpz_poly_q_init(op1);
    fmpz_poly_q_set_str(op1, in1);
    
    fmpz_poly_q_init(op2);
    fmpz_poly_q_set_str(op2, in2);
    
    fmpz_poly_q_sub(op1, op1, op2);
    
    res = fmpz_poly_q_get_str(op1);
    ans = !strcmp(out, res);
    
    if (!ans)
    {
        printf("test_sub_in_place1: failed\n");
        printf("    Expected \"%s\", got \"%s\"\n", out, res);
        abort();
    }
    
    fmpz_poly_q_clear(op1);
    fmpz_poly_q_clear(op2);
    flint_free(res);
}

/* in2 = in1 + in2 */
void test_sub_in_place2(char * in1, char * in2, char * out)
{
    int ans;
    fmpz_poly_q_t op1, op2;
    char * res;
    
    fmpz_poly_q_init(op1);
    fmpz_poly_q_set_str(op1, in1);
    
    fmpz_poly_q_init(op2);
    fmpz_poly_q_set_str(op2, in2);
    
    fmpz_poly_q_sub(op2, op1, op2);
    
    res = fmpz_poly_q_get_str(op2);
    ans = !strcmp(out, res);
    
    if (!ans)
    {
        printf("test_sub_in_place2: failed\n");
        printf("    Expected \"%s\", got \"%s\"\n", out, res);
        abort();
    }
    
    fmpz_poly_q_clear(op1);
    fmpz_poly_q_clear(op2);
    flint_free(res);
}

/* Runs out = in - in */
void test_sub_in_place3(char * in, char * out)
{
    int ans;
    fmpz_poly_q_t rop, op;
    char * res;
    
    fmpz_poly_q_init(op);
    fmpz_poly_q_set_str(op, in);
    
    fmpz_poly_q_init(rop);
    fmpz_poly_q_sub(rop, op, op);
    
    res = fmpz_poly_q_get_str(rop);
    ans = !strcmp(out, res);

    if (!ans)
    {
        printf("test_sub_in_place3: failed\n");
        printf("    Expected \"%s\", got \"%s\"\n", out, res);
        abort();
    }
    
    fmpz_poly_q_clear(rop);
    fmpz_poly_q_clear(op);
    flint_free(res);
}

void test_scalar_mul_si(char * in, len_t x, char * out)
{
    int ans;
    fmpz_poly_q_t rop, op;
    char * res;
    
    fmpz_poly_q_init(op);
    fmpz_poly_q_set_str(op, in);
    
    fmpz_poly_q_init(rop);
    fmpz_poly_q_scalar_mul_si(rop, op, x);
    
    res = fmpz_poly_q_get_str(rop);
    ans = !strcmp(out, res);
    
    if (!ans)
    {
        printf("test_scalar_mul_si: failed\n");
        printf("    Expected \"%s\", got \"%s\"\n", out, res);
        abort();
    }
    
    fmpz_poly_q_clear(op);
    fmpz_poly_q_clear(rop);
    flint_free(res);
}

void test_scalar_mul_mpz(char * in, mpz_t x, char * out)
{
    int ans;
    fmpz_poly_q_t rop, op;
    char * res;
    
    fmpz_poly_q_init(op);
    fmpz_poly_q_set_str(op, in);
    
    fmpz_poly_q_init(rop);
    fmpz_poly_q_scalar_mul_mpz(rop, op, x);
    
    res = fmpz_poly_q_get_str(rop);
    ans = !strcmp(out, res);
    
    if (!ans)
    {
        printf("test_scalar_mul_mpz: failed\n");
        printf("    Expected \"%s\", got \"%s\"\n", out, res);
        abort();
    }
    
    fmpz_poly_q_clear(op);
    fmpz_poly_q_clear(rop);
    flint_free(res);
}

void test_scalar_mul_mpq(char * in, mpq_t x, char * out)
{
    int ans;
    fmpz_poly_q_t rop, op;
    char * res;
    
    fmpz_poly_q_init(op);
    fmpz_poly_q_set_str(op, in);
    
    fmpz_poly_q_init(rop);
    fmpz_poly_q_scalar_mul_mpq(rop, op, x);
    
    res = fmpz_poly_q_get_str(rop);
    ans = !strcmp(out, res);
    
    if (!ans)
    {
        printf("test_scalar_mul_mpq: failed\n");
        printf("    Expected \"%s\", got \"%s\"\n", out, res);
        abort();
    }
    
    fmpz_poly_q_clear(op);
    fmpz_poly_q_clear(rop);
    flint_free(res);
}

void test_scalar_div_si(char * in, len_t x, char * out)
{
    int ans;
    fmpz_poly_q_t rop, op;
    char * res;
    
    fmpz_poly_q_init(op);
    fmpz_poly_q_set_str(op, in);
    
    fmpz_poly_q_init(rop);
    fmpz_poly_q_scalar_div_si(rop, op, x);
    
    res = fmpz_poly_q_get_str(rop);
    ans = !strcmp(out, res);
    
    if (!ans)
    {
        printf("test_scalar_div_si: failed\n");
        printf("    Expected \"%s\", got \"%s\"\n", out, res);
        abort();
    }
    
    fmpz_poly_q_clear(op);
    fmpz_poly_q_clear(rop);
    flint_free(res);
}

void test_scalar_div_mpz(char * in, mpz_t x, char * out)
{
    int ans;
    fmpz_poly_q_t rop, op;
    char * res;
    
    fmpz_poly_q_init(op);
    fmpz_poly_q_set_str(op, in);
    
    fmpz_poly_q_init(rop);
    fmpz_poly_q_scalar_div_mpz(rop, op, x);
    
    res = fmpz_poly_q_get_str(rop);
    ans = !strcmp(out, res);
    
    if (!ans)
    {
        printf("test_scalar_div_mpz: failed\n");
        printf("    Expected \"%s\", got \"%s\"\n", out, res);
        abort();
    }
    
    fmpz_poly_q_clear(op);
    fmpz_poly_q_clear(rop);
    flint_free(res);
}

void test_scalar_div_mpq(char * in, mpq_t x, char * out)
{
    int ans;
    fmpz_poly_q_t rop, op;
    char * res;
    
    fmpz_poly_q_init(op);
    fmpz_poly_q_set_str(op, in);
    
    fmpz_poly_q_init(rop);
    fmpz_poly_q_scalar_div_mpq(rop, op, x);
    
    res = fmpz_poly_q_get_str(rop);
    ans = !strcmp(out, res);
    
    if (!ans)
    {
        printf("test_scalar_div_mpq: failed\n");
        printf("    Expected \"%s\", got \"%s\"\n", out, res);
        abort();
    }
    
    fmpz_poly_q_clear(op);
    fmpz_poly_q_clear(rop);
    flint_free(res);
}

void test_mul(char * in1, char * in2, char * out)
{
    int ans;
    fmpz_poly_q_t rop, op1, op2;
    char * res;
    
    fmpz_poly_q_init(op1);
    fmpz_poly_q_set_str(op1, in1);
    
    fmpz_poly_q_init(op2);
    fmpz_poly_q_set_str(op2, in2);
    
    fmpz_poly_q_init(rop);
    fmpz_poly_q_mul(rop, op1, op2);
    
    res = fmpz_poly_q_get_str(rop);
    ans = !strcmp(out, res);
    
    if (!ans)
    {
        printf("test_mul: failed\n");
        printf("    Expected \"%s\", got \"%s\"\n", out, res);
        abort();
    }
    
    fmpz_poly_q_clear(op1);
    fmpz_poly_q_clear(op2);
    fmpz_poly_q_clear(rop);
    flint_free(res);
}

/* in1 = in1 * in2 */
void test_mul_in_place1(char * in1, char * in2, char * out)
{
    int ans;
    fmpz_poly_q_t op1, op2;
    char * res;
    
    fmpz_poly_q_init(op1);
    fmpz_poly_q_set_str(op1, in1);
    
    fmpz_poly_q_init(op2);
    fmpz_poly_q_set_str(op2, in2);
    
    fmpz_poly_q_mul(op1, op1, op2);
    
    res = fmpz_poly_q_get_str(op1);
    ans = !strcmp(out, res);
    
    if (!ans)
    {
        printf("test_mul_in_place1: failed\n");
        printf("    Expected \"%s\", got \"%s\"\n", out, res);
        abort();
    }
    
    fmpz_poly_q_clear(op1);
    fmpz_poly_q_clear(op2);
    flint_free(res);
}

/* in2 = in1 * in2 */
void test_mul_in_place2(char * in1, char * in2, char * out)
{
    int ans;
    fmpz_poly_q_t op1, op2;
    char * res;
    
    fmpz_poly_q_init(op1);
    fmpz_poly_q_set_str(op1, in1);
    
    fmpz_poly_q_init(op2);
    fmpz_poly_q_set_str(op2, in2);
    
    fmpz_poly_q_mul(op2, op1, op2);
    
    res = fmpz_poly_q_get_str(op2);
    ans = !strcmp(out, res);
    
    if (!ans)
    {
        printf("test_mul_in_place2: failed\n");
        printf("    Expected \"%s\", got \"%s\"\n", out, res);
        abort();
    }
    
    fmpz_poly_q_clear(op1);
    fmpz_poly_q_clear(op2);
    flint_free(res);
}

/* Runs out = in * in */
void test_mul_in_place3(char * in, char * out)
{
    int ans;
    fmpz_poly_q_t rop, op;
    char * res;
    
    fmpz_poly_q_init(op);
    fmpz_poly_q_set_str(op, in);
    
    fmpz_poly_q_init(rop);
    fmpz_poly_q_mul(rop, op, op);
    
    res = fmpz_poly_q_get_str(rop);
    ans = !strcmp(out, res);

    if (!ans)
    {
        printf("test_mul_in_place3: failed\n");
        printf("    Expected \"%s\", got \"%s\"\n", out, res);
        abort();
    }
    
    fmpz_poly_q_clear(rop);
    fmpz_poly_q_clear(op);
    flint_free(res);
}

void test_div(char * in1, char * in2, char * out)
{
    int ans;
    fmpz_poly_q_t rop, op1, op2;
    char * res;
    
    fmpz_poly_q_init(op1);
    fmpz_poly_q_set_str(op1, in1);
    
    fmpz_poly_q_init(op2);
    fmpz_poly_q_set_str(op2, in2);
    
    fmpz_poly_q_init(rop);
    fmpz_poly_q_div(rop, op1, op2);
    
    res = fmpz_poly_q_get_str(rop);
    ans = !strcmp(out, res);
    
    if (!ans)
    {
        printf("test_div: failed\n");
        printf("    Expected \"%s\", got \"%s\"\n", out, res);
        abort();
    }
    
    fmpz_poly_q_clear(op1);
    fmpz_poly_q_clear(op2);
    fmpz_poly_q_clear(rop);
    flint_free(res);
}

/* in1 = in1 / in2 */
void test_div_in_place1(char * in1, char * in2, char * out)
{
    int ans;
    fmpz_poly_q_t op1, op2;
    char * res;
    
    fmpz_poly_q_init(op1);
    fmpz_poly_q_set_str(op1, in1);
    
    fmpz_poly_q_init(op2);
    fmpz_poly_q_set_str(op2, in2);
    
    fmpz_poly_q_div(op1, op1, op2);
    
    res = fmpz_poly_q_get_str(op1);
    ans = !strcmp(out, res);
    
    if (!ans)
    {
        printf("test_div_in_place1: failed\n");
        printf("    Expected \"%s\", got \"%s\"\n", out, res);
        abort();
    }
    
    fmpz_poly_q_clear(op1);
    fmpz_poly_q_clear(op2);
    flint_free(res);
}

/* in2 = in1 / in2 */
void test_div_in_place2(char * in1, char * in2, char * out)
{
    int ans;
    fmpz_poly_q_t op1, op2;
    char * res;
    
    fmpz_poly_q_init(op1);
    fmpz_poly_q_set_str(op1, in1);
    
    fmpz_poly_q_init(op2);
    fmpz_poly_q_set_str(op2, in2);
    
    fmpz_poly_q_div(op2, op1, op2);
    
    res = fmpz_poly_q_get_str(op2);
    ans = !strcmp(out, res);
    
    if (!ans)
    {
        printf("test_div_in_place2: failed\n");
        printf("    Expected \"%s\", got \"%s\"\n", out, res);
        abort();
    }
    
    fmpz_poly_q_clear(op1);
    fmpz_poly_q_clear(op2);
    flint_free(res);
}

/* Runs out = in / in */
void test_div_in_place3(char * in, char * out)
{
    int ans;
    fmpz_poly_q_t rop, op;
    char * res;
    
    fmpz_poly_q_init(op);
    fmpz_poly_q_set_str(op, in);
    
    fmpz_poly_q_init(rop);
    fmpz_poly_q_div(rop, op, op);
    
    res = fmpz_poly_q_get_str(rop);
    ans = !strcmp(out, res);

    if (!ans)
    {
        printf("test_div_in_place3: failed\n");
        printf("    Expected \"%s\", got \"%s\"\n", out, res);
        abort();
    }
    
    fmpz_poly_q_clear(rop);
    fmpz_poly_q_clear(op);
    flint_free(res);
}

void test_pow(char * in, unsigned long x, char * out)
{
    int ans;
    fmpz_poly_q_t rop, op;
    char * res;
    
    fmpz_poly_q_init(op);
    fmpz_poly_q_set_str(op, in);
    
    fmpz_poly_q_init(rop);
    fmpz_poly_q_pow(rop, op, x);
    
    res = fmpz_poly_q_get_str(rop);
    ans = !strcmp(out, res);
    
    if (!ans)
    {
        printf("test_pow: failed\n");
        printf("    Expected \"%s\", got \"%s\"\n", out, res);
        abort();
    }
    
    fmpz_poly_q_clear(op);
    fmpz_poly_q_clear(rop);
    flint_free(res);
}

void test_derivative(char * in, char * out)
{
    int ans;
    fmpz_poly_q_t rop, op;
    char * res;
    
    fmpz_poly_q_init(op);
    fmpz_poly_q_set_str(op, in);
    
    fmpz_poly_q_init(rop);
    fmpz_poly_q_derivative(rop, op);
    
    res = fmpz_poly_q_get_str(rop);
    ans = !strcmp(out, res);
    
    if (!ans)
    {
        printf("test_derivative: failed\n");
        printf("    Expected \"%s\", got \"%s\"\n", out, res);
        abort();
    }
    
    fmpz_poly_q_clear(op);
    fmpz_poly_q_clear(rop);
    flint_free(res);
}

void test_evaluate(char * in, int numa, int numb, char * out)
{
    int ans, pole;
    fmpz_poly_q_t op;
    mpq_t rop, a;
    char *res = NULL;
    
    fmpz_poly_q_init(op);
    fmpz_poly_q_set_str(op, in);
    
    mpq_init(a);
    mpq_set_si(a, numa, numb);
    mpq_init(rop);
    pole = fmpz_poly_q_evaluate(rop, op, a);
    
    if (pole && strcmp(out, "P"))
    {
        printf("test_evaluate: failed\n");
        printf("    Expected \"%s\", got a pole\n", out);
        abort();
    }
    if (!pole && !strcmp(out, "P"))
    {
        res = mpq_get_str(NULL, 10, rop);
        printf("test_evaluate: failed\n");
        printf("    Expected a pole, got \"%s\"\n", res);
        abort();
    }
    if (!pole)
    {
        res = mpq_get_str(NULL, 10, rop);
        ans = (strcmp(out, res) == 0);
        
        if (!ans)
        {
            printf("test_evaluate: failed\n");
            printf("    Expected \"%s\", got \"%s\"\n", out, res);
            abort();
        }
    }
    
    fmpz_poly_q_clear(op);
    mpq_clear(rop);
    mpq_clear(a);
    flint_free(res);
}

void test_get_str_pretty(char * in, char * out)
{
    int ans;
    fmpz_poly_q_t rop;
    char * res;
    
    fmpz_poly_q_init(rop);
    fmpz_poly_q_set_str(rop, in);
    res = fmpz_poly_q_get_str_pretty(rop, "t");
    ans = !strcmp(out, res);
    
    if (!ans)
    {
        printf("test_get_str_pretty: failed\n");
        printf("    Expected \"%s\", got \"%s\"\n", out, res);
        abort();
    }
    
    fmpz_poly_q_clear(rop);
    flint_free(res);
}

void test_addmul(char * in1, char * in2, char * in3, char * out)
{
    int ans;
    fmpz_poly_q_t rop, op1, op2;
    char * res;
    
    fmpz_poly_q_init(rop);
    fmpz_poly_q_set_str(rop, in1);
    fmpz_poly_q_init(op1);
    fmpz_poly_q_set_str(op1, in2);
    fmpz_poly_q_init(op2);
    fmpz_poly_q_set_str(op2, in3);
    
    fmpz_poly_q_addmul(rop, op1, op2);
    
    res = fmpz_poly_q_get_str(rop);
    ans = !strcmp(out, res);
    
    if (!ans)
    {
        printf("test_addmul: failed\n");
        printf("    Expected \"%s\", got \"%s\"\n", out, res);
        abort();
    }
    
    fmpz_poly_q_clear(rop);
    fmpz_poly_q_clear(op1);
    fmpz_poly_q_clear(op2);
    flint_free(res);
}

void test_submul(char * in1, char * in2, char * in3, char * out)
{
    int ans;
    fmpz_poly_q_t rop, op1, op2;
    char * res;
    
    fmpz_poly_q_init(rop);
    fmpz_poly_q_set_str(rop, in1);
    fmpz_poly_q_init(op1);
    fmpz_poly_q_set_str(op1, in2);
    fmpz_poly_q_init(op2);
    fmpz_poly_q_set_str(op2, in3);
    
    fmpz_poly_q_submul(rop, op1, op2);
    
    res = fmpz_poly_q_get_str(rop);
    ans = !strcmp(out, res);
    
    if (!ans)
    {
        printf("test_submul: failed\n");
        printf("    Expected \"%s\", got \"%s\"\n", out, res);
        abort();
    }
    
    fmpz_poly_q_clear(rop);
    fmpz_poly_q_clear(op1);
    fmpz_poly_q_clear(op2);
    flint_free(res);
}

int main(int argc, char *argv[])
{
    int ans;
    char *str, *strout;
    
    fmpz_poly_t zpoly;
    fmpz_poly_q_t qpoly1;
    
    mpz_t mpzzero, mpzone, mpztwo;
    mpq_t mpqzero, mpqone, mpqtwo, mpqtwoinv;
    
    printf("all... ");
    fflush(stdout);
    
    /* Accessing numerator and denominator ***********************************/
    
    fmpz_poly_q_init(qpoly1);
    fmpz_poly_q_set_str(qpoly1, "2  -1 1/2  0 1");
    str = "2  -1 1";
    strout = fmpz_poly_get_str(fmpz_poly_q_numref(qpoly1));
    ans = !strcmp(str, strout);
    if (!ans)
    {
        printf("test_numref: failed\n");
        printf("    Expected \"%s\", got \"%s\"\n", str, strout);
        printf("    qpoly1 = \""), fmpz_poly_q_print(qpoly1), printf("\"\n");
        abort();
    }
    fmpz_poly_q_clear(qpoly1);
    flint_free(strout);
    
    fmpz_poly_q_init(qpoly1);
    fmpz_poly_q_set_str(qpoly1, "2  -1 1/2  0 1");
    str = "2  0 1";
    strout = fmpz_poly_get_str(fmpz_poly_q_denref(qpoly1));
    ans = !strcmp(str, strout);
    if (!ans)
    {
        printf("test_denref: failed\n");
        printf("    Expected \"%s\", got \"%s\"\n", str, strout);
        abort();
    }
    fmpz_poly_q_clear(qpoly1);
    flint_free(strout);
    
    fmpz_poly_q_init(qpoly1);
    fmpz_poly_init(zpoly);
    fmpz_poly_q_set_str(qpoly1, "2  -1 1/2  0 1");
    fmpz_poly_set(zpoly, fmpz_poly_q_numref(qpoly1));
    str = "2  -1 1";
    strout = fmpz_poly_get_str(zpoly);
    ans = !strcmp(str, strout);
    if (!ans)
    {
        printf("test_get_num: failed\n");
        printf("    Expected \"%s\", got \"%s\"\n", str, strout);
        abort();
    }
    fmpz_poly_q_clear(qpoly1);
    fmpz_poly_clear(zpoly);
    flint_free(strout);
    
    fmpz_poly_q_init(qpoly1);
    fmpz_poly_init(zpoly);
    fmpz_poly_q_set_str(qpoly1, "2  -1 1/2  0 1");
    fmpz_poly_set(zpoly, fmpz_poly_q_denref(qpoly1));
    
    str = "2  0 1";
    strout = fmpz_poly_get_str(zpoly);
    ans = !strcmp(str, strout);
    if (!ans)
    {
        printf("test_get_den: failed\n");
        printf("    Expected \"%s\", got \"%s\"\n", str, strout);
        abort();
    }
    fmpz_poly_q_clear(qpoly1);
    fmpz_poly_clear(zpoly);
    flint_free(strout);
    
    fmpz_poly_q_init(qpoly1);
    fmpz_poly_init(zpoly);
    fmpz_poly_q_set_str(qpoly1, "1  1/1  1");
    fmpz_poly_set_str(zpoly, "2  0 1");
    fmpz_poly_set(fmpz_poly_q_numref(qpoly1), zpoly);
    str = "2  0 1";
    strout = fmpz_poly_get_str(fmpz_poly_q_numref(qpoly1));
    ans = !strcmp(str, strout);
    if (!ans)
    {
        printf("test_set_num: failed\n");
        printf("    Expected \"%s\", got \"%s\"\n", str, strout);
        abort();
    }
    fmpz_poly_q_clear(qpoly1);
    fmpz_poly_clear(zpoly);
    flint_free(strout);
    
    fmpz_poly_q_init(qpoly1);
    fmpz_poly_init(zpoly);
    fmpz_poly_q_set_str(qpoly1, "1  1/1  1");
    fmpz_poly_set_str(zpoly, "2  0 1");
    fmpz_poly_set(fmpz_poly_q_denref(qpoly1), zpoly);
    str = "2  0 1";
    strout = fmpz_poly_get_str(fmpz_poly_q_denref(qpoly1));
    ans = !strcmp(str, strout);
    if (!ans)
    {
        printf("test_set_den: failed\n");
        printf("    Expected \"%s\", got \"%s\"\n", str, strout);
        abort();
    }
    fmpz_poly_q_clear(qpoly1);
    fmpz_poly_clear(zpoly);
    flint_free(strout);
    
    /* Canonicalise **********************************************************/
    
    fmpz_poly_q_init(qpoly1);
    str = "2  -1 1/2  0 1";
    fmpz_poly_q_set_str(qpoly1, str);
    strout = fmpz_poly_q_get_str(qpoly1);
    ans = !strcmp(str, strout);
    if (!ans)
    {
        printf("test_canonicalize: failed\n");
        printf("    Expected \"%s\", got \"%s\"\n", str, strout);
        abort();
    }
    fmpz_poly_q_clear(qpoly1);
    flint_free(strout);
    
    fmpz_poly_q_init(qpoly1);
    str = "2  -1 -1/2  0 1";
    fmpz_poly_q_set_str(qpoly1, "2  1 1/2  0 -1");
    strout = fmpz_poly_q_get_str(qpoly1);
    ans = !strcmp("2  -1 -1/2  0 1", strout);
    if (!ans)
    {
        printf("test_canonicalize: failed\n");
        printf("    Expected \"%s\", got \"%s\"\n", str, strout);
        abort();
    }
    flint_free(strout);
    fmpz_poly_q_clear(qpoly1);
    
    /* Initialization, memory management and basic operations ****************/
    
    test_set("0", "0");
    test_set("0/1  1", "0");
    test_set("3  -1 0 1/2  0 1", "3  -1 0 1/2  0 1");
    test_set("3  -1 0 1/2  1 1", "2  -1 1");
    
    test_set_si(-1, "1  -1");
    test_set_si(13, "1  13");
    test_set_si(0, "0");
    
    test_swap("3  -1 0 1/2  0 1", "1  2/1  3", "1  2/1  3", "3  -1 0 1/2  0 1");
    
    test_zero("0", "0");
    test_zero("0/1  1", "0");
    test_zero("3  -1 0 1/2  0 1", "0");
    
    test_neg("0", "0");
    test_neg("1  1/1  2", "1  -1/1  2");
    test_neg("3  -1 0 1/2  0 1", "3  1 0 -1/2  0 1");
    
    test_inv("1  1/1  2", "1  2");
    test_inv("3  -1 0 1/2  0 1", "2  0 1/3  -1 0 1");
    test_inv("3  -1 0 -1/2  0 1", "2  0 -1/3  1 0 1");
    
    test_inv_inplace("1  1/1  2", "1  2");
    test_inv_inplace("3  -1 0 1/2  0 1", "2  0 1/3  -1 0 1");
    test_inv_inplace("3  -1 0 -1/2  0 1", "2  0 -1/3  1 0 1");
    
    test_is_zero("0", 1);
    test_is_zero("0/1  1", 1);
    test_is_zero("3  -1 0 1/2  0 1", 0);
    test_is_zero("3  -1 0 1/2  1 1", 0);
    
    test_is_one("0", 0);
    test_is_one("0/1  1", 0);
    test_is_one("1  1/1  1", 1);
    test_is_one("2  1 1/2  1 1", 1);
    test_is_one("3  -1 0 1/2  0 1", 0);
    
    test_equal("1  1/1  2", "1  1/1  2", 1);
    test_equal("1  1/1  2", "1  1/1  2", 1);
    test_equal("3  -1 0 1/2  1 1", "2  -1 1", 1);
    test_equal("3  -1 0 1/2  -1 1", "2  -1 1", 0);
    
    /* Addition and subtraction **********************************************/
    
    test_add("3  1 0 1/2  0 1", "2  0 -1/3  1 0 1", "5  1 0 1 0 1/4  0 1 0 1");
    test_add("3  -1 0 1/2  1 1", "1  2/2  -1 1", "3  3 -2 1/2  -1 1");
    test_add("0/2  1 1", "1  2/1  1", "1  2");
    test_add("1  -3/1  4", "0/3  1 0 1", "1  -3/1  4");
    test_add("2  1 1/1  1", "2  -1 1/1  1", "2  0 2");
    test_add("2  1 1/2  0 1", "2  2 1/2  -1 1", "3  -1 2 2/3  0 -1 1");
    test_add("2  -1 1/2  2 1", "3  4 4 1/2  1 1", "4  7 12 7 1/3  2 3 1");
    test_add("2  1 1/2  -1 1", "2  1 1", "3  0 1 1/2  -1 1");
    test_add("1  1/2  1 1", "2  0 1/2  1 1", "1  1");
    test_add("2  1 1/3  4 -4 1", "1  1/2  -2 1", "2  -1 2/3  4 -4 1");
    test_add("3  0 1 1/3  1 2 1", "2  0 -1/2  1 1", "0");
    test_add("2  1 1/2  0 1", "2  -1 1/2  0 1", "1  2");
    test_add("1  1/3  3 5 2", "1  1/3  6 7 2", "1  1/3  2 3 1");
    
    test_add_in_place1("3  1 0 1/2  0 1", "2  0 -1/3  1 0 1", "5  1 0 1 0 1/4  0 1 0 1");
    test_add_in_place1("3  -1 0 1/2  1 1", "1  2/2  -1 1", "3  3 -2 1/2  -1 1");
    test_add_in_place1("0/2  1 1", "1  2/1  1", "1  2");
    test_add_in_place1("1  -3/1  4", "0/3  1 0 1", "1  -3/1  4");
    test_add_in_place1("2  1 1/1  1", "2  -1 1/1  1", "2  0 2");
    test_add_in_place1("2  1 1/2  0 1", "2  2 1/2  -1 1", "3  -1 2 2/3  0 -1 1");
    test_add_in_place1("2  -1 1/2  2 1", "3  4 4 1/2  1 1", "4  7 12 7 1/3  2 3 1");
    test_add_in_place1("2  1 1/2  -1 1", "2  1 1", "3  0 1 1/2  -1 1");
    test_add_in_place1("1  1/2  1 1", "2  0 1/2  1 1", "1  1");
    test_add_in_place1("2  1 1/3  4 -4 1", "1  1/2  -2 1", "2  -1 2/3  4 -4 1");
    test_add_in_place1("3  0 1 1/3  1 2 1", "2  0 -1/2  1 1", "0");
    test_add_in_place1("2  1 1/2  0 1", "2  -1 1/2  0 1", "1  2");
    
    test_add_in_place2("3  1 0 1/2  0 1", "2  0 -1/3  1 0 1", "5  1 0 1 0 1/4  0 1 0 1");
    test_add_in_place2("3  -1 0 1/2  1 1", "1  2/2  -1 1", "3  3 -2 1/2  -1 1");
    test_add_in_place2("0/2  1 1", "1  2/1  1", "1  2");
    test_add_in_place2("1  -3/1  4", "0/3  1 0 1", "1  -3/1  4");
    test_add_in_place2("2  1 1/1  1", "2  -1 1/1  1", "2  0 2");
    test_add_in_place2("2  1 1/2  0 1", "2  2 1/2  -1 1", "3  -1 2 2/3  0 -1 1");
    test_add_in_place2("2  -1 1/2  2 1", "3  4 4 1/2  1 1", "4  7 12 7 1/3  2 3 1");
    test_add_in_place2("2  1 1/2  -1 1", "2  1 1", "3  0 1 1/2  -1 1");
    test_add_in_place2("1  1/2  1 1", "2  0 1/2  1 1", "1  1");
    test_add_in_place2("2  1 1/3  4 -4 1", "1  1/2  -2 1", "2  -1 2/3  4 -4 1");
    test_add_in_place2("3  0 1 1/3  1 2 1", "2  0 -1/2  1 1", "0");
    test_add_in_place2("2  1 1/2  0 1", "2  -1 1/2  0 1", "1  2");
    
    test_add_in_place3("2  1 1", "2  2 2");
    test_add_in_place3("2  1 1/1  2", "2  1 1");
    
    test_sub("3  1 0 1/2  0 1", "2  0 -1/3  1 0 1", "5  1 0 3 0 1/4  0 1 0 1");
    test_sub("3  -1 0 1/2  1 1", "1  2/2  -1 1", "3  -1 -2 1/2  -1 1");
    test_sub("0/2  1 1", "1  2/1  1", "1  -2");
    test_sub("1  -3/1  4", "0/3  1 0 1", "1  -3/1  4");
    test_sub("2  1 1/1  1", "2  -1 1/1  1", "1  2");
    test_sub("2  1 1/2  0 1", "2  2 1/2  -1 1", "2  -1 -2/3  0 -1 1");
    test_sub("2  -1 1/2  2 1", "3  4 4 1/2  1 1", "4  -9 -12 -5 -1/3  2 3 1");
    test_sub("2  -1 1/2  0 1", "1  1", "1  -1/2  0 1");
    test_sub("3  1 0 1/2  0 1", "2  0 -1/3  1 0 1", "5  1 0 3 0 1/4  0 1 0 1");
    test_sub("3  -1 0 1/2  1 1", "1  2/2  -1 1", "3  -1 -2 1/2  -1 1");
    test_sub("0/2  1 1", "1  2/1  1", "1  -2");
    test_sub("1  -3/1  4", "0/3  1 0 1", "1  -3/1  4");
    test_sub("2  1 1/1  1", "2  -1 1/1  1", "1  2");
    test_sub("2  1 1/2  0 1", "2  2 1/2  -1 1", "2  -1 -2/3  0 -1 1");
    test_sub("2  -1 1/2  2 1", "3  4 4 1/2  1 1", "4  -9 -12 -5 -1/3  2 3 1");
    test_sub("2  1 1/2  -1 1", "2  1 1", "3  2 1 -1/2  -1 1");
    test_sub("1  1/2  1 1", "2  0 1/2  1 1", "2  1 -1/2  1 1");
    test_sub("2  1 1/3  4 -4 1", "1  1/2  -2 1", "1  3/3  4 -4 1");
    test_sub("3  0 1 1/3  1 2 1", "2  0 -1/2  1 1", "2  0 2/2  1 1");
    test_sub("2  1 1/2  0 1", "2  -1 1/2  0 1", "1  2/2  0 1");
    test_sub("1  1/3  3 5 2", "1  1/3  6 7 2", "1  1/4  6 13 9 2");
    test_sub("2  1 1/2  0 2", "2  1 1/2  0 2", "0");
    test_sub("2  -1 2/2  0 1", "2  -1 1/2  0 1", "1  1");
    
    test_sub_in_place1("3  1 0 1/2  0 1", "2  0 -1/3  1 0 1", "5  1 0 3 0 1/4  0 1 0 1");
    test_sub_in_place1("3  -1 0 1/2  1 1", "1  2/2  -1 1", "3  -1 -2 1/2  -1 1");
    test_sub_in_place1("0/2  1 1", "1  2/1  1", "1  -2");
    test_sub_in_place1("1  -3/1  4", "0/3  1 0 1", "1  -3/1  4");
    test_sub_in_place1("2  1 1/1  1", "2  -1 1/1  1", "1  2");
    test_sub_in_place1("2  1 1/2  0 1", "2  2 1/2  -1 1", "2  -1 -2/3  0 -1 1");
    test_sub_in_place1("2  -1 1/2  2 1", "3  4 4 1/2  1 1", "4  -9 -12 -5 -1/3  2 3 1");
    
    test_sub_in_place2("3  1 0 1/2  0 1", "2  0 -1/3  1 0 1", "5  1 0 3 0 1/4  0 1 0 1");
    test_sub_in_place2("3  -1 0 1/2  1 1", "1  2/2  -1 1", "3  -1 -2 1/2  -1 1");
    test_sub_in_place2("0/2  1 1", "1  2/1  1", "1  -2");
    test_sub_in_place2("1  -3/1  4", "0/3  1 0 1", "1  -3/1  4");
    test_sub_in_place2("2  1 1/1  1", "2  -1 1/1  1", "1  2");
    test_sub_in_place2("2  1 1/2  0 1", "2  2 1/2  -1 1", "2  -1 -2/3  0 -1 1");
    test_sub_in_place2("2  -1 1/2  2 1", "3  4 4 1/2  1 1", "4  -9 -12 -5 -1/3  2 3 1");
    
    test_sub_in_place3("2  -1 1/2  2 1", "0");
    
    test_addmul("1  1/2  0 2", "2  3 1/1  4", "3  1 0 1/4  -2 0 0 1", "5  -4 3 1 5 1/5  0 -8 0 0 4");
    
    test_submul("1  1/2  0 2", "2  3 1/1  4", "3  1 0 1/4  -2 0 0 1", "5  -4 -3 -1 -1 -1/5  0 -8 0 0 4");
    
    /* Scalar multiplication and devision ************************************/
    
    mpz_init_set_si(mpzzero, 0);
    mpz_init_set_si(mpzone, 1);
    mpz_init_set_si(mpztwo, 2);
    
    mpq_init(mpqzero); mpq_set_si(mpqzero, 0, 1);
    mpq_init(mpqone); mpq_set_si(mpqone, 1, 1);
    mpq_init(mpqtwo); mpq_set_si(mpqtwo, 2, 1);
    mpq_init(mpqtwoinv); mpq_set_si(mpqtwoinv, 1, 2);
    
    test_scalar_mul_si("0", 1, "0");
    test_scalar_mul_si("0", 0, "0");
    test_scalar_mul_si("1  2", 0, "0");
    test_scalar_mul_si("1  1/1  2", -2, "1  -1");
    test_scalar_mul_si("2  1 1/2  -2 3", 5, "2  5 5/2  -2 3");
    test_scalar_mul_si("2  1 1/2  -2 2", 3, "2  3 3/2  -2 2");
    
    test_scalar_mul_mpz("0", mpzone, "0");
    test_scalar_mul_mpz("0", mpzzero, "0");
    test_scalar_mul_mpz("1  2", mpzzero, "0");
    test_scalar_mul_mpz("1  1/1  2", mpztwo, "1  1");
    
    test_scalar_mul_mpq("0", mpqone, "0");
    test_scalar_mul_mpq("0", mpqzero, "0");
    test_scalar_mul_mpq("1  2", mpqzero, "0");
    test_scalar_mul_mpq("1  1/1  2", mpqtwo, "1  1");
    test_scalar_mul_mpq("1  -2/1  1", mpqtwoinv, "1  -1");
    
    test_scalar_div_si("0", 1, "0");
    test_scalar_div_si("1  2", 2, "1  1");
    test_scalar_div_si("1  1/1  2", -2, "1  -1/1  4");
    test_scalar_div_si("3  -5 0 3/2  1 1", 2, "3  -5 0 3/2  2 2");
    test_scalar_div_si("3  2 8 4/2  0 1", 3, "3  2 8 4/2  0 3");
    test_scalar_div_si("3  2 8 4/2  0 1", -3, "3  -2 -8 -4/2  0 3");
    test_scalar_div_si("3  -27 0 9/2  0 1", -3, "3  9 0 -3/2  0 1");
    
    test_scalar_div_mpz("0", mpzone, "0");
    test_scalar_div_mpz("1  2", mpztwo, "1  1");
    test_scalar_div_mpz("1  1/1  2", mpztwo, "1  1/1  4");
    
    test_scalar_div_mpq("0", mpqone, "0");
    test_scalar_div_mpq("1  2", mpqone, "1  2");
    test_scalar_div_mpq("1  1/1  2", mpqtwo, "1  1/1  4");
    test_scalar_div_mpq("1  -2/1  1", mpqtwoinv, "1  -4");
    
    mpz_clear(mpzzero);
    mpz_clear(mpzone);
    mpz_clear(mpztwo);
    mpq_clear(mpqzero);
    mpq_clear(mpqone);
    mpq_clear(mpqtwo);
    mpq_clear(mpqtwoinv);
    
    /* Multiplication, division and powing *********************************/
    
    test_mul("3  1 0 1/2  0 1", "2  0 -1/3  1 0 1", "1  -1");
    test_mul("3  -1 0 1/2  1 1", "1  2/2  -1 1", "1  2");
    test_mul("0/2  1 1", "1  2/1  1", "0");
    test_mul("1  -3/1  4", "0/3  1 0 1", "0");
    test_mul("2  1 1/1  1", "2  -1 1/1  1", "3  -1 0 1");
    test_mul("2  1 1/2  0 1", "2  2 1/2  -1 1", "3  2 3 1/3  0 -1 1");
    test_mul("2  -1 1/2  2 1", "3  4 4 1/2  1 1", "3  -2 1 1/2  1 1");
    
    test_mul_in_place1("3  1 0 1/2  0 1", "2  0 -1/3  1 0 1", "1  -1");
    test_mul_in_place1("3  -1 0 1/2  1 1", "1  2/2  -1 1", "1  2");
    test_mul_in_place1("0/2  1 1", "1  2/1  1", "0");
    test_mul_in_place1("1  -3/1  4", "0/3  1 0 1", "0");
    test_mul_in_place1("2  1 1/1  1", "2  -1 1/1  1", "3  -1 0 1");
    test_mul_in_place1("2  1 1/2  0 1", "2  2 1/2  -1 1", "3  2 3 1/3  0 -1 1");
    test_mul_in_place1("2  -1 1/2  2 1", "3  4 4 1/2  1 1", "3  -2 1 1/2  1 1");
    
    test_mul_in_place2("3  1 0 1/2  0 1", "2  0 -1/3  1 0 1", "1  -1");
    test_mul_in_place2("3  -1 0 1/2  1 1", "1  2/2  -1 1", "1  2");
    test_mul_in_place2("0/2  1 1", "1  2/1  1", "0");
    test_mul_in_place2("1  -3/1  4", "0/3  1 0 1", "0");
    test_mul_in_place2("2  1 1/1  1", "2  -1 1/1  1", "3  -1 0 1");
    test_mul_in_place2("2  1 1/2  0 1", "2  2 1/2  -1 1", "3  2 3 1/3  0 -1 1");
    test_mul_in_place2("2  -1 1/2  2 1", "3  4 4 1/2  1 1", "3  -2 1 1/2  1 1");
    
    test_mul_in_place3("2  0 1/2  1 1", "3  0 0 1/3  1 2 1");
    
    test_div("3  -1 0 1/1  2", "2  1 1/1  1", "2  -1 1/1  2");
    test_div("0/2  1 1", "2  1 1/1  1", "0");
    test_div("3  -1 0 1/1  4", "2  -1 -1/1  2", "2  1 -1/1  2");
    test_div("2  1 1", "2  1 -1/2  1 -1", "2  1 1");
    test_div("2  1 1/3  4 4 1", "2  -1 1/3  6 5 1", "3  3 4 1/3  -2 1 1");
    
    test_div_in_place1("3  -1 0 1/1  2", "2  1 1/1  1", "2  -1 1/1  2");
    test_div_in_place1("0/2  1 1", "2  1 1/1  1", "0");
    test_div_in_place1("3  -1 0 1/1  4", "2  -1 -1/1  2", "2  1 -1/1  2");
    test_div_in_place1("2  1 1", "2  1 -1/2  1 -1", "2  1 1");
    test_div_in_place1("2  1 1/3  4 4 1", "2  -1 1/3  6 5 1", "3  3 4 1/3  -2 1 1");
    test_div_in_place1("0", "1  2/2  3 5", "0");
    
    test_div_in_place2("3  -1 0 1/1  2", "2  1 1/1  1", "2  -1 1/1  2");
    test_div_in_place2("0/2  1 1", "2  1 1/1  1", "0");
    test_div_in_place2("3  -1 0 1/1  4", "2  -1 -1/1  2", "2  1 -1/1  2");
    test_div_in_place2("2  1 1", "2  1 -1/2  1 -1", "2  1 1");
    test_div_in_place2("2  1 1/3  4 4 1", "2  -1 1/3  6 5 1", "3  3 4 1/3  -2 1 1");
    
    test_div_in_place3("3  -1 0 1/1  2", "1  1");
    
    test_pow("2  0 -1/1  2", 3, "4  0 0 0 -1/1  8");
    test_pow("0", 0, "1  1");
    test_pow("2  1 -1", 0, "1  1");
    test_pow("2  1 1/2  0 1", 0, "1  1");
    
    /* Derivative ************************************************************/
    
    test_derivative("0", "0");
    test_derivative("1  2", "0");
    test_derivative("1  -1/1  2", "0");
    test_derivative("2  0 1", "1  1");
    test_derivative("3  1 0 1", "2  0 2");
    test_derivative("1  1/2  0 1", "1  -1/3  0 0 1");
    test_derivative("2  2 1/2  -1 1", "1  -3/3  1 -2 1");
    
    test_derivative("2  0 1/3  1 2 1", "2  1 -1/4  1 3 3 1");

    /* Bug which allowed constant factors */
    test_derivative("3  5 1 -2/2  10 2", "3  0 -10 -1/3  25 10 1");
    
    /* Evaluation ************************************************************/
    
    test_evaluate("1  1/1  2", -2, 3, "1/2");
    test_evaluate("3  1 0 1/2  0 1", -1, 2, "-5/2");
    test_evaluate("2  3 1/2  -1 1", 1, 1, "P");
    test_evaluate("2  3 1/2  -1 1", 2, 3, "-11");
    test_evaluate("2  3 1/2  -1 2", 1, 2, "P");
    test_evaluate("2  1 1/2  -1 1", 2, 1, "3");
    
    /* String methods ********************************************************/
    
    fmpz_poly_q_init(qpoly1);
    ans = fmpz_poly_q_set_str(qpoly1, "1  3/xyz");
    if ((ans == 0) || !fmpz_poly_q_is_zero(qpoly1))
    {
        printf("test_set_str: failed\n");
        abort();
    }
    fmpz_poly_q_clear(qpoly1);
    
    fmpz_poly_q_init(qpoly1);
    ans = fmpz_poly_q_set_str(qpoly1, "abc/1  3");
    if ((ans == 0) || !fmpz_poly_q_is_zero(qpoly1))
    {
        printf("test_set_str: failed\n");
        abort();
    }
    fmpz_poly_q_clear(qpoly1);
    
    fmpz_poly_q_init(qpoly1);
    ans = fmpz_poly_q_set_str(qpoly1, "abc/xyz");
    if ((ans == 0) || !fmpz_poly_q_is_zero(qpoly1))
    {
        printf("test_set_str: failed\n");
        abort();
    }
    fmpz_poly_q_clear(qpoly1);
    
    test_get_str_pretty("1  -3", "-3");
    test_get_str_pretty("3  1 2 1", "t^2+2*t+1");
    test_get_str_pretty("1  -2/2  1 1", "-2/(t+1)");
    test_get_str_pretty("2  1 1/2  -1 1", "(t+1)/(t-1)");
    test_get_str_pretty("2  1 1/1  2", "(t+1)/2");
    test_get_str_pretty("1  1/1  2", "1/2");

    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

