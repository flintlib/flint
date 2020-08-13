/*
    Copyright (C) 2013 Tom Bachmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <sstream>
#include <string>

#include "fmpz_poly_matxx.h"
#include "flintxx/test/helpers.h"

using namespace flint;

void
test_init()
{
    fmpz_poly_matxx A(3, 4);
    tassert(A.rows() == 3 && A.cols() == 4);
    tassert(A.at(0, 0) == fmpz_polyxx::from_ground(0));
    A.at(0, 0) = fmpz_polyxx::from_ground(1);

    fmpz_poly_matxx B(A);
    tassert(B.rows() == 3 && B.cols() == 4);
    tassert(B.at(0, 0) == fmpz_polyxx::from_ground(1));
    B.at(0, 0) = fmpz_polyxx::from_ground(0);
    tassert(A.at(0, 0) == fmpz_polyxx::from_ground(1));

    tassert(fmpz_poly_matxx::zero(3, 3).is_zero());
    tassert(fmpz_poly_matxx::one(3, 3).is_one());
}

template<class Expr>
bool has_explicit_temporaries(const Expr&)
{
    return Expr::ev_traits_t::rule_t::temporaries_t::len != 0;
}
void
test_arithmetic()
{
    fmpz_poly_matxx A(10, 10);
    fmpz_poly_matxx v(10, 1);
    for(unsigned i = 0;i < 10;++i)
        v.at(i, 0) = fmpz_polyxx::from_ground(i);

    tassert(transpose(v).rows() == 1);
    tassert(v.transpose().cols() == 10);
    tassert((fmpzxx(2)*v).rows() == 10);
    tassert((v*fmpzxx(2)).rows() == 10);
    tassert((v*transpose(v)).rows() == 10
            && (v*transpose(v)).cols() == 10);
    tassert(mul_classical(v, transpose(v)).rows() == 10);
    tassert(mul_KS(v, transpose(v)).cols() == 10);

    tassert(!has_explicit_temporaries(trace(transpose(v))));
    tassert(!has_explicit_temporaries(trace(A + v*transpose(v))));
    tassert(!has_explicit_temporaries(A + v*transpose(v)));
    tassert(!has_explicit_temporaries(trace((v*transpose(v) + A))));
    tassert(!has_explicit_temporaries(trace(v*transpose(v) + v*transpose(v))));
    tassert(!has_explicit_temporaries(v*transpose(v) + v*transpose(v)));

    tassert((transpose(v)).trace() == fmpz_polyxx::from_ground(0));
    tassert(trace(A + v*transpose(v)) == fmpz_polyxx::from_ground(285));
    tassert(trace(v*transpose(v) + A) == fmpz_polyxx::from_ground(285));
    tassert(trace(v*transpose(v) + v*transpose(v))
            == fmpz_polyxx::from_ground(2*285));
    tassert(trace((A+A)*(fmpz_polyxx(1) + fmpz_polyxx(1))).is_zero());

    for(unsigned i = 0;i < 10; ++i)
        for(unsigned j = 0; j < 10; ++j)
            A.at(i, j) = fmpz_polyxx::from_ground(i*j);
    tassert(A == v*transpose(v));
    tassert(A != transpose(v)*v);
    A.at(0, 0) = fmpz_polyxx::from_ground(15);
    tassert(A != v*transpose(v));

    A.at(0, 0) = fmpz_polyxx::from_ground(0);
    for(unsigned i = 0;i < 10; ++i)
        for(unsigned j = 0; j < 10; ++j)
            A.at(i, j) *= 2;
    tassert(A == v*transpose(v) + v*transpose(v));
    tassert(A - v*transpose(v) == v*transpose(v));
    tassert(((-A) + A).is_zero());
    tassert((A + A).at(0, 0) == A.at(0, 0) + A.at(0, 0));

    tassert((A + A) == fmpzxx(2)*A && fmpz_polyxx::from_ground(2)*A == (A + A));

    frandxx rand;
    fmpzxx x(17);
    A.set_randtest(rand, 3, 5);
    fmpz_matxx B(A.rows(), A.cols());
    B = A(x);
    for(slong i = 0;i < A.rows();++i)
        for(slong j = 0;j < A.cols();++j)
            tassert(B.at(i, j) == A.at(i, j)(x));
    tassert(A(x) == evaluate(A, x));
}

void
test_functions()
{
    fmpz_poly_matxx A(2, 3), B(2, 2), empty(0, 15);
    B.at(0, 0) = fmpz_polyxx::from_ground(1);
    tassert(A.is_zero() && !A.is_empty() && !A.is_square() && !A.is_one());
    tassert(!B.is_zero() == B.is_square());
    tassert(empty.is_zero() && empty.is_empty());

    // transpose tested in arithmetic
    // mul tested in arithmetic
    // trace tested in arithmetic

    A.at(0, 0).set_coeff(35, 17);
    tassert(A.max_length() == 36);
    tassert(A.max_bits() == 5);

    frandxx rand;
    A.set_randtest(rand, 4, 10);
    B.set_randtest(rand, 4, 10);
    tassert(B*A == B.mul_classical(A));
    tassert(B*A == B.mul_KS(A));

    fmpz_poly_matxx tmp(B*A);
    tmp.truncate(3);
    tassert(tmp == B.mullow(A, 3));

    tassert(sqr(B) == B*B);
    tassert(B.sqr().sqr() == pow(B, 4u));
    tassert(B.sqrlow(3) == B.mullow(B, 3));
    tmp = pow(B, 5u);
    tmp.truncate(7);
    tassert(tmp == B.pow_trunc(5u, 7));

    B.set_randtest(rand, 4, 10);
    tassert(B.det() == B.det_fflu());
    tassert(B.det()(fmpzxx(123)) == B(fmpzxx(123)).det());
    tassert(B.det() == B.det_interpolate());

    fmpz_matxx Bp(2, 2);
    Bp.set_randdet(rand, fmpzxx(2*3*5));
    tassert(fmpz_poly_matxx::from_ground(Bp).det()
            == fmpz_polyxx::from_ground(2*3*5));

    fmpz_poly_matxx C(3, 3);
    C.at(0, 0).set_coeff(0, 1);
    C.at(1, 1).set_coeff(0, 1);
    tassert(rank(C) == 2);

    Bp.set_randrank(rand, 1, 10);
    B = fmpz_poly_matxx::from_ground(Bp);
    tassert(!inv(B).get<0>());
    Bp.set_randrank(rand, 2, 10);
    B = fmpz_poly_matxx::from_ground(Bp);
    fmpz_poly_matxx Binv(2, 2); bool worked; fmpz_polyxx d;
    ltupleref(worked, Binv, d) = inv(B);
    tassert(worked);
    fmpz_poly_matxx eye(2, 2);
    eye.at(0, 0).set_coeff(0, 1);eye.at(1, 1).set_coeff(0, 1);
    tassert(eye.is_one());
    tassert(Binv * B == d*eye);

    fmpz_poly_matxx X(2, 3);
    ltupleref(worked, X, d) = solve(B, A);
    tassert(worked == true && (B*X) == A*d);
    ltupleref(worked, X, d) = B.solve_fflu(A);
    tassert(worked == true && (B*X) == A*d);
    tassert(solve(B, A).get<1>() == X);

    permxx perm(B.rows());
    fmpz_poly_matxx F(B.rows(), B.rows());
    fmpz_polyxx Fd, Xd;
    slong rk;
    // Note: fflu has a (false?) dependency on perm - need perm = id before call
    ltupleref(rk, F, Fd) = B.fflu(&perm, false);
    ltupleref(worked, X, Xd) = B.solve_fflu(A);
    tassert(worked && Xd*solve_fflu_precomp(perm, F, A) == Fd * X);

    slong nullity;
    tassert(nullspace(A).get<1>().rows() == 3);
    tassert(nullspace(A).get<1>().cols() == 3);
    ltupleref(nullity, C) = nullspace(A);
    tassert(nullity == 3 - rank(A));
    tassert(C.rank() == nullity);
    tassert((A*C).is_zero());

    if(0)
        print_pretty(A, "x"); // make sure this compiles
}

void
test_extras()
{
    fmpz_poly_matxx A(2, 2);
    A.at(0, 0).set_coeff(0, 1);

    fmpz_poly_matxx_srcref Asr(A);
    const fmpz_poly_matxx& Acr = A;
    tassert(A.at(0, 0) == Acr.at(0, 0));
    tassert(A.at(0, 0) == Asr.at(0, 0));
}

void
test_randomisation()
{
    frandxx rand, rand2;
    fmpz_poly_matxx A(2, 2);
    A.set_randtest(rand, 4, 5);
    tassert(abs(A.at(0, 0).get_coeff(0)) <= 31);
    tassert(A == fmpz_poly_matxx::randtest(2, 2, rand2, 4, 5));
    A.set_randtest_unsigned(rand, 4, 5);
    tassert(A.at(0, 0).get_coeff(0) >= 0);
    tassert(A == fmpz_poly_matxx::randtest_unsigned(2, 2, rand2, 4, 5));
    A.set_randtest_sparse(rand, 4, 5, 0.5);
    tassert(abs(fmpz_polyxx_get_coeff(A.at(0, 0), 0)) <= 31);
    tassert(A == fmpz_poly_matxx::randtest_sparse(2, 2, rand2, 4, 5, 0.5));
}

void
test_row_reduction()
{
    frandxx state;
    fmpz_poly_matxx A = fmpz_poly_matxx::randtest(5, 5, state, 7, 15);
    slong rank1, rank2;
    fmpz_polyxx den1, den2;
    fmpz_poly_matxx res1(A.rows(), A.cols()), res2(A.rows(), A.cols());

    tassert(find_pivot_any(A, 2, 4, 1)
            == fmpz_poly_mat_find_pivot_any(A._mat(), 2, 4, 1));
    tassert(find_pivot_partial(A, 2, 4, 1)
            == fmpz_poly_mat_find_pivot_partial(A._mat(), 2, 4, 1));
    tassert(A.fflu(0, false).get<1>().rows() == A.rows());
    permxx p1(5), p2(5);
    ltupleref(rank1, res1, den1) = fflu(A, &p1);
    rank2 = fmpz_poly_mat_fflu(res2._mat(), den2._poly(), p2._data(),
            A._mat(), false);
    tassert(rank1 == rank2 && res1 == res2 && p1 == p2 && den1 == den2);
    tassert(rank1 == A.fflu(0, false).get<0>());

    ltupleref(rank1, res1, den1) = rref(A);
    rank2 = fmpz_poly_mat_rref(res2._mat(), den2._poly(), A._mat());
    tassert(rank1 == rank2 && res1 == res2 && p1 == p2 && den1 == den2);
}

void
test_prod()
{
    fmpz_poly_mat_vecxx v1(10, 3, 3), v2(10, 3, 3), v3(9, 3, 3), v4(v1);
    tassert(v1 == v2);
    tassert(v1 != v3);
    v1[0].at(0, 0).set_coeff(0, 7u);
    tassert(v1 != v4);

    frandxx rand;
    fmpz_poly_matxx prod = fmpz_poly_matxx::one(3, 3);
    for(slong i = 0;i < v1.size();++i)
    {
        v1[i].set_randtest(rand, 4, 17);
        prod *= v1[i];
    }
    tassert(flint::prod(v1) == prod);
}

int
main()
{
    std::cout << "fmpz_poly_matxx....";

    test_init();
    test_arithmetic();
    test_functions();
    test_extras();
    test_randomisation();
    test_row_reduction();
    test_prod();

    std::cout << "PASS" << std::endl;
    return 0;
}
