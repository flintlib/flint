/*
    Copyright (C) 2013 Tom Bachmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <iostream>

#include "nmod_poly_matxx.h"

#include "flintxx/test/helpers.h"

using namespace flint;

void
test_init()
{
    mp_limb_t M = 1039;
    nmod_poly_matxx A(3, 4, M);
    nmodxx_ctx_srcref ctx = A.estimate_ctx();
    tassert(ctx.n() == M);
    tassert((A + A).modulus() == M);
    tassert(A.rows() == 3 && A.cols() == 4);
    tassert(A.at(0, 0) == nmod_polyxx::from_ground(0, ctx));
    A.at(0, 0) = nmod_polyxx::from_ground(1, ctx);

    nmod_poly_matxx B(A);
    tassert(A == B);
    tassert(B.rows() == 3 && B.cols() == 4);
    tassert(B.at(0, 0) == nmod_polyxx::from_ground(1, ctx));
    B.at(0, 0) = nmod_polyxx::from_ground(0, ctx);
    tassert(A.at(0, 0) == nmod_polyxx::from_ground(1, ctx));
    tassert(A != B);

    B = A;
    tassert(A == B);

    frandxx state;
    nmod_matxx C(A.rows(), A.cols(), A.modulus());
    C.set_randtest(state);
    A = nmod_poly_matxx::from_ground(C);
    for(slong i = 0;i < A.rows();++i)
        for(slong j = 0;j < A.cols();++j)
            tassert(A.at(i, j) == nmod_polyxx::from_ground(C.at(i, j)));

    tassert(nmod_poly_matxx::zero(2, 2, M).is_zero()
            && nmod_poly_matxx::one(2, 2, M).is_one());
}

template<class Expr>
bool has_explicit_temporaries(const Expr&)
{
    return Expr::ev_traits_t::rule_t::temporaries_t::len != 0;
}
void
test_arithmetic()
{
    mp_limb_t M = 1039;
    nmod_poly_matxx A(10, 10, M);
    nmod_poly_matxx v(10, 1, M);
    nmodxx_ctx_srcref ctx = A.estimate_ctx();
    for(unsigned i = 0;i < 10;++i)
        v.at(i, 0) = nmod_polyxx::from_ground(i, ctx);
    nmod_polyxx two = nmod_polyxx::from_ground(2, ctx);

    tassert(transpose(v).rows() == 1);
    tassert(v.transpose().cols() == 10);
    tassert((two*v).rows() == 10);
    tassert((v*two).rows() == 10);
    tassert((v*transpose(v)).rows() == 10 && (v*transpose(v)).cols() == 10);

    tassert(!has_explicit_temporaries(trace(transpose(v))));
    tassert(!has_explicit_temporaries(trace(A + v*transpose(v))));
    tassert(!has_explicit_temporaries(A + v*transpose(v)));
    tassert(!has_explicit_temporaries(trace((v*transpose(v) + A))));
    tassert(!has_explicit_temporaries(trace(v*transpose(v) + v*transpose(v))));
    tassert(!has_explicit_temporaries(v*transpose(v) + v*transpose(v)));

    tassert(trace(transpose(v)) == nmod_polyxx::from_ground(0, ctx));
    tassert(trace(A + v*transpose(v)) == nmod_polyxx::from_ground(285, ctx));
    tassert(trace(v*transpose(v) + A) == nmod_polyxx::from_ground(285, ctx));
    tassert(trace(v*transpose(v) + v*transpose(v))
            == nmod_polyxx::from_ground(2*285, ctx));
    tassert(trace((A+A)*(nmod_polyxx::from_ground(1, ctx)
                    + nmod_polyxx::from_ground(1, ctx)))
            == nmod_polyxx::from_ground(0, ctx));

    for(unsigned i = 0;i < 10; ++i)
        for(unsigned j = 0; j < 10; ++j)
            A.at(i, j) = nmod_polyxx::from_ground(i*j, ctx);
    tassert(A == v*transpose(v));
    tassert(A != transpose(v)*v);
    A.at(0, 0) = nmod_polyxx::from_ground(15, ctx);
    tassert(A != v*transpose(v));

    A.at(0, 0) = nmod_polyxx::from_ground(0, ctx);
    for(unsigned i = 0;i < 10; ++i)
        for(unsigned j = 0; j < 10; ++j)
            A.at(i, j) *= two;
    tassert(A == v*transpose(v) + v*transpose(v));
    tassert(A - v*transpose(v) == v*transpose(v));
    tassert(((-A) + A).is_zero());
    tassert((A + A).at(0, 0) == A.at(0, 0) + A.at(0, 0));

    tassert(A * nmodxx::red(17, ctx) == A * nmod_polyxx::from_ground(17, ctx));

    frandxx rand;
    nmodxx x = nmodxx::red(17, ctx);
    A.set_randtest(rand, 5);
    nmod_matxx B(A.rows(), A.cols(), M);
    B = A(x);
    for(slong i = 0;i < A.rows();++i)
        for(slong j = 0;j < A.cols();++j)
            tassert(B.at(i, j) == A.at(i, j)(x));
    tassert(A(x) == evaluate(A, x));
}

void
test_functions()
{
    mp_limb_t M = 1031;
    nmod_poly_matxx A(2, 3, M), B(2, 2, M), empty(0, 15, M);
    nmodxx_ctx_srcref ctx = A.estimate_ctx();
    B.at(0, 0) = nmod_polyxx::from_ground(1, ctx);
    tassert(A.is_zero() && !A.is_empty() && !A.is_square() && !A.is_one());
    tassert(!B.is_zero() == B.is_square());
    tassert(empty.is_zero() && empty.is_empty());
    B.at(1, 1) = B.at(0, 0);
    tassert(B.is_one());

    // transpose tested in arithmetic
    // mul tested in arithmetic
    // trace tested in arithmetic

    A.at(0, 0).set_coeff(35, 1);
    tassert(A.max_length() == 36);

    frandxx rand;
    A.set_randtest(rand, 5);
    B.set_randtest(rand, 5);
    tassert(B*A == B.mul_classical(A));
    tassert(B*A == B.mul_KS(A));
    tassert(B*A == B.mul_interpolate(A));

    tassert(B.sqr() == B*B);
    tassert(B.sqr_classical() == B*B);
    tassert(B.sqr_KS() == B*B);
    tassert(B.sqr_interpolate() == B*B);

    tassert(B.pow(5u) == B*B.sqr().sqr());

    nmod_matxx Bp(B.rows(), B.cols(), B.modulus());
    Bp.set_randrank(rand, 1);
    tassert(nmod_poly_matxx::from_ground(Bp).det().is_zero());
    Bp.set_randrank(rand, 2);
    tassert(nmod_poly_matxx::from_ground(Bp).det()
            == nmod_polyxx::from_ground(Bp.det()));

    Bp.set_randrank(rand, 1);
    tassert(inv(nmod_poly_matxx::from_ground(Bp)).get<0>() == false);

    Bp.set_randrank(rand, 2);
    bool worked;nmod_polyxx den(B.modulus());
    ltupleref(worked, B, den) = inv(nmod_poly_matxx::from_ground(Bp));
    tassert(worked && B*nmod_poly_matxx::from_ground(Bp)*A == A*den);

    tassert(rank(B) == 2);

    Bp.set_randrank(rand, 1);
    tassert(nmod_poly_matxx::from_ground(Bp).solve(A).get<0>() == false);
    Bp.set_randrank(rand, 2);
    nmod_poly_matxx P(A.rows(), A.cols(), A.modulus());
    ltupleref(worked, P, den) = nmod_poly_matxx::from_ground(Bp).solve(A);
    tassert(worked && nmod_poly_matxx::from_ground(Bp)*P == A*den);
    B = nmod_poly_matxx::from_ground(Bp);
    tassert(B.solve(A) == B.solve_fflu(A));

    permxx perm(B.rows());
    tassert(solve_fflu_precomp(perm, B.fflu(&perm, false).get<1>().evaluate(), A)
            == B.solve_fflu(A).get<1>());

    Bp.set_randtest(rand);
    B = nmod_poly_matxx::from_ground(Bp);
    slong nullity;nmod_poly_matxx C(2, 2, M);
    tassert(nullspace(B).get<1>().rows() == 2);
    tassert(nullspace(B).get<1>().cols() == 2);
    ltupleref(nullity, C) = nullspace(B);
    tassert(nullity == 2 - rank(B));
    tassert(C.rank() == nullity);
    tassert((B*C).is_zero());

    B.set_zero();tassert(B.is_zero());
}

void
test_randomisation()
{
    frandxx rand, rand2;
    mp_limb_t M = 1031;
    nmod_poly_matxx A(2, 2, M);

    A.set_randtest(rand, 17);
    tassert(A.at(0, 0).length() <= 17);
    tassert(A == nmod_poly_matxx::randtest(2, 2, M, rand2, 17));
    A.set_randtest_sparse(rand, 17, 0.5);
    tassert(A.at(0, 0).length() <= 17);
    tassert(A == nmod_poly_matxx::randtest_sparse(2, 2, M, rand2, 17, 0.5));
}

void
test_row_reduction()
{
    frandxx state;
    nmod_poly_matxx A = nmod_poly_matxx::randtest(5, 5, 1031, state, 7);
    slong rank1, rank2;
    nmod_polyxx den1(A.modulus()), den2(A.modulus());
    nmod_poly_matxx res1(A.rows(), A.cols(), A.modulus());
    nmod_poly_matxx res2(res1);

    tassert(find_pivot_any(A, 2, 4, 1)
            == nmod_poly_mat_find_pivot_any(A._mat(), 2, 4, 1));
    tassert(find_pivot_partial(A, 2, 4, 1)
            == nmod_poly_mat_find_pivot_partial(A._mat(), 2, 4, 1));
    tassert(A.fflu(0, false).get<1>().rows() == A.rows());
    permxx p1(5), p2(5);
    ltupleref(rank1, res1, den1) = fflu(A, &p1);
    rank2 = nmod_poly_mat_fflu(res2._mat(), den2._poly(), p2._data(),
            A._mat(), false);
    tassert(rank1 == rank2 && res1 == res2 && p1 == p2 && den1 == den2);
    tassert(rank1 == A.fflu(0, false).get<0>());

    ltupleref(rank1, res1, den1) = rref(A);
    rank2 = nmod_poly_mat_rref(res2._mat(), den2._poly(), A._mat());
    tassert(rank1 == rank2 && res1 == res2 && p1 == p2 && den1 == den2);
}


void
test_printing()
{
    if(0)
        print_pretty(nmod_poly_matxx::zero(2, 2, 7), "x"); // make sure this compiles
}

int
main()
{
    std::cout << "nmod_poly_matxx....";

    test_init();
    test_arithmetic();
    test_functions();
    test_randomisation();
    test_row_reduction();
    test_printing();

    std::cout << "PASS" << std::endl;
    return 0;
}

