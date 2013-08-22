/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2013 Tom Bachmann

******************************************************************************/

#include <iostream>

#include "nmod_matxx.h"

#include "flintxx/test/helpers.h"

using namespace flint;

void
test_init()
{
    mp_limb_t M = 1039;
    nmod_matxx A(3, 4, M);
    nmodxx_ctx_srcref ctx = A.estimate_ctx();
    tassert(ctx.n() == M);
    tassert((A + A).modulus() == M);
    tassert(A.rows() == 3 && A.cols() == 4);
    tassert(A.at(0, 0) == nmodxx::red(0, ctx));
    A.at(0, 0) = nmodxx::red(1, ctx);

    nmod_matxx B(A);
    tassert(A == B);
    tassert(B.rows() == 3 && B.cols() == 4);
    tassert(B.at(0, 0) == nmodxx::red(1, ctx));
    B.at(0, 0) = nmodxx::red(0, ctx);
    tassert(A.at(0, 0) == nmodxx::red(1, ctx));
    tassert(A != B);

    B = A;
    tassert(A == B);
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
    nmod_matxx A(10, 10, M);
    nmod_matxx v(10, 1, M);
    nmodxx_ctx_srcref ctx = A.estimate_ctx();
    for(unsigned i = 0;i < 10;++i)
        v.at(i, 0) = nmodxx::red(i, ctx);
    nmodxx two = nmodxx::red(2, ctx);

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

    tassert(trace(transpose(v)) == nmodxx::red(0, ctx));
    tassert(trace(A + v*transpose(v)) == nmodxx::red(285, ctx));
    tassert(trace(v*transpose(v) + A) == nmodxx::red(285, ctx));
    tassert(trace(v*transpose(v) + v*transpose(v)) == nmodxx::red(2*285, ctx));
    tassert(trace((A+A)*(nmodxx::red(1, ctx) + nmodxx::red(1, ctx)))
            == nmodxx::red(0, ctx));

    for(unsigned i = 0;i < 10; ++i)
        for(unsigned j = 0; j < 10; ++j)
            A.at(i, j) = nmodxx::red(i*j, ctx);
    tassert(A == v*transpose(v));
    tassert(A != transpose(v)*v);
    A.at(0, 0) = nmodxx::red(15, ctx);
    tassert(A != v*transpose(v));

    A.at(0, 0) = nmodxx::red(0, ctx);
    for(unsigned i = 0;i < 10; ++i)
        for(unsigned j = 0; j < 10; ++j)
            A.at(i, j) *= two;
    tassert(A == v*transpose(v) + v*transpose(v));
    tassert(A - v*transpose(v) == v*transpose(v));
    tassert(((-A) + A).is_zero());
    tassert((A + A).at(0, 0) == A.at(0, 0) + A.at(0, 0));
}

void
test_functions()
{
    mp_limb_t M = 1031;
    nmod_matxx A(2, 3, M), B(2, 2, M), empty(0, 15, M);
    nmodxx_ctx_srcref ctx = A.estimate_ctx();
    B.at(0, 0) = nmodxx::red(1, ctx);
    tassert(A.is_zero() && !A.is_empty() && !A.is_square());
    tassert(!B.is_zero() == B.is_square());
    tassert(empty.is_zero() && empty.is_empty());

    // transpose tested in arithmetic
    // mul tested in arithmetic
    // trace tested in arithmetic

    frandxx rand;
    A.set_randtest(rand);
    B.set_randtest(rand);
    tassert(B*A == B.mul_classical(A));
    tassert(B*A == B.mul_strassen(A));

    B.set_randrank(rand, 1);
    tassert(B.det() == nmodxx::red(0, ctx));
    B.set_randrank(rand, 2);
    tassert(B.det() != nmodxx::red(0, ctx));

    B.set_randrank(rand, 1);
    assert_exception(B.inv().evaluate());

    B.set_randrank(rand, 2);
    nmod_matxx eye(2, 2, M);
    eye.at(0, 0) = nmodxx::red(1, ctx);eye.at(1, 1) = nmodxx::red(1, ctx);
    tassert(B.inv() * B == eye);

    A.set_randrank(rand, 2);
    tassert(rank(A) == 2);

#if 0
    tassert(B.charpoly().get_coeff(0) == B.det());
    tassert(charpoly(B).get_coeff(1) == -B.trace());
    tassert(charpoly(B).lead() == 1);
#endif

    B.set_randtril(rand, false);
    tassert(B*B.solve_tril(A, false) == A);
    tassert(B.solve_tril_classical(A, false) == B.solve_tril(A, false));
    tassert(B.solve_tril_recursive(A, false) == B.solve_tril(A, false));
    B.set_randtriu(rand, true);
    tassert(B*B.solve_triu(A, true) == A);
    tassert(B.solve_triu_classical(A, true) == B.solve_triu(A, true));
    tassert(B.solve_triu_recursive(A, true) == B.solve_triu(A, true));

    B.set_randrank(rand, 2);
    tassert(B*B.solve(A) == A);
    nmod_vecxx X(2, ctx); X[0] = nmodxx::red(1, ctx); X[1] = nmodxx::red(2, ctx);
    X = B.solve(X);
    tassert(B.at(0, 0)*X[0] + B.at(0, 1) * X[1] == nmodxx::red(1, ctx));
    tassert(B.at(1, 0)*X[0] + B.at(1, 1) * X[1] == nmodxx::red(2, ctx));

    B.set_randrank(rand, 1);
    assert_exception(B.solve(A).evaluate());
    assert_exception(B.solve(X).evaluate());

    slong nullity;nmod_matxx C(3, 3, M);
    tassert(nullspace(A).get<1>().rows() == 3);
    tassert(nullspace(A).get<1>().cols() == 3);
    ltupleref(nullity, C) = nullspace(A);
    tassert(nullity == 3 - rank(A));
    tassert(C.rank() == nullity);
    tassert((A*C).is_zero());

    A.set_rref();
    tassert(A.at(1, 0) == nmodxx::red(0, ctx));
}

void
test_extras()
{
    // TODO
}

void
test_randomisation()
{
    frandxx rand;
    mp_limb_t M = 1031;
    nmod_matxx A(2, 2, M);
    nmodxx_ctx_srcref ctx = A.estimate_ctx();

    // not really anything we can test about these ...
    // just make sure the call works
    A.set_randtest(rand);
    A.set_randfull(rand);


    nmod_vecxx v(2, ctx);v[0] = nmodxx::red(5, ctx);v[1] = nmodxx::red(7, ctx);
    A.set_randpermdiag(rand, v);
    tassert(A.at(0, 0) + A.at(0, 1) + A.at(1, 0) + A.at(1, 1)
            == nmodxx::red(5 + 7, ctx));

    A.set_randrank(rand, 1);
    tassert(A.rank() == 1);
    A.apply_randops(rand, 17);
    tassert(A.rank() == 1);

    A.set_randtril(rand, true);
    tassert(A.at(0, 0) == nmodxx::red(1, ctx));
    tassert(A.at(1, 1) == nmodxx::red(1, ctx));
    tassert(A.at(0, 1) == nmodxx::red(0, ctx));

    A.set_randtriu(rand, false);
    tassert(A.at(1, 0) == nmodxx::red(0, ctx));
}

int
main()
{
    std::cout << "nmod_matxx....";

    test_init();
    test_arithmetic();
    test_functions();
    test_extras();
    test_randomisation();

    std::cout << "PASS" << std::endl;
    return 0;
}

