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

#include<sstream>

#include "flintxx/expression.h"
#include "flintxx/tuple.h"

#include "flintxx/test/helpers.h"
#include "flintxx/test/myint.h"

using namespace flint;
using namespace mp;
using namespace traits;

void
test_initialization()
{
    myint a;
    tassert(a._data().payload == -1 && a._data().extra == 42);

    a._data().payload = 17;
    myint b = a;
    tassert(b._data().payload == 17);

    myint c(8);
    tassert(c._data().payload == 8 && c._data().extra == 2);

    myint d('a');
    tassert(d._data().payload == 'a' && d._data().extra == 3);

    myint e(c + c);
    tassert(e._data().payload == 16);
}

void
test_destruction()
{
    bool destroyed = false;
    {
        myint a;
        a._data().destroyed = &destroyed;
    }
    tassert(destroyed);
}

void
test_printing()
{
    std::ostringstream o;
    myint a(4);
    o << a;
    tassert(o.str() == "4");
}

void
test_assignment()
{
    myint a;
    myint b(43);
    a = b;
    tassert(a._data().payload == 43);
    tassert(a._data().extra == 4);

    a = 44l;
    tassert(a._data().payload == 44 && a._data().extra == 5);

    mylong c(0);
    c = a + b;
    tassert(c == 87l);
}

void
test_swap()
{
    myint a(2), b(3);
    swap(a, b);
    tassert(a == 3 && b == 2 && a._data().extra == 1234);
}

void
test_traits()
{
    typedef myint immediate_expression;
    typedef int immediate_nonexpression;
    typedef my_expression<
        operations::plus,
        make_tuple<const myint&, const myint&>::type
      > lazy_expression;

    tassert(is_expression<immediate_expression>::val == true);
    tassert(is_expression<immediate_nonexpression>::val == false);
    tassert(is_expression<lazy_expression>::val == true);

    tassert(is_immediate_expr<immediate_expression>::val == true);
    tassert(is_immediate_expr<immediate_nonexpression>::val == false);
    tassert(is_immediate_expr<lazy_expression>::val == false);

    tassert(is_immediate<immediate_expression>::val == true);
    tassert(is_immediate<immediate_nonexpression>::val == true);
    tassert(is_immediate<lazy_expression>::val == false);

    tassert(is_lazy_expr<immediate_expression>::val == false);
    tassert(is_lazy_expr<immediate_nonexpression>::val == false);
    tassert(is_lazy_expr<lazy_expression>::val == true);
    //tassert(is_lazy_expr<lazy_expression&>::val == false);
    tassert(is_lazy_expr<const immediate_expression&>::val == false);
}

void
test_equals()
{
    myint a(3);
    myint b(4);
    myint c(3);

    tassert(a != b);
    tassert(a == c);

    tassert(a == 3);
    tassert(3 == a);

    tassert(a != 4);
    tassert(4 != a);

    tassert((a + b) == (a + 4));
}

void
test_arithmetic()
{
    myint a(3);
    myint b(4);
    myint c(7);
    myint d(1);
    myint e(2);

    tassert((a + b).evaluate() == 7);
    tassert(a + b == c);
    tassert(a + b == 7);

    tassert(a - b == -1);
    tassert(a * b == 12);
    tassert(c / e == 3);
    tassert(c % e == 1);
    tassert(-a == -3);

    // Test arithmetic with immediates
    tassert(a + 4 == 7);
    tassert(4 + a == 7);
    tassert(a + 4l == 7);
    tassert(4u + a == 7);

    // Test composite arithmetic
    tassert((a + 1) + (b + 2) == 10);
    tassert((a + d) + (b + 2) == 10);
    tassert((a + d) + (b + e) == 10);
    tassert((3 + d) + (b + e) == 10);
    tassert((3 + d) + (4 + e) == 10);

    tassert(a + (b + c) == 14);
    tassert(3 + (b + c) == 14);
    tassert(3 + (4 + c) == 14);
    tassert(3 + (b + 7) == 14);
    tassert(a + (b + 7) == 14);

    tassert((b + c) + a == 14);
    tassert((b + c) + 3 == 14);
    tassert((4 + c) + 3== 14);
    tassert((b + 7) + 3== 14);
    tassert((b + 7) + a== 14);

    // test combining unary and binary arithmetic
    tassert(-(-a) == 3);
    tassert(-(a - b) == b - a);
    tassert(-((-a) + (-b)) == a + b);

    // Test mixed arithmetic
    mylong al(3l);
    mylong bl(4l);
    mylong cl(7l);

    tassert(al + bl == cl);
    tassert(al + bl == 7l);
    tassert(al + b  == 7l);
    tassert(a + bl  == 7l);
    tassert((a + bl) + cl  == 14l);
    tassert((a + b) + cl  == 14l);
    tassert((al + b) + c  == 14l);
    tassert(cl + (a + bl)  == 14l);
    tassert(cl + (a + b)  == 14l);
    tassert(c + (al + b)  == 14l);
    tassert((a + bl) + (cl + d) == 15l);
    tassert((a + bl) + (c + d) == 15l);

    tassert((a << 2) == 12);
    tassert(((a << 2) >> 2) == a);
}

void
test_conversion()
{
    myint a(4);
    mylong b(4l);
    tassert(typed_equals(a.to<int>(), 4));
    tassert(typed_equals(a.to<mylong>(), b));
    tassert(typed_equals((a + a).to<int>(), 8));
}

void
test_assignment_arith()
{
    myint a(1);
    myint b(2);
    a += b;
    tassert(a == 3 && b == 2);
    a += a*b;
    tassert(a == 9);
    a += 1;
    tassert(a == 10);
}

template<class T, class Expr>
bool is_subexpr(const Expr& e)
{
    return tools::has_subexpr<tools::equal_types_pred<T>, Expr>::val;
}
void
test_tools()
{
    typedef tools::equal_types_pred<myint> intpred;
    myint a(1);
    mylong b(2l);
    tassert(tools::find_subexpr<intpred>(a) == 1);
    tassert(tools::find_subexpr<intpred>(a + b) == 1);
    tassert(tools::find_subexpr_T<mylong>(a + b) == 2l);
    tassert(tools::find_subexpr<intpred>(b + (a + 2) + b) == 1);
    tassert(is_subexpr<myint>(a+b));
    tassert(!is_subexpr<mylong>(a+a));
}

void
test_temporaries()
{
    myint a(1);
#define T2 ((a + a) + (a + a))
#define T3 (T2 + T2)
    tassert(count_temporaries(T2) == 2);
    tassert(count_temporaries(T3) == 3);
    tassert(count_temporaries(T3 + T2) == 3);
    tassert(count_temporaries(T2 + T3) == 3);
}

void
test_references()
{
    mylong a(4);
    mylong_ref ar(a);
    mylong_cref acr(a);
    tassert(a == ar);
    tassert(a == acr);
    tassert(ar == acr);
    tassert(ar == 4l && acr == 4l);

    ar = 5l;
    tassert(a == 5l && acr == 5l && ar == 5l);

    mylong b(6);
    ar = b;
    tassert(a == b);

    mylong_cref bcr(b);
    b = 7l;
    ar = bcr;
    tassert(a == b);

    a = 4l;
    b = 5l;
    tassert(a + bcr == 9l);
    tassert(ar + bcr == 9l);
    a = acr + b;
    tassert(a == 9l);
    ar = acr + bcr;
    tassert(a == 14l);

    a = 4l;
    tassert((a + a) + bcr == 13l);
    tassert((acr + acr) + b == 13l);
    tassert(((a + bcr) + acr + (ar + bcr)) + ((a + a) + (bcr + bcr)) == 40l);
    a = ((a + bcr) + acr + (ar + bcr)) + ((a + a) + (bcr + bcr));
    tassert(a == 40l);
    a = 4l;
    ar = ((a + bcr) + acr + (ar + bcr)) + ((a + a) + (bcr + bcr));
    tassert(a == 40l);
}

struct S { };
S operator+(const S& s, const S&) {return s;}
S operator-(const S&s) {return s;}
template<class T>
S operator*(const S& s, const T&) {return s;}
void
test_unrelated_overload()
{
    // Test a bug in our operator overloading logic which used to not allow
    // this to compile.
    S s;
    s = s + s;
    s = -s;
    s = s * (myint(5) + myint(5));
}

void
test_multiary()
{
    tassert(fourary_test(1, 2, 3, 4) == myint(1 + 2 + 3 + 4));
    tassert(fiveary_test(1, 2, 3, 4, 5u) == 1 + 2 + 3 + 4 + 5);
}

int
main()
{
    std::cout << "expression....";

    test_traits();
    test_initialization();
    test_destruction();
    test_printing();
    test_assignment();
    test_swap();
    test_equals();
    test_arithmetic();
    test_conversion();
    test_assignment_arith();
    test_tools();
    test_temporaries();
    test_references();
    test_multiary();

    test_unrelated_overload();

    std::cout << "PASS" << std::endl;

    // TODO test that certain things *don't* compile?

    return 0;
}
