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

#include "cxx/expression.h"
#include "cxx/tuple.h"

#include "cxx/test/helpers.h"

using namespace flint;
using namespace mp;
using namespace traits;

template<class Operation, class Data>
class my_expression
    : public expression<derived_wrapper<my_expression>, Operation, Data>
{
public:
    my_expression() {};
    template<class T>
    explicit my_expression(const T& t)
        : expression<derived_wrapper< ::my_expression>,
              Operation, Data>(t) {}

    template<class T>
    my_expression& operator=(const T& t)
    {
        this->set(t);
        return *this;
    }

protected:
    explicit my_expression(const Data& d)
        : expression<derived_wrapper< ::my_expression>,
              Operation, Data>(d) {}

    template<class D, class O, class Da>
    friend class flint::expression;
};

struct data
{
    int payload;
    bool* destroyed;
    int extra;

    data() : payload(-1), destroyed(0), extra(42) {}

    ~data()
    {
        if(destroyed)
            *destroyed = true;
    }

    data(const data& d)
        : payload(d.payload), destroyed(0), extra(1) {}
    data(int i)
        : payload(i), destroyed(0), extra(2) {}
    data(char i)
        : payload(i), destroyed(0), extra(3) {}
};

typedef my_expression<operations::immediate, data> myint;

template<class Operation, class Data>
class my_expression2
    : public expression<derived_wrapper<my_expression2>, Operation, Data>
{
public:
    // cannot have a default constructor

    template<class T>
    explicit my_expression2(const T& t)
        : expression<derived_wrapper< ::my_expression2>,
              Operation, Data>(t) {}

    template<class T>
    my_expression2& operator=(const T& t)
    {
        this->set(t);
        return *this;
    }

    my_expression2 create_temporary() const
    {
        return my_expression2(0l);
    }

protected:
    explicit my_expression2(const Data& d)
        : expression<derived_wrapper< ::my_expression2>,
              Operation, Data>(d) {}

    template<class D, class O, class Da>
    friend class flint::expression;
};
struct long_data
{
    long payload;
    // no default constructor
    long_data(long d) : payload(d) {}
    long_data(const myint& m) : payload(m._data().payload) {}
};
typedef my_expression2<operations::immediate, long_data> mylong;

namespace flint {
namespace rules {

template<>
struct print<myint>
{
    static void doit(const myint& v, std::ostream& o)
    {
        o << v._data().payload;
    }
};

template<>
struct assignment<myint, myint>
{
    static void doit(myint& to, const myint& from)
    {
        to._data().payload = from._data().payload;
        to._data().extra = 4;
    }
};

template<>
struct assignment<myint, long>
{
    static void doit(myint& to, long from)
    {
        to._data().payload = from;
        to._data().extra = 5;
    }
};

template<>
struct equals<myint, myint>
{
    static bool get(const myint& i1, const myint& i2)
    {
        return i1._data().payload == i2._data().payload;
    }
};

template<>
struct equals<myint, int>
{
    static bool get(const myint& i1, int i2)
    {
        return i1._data().payload == i2;
    }
};

template<>
struct conversion<int, myint>
{
    static int get(const myint& from)
    {
        return from._data().payload;
    }
};

template<>
struct commutative_binary_expression<myint, operations::plus, myint>
{
    typedef myint return_t;
    static void doit(myint& to, const myint& a1, const myint& a2)
    {
        to._data().payload = a1._data().payload + a2._data().payload;
    }
};

template<class T>
struct commutative_binary_expression<myint,
    typename mp::enable_if<traits::is_integer<T>, operations::plus>::type,
    T>
{
    typedef myint return_t;
    static void doit(myint& to, const myint& a1, T a2)
    {
        to._data().payload = a1._data().payload + a2;
    }
};

template<>
struct commutative_binary_expression<myint, operations::times, myint>
{
    typedef myint return_t;
    static void doit(myint& to, const myint& a1, const myint& a2)
    {
        to._data().payload = a1._data().payload * a2._data().payload;
    }
};

template<>
struct binary_expression<myint, operations::minus, myint>
{
    typedef myint return_t;
    static void doit(myint& to, const myint& a1, const myint& a2)
    {
        to._data().payload = a1._data().payload - a2._data().payload;
    }
};

template<>
struct binary_expression<myint, operations::divided_by, myint>
{
    typedef myint return_t;
    static void doit(myint& to, const myint& a1, const myint& a2)
    {
        to._data().payload = a1._data().payload / a2._data().payload;
    }
};

template<>
struct binary_expression<myint, operations::modulo, myint>
{
    typedef myint return_t;
    static void doit(myint& to, const myint& a1, const myint& a2)
    {
        to._data().payload = a1._data().payload % a2._data().payload;
    }
};

template<>
struct unary_expression<operations::negate, myint>
{
    typedef myint return_t;
    static void doit(myint& to, const myint& from)
    {
        to._data().payload = - from._data().payload;
    }
};


/////////////////////////////////////////////////////////////////////////////
// Minimal rules for mylong
/////////////////////////////////////////////////////////////////////////////

template<>
struct equals<mylong, mylong>
{
    static bool get(const mylong& i1, const mylong& i2)
    {
        return i1._data().payload == i2._data().payload;
    }
};

template<>
struct equals<mylong, long>
{
    static bool get(const mylong& i1, long i2)
    {
        return i1._data().payload == i2;
    }
};

#if 0
template<bool c, class Op, class Data>
struct evaluation<
    my_expression<
        operations::plus,
        make_tuple<const mylong&, const mylong&>::type>,
    Op, Data,
    c, 0>
{
    typedef mylong return_t;
    typedef empty_tuple temporaries_t;
    typedef my_expression<
        operations::plus,
        make_tuple<const mylong&, const mylong&>::type> expr_t;
    static void doit(const expr_t& input, temporaries_t temps, return_t* output)
    {
        output->_data().payload = input._data().first()._data().payload
                                + input._data().second()._data().payload;
    }
};
#endif

template<>
struct commutative_binary_expression<mylong, operations::plus, mylong>
{
    typedef mylong return_t;
    static void doit(mylong& to, const mylong& a1, const mylong& a2)
    {
        to._data().payload = a1._data().payload + a2._data().payload;
    }
};

template<>
struct commutative_binary_expression<myint, operations::plus, mylong>
{
    typedef mylong return_t;
    static void doit(mylong& to, const myint& a1, const mylong& a2)
    {
        to._data().payload = a1._data().payload + a2._data().payload;
    }
};
} // rules
} // flint

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
}

template<class T, class U>
bool typed_equals(const T&, const U&)
{
    return false;
}

template<class T>
bool typed_equals(const T& a, const T& b)
{
    return a == b;
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
    return tools::is_super_sub_expr<Expr, T>::val;
}
void
test_tools()
{
    myint a(1);
    mylong b(2l);
    tassert(tools::find_subexpr<myint>(a) == 1);
    tassert(tools::find_subexpr<myint>(a + b) == 1);
    tassert(tools::find_subexpr<mylong>(a + b) == 2l);
    tassert(tools::find_subexpr<myint>(b + (a + 2) + b) == 1);
    tassert(is_subexpr<myint>(a+b));
    tassert(!is_subexpr<mylong>(a+a));
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
    test_equals();
    test_arithmetic();
    test_conversion();
    test_assignment_arith();
    test_tools();

    std::cout << "PASS" << std::endl;

    // TODO test that certain things *don't* compile?

    return 0;
}
