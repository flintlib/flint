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
    explicit my_expression(const T& t) : my_expression::expression(t) {}

    template<class T>
    my_expression& operator=(const T& t)
    {
        this->set(t);
        return *this;
    }

protected:
    explicit my_expression(const Data& d) : my_expression::expression(d) {}

    template<class D, class O, class Da>
    friend class flint::expression;
};

struct data
{
    int payload;
    bool* destroyed;
    int extra;

    data() : payload(-1), destroyed(0), extra(0) {}

    ~data()
    {
        if(destroyed)
            *destroyed = true;
    }

    data(const data& d)
    {
        payload = d.payload;
        destroyed = 0;
        extra = d.extra | 1 << 30;
    }
};

typedef my_expression<operations::immediate, data> myint;

namespace flint {
namespace rules {
template<>
struct empty_initialization<myint>
{
    static void doit(myint& v)
    {
        v._data().extra = 42;
    }
};

template<>
struct initialization<myint, myint>
{
    static void doit(myint& to, const myint& from)
    {
        to._data().payload = from._data().payload;
        to._data().extra = 1;
    }
};

template<>
struct initialization<myint, int>
{
    static void doit(myint& to, int i)
    {
        to._data().payload = i;
        to._data().extra = 2;
    }
};

template<>
struct initialization<myint, char>
{
    static void doit(myint& to, char i)
    {
        to._data().payload = i;
        to._data().extra = 3;
    }
};

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

// TODO
template<bool c>
struct evaluation<
    my_expression<
        operations::plus,
        make_tuple<const myint&, const myint&>::type>,
    c, 0>
{
    typedef myint return_t;
    typedef empty_tuple temporaries_t;
    typedef my_expression<
        operations::plus,
        make_tuple<const myint&, const myint&>::type> expr_t;
    static void doit(const expr_t& input, temporaries_t temps, return_t* output)
    {
        output->_data().payload = input._data().first()._data().payload
                                + input._data().second()._data().payload;
    }
};

template<class T, bool c>
struct evaluation<
    my_expression<
        operations::plus,
        tuple<const myint&, tuple<T, empty_tuple> > >,
    c, 0,
    typename mp::enable_if<traits::is_integer<T> >::type>
{
    typedef myint return_t;
    typedef empty_tuple temporaries_t;
    typedef my_expression<
        operations::plus,
        typename make_tuple<const myint&, T>::type> expr_t;
    static void doit(const expr_t& input, temporaries_t temps, return_t* output)
    {
        output->_data().payload = input._data().first()._data().payload
                                + input._data().second();
    }
};

template<class T, bool c>
struct evaluation<
    my_expression<
        operations::plus,
        tuple<T, tuple<const myint&, empty_tuple> > >,
    c, 0,
    typename mp::enable_if<traits::is_integer<T> >::type>
{
    typedef myint return_t;
    typedef empty_tuple temporaries_t;
    typedef my_expression<
        operations::plus,
        typename make_tuple<T, const myint&>::type> expr_t;
    static void doit(const expr_t& input, temporaries_t temps, return_t* output)
    {
        output->_data().payload = input._data().second()._data().payload
                                + input._data().first();
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
    tassert(b._data().payload == 17 && b._data().extra == (1 | 1 << 30));

    myint c(8);
    tassert(c._data().payload == 8 && c._data().extra == 2);

    myint d('a');
    tassert(d._data().payload == 'a' && d._data().extra == 3);
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
    // tassert(a._data().extra == 4); // this won't happen because of copy init

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
}

void
test_arithmetic()
{
    myint a(3);
    myint b(4);
    myint c(7);

    tassert((a + b).evaluate() == 7);
    tassert(a + b == c);
    tassert(a + b == 7);

    tassert(a + 4 == 7);
    tassert(4 + a == 7);
    tassert(a + 4l == 7);
    tassert(4u + a == 7);
}

int
main()
{
    test_traits();
    test_initialization();
    test_destruction();
    test_printing();
    test_assignment();
    test_equals();
    test_arithmetic();

    // TODO test that certain things *don't* compile?

    return 0;
}
