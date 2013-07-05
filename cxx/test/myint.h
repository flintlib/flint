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

#ifndef CXX_TEST_MYINT_H
#define CXX_TEST_MYINT_H CXX_TEST_MYINT_H

#include "cxx/expression.h"

namespace flint {
template<class Operation, class Data>
class my_expression
    : public expression<derived_wrapper<my_expression>, Operation, Data>
{
public:
    my_expression() {};
    template<class T>
    explicit my_expression(const T& t)
        : expression<derived_wrapper< ::flint::my_expression>,
              Operation, Data>(t) {}

    template<class T>
    my_expression& operator=(const T& t)
    {
        this->set(t);
        return *this;
    }

protected:
    explicit my_expression(const Data& d)
        : expression<derived_wrapper< ::flint::my_expression>,
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
        : expression<derived_wrapper< ::flint::my_expression2>,
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
        : expression<derived_wrapper< ::flint::my_expression2>,
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

#endif
