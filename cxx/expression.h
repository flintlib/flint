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

#ifndef CXX_EXPRESSION_H
#define CXX_EXPRESSION_H

// TODO
// * static asserts

#include <iosfwd>
#include <string>

#include "cxx/evaluation_tools.h"
#include "cxx/expression_traits.h"
#include "cxx/mp.h"
#include "cxx/rules.h"
#include "cxx/traits.h"
#include "cxx/tuple.h"

namespace flint {
namespace detail {
template<class Operation, class Expr, class Data>
struct evaluation_traits
{
    typedef typename Expr::derived_t derived_t;
    typedef typename mp::find_evaluation<
        Operation, Data, false>::type rule_t;
    typedef typename mp::find_evaluation<
        Operation, Data, true>::type temp_rule_t;
    typedef typename rule_t::return_t evaluation_return_t;
    typedef evaluation_return_t evaluated_t;

    static evaluation_return_t evaluate(const derived_t& from)
    {
        evaluation_return_t res;
        evaluate_into_fresh(res, from);
        return res;
    }

    static void evaluate_into(evaluation_return_t& to, const derived_t& from)
    {
        typedef mp::back_tuple<typename rule_t::temporaries_t> back_t;
        typename back_t::type temps_backing;
        typename rule_t::temporaries_t temps;
        back_t::init(temps, temps_backing, 0);
        rule_t::doit(from._data(), temps, &to);
    }

    static void evaluate_into_fresh(evaluation_return_t& to, const derived_t& from)
    {
        typedef mp::back_tuple<
            typename temp_rule_t::temporaries_t,
            evaluation_return_t
          > back_t;
        typename back_t::type temps_backing;
        typename temp_rule_t::temporaries_t temps;
        back_t::init(temps, temps_backing, &to);
        temp_rule_t::doit(from._data(), temps, &to);
    }
};

template<class Expr, class Data>
struct evaluation_traits<operations::immediate, Expr, Data>
{
    typedef typename Expr::derived_t derived_t;
    typedef typename Expr::derived_t evaluated_t;
    typedef evaluated_t& evaluation_return_t;

    static evaluated_t& evaluate(derived_t& d) {return d;}
    static const evaluated_t& evaluate(const derived_t& d) {return d;}

    static void evaluate_into(derived_t& to, const derived_t& from)
    {
        rules::assignment<derived_t, derived_t>::doit(to, from);
    }

    static void evaluate_into_fresh(derived_t& to, const derived_t& from)
    {
        evaluate_into(to, from);
    }
};

// See copy constructor of expression for explanation
template<class Data, class Enable = void>
struct copy_initialization_helper
{
    Data data;

    copy_initialization_helper(const Data& d) : data(d) {}
    copy_initialization_helper() {}

    // This does nothing, by design!
    copy_initialization_helper(const copy_initialization_helper& o) {}

    template<class T, class U>
    void doit(T& to, const U& from)
    {
        rules::initialization<T, U>::doit(to, from);
    }
};

template<class Data>
struct copy_initialization_helper<Data,
    typename mp::enable_if<traits::is_tuple<Data> >::type>
{
    Data data;

    copy_initialization_helper(const Data& d)
        : data(d) {}

    copy_initialization_helper(const copy_initialization_helper& o)
        : data(o.data)
    {
    }

    template<class T, class U>
    void doit(T& to, const U& from)
    {
    }
};
} // detail

template<class Derived, class Operation, class Data>
class expression
    : public detail::EXPRESSION
{
private:
    detail::copy_initialization_helper<Data> data;

protected:
    explicit expression(const Data& d) : data(d) {}

public:
    typedef detail::evaluation_traits<Operation, expression, Data> ev_traits_t;
    typedef typename Derived::template type<Operation, Data>::result derived_t;
    typedef typename ev_traits_t::evaluated_t evaluated_t;
    typedef typename ev_traits_t::evaluation_return_t evaluation_return_t;
    typedef Data data_t;
    typedef Operation operation_t;

private:
    derived_t& downcast() {return *static_cast<derived_t*>(this);}
    const derived_t& downcast() const
    {
        return *static_cast<const derived_t*>(this);
    }

public:
    // TODO strip qualifiers?
    template<class T>
    explicit expression(const T& t,
            typename mp::enable_if<traits::is_implemented<
                rules::initialization<derived_t, T> > >::type* = 0)
    {
        rules::initialization<derived_t, T>::doit(downcast(), t);
    }

    template<class T>
    explicit expression(const T& t,
            typename mp::disable_if<traits::is_implemented<
                rules::initialization<derived_t, T> > >::type* = 0,
            typename mp::enable_if<traits::is_expression<T> >::type* = 0)
    {
        rules::empty_initialization<derived_t>::doit(downcast());
        T::ev_traits_t::evaluate_into_fresh(downcast(), t);
    }

    // We pay here for using the same class for actual data wrappers,
    // and expression templates. The problem is that Data may not be
    // default-initialisable (contains reference types), but it may also not
    // be copy-initialisable (arrays...).
    // The solution here is to wrap Data into copy_initialization_helper.
    // If Data is a reference or a tuple (the expression template case),
    // then data.doit is no-op the copy constructor of data does the work.
    // If Data is not a reference or tuple, then the copy constructor is no-op,
    // and doit defrers to rules.
    expression(const expression& e)
        : data(e.data)
    {
        data.doit(downcast(), e.downcast());
    }

    // NB: the compiler is not allowed to eagerly instantiate!
    expression() {rules::empty_initialization<derived_t>::doit(downcast());}

    ~expression() {rules::destruction<derived_t>::doit(downcast());}

    expression& operator=(const expression& o)
    {
        this->set(o.downcast());
        return *this;
    }

    Data& _data() {return data.data;}
    const Data& _data() const {return data.data;}

    void print(std::ostream& o) const
    {
        rules::print<evaluated_t>::doit(evaluate(), o);
    }

    std::string to_string(int base = 10) const
    {
        return rules::to_string<evaluated_t>::get(evaluate(), base);
    }

    template<class T>
    T to() const
    {
        return rules::conversion<T, evaluated_t>::get(evaluate());
    }

    typename traits::make_const<evaluation_return_t>::type evaluate() const
    {
        return ev_traits_t::evaluate(downcast());
    }
    evaluation_return_t evaluate() {return ev_traits_t::evaluate(downcast());}

    template<class T>
    void set(const T& t,
            typename mp::enable_if<traits::is_expression<T> >::type* = 0)
    {
        T::ev_traits_t::evaluate_into(downcast(), t);
    }
    template<class T>
    void set(const T& t,
            typename mp::disable_if<traits::is_expression<T> >::type* = 0)
    {
        rules::assignment<derived_t, T>::doit(downcast(), t);
    }

    template<class T>
    bool equals(const T& t,
            typename mp::enable_if<traits::is_lazy_expr<T> >::type* = 0) const
    {
        return equals(t.evaluate());
    }
    template<class T>
    bool equals(const T& t,
            typename mp::disable_if<traits::is_lazy_expr<T> >::type* = 0) const
    {
        return rules::equals<evaluated_t, T>::get(evaluate(), t);
    }

    template<class Op, class NData>
    struct make_helper
    {
        typedef typename Derived::template type<Op, NData>::result type;
        static type make(const NData& ndata)
        {
            return type(ndata);
        }
    };
};

template<template<class O, class D> class Derived>
struct derived_wrapper
{
    template<class Operation, class Data>
    struct type
    {
        typedef Derived<Operation, Data> result;
    };
};


// operators

namespace detail {
template<class Expr>
struct storage_traits
    : mp::if_<
          traits::is_immediate<Expr>,
          typename traits::forwarding<Expr>::type,
          Expr
        > { };

template<class Expr1, class Op, class Expr2>
struct binary_op_helper
{
    typedef mp::make_tuple<
        typename storage_traits<Expr1>::type,
        typename storage_traits<Expr2>::type
      > maker;
    typedef typename maker::type type;
    typedef mp::find_evaluation<Op, type, true> ev_t;
    typedef typename ev_t::type::return_t Expr;
    typedef typename Expr::template make_helper<Op, type> make_helper;
    typedef typename make_helper::type return_t;

    typedef mp::enable_if<traits::is_implemented<ev_t>, return_t> enable;

    static return_t make(const Expr1& left, const Expr2& right)
    {
        return make_helper::make(maker::make(left, right));
    }
};

template<class Op, class Expr>
struct unary_op_helper
{
    typedef tuple<typename storage_traits<Expr>::type, empty_tuple> type;
    typedef mp::find_evaluation<Op, type, true> ev_t;
    typedef typename ev_t::type::return_t Rexpr;
    typedef typename Rexpr::template make_helper<Op, type> make_helper;
    typedef typename make_helper::type return_t;

    typedef mp::enable_if<traits::is_implemented<ev_t>, return_t> enable;

    static return_t make(const Expr& e)
    {
        return make_helper::make(type(e, empty_tuple()));
    }
};

template<class Expr1, class Expr2>
struct order_op_helper
{
    typedef typename tools::evaluation_helper<Expr1>::type ev1_t;
    typedef typename tools::evaluation_helper<Expr2>::type ev2_t;
    typedef tools::symmetric_cmp<ev1_t, ev2_t> scmp;

    typedef mp::enable_if<
          mp::and_<
            traits::is_implemented<scmp>,
            mp::or_<
                traits::is_expression<Expr1>,
                traits::is_expression<Expr2>
              >
          >,
          bool> enable;

    static int get(const Expr1& e1, const Expr2& e2)
    {
        return scmp::get(tools::evaluation_helper<Expr1>::get(e1),
            tools::evaluation_helper<Expr2>::get(e2));
    }
};
}

template<class Expr>
inline typename mp::enable_if<traits::is_expression<Expr>, std::ostream&>::type
operator<<(std::ostream& o, const Expr& e)
{
    e.print(o);
    return o;
}

template<class Expr1, class Expr2>
inline typename mp::enable_if<traits::is_expression<Expr1>, bool>::type
operator==(const Expr1& e1, const Expr2& e2)
{
  return e1.equals(e2);
}

template<class Expr1, class Expr2>
inline typename mp::enable_if<mp::and_<
        mp::not_<traits::is_expression<Expr1> >,
        traits::is_expression<Expr2> >,
    bool>::type
operator==(const Expr1& e1, const Expr2& e2)
{
  return e2.equals(e1);
}

template<class Expr1, class Expr2>
inline typename mp::enable_if<mp::or_<
        traits::is_expression<Expr1>,
        traits::is_expression<Expr2> >,
    bool>::type
operator!=(const Expr1& e1, const Expr2& e2)
{
  return !(e1 == e2);
}

template<class Expr1, class Expr2>
inline typename detail::order_op_helper<Expr1, Expr2>::enable::type
operator<(const Expr1& e1, const Expr2& e2)
{
    return detail::order_op_helper<Expr1, Expr2>::get(e1, e2) < 0;
}

template<class Expr1, class Expr2>
inline typename detail::order_op_helper<Expr1, Expr2>::enable::type
operator<=(const Expr1& e1, const Expr2& e2)
{
    return detail::order_op_helper<Expr1, Expr2>::get(e1, e2) <= 0;
}

template<class Expr1, class Expr2>
inline typename detail::order_op_helper<Expr1, Expr2>::enable::type
operator>(const Expr1& e1, const Expr2& e2)
{
    return detail::order_op_helper<Expr1, Expr2>::get(e1, e2) > 0;
}

template<class Expr1, class Expr2>
inline typename detail::order_op_helper<Expr1, Expr2>::enable::type
operator>=(const Expr1& e1, const Expr2& e2)
{
    return detail::order_op_helper<Expr1, Expr2>::get(e1, e2) >= 0;
}

template<class Expr1, class Expr2>
inline typename detail::binary_op_helper<
    Expr1, operations::plus, Expr2>::enable::type
operator+(const Expr1& e1, const Expr2& e2)
{
    return detail::binary_op_helper<Expr1, operations::plus, Expr2>::make(e1, e2);
}

template<class Expr1, class Expr2>
inline typename detail::binary_op_helper<
    Expr1, operations::minus, Expr2>::enable::type
operator-(const Expr1& e1, const Expr2& e2)
{
    return detail::binary_op_helper<Expr1, operations::minus, Expr2>::make(e1, e2);
}

template<class Expr1, class Expr2>
inline typename detail::binary_op_helper<
    Expr1, operations::times, Expr2>::enable::type
operator*(const Expr1& e1, const Expr2& e2)
{
    return detail::binary_op_helper<Expr1, operations::times, Expr2>::make(e1, e2);
}

template<class Expr1, class Expr2>
inline typename detail::binary_op_helper<
    Expr1, operations::divided_by, Expr2>::enable::type
operator/(const Expr1& e1, const Expr2& e2)
{
    return detail::binary_op_helper<Expr1, operations::divided_by, Expr2>::make(e1, e2);
}

template<class Expr1, class Expr2>
inline typename detail::binary_op_helper<
    Expr1, operations::modulo, Expr2>::enable::type
operator%(const Expr1& e1, const Expr2& e2)
{
    return detail::binary_op_helper<Expr1, operations::modulo, Expr2>::make(e1, e2);
}

template<class Expr>
inline typename detail::unary_op_helper<operations::negate, Expr>::enable::type
operator-(const Expr& e)
{
    return detail::unary_op_helper<operations::negate, Expr>::make(e);
}

template<class Expr1, class Expr2>
inline typename mp::enable_if<traits::is_immediate_expr<Expr1>, Expr1&>::type
operator+=(Expr1& e1, const Expr2& e2)
{
    e1.set(e1 + e2);
    return e1;
}

template<class Expr1, class Expr2>
inline typename mp::enable_if<traits::is_immediate_expr<Expr1>, Expr1&>::type
operator-=(Expr1& e1, const Expr2& e2)
{
    e1.set(e1 - e2);
    return e1;
}

template<class Expr1, class Expr2>
inline typename mp::enable_if<traits::is_immediate_expr<Expr1>, Expr1&>::type
operator*=(Expr1& e1, const Expr2& e2)
{
    e1.set(e1 * e2);
    return e1;
}

template<class Expr1, class Expr2>
inline typename mp::enable_if<traits::is_immediate_expr<Expr1>, Expr1&>::type
operator/=(Expr1& e1, const Expr2& e2)
{
    e1.set(e1 / e2);
    return e1;
}

template<class Expr1, class Expr2>
inline typename mp::enable_if<traits::is_immediate_expr<Expr1>, Expr1&>::type
operator%=(Expr1& e1, const Expr2& e2)
{
    e1.set(e1 % e2);
    return e1;
}
}

// TODO remove this?
#include "cxx/default_rules.h"


////////////////////////////////////////////////////////////////////////
// HELPER MACROS
////////////////////////////////////////////////////////////////////////

// This set of macros should be called in namespace flint.

// Introduce a new binary operation called "name"
// NB: because of ADL bugs in g++ <= 4.4, the operation tag is called "name_op",
// whereas the function corresponding to it is just called "name"
#define FLINT_DEFINE_BINOP(name) \
namespace operations { \
struct name##_op { }; \
} \
template<class T1, class T2> \
inline typename detail::binary_op_helper<T1, operations::name##_op, T2>::enable::type \
name(const T1& t1, const T2& t2) \
{ \
    return detail::binary_op_helper<T1, operations::name##_op, T2>::make(t1, t2); \
}

// Introduce a new unary operation called "name"
#define FLINT_DEFINE_UNOP(name) \
namespace operations { \
struct name##_op { }; \
} \
template<class T1> \
inline typename detail::unary_op_helper<operations::name##_op, T1>::enable::type \
name(const T1& t1) \
{ \
    return detail::unary_op_helper<operations::name##_op, T1>::make(t1); \
}

// To be called in any namespace

// Make the binary operation "name" available in current namespace
#define FLINT_DEFINE_BINOP_HERE(name) \
template<class T1, class T2> \
inline typename ::flint::detail::binary_op_helper<\
    T1, ::flint::operations::name##_op, T2>::enable::type \
name(const T1& t1, const T2& t2) \
{ \
  return ::flint::detail::binary_op_helper< \
      T1, ::flint::operations::name##_op, T2>::make(t1, t2); \
}

// Make the unary operation "name" available in current namespace
#define FLINT_DEFINE_UNOP_HERE(name) \
template<class T1> \
inline typename ::flint::detail::unary_op_helper<\
    ::flint::operations::name##_op, T1>::enable::type \
name(const T1& t1) \
{ \
  return ::flint::detail::unary_op_helper< ::flint::operations::name##_op, T1>::make(t1); \
}

#endif
