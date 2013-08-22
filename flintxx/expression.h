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

#include "evaluation_tools.h"
#include "expression_traits.h"
#include "mp.h"
#include "rules.h"
#include "traits.h"
#include "tuple.h"

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
        evaluated_t res = 
            rules::instantiate_temporaries<derived_t, evaluated_t>::get(from);
        evaluate_into_fresh(res, from);
        return res;
    }

    template<class T>
    static void evaluate_into(T& to, const derived_t& from)
    {
        typedef mp::back_tuple<typename rule_t::temporaries_t> back_t;
        typename back_t::type temps_backing =
            mp::htuples::fill<typename back_t::type>(
                tools::temporaries_filler(from));
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
        typename back_t::type temps_backing =
            mp::htuples::fill<typename back_t::type>(
                tools::temporaries_filler(from));
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

    template<class T>
    static void evaluate_into(T& to, const derived_t& from)
    {
        rules::assignment<T, derived_t>::doit(to, from);
    }

    static void evaluate_into_fresh(derived_t& to, const derived_t& from)
    {
        evaluate_into(to, from);
    }
};
} // detail

// Note that, wihle Data does not have to be default constructible,
// it *does* need to be copy-constructible, and have a working destructor.
template<class Derived, class Operation, class Data>
class expression
{
private:
    Data data;

protected:
    explicit expression(const Data& d) : data(d) {}

public:
    // internal -- see is_expression implementation.
    typedef void IS_EXPRESSION_MARKER;

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

    // Some helpers for initialization, since it is not possible to
    // conditionally enable constructors in C++98
    template<class T>
    static data_t get_data(const T& t,
        typename mp::disable_if<traits::is_lazy_expr<T> >::type* = 0)
    {
        return data_t(t);
    }
    template<class T>
    static data_t get_data(T& t,
        typename mp::disable_if<traits::is_lazy_expr<T> >::type* = 0)
    {
        return data_t(t);
    }
    template<class T>
    static data_t get_data(const T& t,
        typename mp::enable_if<traits::is_lazy_expr<T> >::type* = 0,
        typename mp::disable_if<
            mp::equal_types<typename T::evaluated_t, derived_t> >::type* = 0)
    {
        return data_t(t.evaluate());
    }
    template<class T>
    static data_t get_data(const T& t,
        typename mp::enable_if<traits::is_lazy_expr<T> >::type* = 0,
        typename mp::enable_if<
            mp::equal_types<typename T::evaluated_t, derived_t> >::type* = 0)
    {
        return data_t(t.evaluate()._data());
    }

    // Invoke the data copy constructor when appropriate
    static data_t get_data(const derived_t& o)
    {
        return data_t(o._data());
    }
    static data_t get_data(derived_t& o)
    {
        return data_t(o._data());
    }

    // Having the empty constructor here delays its instantiation, and allows
    // compiling even if data is *not* default constructible.
    static data_t get_data() {return data_t();}

public:
    template<class T>
    explicit expression(const T& t)
        : data(get_data(t)) {}

    template<class T>
    explicit expression(T& t)
        : data(get_data(t)) {}

    template<class T, class U>
    expression(const T& t, const U& u)
        : data(t, u) {}

    template<class T, class U, class V>
    expression(const T& t, const U& u, const V& v)
        : data(t, u, v) {}

    expression() : data(get_data()) {}

    expression& operator=(const expression& o)
    {
        this->set(o.downcast());
        return *this;
    }

    // See rules::instantiate_temporaries for explanation.
    evaluated_t create_temporary() const
    {
        return evaluated_t();
    }

    Data& _data() {return data;}
    const Data& _data() const {return data;}

    void print(std::ostream& o) const
    {
        tools::print_using_str<evaluated_t>::doit(evaluate(), o);
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
            typename mp::enable_if<traits::is_expression<T> >::type* = 0,
            typename mp::enable_if<traits::can_evaluate_into<
                derived_t, typename T::evaluated_t> >::type* = 0)
    {
        T::ev_traits_t::evaluate_into(downcast(), t);
    }
    template<class T>
    void set(const T& t,
            typename mp::enable_if<traits::is_expression<T> >::type* = 0,
            typename mp::disable_if<traits::can_evaluate_into<
                derived_t, typename T::evaluated_t> >::type* = 0)
    {
        rules::assignment<derived_t, typename T::evaluated_t>::doit(
            downcast(), t.evaluate());
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
        return tools::equals_using_cmp<evaluated_t, T>::get(evaluate(), t);
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

template<template<class E, class O, class D> class Derived, class Extra>
struct derived_wrapper2
{
    template<class Operation, class Data>
    struct type
    {
        typedef Derived<Extra, Operation, Data> result;
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
template<>
struct storage_traits<detail::UNUSED> {typedef detail::UNUSED type;};

template<class ev_t, class Op, class type>
struct nary_op_helper_step2
{
    typedef typename ev_t::return_t Expr;
    typedef typename Expr::template make_helper<Op, type> make_helper;
    typedef typename make_helper::type return_t;
};
template<class Op, class type>
struct nary_op_helper_step2<rules::UNIMPLEMENTED, Op, type>
{
    struct return_t { };
    struct make_helper { };
};

template<class Op, class Data>
struct nary_op_helper
{
    typedef typename mp::find_evaluation<Op, Data, true>::type ev_t;
    typedef nary_op_helper_step2<ev_t, Op, Data> nohs2;
    typedef typename nohs2::return_t return_t;
    typedef typename nohs2::make_helper make_helper;

    typedef traits::is_implemented<ev_t> cond;
    typedef mp::enable_if<cond, return_t> enable;
};

template<class Op, class Maker>
struct nary_op_helper_maker
    : nary_op_helper<Op, typename Maker::type>
{
    typedef Maker maker;
};

#define FLINTXX_NARY_OP_HELPER_MACRO(arg) typename storage_traits< arg >::type
template<class Op, FLINTXX_MAKE_TUPLE_TEMPLATE_ARGS>
struct nary_op_helper2
    : nary_op_helper_maker<Op, mp::make_tuple<
          FLINTXX_MAKE_TUPLE_TYPES_APPLYMACRO(FLINTXX_NARY_OP_HELPER_MACRO) > >
{
    typedef nary_op_helper2 noh2;
    static typename noh2::return_t make(FLINTXX_MAKE_TUPLE_FUNC_ARGS)
    {
        return noh2::make_helper::make(noh2::maker::make(
              FLINTXX_MAKE_TUPLE_FUNC_ARG_NAMES));
    }
};

template<class Expr1, class Op, class Expr2>
struct binary_op_helper
    : nary_op_helper2<Op, Expr1, Expr2>
{ };

template<class Op, class Expr>
struct unary_op_helper : nary_op_helper2<Op, Expr> { };

template<class Ret, class Op, class Expr>
struct unary_op_helper_with_rettype
{
    typedef mp::make_tuple<typename storage_traits<Expr>::type> maker;
    typedef typename Ret::template make_helper<
        Op, typename maker::type>::type return_t;
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

template<class Expr1, class Expr2>
inline typename detail::binary_op_helper<
    Expr1, operations::binary_and, Expr2>::enable::type
operator&(const Expr1& e1, const Expr2& e2)
{
    return detail::binary_op_helper<Expr1, operations::binary_and, Expr2>::make(e1, e2);
}

template<class Expr1, class Expr2>
inline typename detail::binary_op_helper<
    Expr1, operations::binary_or, Expr2>::enable::type
operator|(const Expr1& e1, const Expr2& e2)
{
    return detail::binary_op_helper<Expr1, operations::binary_or, Expr2>::make(e1, e2);
}

template<class Expr1, class Expr2>
inline typename detail::binary_op_helper<
    Expr1, operations::binary_xor, Expr2>::enable::type
operator^(const Expr1& e1, const Expr2& e2)
{
    return detail::binary_op_helper<Expr1, operations::binary_xor, Expr2>::make(e1, e2);
}

template<class Expr1, class Expr2>
inline typename detail::binary_op_helper<
    Expr1, operations::shift, Expr2>::enable::type
operator<<(const Expr1& e1, const Expr2& e2)
{
    return detail::binary_op_helper<Expr1, operations::shift, Expr2>::make(e1, e2);
}

template<class Expr1, class Expr2>
inline typename detail::binary_op_helper<
    Expr1, operations::shift, Expr2>::enable::type
operator>>(const Expr1& e1, const Expr2& e2)
{
    return detail::binary_op_helper<Expr1, operations::shift, Expr2>::make(e1, -e2);
}

template<class Expr>
inline typename detail::unary_op_helper<operations::negate, Expr>::enable::type
operator-(const Expr& e)
{
    return detail::unary_op_helper<operations::negate, Expr>::make(e);
}

template<class Expr>
inline typename detail::unary_op_helper<operations::complement, Expr>::enable::type
operator~(const Expr& e)
{
    return detail::unary_op_helper<operations::complement, Expr>::make(e);
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

template<class Expr1, class Expr2>
inline typename mp::enable_if<traits::is_immediate_expr<Expr1>, Expr1&>::type
operator|=(Expr1& e1, const Expr2& e2)
{
    e1.set(e1 | e2);
    return e1;
}

template<class Expr1, class Expr2>
inline typename mp::enable_if<traits::is_immediate_expr<Expr1>, Expr1&>::type
operator&=(Expr1& e1, const Expr2& e2)
{
    e1.set(e1 & e2);
    return e1;
}

template<class Expr1, class Expr2>
inline typename mp::enable_if<traits::is_immediate_expr<Expr1>, Expr1&>::type
operator^=(Expr1& e1, const Expr2& e2)
{
    e1.set(e1 ^ e2);
    return e1;
}

// TODO move to std?
template<class Expr1, class Expr2>
inline typename mp::enable_if<typename traits::is_implemented<
    rules::swap<Expr1, Expr2> > >::type swap(Expr1& e1, Expr2& e2)
{
    rules::swap<Expr1, Expr2>::doit(e1, e2);
}
}

// TODO remove this?
#include "default_rules.h"


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

#define FLINT_BINOP_ENABLE_RETTYPE(name, T1, T2) \
    typename detail::binary_op_helper<T1, operations::name##_op, T2>::enable::type

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

#define FLINT_UNOP_ENABLE_RETTYPE(name, T) \
    typename detail::unary_op_helper<operations::name##_op, T>::return_t
#define FLINT_UNOP_BUILD_RETTYPE(name, rettype, T) \
    typename detail::unary_op_helper_with_rettype<rettype, \
        operations::name##_op, T>::return_t

#define FLINT_DEFINE_THREEARY(name) \
namespace operations { \
struct name##_op { }; \
} \
template<class T1, class T2, class T3> \
inline typename detail::nary_op_helper2<operations::name##_op, T1, T2, T3>::enable::type \
name(const T1& t1, const T2& t2, const T3& t3) \
{ \
    return detail::nary_op_helper2<operations::name##_op, T1, T2, T3>::make( \
        t1, t2, t3); \
}

#define FLINT_THREEARY_ENABLE_RETTYPE(name, T1, T2, T3) \
    typename detail::nary_op_helper2<operations::name##_op, T1, T2, T3>::enable::type

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
