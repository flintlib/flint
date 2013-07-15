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

// Helpers to define concrete subclasses of expression.
// Contrary to other parts of this library, they are tailored very
// specifically towards FLINT.

#include "cxx/mp.h"
#include "cxx/expression_traits.h"

#ifndef CXX_FLINT_CLASSES_H
#define CXX_FLINT_CLASSES_H

namespace flint {
namespace flint_classes {
template<class Wrapped, class Inner>
struct ref_data
{
    typedef void IS_REF_OR_CREF;
    typedef Wrapped wrapped_t;

    typedef Inner* data_ref_t;
    typedef const Inner* data_srcref_t;

    Inner* inner;

    ref_data(Wrapped& o) : inner(o._data().inner) {}

    static ref_data make(Inner* f) {return ref_data(f);}

private:
    ref_data(Inner* fp) : inner(fp) {}
};

template<class Wrapped, class Ref, class Inner>
struct srcref_data
{
    typedef void IS_REF_OR_CREF;
    typedef Wrapped wrapped_t;

    typedef const Inner* data_ref_t;
    typedef const Inner* data_srcref_t;

    const Inner* inner;

    srcref_data(const Wrapped& o) : inner(o._data().inner) {}
    srcref_data(Ref o) : inner(o._data().inner) {}

    static srcref_data make(const Inner* f) {return srcref_data(f);}

private:
    srcref_data(const Inner* fp) : inner(fp) {}
};

template<class T, class Enable = void> struct is_flint_class : mp::false_ { };
template<class T>
struct is_flint_class<T, typename T::IS_FLINT_CLASS> : mp::true_ { };

template<class T, class Enable = void>
struct to_nonref {typedef typename T::evaluated_t type;};

template<class T>
struct to_nonref<T, typename T::data_t::IS_REF_OR_CREF>
{
    typedef typename T::data_t::wrapped_t type;
};

template<class T>
struct c_base_t
{
    typedef typename T::c_base_t type;
};

template<class T>
struct to_ref
{
    typedef typename T::template make_helper<operations::immediate, ref_data<
        typename to_nonref<T>::type, typename c_base_t<T>::type> >::type type;
};
template<class T>
struct to_srcref
{
    typedef typename T::template make_helper<operations::immediate, srcref_data<
        typename to_nonref<T>::type, typename to_ref<T>::type,
        typename c_base_t<T>::type> >::type type;
};

template<class Compare, class Ref>
struct is_ref : mp::equal_types<Ref, typename to_ref<Compare>::type> { };

template<class Compare, class Ref>
struct is_srcref : mp::equal_types<Ref, typename to_srcref<Compare>::type> { };

template<class T>
struct is_nonref : mp::equal_types<T, typename to_nonref<T>::type > { };

template<class T, class U, class Enable = void>
struct enableimplicit : mp::false_ { };

template<class T, class U>
struct enableimplicit<T, U, typename mp::enable_if<mp::and_<
        is_flint_class<T>, is_flint_class<U>
      > >::type>
    : mp::and_<
        traits::is_immediate_expr<T>,
        traits::is_immediate_expr<U>,
        mp::or_<
            mp::and_<is_ref<U, T>, is_nonref<U> >,
            mp::and_<is_srcref<U, T>, is_nonref<U> >,
            mp::and_<is_srcref<U, T>, is_ref<U, U> >
      > > { };

template<class Expr, class Enable = void>
struct maybe_data_ref
{
    typedef void data_ref_t;
    typedef void data_srcref_t;
};
template<class Expr>
struct maybe_data_ref<Expr,
    typename mp::enable_if<
        // NB: cannot use is_immediate, since Expr may be incomplete type!
        mp::equal_types<typename Expr::operation_t, operations::immediate> >::type>
{
    typedef typename Expr::data_t::data_ref_t data_ref_t;
    typedef typename Expr::data_t::data_srcref_t data_srcref_t;
};

template<class Base, class T, class Enable = void>
struct is_source : mp::false_ { };
template<class Base, class T>
struct is_source<Base, T, typename mp::enable_if<
    mp::and_<is_flint_class<T>, is_flint_class<Base> > >::type>
    : mp::or_<mp::equal_types<T, Base>, is_ref<Base, T>, is_srcref<Base, T> > { };

template<class Base, class T, class Enable = void>
struct is_target : mp::false_ { };
template<class Base, class T>
struct is_target<Base, T, typename mp::enable_if<
    mp::and_<is_flint_class<T>, is_flint_class<Base> > >::type>
    : mp::or_<mp::equal_types<T, Base>, is_ref<Base, T> > { };

template<class Base>
struct is_source_base
{
    template<class T> struct type : is_source<Base, T> { };
};
template<class Base>
struct is_target_base
{
    template<class T> struct type : is_target<Base, T> { };
};
} // flint_classes

namespace traits {
// XXX why do we need to disable the case where T, U are equal?
// Is <T, T> not more special?
template<class T, class U>
struct can_evaluate_into<T, U,
    typename mp::enable_if<mp::and_<flint_classes::is_flint_class<T>,
        flint_classes::is_flint_class<U>,
        mp::not_<mp::equal_types<T, U> > > >::type>
    : flint_classes::is_ref<U, T> { };
} // traits
} // flint

#define FLINTXX_DEFINE_BASICS(name)                                            \
public:                                                                       \
    typedef void IS_FLINT_CLASS;                                              \
                                                                              \
    template<class T>                                                         \
    struct doimplicit                                                         \
        : flint_classes::enableimplicit<name, T> { };                         \
                                                                              \
    template<class T>                                                         \
    name& operator=(const T& t)                                               \
    {                                                                         \
        this->set(t);                                                         \
        return *this;                                                         \
    }                                                                         \
                                                                              \
protected:                                                                    \
    explicit name(const Data& d) : base_t(d) {}                               \
                                                                              \
    template<class D, class O, class Da>                                      \
    friend class expression;

#define _FLINTXX_DEFINE_CTORS(name, ctx, ctxarg, MAYBE_COMMA)                  \
public:                                                                       \
    name(ctx) : base_t(ctxarg) {}                                             \
    template<class T>                                                         \
    explicit name(ctx MAYBE_COMMA const T& t,                                 \
            typename mp::disable_if<doimplicit<T> >::type* = 0)               \
        : base_t(ctxarg MAYBE_COMMA t) {}                                     \
    template<class T>                                                         \
    explicit name(ctx MAYBE_COMMA T& t,                                       \
            typename mp::disable_if<doimplicit<T> >::type* = 0)               \
        : base_t(ctxarg MAYBE_COMMA t) {}                                     \
    template<class T>                                                         \
    name(ctx MAYBE_COMMA const T& t,                                          \
            typename mp::enable_if<doimplicit<T> >::type* = 0)                \
        : base_t(ctxarg MAYBE_COMMA t) {}                                     \
    template<class T>                                                         \
    name(ctx MAYBE_COMMA T& t,                                                \
            typename mp::enable_if<doimplicit<T> >::type* = 0)                \
        : base_t(ctxarg MAYBE_COMMA t) {}

#define FLINTXX_EMPTY
#define FLINTXX_DEFINE_CTORS(name) \
    _FLINTXX_DEFINE_CTORS(name, FLINTXX_EMPTY, FLINTXX_EMPTY, FLINTXX_EMPTY)  \
    template<class T, class U>                                                \
    name(const T& t, const U& u) : base_t(t, u) {}

#define FLINTXX_COMMA ,

#define FLINTXX_DEFINE_C_REF(name, ctype, accessname)                          \
public:                                                                       \
    typedef ctype c_base_t;                                                   \
    typedef flint_classes::maybe_data_ref<name> wrapped_traits;               \
    typename wrapped_traits::data_ref_t accessname()                          \
    {                                                                         \
        return this->_data().inner;                                           \
    }                                                                         \
    typename wrapped_traits::data_srcref_t accessname() const                 \
    {                                                                         \
        return this->_data().inner;                                           \
    }                                                                         \
                                                                              \
    /* These only make sense with the reference types */                      \
    template<class T>                                                         \
    static name make(T& f)                                                    \
    {                                                                         \
        return name(Data::make(f));                                           \
    }                                                                         \
    template<class T>                                                         \
    static name make(const T& f)                                              \
    {                                                                         \
        return name(Data::make(f));                                           \
    }

#define FLINTXX_COND_S(Base) flint_classes::is_source_base<Base>::template type
#define FLINTXX_COND_T(Base) flint_classes::is_target_base<Base>::template type

#define FLINTXX_DEFINE_TO_STR(Base, eval) \
template<class T> \
struct to_string<T, \
    typename mp::enable_if<FLINTXX_COND_S(Base)<T> >::type> \
{ \
    static std::string get(const T& from, int base) \
    { \
        char* str = eval; \
        std::string res(str); \
        std::free(str); \
        return res; \
    } \
};

#define FLINTXX_DEFINE_CMP(Base, eval) \
template<class T, class U> \
struct cmp<T, U, \
    typename mp::enable_if<mp::and_<FLINTXX_COND_S(Base)<T>, \
        FLINTXX_COND_S(Base)<U> > >::type> \
{ \
    static int get(const T& e1, const U& e2) \
    { \
        return eval; \
    } \
};

#endif
