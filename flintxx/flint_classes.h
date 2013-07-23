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

#include "flintxx/mp.h"
#include "flintxx/expression_traits.h"

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

template<class T, class Right1, class Right2>
struct ternary_assign_helper
{
    typedef tools::evaluate_2<Right1, Right2> ev2_t;
    typedef typename ev2_t::temporaries_t temporaries_t;
    typedef mp::back_tuple<temporaries_t> back_t;

    typename back_t::type backing;
    ev2_t ev2;

    static temporaries_t backtemps(typename back_t::type& backing)
    {
        temporaries_t temps;
        back_t::init(temps, backing);
        return temps;
    }

    ternary_assign_helper(typename ev2_t::arg1_t r1, typename ev2_t::arg2_t r2)
        : backing(mp::htuples::fill<typename back_t::type>(
                    tools::temporaries_filler(r1+r2 /* XXX */))),
          ev2(backtemps(backing), r1, r2) {}
    const T& getleft() {return ev2.get1();}
    const T& getright() {return ev2.get2();}
};

template<class T, class Right1, class Right2>
struct enable_ternary_assign
    : mp::enable_if<mp::and_<
          traits::is_T_expr<typename traits::basetype<Right1>::type, T>,
          traits::is_T_expr<typename traits::basetype<Right2>::type, T> >, T&> { };
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

#define FLINTXX_DEFINE_BASICS(name)                                           \
public:                                                                       \
    typedef void IS_FLINT_CLASS;                                              \
    typedef typename base_t::evaluated_t evaluated_t;                         \
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

#define FLINTXX_DEFINE_CTORS(name)                                           \
public:                                                                       \
    name() : base_t() {}                                                      \
    template<class T>                                                         \
    explicit name(const T& t,                                                 \
            typename mp::disable_if<doimplicit<T> >::type* = 0)               \
        : base_t(t) {}                                                        \
    template<class T>                                                         \
    explicit name(T& t,                                                       \
            typename mp::disable_if<doimplicit<T> >::type* = 0)               \
        : base_t(t) {}                                                        \
    template<class T>                                                         \
    name(const T& t,                                                          \
            typename mp::enable_if<doimplicit<T> >::type* = 0)                \
        : base_t(t) {}                                                        \
    template<class T>                                                         \
    name(T& t,                                                                \
            typename mp::enable_if<doimplicit<T> >::type* = 0)                \
        : base_t(t) {}                                                        \
    template<class T, class U>                                                \
    name(const T& t, const U& u) : base_t(t, u) {}

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
    }                                                                         \
    template<class T, class U>                                                \
    static name make(T& f, const U& u)                                        \
    {                                                                         \
        return name(Data::make(f, u));                                        \
    }                                                                         \
    template<class T, class U>                                                \
    static name make(const T& f, const U& u)                                  \
    {                                                                         \
        return name(Data::make(f, u));                                        \
    }

#define FLINTXX_COND_S(Base) flint_classes::is_source_base<Base>::template type
#define FLINTXX_COND_T(Base) flint_classes::is_target_base<Base>::template type

#define FLINTXX_DEFINE_TO_STR(Base, eval) \
template<class T> \
struct to_string<T, \
    typename mp::enable_if< FLINTXX_COND_S(Base)<T> >::type> \
{ \
    static std::string get(const T& from, int base) \
    { \
        char* str = eval; \
        std::string res(str); \
        std::free(str); \
        return res; \
    } \
};

#define FLINTXX_DEFINE_SWAP(Base, eval) \
template<class T, class U> \
struct swap<T, U, typename mp::enable_if< mp::and_< \
    FLINTXX_COND_T(Base)<T>, FLINTXX_COND_T(Base)<U> > >::type> \
{ \
    static void doit(T& e1, U& e2) \
    { \
        eval; \
    } \
};

#define FLINTXX_DEFINE_CONVERSION_TMP(totype, Base, eval) \
template<class T> \
struct conversion<totype, T, \
    typename mp::enable_if< FLINTXX_COND_S(Base)<T> >::type> \
{ \
    static totype get(const T& from) \
    { \
        totype to; \
        eval; \
        return to; \
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

#define FLINTXX_DEFINE_EQUALS(Base, eval) \
template<class T, class U> \
struct equals<T, U, typename mp::enable_if<mp::and_< \
    FLINTXX_COND_S(Base)<T>, FLINTXX_COND_S(Base)<U> > >::type> \
{ \
    static bool get(const T& e1, const U& e2) \
    { \
        return eval; \
    } \
};

#define FLINTXX_DEFINE_RANDFUNC(CBase, name) \
static CBase##xx_expression name(frandxx& state, mp_bitcnt_t bits) \
{ \
    CBase##xx_expression res; \
    CBase##_##name(res._data().inner, state._data(), bits); \
    return res; \
}

// XXX due to circular definition problems, this cannot use the usual
// unary_op_helper type approach, and the unop must return the same type
// of expression
#define FLINTXX_DEFINE_MEMBER_UNOP(name, funcname) \
typename base_t::template make_helper<operations::funcname##_op, \
    tuple<typename detail::storage_traits<typename base_t::derived_t>::type, \
            empty_tuple> >::type \
name() const \
{ \
    return funcname(*this); \
}


#define FLINTXX_UNADORNED_MAKETYPES(Base, left, op, right) \
    Base##_expression< op, tuple< left, tuple< right, empty_tuple> > >

// Optimized evaluation rules using ternary arithmetic (addmul, submul)
// NB: this has to be called in namespace flint, not flint::rules!
#define FLINTXX_DEFINE_TERNARY(Base, addmuleval, submuleval, maketypes)       \
namespace rules {                                                             \
/* a +- b*c */                                                                \
template<class Op, class Left, class Right1, class Right2>                    \
struct evaluation<Op,                                                         \
    tuple<Left, tuple<                                                        \
        maketypes(Base, Right1, operations::times, Right2),                   \
        /* NB: there is no particular reason to have the enable_if here,      \
               many other similar places would do */                          \
        typename mp::enable_if<mp::or_<                                       \
                mp::equal_types<Op, operations::plus>,                        \
                mp::equal_types<Op, operations::minus> >,                     \
            empty_tuple>::type> >,                                            \
    true, 1,                                                                  \
    typename tools::ternary_helper<Base, Left, Right1, Right2>::enable::type> \
{                                                                             \
    /* Helpful for testing. */                                                \
    static const unsigned TERNARY_OP_MARKER = 0;                              \
                                                                              \
    typedef Base return_t;                                                    \
    typedef tools::ternary_helper<Base, Left, Right1, Right2> th;             \
    typedef typename th::temporaries_t temporaries_t;                         \
    typedef tuple<Left, tuple<                                                \
        maketypes(Base, Right1, operations::times, Right2),                   \
        empty_tuple> > data_t;                                                \
    static const bool is_add = mp::equal_types<Op, operations::plus>::val;    \
                                                                              \
    static void doit(const data_t& input, temporaries_t temps, return_t* res) \
    {                                                                         \
        const Base* left = 0;                                                 \
        const Base* right = 0;                                                \
        th::doit(input.first(), input.second()._data().first(),               \
                input.second()._data().second(), temps, res, right, left);    \
        const Base& e1 = *left;                                               \
        const Base& e2 = *right;                                              \
        Base& to = *res;                                                      \
        if(is_add)                                                            \
        {                                                                     \
            addmuleval;                                                       \
        }                                                                     \
        else                                                                  \
        {                                                                     \
            submuleval;                                                       \
        }                                                                     \
    }                                                                         \
};                                                                            \
                                                                              \
/* b*c + a */                                                                 \
template<class Right, class Left1, class Left2>                               \
struct evaluation<operations::plus,                                           \
    tuple<maketypes(Base, Left1, operations::times, Left2),                   \
        tuple<Right, empty_tuple> >,                                          \
    true, 1,                                                                  \
    typename tools::ternary_helper<Base,                                      \
        Right, Left1, Left2, operations::times>::enable::type>                \
{                                                                             \
    /* Helpful for testing. */                                                \
    static const unsigned TERNARY_OP_MARKER = 0;                              \
                                                                              \
    typedef Base return_t;                                                    \
    typedef tools::ternary_helper<Base, Right, Left1, Left2> th;              \
    typedef typename th::temporaries_t temporaries_t;                         \
    typedef tuple<maketypes(Base, Left1, operations::times, Left2),           \
        tuple<Right, empty_tuple> > data_t;                                   \
                                                                              \
    static void doit(const data_t& input, temporaries_t temps, return_t* res) \
    {                                                                         \
        const Base* left = 0;                                                 \
        const Base* right = 0;                                                \
        th::doit(input.second(), input.first()._data().first(),               \
                input.first()._data().second(), temps, res, right, left);     \
        const Base& e1 = *left;                                               \
        const Base& e2 = *right;                                              \
        Base& to = *res;                                                      \
        addmuleval;                                                           \
    }                                                                         \
};                                                                            \
} /* rules */                                                                 \
                                                                              \
/* TODO enable these with references on left hand side(?) */                  \
/* a += b*c */                                                                \
template<class Right1, class Right2>                                          \
inline typename flint_classes::enable_ternary_assign<Base, Right1, Right2>::type \
operator+=(Base& to,                                                          \
        const maketypes(Base, Right1, operations::times, Right2)& other)      \
{                                                                             \
    flint_classes::ternary_assign_helper<Base, Right1, Right2> tah(           \
            other._data().first(), other._data().second());                   \
    const Base& e1 = tah.getleft();                                           \
    const Base& e2 = tah.getright();                                          \
    addmuleval;                                                               \
    return to;                                                                \
}                                                                             \
                                                                              \
/* a -= b*c */                                                                \
template<class Right1, class Right2>                                          \
inline typename flint_classes::enable_ternary_assign<Base, Right1, Right2>::type \
operator-=(Base& to,                                                          \
        const maketypes(Base, Right1, operations::times, Right2)& other)      \
{                                                                             \
    flint_classes::ternary_assign_helper<Base, Right1, Right2> tah(           \
            other._data().first(), other._data().second());                   \
    const Base& e1 = tah.getleft();                                           \
    const Base& e2 = tah.getright();                                          \
    submuleval;                                                               \
    return to;                                                                \
}

#endif
