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

#ifndef NMOD_POLYXX_H
#define NMOD_POLYXX_H

#include <cstdlib>

#include "nmod_poly.h"

#include "fmpzxx.h"
#include "nmod_vecxx.h"

#include "flintxx/expression.h"
#include "flintxx/flint_classes.h"
#include "flintxx/flint_exception.h"
#include "flintxx/frandxx.h"
#include "flintxx/stdmath.h"
#include "flintxx/traits.h"

// TODO exhibit this as a specialisation of a generic poly<nmodxx>
// TODO input

namespace flint {
FLINT_DEFINE_BINOP(nmod_polyxx_get_coeff)

template<class Operation, class Data>
class nmod_polyxx_expression
    : public expression<derived_wrapper<nmod_polyxx_expression>,
                            Operation, Data>
{
public:
    typedef expression<derived_wrapper< ::flint::nmod_polyxx_expression>,
              Operation, Data> base_t;

    FLINTXX_DEFINE_BASICS(nmod_polyxx_expression)
    FLINTXX_DEFINE_CTORS(nmod_polyxx_expression)
    FLINTXX_DEFINE_C_REF(nmod_polyxx_expression, nmod_poly_struct, _poly)

    // these only make sense with immediates
    void realloc(slong alloc) {nmod_poly_realloc(_poly(), alloc);}
    void fit_length(slong len) {nmod_poly_fit_length(_poly(), len);}
    void _normalise() {_nmod_poly_normalise(_poly());}
    nmodxx_ctx_srcref _ctx() const
        {return nmodxx_ctx_srcref::make(_poly()->mod);}

    // These only make sense with target immediates
    void set_coeff(slong n, ulong c) {nmod_poly_set_coeff_ui(_poly(), n, c);}
    template<class Nmod>
    typename mp::enable_if<traits::is_nmodxx<Nmod> >::type
    set_coeff(slong j, const Nmod& c)
    {
        // TODO this does not need reduction
        nmod_poly_set_coeff_ui(_poly(), j, c.template to<mp_limb_t>());
    }
    void truncate(slong n) {nmod_poly_truncate(_poly(), n);}

    void set_randtest(frandxx& state, slong len)
        {nmod_poly_randtest(_poly(), state._data(), len);}
    void set_randtest_irreducible(frandxx& state, slong len)
        {nmod_poly_randtest_irreducible(_poly(), state._data(), len);}

    // These work on any expression without evaluation
    nmodxx_ctx_srcref estimate_ctx() const;
    mp_limb_t modulus() const {return estimate_ctx().n();}

    evaluated_t create_temporary() const
    {
        return evaluated_t(modulus());
    }

    // These cause evaluation
    slong length() const {return nmod_poly_length(this->evaluate()._poly());}
    slong degree() const {return nmod_poly_degree(this->evaluate()._poly());}
    bool is_one() const {return nmod_poly_is_one(this->evaluate()._poly());}
    bool is_zero() const {return nmod_poly_is_zero(this->evaluate()._poly());}
    bool is_squarefree() const
        {return nmod_poly_is_squarefree(this->evaluate()._poly());}
    bool is_irreducible() const
        {return nmod_poly_is_irreducible(this->evaluate()._poly());}
    slong max_bits() const {return nmod_poly_max_bits(this->evaluate()._poly());}

    // Lazy members
    FLINTXX_DEFINE_MEMBER_BINOP_(get_coeff, nmod_polyxx_get_coeff)

    FLINTXX_DEFINE_MEMBER_BINOP(poly_shift_left)
    FLINTXX_DEFINE_MEMBER_BINOP(poly_shift_right)
    FLINTXX_DEFINE_MEMBER_BINOP(reverse)

    FLINTXX_DEFINE_MEMBER_BINOP_(bit_pack, poly_bit_pack)

    FLINTXX_DEFINE_MEMBER_UNOP(make_monic)
};

namespace detail {
struct nmod_poly_data;
}

typedef nmod_polyxx_expression<operations::immediate, detail::nmod_poly_data>
           nmod_polyxx;
typedef nmod_polyxx_expression<operations::immediate,
            flint_classes::ref_data<nmod_polyxx, nmod_poly_struct> >
           nmod_polyxx_ref;
typedef nmod_polyxx_expression<operations::immediate,
            flint_classes::srcref_data<
                nmod_polyxx, nmod_polyxx_ref, nmod_poly_struct> >
           nmod_polyxx_srcref;

#define NMOD_POLYXX_COND_S FLINTXX_COND_S(nmod_polyxx)
#define NMOD_POLYXX_COND_T FLINTXX_COND_T(nmod_polyxx)

namespace detail {
struct nmod_poly_data
{
    nmod_poly_t inner;
    typedef nmod_poly_t& data_ref_t;
    typedef const nmod_poly_t& data_srcref_t;

    nmod_poly_data(mp_limb_t n) {nmod_poly_init(inner, n);}
    nmod_poly_data(mp_limb_t n, slong alloc) {nmod_poly_init2(inner, n, alloc);}
    ~nmod_poly_data() {nmod_poly_clear(inner);}

    nmod_poly_data(const nmod_poly_data& o)
    {
        nmod_poly_init2_preinv(inner, o.inner->mod.n,
                o.inner->mod.ninv, o.inner->length);
        nmod_poly_set(inner, o.inner);
    }

    nmod_poly_data(nmod_polyxx_srcref r)
    {
        nmod_poly_init2_preinv(inner, r.modulus(),
                r._poly()->mod.ninv, r.length());
        nmod_poly_set(inner, r._poly());
    }
};
} // detail
namespace traits {
template<> struct has_nmodxx_ctx<nmod_polyxx> : mp::true_ { };
template<> struct has_nmodxx_ctx<nmod_polyxx_ref> : mp::true_ { };
template<> struct has_nmodxx_ctx<nmod_polyxx_srcref> : mp::true_ { };
} // traits
template<class Operation, class Data>
inline nmodxx_ctx_srcref
nmod_polyxx_expression<Operation, Data>::estimate_ctx() const
{
    return tools::find_nmodxx_ctx(*this);
}

namespace rules {
#define NMOD_POLYXX_COND_S FLINTXX_COND_S(nmod_polyxx)
#define NMOD_POLYXX_COND_T FLINTXX_COND_T(nmod_polyxx)

FLINT_DEFINE_DOIT_COND2(assignment, NMOD_POLYXX_COND_T, NMOD_POLYXX_COND_S,
        nmod_poly_set(to._poly(), from._poly()))

FLINTXX_DEFINE_ASSIGN_STR(nmod_polyxx, execution_check(
            nmod_poly_set_str(to._poly(), from), "assign string", "nmod_polyxx"))

FLINTXX_DEFINE_EQUALS(nmod_polyxx, nmod_poly_equal(e1._poly(), e2._poly()))
FLINTXX_DEFINE_TO_STR(nmod_polyxx, nmod_poly_get_str(from._poly()))
FLINTXX_DEFINE_SWAP(nmod_polyxx, nmod_poly_swap(e1._poly(), e2._poly()))

FLINT_DEFINE_BINARY_EXPR_COND2(nmod_polyxx_get_coeff_op, nmodxx,
        NMOD_POLYXX_COND_S, traits::fits_into_slong,
        to.set_nored(nmod_poly_get_coeff_ui(e1._poly(), e2)))

FLINT_DEFINE_BINARY_EXPR_COND2(plus, nmod_polyxx,
        NMOD_POLYXX_COND_S, NMOD_POLYXX_COND_S,
        nmod_poly_add(to._poly(), e1._poly(), e2._poly()))
FLINT_DEFINE_BINARY_EXPR_COND2(minus, nmod_polyxx,
        NMOD_POLYXX_COND_S, NMOD_POLYXX_COND_S,
        nmod_poly_sub(to._poly(), e1._poly(), e2._poly()))

FLINT_DEFINE_UNARY_EXPR_COND(negate, nmod_polyxx, NMOD_POLYXX_COND_S,
        nmod_poly_neg(to._poly(), from._poly()))

FLINT_DEFINE_BINARY_EXPR_COND2(reverse_op, nmod_polyxx,
        NMOD_POLYXX_COND_S, traits::fits_into_slong,
        nmod_poly_reverse(to._poly(), e1._poly(), e2))

FLINT_DEFINE_BINARY_EXPR_COND2(poly_shift_left_op, nmod_polyxx,
        NMOD_POLYXX_COND_S, traits::fits_into_slong,
        nmod_poly_shift_left(to._poly(), e1._poly(), e2))
FLINT_DEFINE_BINARY_EXPR_COND2(poly_shift_right_op, nmod_polyxx,
        NMOD_POLYXX_COND_S, traits::fits_into_slong,
        nmod_poly_shift_right(to._poly(), e1._poly(), e2))

FLINT_DEFINE_CBINARY_EXPR_COND2(times, nmod_polyxx,
        NMOD_POLYXX_COND_S, NMODXX_COND_S,
        nmod_poly_scalar_mul_nmod(to._poly(), e1._poly(), e2._limb()))
FLINT_DEFINE_UNARY_EXPR_COND(make_monic_op, nmod_polyxx, NMOD_POLYXX_COND_S,
        nmod_poly_make_monic(to._poly(), from._poly()))

FLINT_DEFINE_BINARY_EXPR_COND2(poly_bit_pack_op, fmpzxx,
        NMOD_POLYXX_COND_S, traits::fits_into_mp_bitcnt_t,
        nmod_poly_bit_pack(to._fmpz(), e1._poly(), e2))
FLINT_DEFINE_BINARY_EXPR_COND2(poly_bit_unpack_op, nmod_polyxx,
        FMPZXX_COND_S, traits::fits_into_mp_bitcnt_t,
        nmod_poly_bit_unpack(to._poly(), e1._fmpz(), e2))
} // rules
} // flint

#endif
