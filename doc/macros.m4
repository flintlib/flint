dnl############################################################################
dnl Change quotation marks to avoid conflicts
dnl############################################################################
changequote({{{,}}})dnl
dnl############################################################################
dnl notation
dnl############################################################################
define({{{_neg_}}},dnl
{{{-$1}}})dnl
define({{{_add_}}},dnl
{{{$1 + $2}}})dnl
define({{{_sub_}}},dnl
{{{$1 - $2}}})dnl
define({{{_mul_}}},dnl
{{{$1 \cdot $2}}})dnl
define({{{_div_}}},dnl
{{{$1 / $2}}})dnl
define({{{_pow_}}},dnl
{{{{$1}^{$2}}}})dnl
define({{{_addmul_}}},dnl
{{{_add_($1, _mul_($2, $3))}}})dnl
define({{{_submul_}}},dnl
{{{_sub_($1, _mul_($2, $3))}}})dnl
define({{{_fmma_}}},dnl
{{{_add_(_mul_($1, $2), _mul_($3, $4))}}})dnl
define({{{_fmms_}}},dnl
{{{_sub_(_mul_($1, $2), _mul_($3, $4))}}})dnl
define({{{_lt_}}},dnl
{{{$1 < $2}}})dnl
define({{{_gt_}}},dnl
{{{$1 > $2}}})dnl
define({{{_eq_}}},dnl
{{{$1 = $2}}})dnl
dnl
define({{{_boolzero_}}},dnl
{{{`0`}}})dnl
define({{{_boolpos_}}},dnl
{{{`1`}}})dnl
define({{{_boolneg}}},dnl
{{{`-1`}}})dnl
define({{{_zero_}}},dnl
{{{`0`}}})dnl
define({{{_one_}}},dnl
{{{`1`}}})dnl
dnl############################################################################
dnl memory management
dnl############################################################################
define({{{desc_init_set}}},{{{
    Initialises `$1` and sets it to `$2`.dnl
}}})dnl
dnl############################################################################
dnl set
dnl############################################################################
define({{{desc_set}}},{{{
    Sets `$1` to `$2`.dnl
}}})dnl
define({{{desc_zero}}},{{{
    Sets `$1` to _zero_.dnl
}}})dnl
define({{{desc_one}}},{{{
    Sets `$1` to _one_.dnl
}}})dnl
dnl############################################################################
dnl negation, absolute value etc.
dnl############################################################################
define({{{desc_neg}}},{{{
    Sets `$1` to `_neg_($2)`.dnl
}}})dnl
define({{{desc_abs}}},{{{
    Sets `$1` to the absolute value of `$2`.dnl
}}})dnl
dnl############################################################################
dnl basic arithmetic operations
dnl############################################################################
define({{{desc_add}}},{{{
    Sets `$1` to `_add_($2, $3)`.dnl
}}})dnl
define({{{desc_sub}}},{{{
    Sets `$1` to `_sub_($2, $3)`.dnl
}}})dnl
define({{{desc_mul}}},{{{
    Sets `$1` to `_mul_($2, $3)`.dnl
}}})dnl
define({{{desc_divexact}}},{{{
    Sets `$1` to `_div_($2, $3)` under the assumption that the division is
    exact.  If `$3` is _zero_, an exception is raised.dnl
}}})dnl
dnl############################################################################
dnl extended basic arithmetic operations
dnl############################################################################
define({{{desc_addmul}}},{{{
    Sets `$1` to `_addmul_($1, $2, $3)`.dnl
}}})dnl
define({{{desc_submul}}},{{{
    Sets `$1` to `_submul_($1, $2, $3)`.dnl
}}})dnl
define({{{desc_fmma}}},{{{
    Sets `$1` to `_fmma_($2, $3, $4, $5)`.dnl
}}})dnl
define({{{desc_fmms}}},{{{
    Sets `$1` to `_fmms_($2, $3, $4, $5)`.dnl
}}})dnl
dnl############################################################################
dnl powering
dnl############################################################################
define({{{desc_pow}}},{{{
    Sets `$1` to `_pow_($2, $3)`.  Defines `_eq_(_pow_(0, 0), 1)`.dnl
}}})dnl
dnl############################################################################
dnl comparisons
dnl############################################################################
define({{{desc_cmp}}},{{{
    Returns _boolneg_ if `_lt_($1, $2)`, _boolpos_ if `_gt_($1, $2)`, otherwise
    returns zero.dnl
}}})dnl
define({{{desc_equal}}},{{{
    Returns _boolpos_ if `_eq_($1, $2)`, otherwise returns _boolzero_.dnl
}}})dnl
define({{{desc_is_zero}}},{{{
    Returns _boolpos_ if `_eq_($1, 0)`, otherwise returns _boolzero_.dnl
}}})dnl
define({{{desc_is_one}}},{{{
    Returns _boolpos_ if `_eq_($1, 1)`, otherwise returns _boolzero_.dnl
}}})dnl
define({{{desc_sgn}}},{{{
    Returns the sign of `$1`.  That is, returns `-1` if `_lt_($1, 0)`, `1` if
    `_gt_($1, 0)` and `0` if `_eq_($1, 0)`.dnl
}}})dnl
dnl############################################################################
dnl factoring
dnl############################################################################
define({{{desc_sqrt_nonordered_ring}}},{{{
    If `$2` is a perfect square, sets `$1` to a square root of `$2`
    and returns _boolpos_.  Otherwise returns _boolzero_.dnl
}}})dnl
define({{{desc_gcd_int}}},{{{
    Sets `$1` to the greatest common divisor of `$2` and `$3`.  The result is
    always non-negative.dnl
}}})dnl
define({{{desc_divisible}}},{{{
    Returns _boolpos_ if there is an `x` such that `_eq_($1, _mul_(x, $2))`,
    and returns _boolzero_ if there is none.dnl
}}})dnl
define({{{desc_divides}}},{{{
    Returns _boolpos_ if there is an `x` such that `_eq_($2, _mul_(x, $3))` and
    sets `_eq_($1, x)`, and returns _boolzero_ if there is none and sets
    `_eq_($1, 0)`.dnl
}}})dnl
