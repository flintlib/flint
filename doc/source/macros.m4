dnl############################################################################
dnl Change quotation marks to avoid conflicts
dnl############################################################################
changequote({{{,}}})dnl
dnl############################################################################
dnl helper stuff
dnl############################################################################
define({{{_neg_}}},dnl
{{{-$1}}})dnl
dnl
define({{{_add_}}},dnl
{{{$1 + $2}}})dnl
dnl
define({{{_sub_}}},dnl
{{{$1 - $2}}})dnl
dnl
define({{{_mul_}}},dnl
{{{$1 \cdot $2}}})dnl
dnl
define({{{_div_}}},dnl
{{{$1 / $2}}})dnl
dnl
define({{{_addmul_}}},dnl
{{{_add_($1, _mul_($2, $3))}}})dnl
dnl
define({{{_submul_}}},dnl
{{{_sub_($1, _mul_($2, $3))}}})dnl
dnl
define({{{_fmma_}}},dnl
{{{_add_(_mul_($1, $2), _mul_($3, $4))}}})dnl
dnl
define({{{_fmms_}}},dnl
{{{_sub_(_mul_($1, $2), _mul_($3, $4))}}})dnl
dnl
define({{{_lt_}}},dnl
{{{$1 < $2}}})dnl
dnl
define({{{_gt_}}},dnl
{{{$1 > $2}}})dnl
dnl
define({{{_equal_}}},dnl
{{{$1 = $2}}})dnl
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
    Sets `$1` to zero.dnl
}}})dnl
define({{{desc_one}}},{{{
    Sets `$1` to one.dnl
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
    exact.  If `$3` is zero, an exception is raised.dnl
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
dnl sqrt
dnl############################################################################
define({{{desc_sqrt_nonordered_ring}}},{{{
    If `$2` is a perfect square, sets `$1` to a square root of `$2`
    and returns nonzero.  Otherwise returns zero.dnl
}}})dnl
dnl############################################################################
dnl comparisons
dnl############################################################################
define({{{desc_cmp}}},{{{
    Returns a negative value if `_lt_($1, $2)`, positive value if
    `_gt_($1, $2)`, otherwise returns zero.dnl
}}})dnl
define({{{desc_equal}}},{{{
    Returns nonzero if `_equal_($1, $2)`, otherwise returns zero.dnl
}}})dnl
define({{{desc_is_zero}}},{{{
    Returns nonzero if `_equal_($1, 0)`, otherwise returns zero.dnl
}}})dnl
define({{{desc_is_one}}},{{{
    Returns nonzero if `_equal_($1, 1)`, otherwise returns zero.dnl
}}})dnl
