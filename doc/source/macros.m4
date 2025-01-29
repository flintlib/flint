dnl Change quotes to something that is not used anywhere
changequote({{{,}}})dnl
dnl
dnl############################################################################
dnl helper stuff
dnl############################################################################
define({{{_prefix}}},{{{.. function::}}})dnl
define({{{_indent}}},{{{             }}})dnl
dnl
define({{{_res_}}},{{{r}}})dnl
define({{{_ip1_}}},{{{a}}})dnl
define({{{_ip2_}}},{{{b}}})dnl
dnl
define({{{_neg_}}},{{{-$1}}})dnl
define({{{_add_}}},{{{$1 + $2}}})dnl
define({{{_sub_}}},{{{$1 - $2}}})dnl
define({{{_mul_}}},{{{$1 \cdot $2}}})dnl
define({{{_div_}}},{{{$1 / $2}}})dnl
dnl
define({{{_addmul_}}},{{{_add_($1, _mul_($2, $3))}}})dnl
define({{{_submul_}}},{{{_sub_($1, _mul_($2, $3))}}})dnl
dnl
define({{{_lt_}}},{{{$1 < $2}}})dnl
define({{{_gt_}}},{{{$1 < $2}}})dnl
define({{{_equal_}}},{{{$1 = $2}}})dnl
dnl############################################################################
dnl set
dnl############################################################################
define({{{func_set}}},dnl
{{{dnl
void $1_set($1_t _res_, const $1_t _ip1_)dnl
}}})dnl
define({{{func_set_si}}},dnl
{{{dnl
void $1_set_si($1_t _res_, slong _ip1_)dnl
}}})dnl
define({{{func_set_ui}}},dnl
{{{dnl
void $1_set_ui($1_t _res_, ulong _ip1_)dnl
}}})dnl
define({{{desc_set}}},{{{
    Sets `_res_` to `_ip1_`.dnl
}}})dnl
dnl############################################################################
dnl set to common constants
dnl############################################################################
define({{{func_zero}}},dnl
{{{dnl
void $1_zero($1_t _res_)dnl
}}})dnl
define({{{desc_zero}}},{{{
    Sets `_res_` to zero.dnl
}}})dnl
define({{{func_one}}},dnl
{{{dnl
void $1_one($1_t _res_)dnl
}}})dnl
define({{{desc_one}}},{{{
    Sets `_res_` to one.dnl
}}})dnl
dnl############################################################################
dnl negation
dnl############################################################################
define({{{func_neg}}},dnl
{{{dnl
void $1_neg($1_t _res_, const $1_t _ip1_)dnl
}}})dnl
define({{{desc_neg}}},{{{
    Sets `_res_` to `_neg_(_ip1_)`.dnl
}}})dnl
dnl############################################################################
dnl absolute value
dnl############################################################################
define({{{func_abs}}},dnl
{{{dnl
void $1_abs($1_t _res_, const $1_t _ip1_)dnl
}}})dnl
define({{{desc_abs}}},{{{
    Sets `_res_` to the absolute value of `_ip1_`.dnl
}}})dnl
dnl############################################################################
dnl addition
dnl############################################################################
define({{{func_add}}},dnl
{{{dnl
void $1_add($1_t _res_, const $1_t _ip1_, const $1_t _ip2_)dnl
}}})dnl
define({{{func_add_si}}},dnl
{{{dnl
void $1_add_si($1_t _res_, const $1_t _ip1_, slong _ip2_)dnl
}}})dnl
define({{{func_add_ui}}},dnl
{{{dnl
void $1_add_ui($1_t _res_, const $1_t _ip1_, ulong _ip2_)dnl
}}})dnl
define({{{desc_add}}},{{{
    Sets `_res_` to `_add_(_ip1_, _ip2_)`.dnl
}}})dnl
dnl############################################################################
dnl subtraction
dnl############################################################################
define({{{func_sub}}},dnl
{{{dnl
void $1_sub($1_t _res_, const $1_t _ip1_, const $1_t _ip2_)dnl
}}})dnl
define({{{func_sub_si}}},dnl
{{{dnl
void $1_sub_si($1_t _res_, const $1_t _ip1_, slong _ip2_)dnl
}}})dnl
define({{{func_sub_ui}}},dnl
{{{dnl
void $1_sub_ui($1_t _res_, const $1_t _ip1_, ulong _ip2_)dnl
}}})dnl
define({{{desc_sub}}},{{{
    Sets `_res_` to `_sub_(_ip1_, _ip2_)`.dnl
}}})dnl
dnl############################################################################
dnl multiplication
dnl############################################################################
define({{{func_mul}}},dnl
{{{dnl
void $1_mul($1_t _res_, const $1_t _ip1_, const $1_t _ip2_)dnl
}}})dnl
define({{{func_mul_si}}},dnl
{{{dnl
void $1_mul_si($1_t _res_, const $1_t _ip1_, slong _ip2_)dnl
}}})dnl
define({{{func_mul_ui}}},dnl
{{{dnl
void $1_mul_ui($1_t _res_, const $1_t _ip1_, ulong _ip2_)dnl
}}})dnl
define({{{desc_mul}}},{{{
    Sets `_res_` to `_mul_(_ip1_, _ip2_)`.dnl
}}})dnl
dnl############################################################################
dnl exact division
dnl############################################################################
define({{{func_divexact}}},dnl
{{{dnl
void $1_divexact($1_t _res_, const $1_t _ip1_, const $1_t _ip2_)dnl
}}})dnl
define({{{func_divexact_si}}},dnl
{{{dnl
void $1_divexact_si($1_t _res_, const $1_t _ip1_, slong _ip2_)dnl
}}})dnl
define({{{func_divexact_ui}}},dnl
{{{dnl
void $1_divexact_ui($1_t _res_, const $1_t _ip1_, ulong _ip2_)dnl
}}})dnl
define({{{desc_divexact}}},{{{
    Sets `_res_` to `_div_(_ip1_, _ip2_)` under the assumption that the
    division is exact.  If `_ip2_` is zero, an exception is raised.dnl
}}})dnl
dnl############################################################################
dnl addmul
dnl############################################################################
define({{{func_addmul}}},dnl
{{{dnl
void $1_addmul($1_t _res_, const $1_t _ip1_, const $1_t _ip2_)dnl
}}})dnl
define({{{func_addmul_si}}},dnl
{{{dnl
void $1_addmul_si($1_t _res_, const $1_t _ip1_, slong _ip2_)dnl
}}})dnl
define({{{func_addmul_ui}}},dnl
{{{dnl
void $1_addmul_ui($1_t _res_, const $1_t _ip1_, ulong _ip2_)dnl
}}})dnl
define({{{desc_addmul}}},{{{
    Sets `_res_` to `_addmul_(_res_, _ip1_, _ip2_)`.dnl
}}})dnl
dnl############################################################################
dnl submul
dnl############################################################################
define({{{func_submul}}},dnl
{{{dnl
void $1_submul($1_t _res_, const $1_t _ip1_, const $1_t _ip2_)dnl
}}})dnl
define({{{func_submul_si}}},dnl
{{{dnl
void $1_submul_si($1_t _res_, const $1_t _ip1_, slong _ip2_)dnl
}}})dnl
define({{{func_submul_ui}}},dnl
{{{dnl
void $1_submul_ui($1_t _res_, const $1_t _ip1_, ulong _ip2_)dnl
}}})dnl
define({{{desc_submul}}},{{{
    Sets `_res_` to `_submul_(_res_, _ip1_, _ip2_)`.dnl
}}})dnl
dnl############################################################################
dnl sqrt
dnl############################################################################
define({{{func_sqrt}}},dnl
{{{dnl
void $1_sqrt($1_t _res_, const $1_t _ip1_)dnl
}}})dnl
define({{{desc_sqrt_nonordered_ring}}},{{{
    If `_ip1_` is a perfect square, sets `_res_` to a square root of `_ip1_`
    and returns nonzero.  Otherwise returns zero.dnl
}}})dnl
dnl############################################################################
dnl comparisons
dnl############################################################################
define({{{func_cmp}}},dnl
{{{dnl
int $1_cmp(const $1_t _ip1_, const $1_t _ip2_)dnl
}}})dnl
define({{{func_cmp_si}}},dnl
{{{dnl
int $1_cmp_si(const $1_t _ip1_, slong _ip2_)dnl
}}})dnl
define({{{func_cmp_ui}}},dnl
{{{dnl
int $1_cmp_ui(const $1_t _ip1_, ulong _ip2_)dnl
}}})dnl
define({{{func_cmp_fmpz}}},dnl
{{{dnl
int $1_cmp_fmpz(const $1_t _ip1_, const fmpz_t _ip2_)dnl
}}})dnl
define({{{desc_cmp}}},{{{
    Returns a negative value if `_lt_(_ip1_, _ip2_)`, positive value if
    `_gt_(_ip1_, _ip2_)`, otherwise returns zero.dnl
}}})dnl
dnl############################################################################
dnl equality
dnl############################################################################
define({{{func_equal}}},dnl
{{{dnl
int $1_equal(const $1_t _ip1_, const $1_t _ip2_)dnl
}}})dnl
define({{{func_equal_si}}},dnl
{{{dnl
int $1_equal_si(const $1_t _ip1_, slong _ip2_)dnl
}}})dnl
define({{{func_equal_ui}}},dnl
{{{dnl
int $1_equal_ui(const $1_t _ip1_, ulong _ip2_)dnl
}}})dnl
define({{{func_equal_fmpz}}},dnl
{{{dnl
int $1_equal_fmpz(const $1_t _ip1_, const fmpz_t _ip2_)dnl
}}})dnl
define({{{desc_equal}}},{{{
    Returns nonzero if `_equal_(_ip1_, _ip2_)`, otherwise returns zero.
}}})dnl
dnl############################################################################
dnl equality to common constants
dnl############################################################################
define({{{func_is_zero}}},dnl
{{{dnl
int $1_is_zero(const $1_t _ip1_)dnl
}}})dnl
define({{{desc_is_zero}}},{{{
    Returns nonzero if `_equal_(_ip1_, 0)`, otherwise returns zero.dnl
}}})dnl
define({{{func_is_one}}},dnl
{{{dnl
int $1_is_one(const $1_t _ip1_)dnl
}}})dnl
define({{{desc_is_one}}},{{{
    Returns nonzero if `_equal_(_ip1_, 1)`, otherwise returns zero.dnl
}}})dnl
