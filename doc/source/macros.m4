dnl Change quotes to something that is not used anywhere
changequote({{{{,}}}})dnl
dnl
dnl############################################################################
dnl function prefix
dnl############################################################################
define({{{__function_prefix}}},{{{.. function::}}})dnl
dnl############################################################################
dnl addition
dnl############################################################################
define({{{func_add}}},dnl
{{{__function_prefix dnl
void $1_add($1_t r, const $1_t a, const $1_t b)dnl
}}})dnl
define({{{func_add_si}}},dnl
{{{__function_prefix dnl
void $1_add_si($1_t r, const $1_t a, slong b)dnl
}}})dnl
define({{{func_add_ui}}},dnl
{{{__function_prefix dnl
void $1_add_ui($1_t r, const $1_t a, ulong b)dnl
}}})dnl
define({{{desc_add}}},{{{dnl
    Sets `r` to `a + b`.dnl
}}})dnl
dnl############################################################################
dnl subtraction
dnl############################################################################
define({{{func_sub}}},dnl
{{{__function_prefix dnl
void $1_sub($1_t r, const $1_t a, const $1_t b)dnl
}}})dnl
define({{{func_sub_si}}},dnl
{{{__function_prefix dnl
void $1_sub_si($1_t r, const $1_t a, slong b)dnl
}}})dnl
define({{{func_sub_ui}}},dnl
{{{__function_prefix dnl
void $1_sub_ui($1_t r, const $1_t a, ulong b)dnl
}}})dnl
define({{{desc_sub}}},{{{dnl
    Sets `r` to `a - b`.dnl
}}})dnl
