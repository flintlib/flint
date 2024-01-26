#
#   Copyright (C) 2024 Albin Ahlbäck
#
#   This file is part of FLINT.
#
#   FLINT is free software: you can redistribute it and/or modify it under
#   the terms of the GNU Lesser General Public License (LGPL) as published
#   by the Free Software Foundation; either version 2.1 of the License, or
#   (at your option) any later version.  See <https://www.gnu.org/licenses/>.
#

include(`config.m4')dnl
dnl
.text

.global	FUNC(flint_mpn_mulhigh_normalised_1)
.p2align	4, 0x90
TYPE(flint_mpn_mulhigh_normalised_1)

FUNC(flint_mpn_mulhigh_normalised_1):
	.cfi_startproc
	mov	0*8(%rdx), %rdx
	mulx	0*8(%rsi), %rax, %rcx
	mov	$0, %rdx
	test	%rcx, %rcx
	setns	%dl
	js	.Lcontinue
	add	%rax, %rax
	adc	%rcx, %rcx
.Lcontinue:
	mov	%rcx, 0*8(%rdi)

	ret
.flint_mpn_mulhigh_normalised_1_end:
SIZE(flint_mpn_mulhigh_normalised_1, .flint_mpn_mulhigh_normalised_1_end)
.cfi_endproc
