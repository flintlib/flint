#
#   Copyright (C) 2023 Albin Ahlbäck
#
#   This file is part of FLINT.
#
#   FLINT is free software: you can redistribute it and/or modify it under
#   the terms of the GNU Lesser General Public License (LGPL) as published
#   by the Free Software Foundation; either version 3 of the License, or
#   (at your option) any later version.  See <https://www.gnu.org/licenses/>.
#

include(`config.m4')dnl
dnl
.text

.global	FUNC(flint_mpn_mul_2_2)
.p2align	4, 0x90
TYPE(flint_mpn_mul_2_2)

FUNC(flint_mpn_mul_2_2):
	.cfi_startproc
	mov	1*8(%rdx), %r10
	mov	0*8(%rdx), %rdx
	xor	%r11d, %r11d
	mulx	0*8(%rsi), %rcx, %r8
	mulx	1*8(%rsi), %rax, %r9
	adcx	%rax, %r8
	mov	%rcx, 0*8(%rdi)
	mov	%r10, %rdx
	mulx	0*8(%rsi), %rcx, %rax
	adox	%rcx, %r8
	adcx	%rax, %r9
	mov	%r8, 1*8(%rdi)
	mulx	1*8(%rsi), %rcx, %rax
	adox	%rcx, %r9
	adcx	%r11, %rax
	adox	%r11, %rax
	mov	%r9, 2*8(%rdi)
	mov	%rax, 3*8(%rdi)

	ret
.flint_mpn_mul_2_2_end:
SIZE(flint_mpn_mul_2_2, .flint_mpn_mul_2_2_end)
.cfi_endproc
