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

.global	FUNC(flint_mpn_sqr_4)
.p2align	4, 0x90
TYPE(flint_mpn_sqr_4)

FUNC(flint_mpn_sqr_4):
	.cfi_startproc
	mov	0*8(%rsi), %rdx
	push	%rbx
	push	%rbp
	xor	%ebp, %ebp
	mulx	1*8(%rsi), %rcx, %r11
	mulx	2*8(%rsi), %r8, %rbx
	mulx	3*8(%rsi), %r9, %r10
	mov	1*8(%rsi), %rdx
	adox	%r11, %r8
	adox	%rbx, %r9
	mulx	2*8(%rsi), %r11, %rbx
	adcx	%r11, %r9
	adcx	%rbx, %r10
	mulx	3*8(%rsi), %rbx, %r11
	mov	2*8(%rsi), %rdx
	adox	%rbx, %r10
	adox	%rbp, %r11
	mulx	3*8(%rsi), %rdx, %rbx
	adcx	%rdx, %r11
	adc	%rbp, %rbx
	mov	0*8(%rsi), %rdx
	add	%rcx, %rcx
	adc	%r8, %r8
	adc	%r9, %r9
	adc	%r10, %r10
	adc	%r11, %r11
	adc	%rbx, %rbx
	setc	%bpl
	mulx	%rdx, %rdx, %rax
	mov	%rdx, 0*8(%rdi)
	add	%rax, %rcx
	mov	1*8(%rsi), %rdx
	mov	%rcx, 1*8(%rdi)
	mulx	%rdx, %rdx, %rax
	adc	%rdx, %r8
	adc	%rax, %r9
	mov	2*8(%rsi), %rdx
	mov	%r8, 2*8(%rdi)
	mov	%r9, 3*8(%rdi)
	mulx	%rdx, %rcx, %rax
	adc	%rcx, %r10
	adc	%rax, %r11
	mov	3*8(%rsi), %rdx
	mov	%r10, 4*8(%rdi)
	mov	%r11, 5*8(%rdi)
	mulx	%rdx, %rcx, %rax
	adc	%rcx, %rbx
	adc	%rbp, %rax
	mov	%rbx, 6*8(%rdi)
	mov	%rax, 7*8(%rdi)
	pop	%rbp
	pop	%rbx

	ret
.flint_mpn_sqr_4_end:
SIZE(flint_mpn_sqr_4, .flint_mpn_sqr_4_end)
.cfi_endproc
