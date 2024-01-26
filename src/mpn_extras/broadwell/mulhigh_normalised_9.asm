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

.global	FUNC(flint_mpn_mulhigh_normalised_9)
.p2align	4, 0x90
TYPE(flint_mpn_mulhigh_normalised_9)

FUNC(flint_mpn_mulhigh_normalised_9):
	.cfi_startproc
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14
	push	%r15
	push	%rdi

	mov	%rdx, %rcx
	mov	0*8(%rdx), %rdx
	xor	%r9d, %r9d

	mulx	7*8(%rsi), %r8, %rax
	mulx	8*8(%rsi), %r8, %r10
	adcx	%r8, %rax
	adcx	%r9, %r10

	mov	1*8(%rcx), %rdx
	mulx	6*8(%rsi), %r8, %rbx
	mulx	7*8(%rsi), %r8, %r11
	adcx	%rbx, %rax
	adox	%r8, %rax
	adcx	%r11, %r10
	mulx	8*8(%rsi), %r8, %r11
	adox	%r8, %r10
	adcx	%r9, %r11
	adox	%r9, %r11

	mov	2*8(%rcx), %rdx
	mulx	5*8(%rsi), %r8, %rbp
	mulx	6*8(%rsi), %r8, %rbx
	adcx	%rbp, %rax
	adox	%r8, %rax
	adcx	%rbx, %r10
	mulx	7*8(%rsi), %r8, %rbx
	adox	%r8, %r10
	adcx	%rbx, %r11
	mulx	8*8(%rsi), %r8, %rbx
	adox	%r8, %r11
	adcx	%r9, %rbx
	adox	%r9, %rbx

	mov	3*8(%rcx), %rdx
	mulx	4*8(%rsi), %r8, %r12
	mulx	5*8(%rsi), %r8, %rbp
	adcx	%r12, %rax
	adox	%r8, %rax
	adcx	%rbp, %r10
	mulx	6*8(%rsi), %r8, %rbp
	adox	%r8, %r10
	adcx	%rbp, %r11
	mulx	7*8(%rsi), %r8, %rbp
	adox	%r8, %r11
	adcx	%rbp, %rbx
	mulx	8*8(%rsi), %r8, %rbp
	adox	%r8, %rbx
	adcx	%r9, %rbp
	adox	%r9, %rbp

	mov	4*8(%rcx), %rdx
	mulx	3*8(%rsi), %r8, %r13
	mulx	4*8(%rsi), %r8, %r12
	adcx	%r13, %rax
	adox	%r8, %rax
	adcx	%r12, %r10
	mulx	5*8(%rsi), %r8, %r12
	adox	%r8, %r10
	adcx	%r12, %r11
	mulx	6*8(%rsi), %r8, %r12
	adox	%r8, %r11
	adcx	%r12, %rbx
	mulx	7*8(%rsi), %r8, %r12
	adox	%r8, %rbx
	adcx	%r12, %rbp
	mulx	8*8(%rsi), %r8, %r12
	adox	%r8, %rbp
	adcx	%r9, %r12
	adox	%r9, %r12

	mov	5*8(%rcx), %rdx
	mulx	2*8(%rsi), %r8, %r14
	mulx	3*8(%rsi), %r8, %r13
	adcx	%r14, %rax
	adox	%r8, %rax
	adcx	%r13, %r10
	mulx	4*8(%rsi), %r8, %r13
	adox	%r8, %r10
	adcx	%r13, %r11
	mulx	5*8(%rsi), %r8, %r13
	adox	%r8, %r11
	adcx	%r13, %rbx
	mulx	6*8(%rsi), %r8, %r13
	adox	%r8, %rbx
	adcx	%r13, %rbp
	mulx	7*8(%rsi), %r8, %r13
	adox	%r8, %rbp
	adcx	%r13, %r12
	mulx	8*8(%rsi), %r8, %r13
	adox	%r8, %r12
	adcx	%r9, %r13
	adox	%r9, %r13

	mov	6*8(%rcx), %rdx
	mulx	1*8(%rsi), %r8, %r15
	mulx	2*8(%rsi), %r8, %r14
	adcx	%r15, %rax
	adox	%r8, %rax
	adcx	%r14, %r10
	mulx	3*8(%rsi), %r8, %r14
	adox	%r8, %r10
	adcx	%r14, %r11
	mulx	4*8(%rsi), %r8, %r14
	adox	%r8, %r11
	adcx	%r14, %rbx
	mulx	5*8(%rsi), %r8, %r14
	adox	%r8, %rbx
	adcx	%r14, %rbp
	mulx	6*8(%rsi), %r8, %r14
	adox	%r8, %rbp
	adcx	%r14, %r12
	mulx	7*8(%rsi), %r8, %r14
	adox	%r8, %r12
	adcx	%r14, %r13
	mulx	8*8(%rsi), %r8, %r14
	adox	%r8, %r13
	adcx	%r9, %r14
	adox	%r9, %r14

	mov	7*8(%rcx), %rdx
	mulx	0*8(%rsi), %r8, %rdi
	mulx	1*8(%rsi), %r8, %r15
	adcx	%rdi, %rax
	adox	%r8, %rax
	adcx	%r15, %r10
	mulx	2*8(%rsi), %r8, %r15
	adox	%r8, %r10
	adcx	%r15, %r11
	mulx	3*8(%rsi), %r8, %r15
	adox	%r8, %r11
	adcx	%r15, %rbx
	mulx	4*8(%rsi), %r8, %r15
	adox	%r8, %rbx
	adcx	%r15, %rbp
	mulx	5*8(%rsi), %r8, %r15
	adox	%r8, %rbp
	adcx	%r15, %r12
	mulx	6*8(%rsi), %r8, %r15
	adox	%r8, %r12
	adcx	%r15, %r13
	mulx	7*8(%rsi), %r8, %r15
	adox	%r8, %r13
	adcx	%r15, %r14
	mulx	8*8(%rsi), %r8, %r15
	adox	%r8, %r14
	adcx	%r9, %r15
	adox	%r9, %r15

	mov	8*8(%rcx), %rdx
	mulx	0*8(%rsi), %r8, %rdi
	adcx	%r8, %rax
	adcx	%rdi, %r10
	mulx	1*8(%rsi), %r8, %rdi
	adox	%r8, %r10
	adox	%rdi, %r11
	mulx	2*8(%rsi), %r8, %rdi
	adcx	%r8, %r11
	adcx	%rdi, %rbx
	mulx	3*8(%rsi), %r8, %rdi
	adox	%r8, %rbx
	adox	%rdi, %rbp
	mulx	4*8(%rsi), %r8, %rdi
	adcx	%r8, %rbp
	adcx	%rdi, %r12
	mulx	5*8(%rsi), %r8, %rdi
	adox	%r8, %r12
	adox	%rdi, %r13
	mulx	6*8(%rsi), %r8, %rdi
	adcx	%r8, %r13
	adcx	%rdi, %r14
	mulx	7*8(%rsi), %r8, %rdi
	adox	%r8, %r14
	adox	%rdi, %r15
	mulx	8*8(%rsi), %r8, %rdi
	adcx	%r8, %r15
	adcx	%r9, %rdi
	pop	%r8
	adox	%r9, %rdi

	mov	$0, %rdx
	test	%rdi, %rdi
	setns	%dl
	js	.Lcontinue
	add	%rax, %rax
	adc	%r10, %r10
	adc	%r11, %r11
	adc	%rbx, %rbx
	adc	%rbp, %rbp
	adc	%r12, %r12
	adc	%r13, %r13
	adc	%r14, %r14
	adc	%r15, %r15
	adc	%rdi, %rdi
.Lcontinue:
	mov	%r10, 0*8(%r8)
	mov	%r11, 1*8(%r8)
	mov	%rbx, 2*8(%r8)
	mov	%rbp, 3*8(%r8)
	mov	%r12, 4*8(%r8)
	mov	%r13, 5*8(%r8)
	mov	%r14, 6*8(%r8)
	mov	%r15, 7*8(%r8)
	mov	%rdi, 8*8(%r8)

	pop	%r15
	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
.flint_mpn_mulhigh_normalised_9_end:
SIZE(flint_mpn_mulhigh_normalised_9, .flint_mpn_mulhigh_normalised_9_end)
.cfi_endproc
