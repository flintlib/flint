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

.global	FUNC(flint_mpn_mulhigh_12)
.p2align	4, 0x90
TYPE(flint_mpn_mulhigh_12)

FUNC(flint_mpn_mulhigh_12):
	.cfi_startproc
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14
	push	%r15
	vmovq	%rsp, %xmm0
	vmovq	%rdi, %xmm1

	mov	%rdx, %rcx
	mov	0*8(%rdx), %rdx

	mulx	10*8(%rsi), %rsp, %rax
	mulx	11*8(%rsi), %rsp, %r8
	add	%rsp, %rax
	adc	$0, %r8
	test	%al, %al

	mov	1*8(%rcx), %rdx
	mulx	9*8(%rsi), %rsp, %r10
	mulx	10*8(%rsi), %rsp, %r9
	adcx	%r10, %rax
	adox	%rsp, %rax
	adcx	%r9, %r8
	mulx	11*8(%rsi), %rsp, %r9
	adox	%rsp, %r8
	mov	$0, %esp
	adcx	%rsp, %r9
	adox	%rsp, %r9

	mov	2*8(%rcx), %rdx
	mulx	8*8(%rsi), %rsp, %r11
	mulx	9*8(%rsi), %rsp, %r10
	adcx	%r11, %rax
	adox	%rsp, %rax
	adcx	%r10, %r8
	mulx	10*8(%rsi), %rsp, %r10
	adox	%rsp, %r8
	adcx	%r10, %r9
	mulx	11*8(%rsi), %rsp, %r10
	adox	%rsp, %r9
	mov	$0, %esp
	adcx	%rsp, %r10
	adox	%rsp, %r10

	mov	3*8(%rcx), %rdx
	mulx	7*8(%rsi), %rsp, %rbx
	mulx	8*8(%rsi), %rsp, %r11
	adcx	%rbx, %rax
	adox	%rsp, %rax
	adcx	%r11, %r8
	mulx	9*8(%rsi), %rsp, %r11
	adox	%rsp, %r8
	adcx	%r11, %r9
	mulx	10*8(%rsi), %rsp, %r11
	adox	%rsp, %r9
	adcx	%r11, %r10
	mulx	11*8(%rsi), %rsp, %r11
	adox	%rsp, %r10
	mov	$0, %esp
	adcx	%rsp, %r11
	adox	%rsp, %r11

	mov	4*8(%rcx), %rdx
	mulx	6*8(%rsi), %rsp, %rbp
	mulx	7*8(%rsi), %rsp, %rbx
	adcx	%rbp, %rax
	adox	%rsp, %rax
	adcx	%rbx, %r8
	mulx	8*8(%rsi), %rsp, %rbx
	adox	%rsp, %r8
	adcx	%rbx, %r9
	mulx	9*8(%rsi), %rsp, %rbx
	adox	%rsp, %r9
	adcx	%rbx, %r10
	mulx	10*8(%rsi), %rsp, %rbx
	adox	%rsp, %r10
	adcx	%rbx, %r11
	mulx	11*8(%rsi), %rsp, %rbx
	adox	%rsp, %r11
	mov	$0, %esp
	adcx	%rsp, %rbx
	adox	%rsp, %rbx

	mov	5*8(%rcx), %rdx
	mulx	5*8(%rsi), %rsp, %r12
	mulx	6*8(%rsi), %rsp, %rbp
	adcx	%r12, %rax
	adox	%rsp, %rax
	adcx	%rbp, %r8
	mulx	7*8(%rsi), %rsp, %rbp
	adox	%rsp, %r8
	adcx	%rbp, %r9
	mulx	8*8(%rsi), %rsp, %rbp
	adox	%rsp, %r9
	adcx	%rbp, %r10
	mulx	9*8(%rsi), %rsp, %rbp
	adox	%rsp, %r10
	adcx	%rbp, %r11
	mulx	10*8(%rsi), %rsp, %rbp
	adox	%rsp, %r11
	adcx	%rbp, %rbx
	mulx	11*8(%rsi), %rsp, %rbp
	adox	%rsp, %rbx
	mov	$0, %esp
	adcx	%rsp, %rbp
	adox	%rsp, %rbp

	mov	6*8(%rcx), %rdx
	mulx	4*8(%rsi), %rsp, %r13
	mulx	5*8(%rsi), %rsp, %r12
	adcx	%r13, %rax
	adox	%rsp, %rax
	adcx	%r12, %r8
	mulx	6*8(%rsi), %rsp, %r12
	adox	%rsp, %r8
	adcx	%r12, %r9
	mulx	7*8(%rsi), %rsp, %r12
	adox	%rsp, %r9
	adcx	%r12, %r10
	mulx	8*8(%rsi), %rsp, %r12
	adox	%rsp, %r10
	adcx	%r12, %r11
	mulx	9*8(%rsi), %rsp, %r12
	adox	%rsp, %r11
	adcx	%r12, %rbx
	mulx	10*8(%rsi), %rsp, %r12
	adox	%rsp, %rbx
	adcx	%r12, %rbp
	mulx	11*8(%rsi), %rsp, %r12
	adox	%rsp, %rbp
	mov	$0, %esp
	adcx	%rsp, %r12
	adox	%rsp, %r12

	mov	7*8(%rcx), %rdx
	mulx	3*8(%rsi), %rsp, %r14
	mulx	4*8(%rsi), %rsp, %r13
	adcx	%r14, %rax
	adox	%rsp, %rax
	adcx	%r13, %r8
	mulx	5*8(%rsi), %rsp, %r13
	adox	%rsp, %r8
	adcx	%r13, %r9
	mulx	6*8(%rsi), %rsp, %r13
	adox	%rsp, %r9
	adcx	%r13, %r10
	mulx	7*8(%rsi), %rsp, %r13
	adox	%rsp, %r10
	adcx	%r13, %r11
	mulx	8*8(%rsi), %rsp, %r13
	adox	%rsp, %r11
	adcx	%r13, %rbx
	mulx	9*8(%rsi), %rsp, %r13
	adox	%rsp, %rbx
	adcx	%r13, %rbp
	mulx	10*8(%rsi), %rsp, %r13
	adox	%rsp, %rbp
	adcx	%r13, %r12
	mulx	11*8(%rsi), %rsp, %r13
	adox	%rsp, %r12
	mov	$0, %esp
	adcx	%rsp, %r13
	adox	%rsp, %r13

	mov	8*8(%rcx), %rdx
	mulx	2*8(%rsi), %rsp, %r15
	mulx	3*8(%rsi), %rsp, %r14
	adcx	%r15, %rax
	adox	%rsp, %rax
	adcx	%r14, %r8
	mulx	4*8(%rsi), %rsp, %r14
	adox	%rsp, %r8
	adcx	%r14, %r9
	mulx	5*8(%rsi), %rsp, %r14
	adox	%rsp, %r9
	adcx	%r14, %r10
	mulx	6*8(%rsi), %rsp, %r14
	adox	%rsp, %r10
	adcx	%r14, %r11
	mulx	7*8(%rsi), %rsp, %r14
	adox	%rsp, %r11
	adcx	%r14, %rbx
	mulx	8*8(%rsi), %rsp, %r14
	adox	%rsp, %rbx
	adcx	%r14, %rbp
	mulx	9*8(%rsi), %rsp, %r14
	adox	%rsp, %rbp
	adcx	%r14, %r12
	mulx	10*8(%rsi), %rsp, %r14
	adox	%rsp, %r12
	adcx	%r14, %r13
	mulx	11*8(%rsi), %rsp, %r14
	adox	%rsp, %r13
	mov	$0, %esp
	adcx	%rsp, %r14
	adox	%rsp, %r14

	mov	9*8(%rcx), %rdx
	mulx	1*8(%rsi), %rsp, %rdi
	mulx	2*8(%rsi), %rsp, %r15
	adcx	%rdi, %rax
	adox	%rsp, %rax
	adcx	%r15, %r8
	mulx	3*8(%rsi), %rsp, %r15
	adox	%rsp, %r8
	adcx	%r15, %r9
	mulx	4*8(%rsi), %rsp, %r15
	adox	%rsp, %r9
	adcx	%r15, %r10
	mulx	5*8(%rsi), %rsp, %r15
	adox	%rsp, %r10
	adcx	%r15, %r11
	mulx	6*8(%rsi), %rsp, %r15
	adox	%rsp, %r11
	adcx	%r15, %rbx
	mulx	7*8(%rsi), %rsp, %r15
	adox	%rsp, %rbx
	adcx	%r15, %rbp
	mulx	8*8(%rsi), %rsp, %r15
	adox	%rsp, %rbp
	adcx	%r15, %r12
	mulx	9*8(%rsi), %rsp, %r15
	adox	%rsp, %r12
	adcx	%r15, %r13
	mulx	10*8(%rsi), %rsp, %r15
	adox	%rsp, %r13
	adcx	%r15, %r14
	mulx	11*8(%rsi), %rsp, %r15
	adox	%rsp, %r14
	mov	$0, %esp
	adcx	%rsp, %r15
	adox	%rsp, %r15

	mov	10*8(%rcx), %rdx
	mulx	0*8(%rsi), %rsp, %rdi
	adcx	%rdi, %rax
	mulx	1*8(%rsi), %rsp, %rdi
	adox	%rsp, %rax
	adcx	%rdi, %r8
	mulx	2*8(%rsi), %rsp, %rdi
	adox	%rsp, %r8
	adcx	%rdi, %r9
	mulx	3*8(%rsi), %rsp, %rdi
	adox	%rsp, %r9
	adcx	%rdi, %r10
	mulx	4*8(%rsi), %rsp, %rdi
	adox	%rsp, %r10
	adcx	%rdi, %r11
	mulx	5*8(%rsi), %rsp, %rdi
	adox	%rsp, %r11
	adcx	%rdi, %rbx
	mulx	6*8(%rsi), %rsp, %rdi
	adox	%rsp, %rbx
	adcx	%rdi, %rbp
	mulx	7*8(%rsi), %rsp, %rdi
	adox	%rsp, %rbp
	adcx	%rdi, %r12
	mulx	8*8(%rsi), %rsp, %rdi
	adox	%rsp, %r12
	adcx	%rdi, %r13
	mulx	9*8(%rsi), %rsp, %rdi
	adox	%rsp, %r13
	adcx	%rdi, %r14
	mulx	10*8(%rsi), %rsp, %rdi
	adox	%rsp, %r14
	adcx	%rdi, %r15
	mulx	11*8(%rsi), %rsp, %rdi
	adox	%rsp, %r15
	mov	$0, %esp
	adcx	%rsp, %rdi
	adox	%rsp, %rdi

	mov	11*8(%rcx), %rdx
	mulx	0*8(%rsi), %rsp, %rcx
	adcx	%rsp, %rax
	adcx	%rcx, %r8
	mulx	1*8(%rsi), %rsp, %rcx
	adox	%rsp, %r8
	adox	%rcx, %r9
	mulx	2*8(%rsi), %rsp, %rcx
	adcx	%rsp, %r9
	adcx	%rcx, %r10
	mulx	3*8(%rsi), %rsp, %rcx
	adox	%rsp, %r10
	adox	%rcx, %r11
	mulx	4*8(%rsi), %rsp, %rcx
	adcx	%rsp, %r11
	adcx	%rcx, %rbx
	mulx	5*8(%rsi), %rsp, %rcx
	adox	%rsp, %rbx
	adox	%rcx, %rbp
	mulx	6*8(%rsi), %rsp, %rcx
	adcx	%rsp, %rbp
	adcx	%rcx, %r12
	mulx	7*8(%rsi), %rsp, %rcx
	adox	%rsp, %r12
	adox	%rcx, %r13
	mulx	8*8(%rsi), %rsp, %rcx
	adcx	%rsp, %r13
	adcx	%rcx, %r14
	mulx	9*8(%rsi), %rsp, %rcx
	adox	%rsp, %r14
	adox	%rcx, %r15
	mulx	10*8(%rsi), %rsp, %rcx
	adcx	%rsp, %r15
	adcx	%rcx, %rdi
	mulx	11*8(%rsi), %rsp, %rcx
	adox	%rsp, %rdi
	mov	$0, %esp
	adcx	%rsp, %rcx
	adox	%rsp, %rcx

	vmovq	%xmm1, %rsp
	mov	%r8, 0*8(%rsp)
	mov	%r9, 1*8(%rsp)
	mov	%r10, 2*8(%rsp)
	mov	%r11, 3*8(%rsp)
	mov	%rbx, 4*8(%rsp)
	mov	%rbp, 5*8(%rsp)
	mov	%r12, 6*8(%rsp)
	mov	%r13, 7*8(%rsp)
	mov	%r14, 8*8(%rsp)
	mov	%r15, 9*8(%rsp)
	mov	%rdi, 10*8(%rsp)
	mov	%rcx, 11*8(%rsp)

	vmovq	%xmm0, %rsp
	pop	%r15
	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
.flint_mpn_mulhigh_12_end:
SIZE(flint_mpn_mulhigh_12, .flint_mpn_mulhigh_12_end)
.cfi_endproc
