#
#   Copyright (C) 2024 Albin Ahlbäck
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

.global	FUNC(flint_mpn_mulhigh_10)
.p2align	4, 0x90
TYPE(flint_mpn_mulhigh_10)

FUNC(flint_mpn_mulhigh_10):
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
	xor	%esp, %esp

	mulx	8*8(%rsi), %r8, %rax
	mulx	9*8(%rsi), %r8, %r9
	adcx	%r8, %rax
	adcx	%rsp, %r9

	mov	1*8(%rcx), %rdx
	mulx	7*8(%rsi), %r8, %r11
	mulx	8*8(%rsi), %r8, %r10
	adcx	%r11, %rax
	adox	%r8, %rax
	adcx	%r10, %r9
	mulx	9*8(%rsi), %r8, %r10
	adox	%r8, %r9
	adcx	%rsp, %r10
	adox	%rsp, %r10

	mov	2*8(%rcx), %rdx
	mulx	6*8(%rsi), %r8, %rbx
	mulx	7*8(%rsi), %r8, %r11
	adcx	%rbx, %rax
	adox	%r8, %rax
	adcx	%r11, %r9
	mulx	8*8(%rsi), %r8, %r11
	adox	%r8, %r9
	adcx	%r11, %r10
	mulx	9*8(%rsi), %r8, %r11
	adox	%r8, %r10
	adcx	%rsp, %r11
	adox	%rsp, %r11

	mov	3*8(%rcx), %rdx
	mulx	5*8(%rsi), %r8, %rbp
	mulx	6*8(%rsi), %r8, %rbx
	adcx	%rbp, %rax
	adox	%r8, %rax
	adcx	%rbx, %r9
	mulx	7*8(%rsi), %r8, %rbx
	adox	%r8, %r9
	adcx	%rbx, %r10
	mulx	8*8(%rsi), %r8, %rbx
	adox	%r8, %r10
	adcx	%rbx, %r11
	mulx	9*8(%rsi), %r8, %rbx
	adox	%r8, %r11
	adcx	%rsp, %rbx
	adox	%rsp, %rbx

	mov	4*8(%rcx), %rdx
	mulx	4*8(%rsi), %r8, %r12
	mulx	5*8(%rsi), %r8, %rbp
	adcx	%r12, %rax
	adox	%r8, %rax
	adcx	%rbp, %r9
	mulx	6*8(%rsi), %r8, %rbp
	adox	%r8, %r9
	adcx	%rbp, %r10
	mulx	7*8(%rsi), %r8, %rbp
	adox	%r8, %r10
	adcx	%rbp, %r11
	mulx	8*8(%rsi), %r8, %rbp
	adox	%r8, %r11
	adcx	%rbp, %rbx
	mulx	9*8(%rsi), %r8, %rbp
	adox	%r8, %rbx
	adcx	%rsp, %rbp
	adox	%rsp, %rbp

	mov	5*8(%rcx), %rdx
	mulx	3*8(%rsi), %r8, %r13
	mulx	4*8(%rsi), %r8, %r12
	adcx	%r13, %rax
	adox	%r8, %rax
	adcx	%r12, %r9
	mulx	5*8(%rsi), %r8, %r12
	adox	%r8, %r9
	adcx	%r12, %r10
	mulx	6*8(%rsi), %r8, %r12
	adox	%r8, %r10
	adcx	%r12, %r11
	mulx	7*8(%rsi), %r8, %r12
	adox	%r8, %r11
	adcx	%r12, %rbx
	mulx	8*8(%rsi), %r8, %r12
	adox	%r8, %rbx
	adcx	%r12, %rbp
	mulx	9*8(%rsi), %r8, %r12
	adox	%r8, %rbp
	adcx	%rsp, %r12
	adox	%rsp, %r12

	mov	6*8(%rcx), %rdx
	mulx	2*8(%rsi), %r8, %r14
	mulx	3*8(%rsi), %r8, %r13
	adcx	%r14, %rax
	adox	%r8, %rax
	adcx	%r13, %r9
	mulx	4*8(%rsi), %r8, %r13
	adox	%r8, %r9
	adcx	%r13, %r10
	mulx	5*8(%rsi), %r8, %r13
	adox	%r8, %r10
	adcx	%r13, %r11
	mulx	6*8(%rsi), %r8, %r13
	adox	%r8, %r11
	adcx	%r13, %rbx
	mulx	7*8(%rsi), %r8, %r13
	adox	%r8, %rbx
	adcx	%r13, %rbp
	mulx	8*8(%rsi), %r8, %r13
	adox	%r8, %rbp
	adcx	%r13, %r12
	mulx	9*8(%rsi), %r8, %r13
	adox	%r8, %r12
	adcx	%rsp, %r13
	adox	%rsp, %r13

	mov	7*8(%rcx), %rdx
	mulx	1*8(%rsi), %r8, %r15
	mulx	2*8(%rsi), %r8, %r14
	adcx	%r15, %rax
	adox	%r8, %rax
	adcx	%r14, %r9
	mulx	3*8(%rsi), %r8, %r14
	adox	%r8, %r9
	adcx	%r14, %r10
	mulx	4*8(%rsi), %r8, %r14
	adox	%r8, %r10
	adcx	%r14, %r11
	mulx	5*8(%rsi), %r8, %r14
	adox	%r8, %r11
	adcx	%r14, %rbx
	mulx	6*8(%rsi), %r8, %r14
	adox	%r8, %rbx
	adcx	%r14, %rbp
	mulx	7*8(%rsi), %r8, %r14
	adox	%r8, %rbp
	adcx	%r14, %r12
	mulx	8*8(%rsi), %r8, %r14
	adox	%r8, %r12
	adcx	%r14, %r13
	mulx	9*8(%rsi), %r8, %r14
	adox	%r8, %r13
	adcx	%rsp, %r14
	adox	%rsp, %r14

	mov	8*8(%rcx), %rdx
	mulx	0*8(%rsi), %r8, %rdi
	mulx	1*8(%rsi), %r8, %r15
	adcx	%rdi, %rax
	adox	%r8, %rax
	adcx	%r15, %r9
	mulx	2*8(%rsi), %r8, %r15
	adox	%r8, %r9
	adcx	%r15, %r10
	mulx	3*8(%rsi), %r8, %r15
	adox	%r8, %r10
	adcx	%r15, %r11
	mulx	4*8(%rsi), %r8, %r15
	adox	%r8, %r11
	adcx	%r15, %rbx
	mulx	5*8(%rsi), %r8, %r15
	adox	%r8, %rbx
	adcx	%r15, %rbp
	mulx	6*8(%rsi), %r8, %r15
	adox	%r8, %rbp
	adcx	%r15, %r12
	mulx	7*8(%rsi), %r8, %r15
	adox	%r8, %r12
	adcx	%r15, %r13
	mulx	8*8(%rsi), %r8, %r15
	adox	%r8, %r13
	adcx	%r15, %r14
	mulx	9*8(%rsi), %r8, %r15
	adox	%r8, %r14
	adcx	%rsp, %r15
	adox	%rsp, %r15

	mov	9*8(%rcx), %rdx
	mulx	0*8(%rsi), %r8, %rdi
	adcx	%r8, %rax
	adcx	%rdi, %r9
	mulx	1*8(%rsi), %r8, %rdi
	adox	%r8, %r9
	adox	%rdi, %r10
	mulx	2*8(%rsi), %r8, %rdi
	adcx	%r8, %r10
	adcx	%rdi, %r11
	mulx	3*8(%rsi), %r8, %rdi
	adox	%r8, %r11
	adox	%rdi, %rbx
	mulx	4*8(%rsi), %r8, %rdi
	adcx	%r8, %rbx
	adcx	%rdi, %rbp
	mulx	5*8(%rsi), %r8, %rdi
	adox	%r8, %rbp
	adox	%rdi, %r12
	mulx	6*8(%rsi), %r8, %rdi
	adcx	%r8, %r12
	adcx	%rdi, %r13
	mulx	7*8(%rsi), %r8, %rdi
	adox	%r8, %r13
	adox	%rdi, %r14
	mulx	8*8(%rsi), %r8, %rdi
	adcx	%r8, %r14
	adcx	%rdi, %r15
	mulx	9*8(%rsi), %r8, %rdi
	adox	%r8, %r15
	adcx	%rsp, %rdi
	adox	%rsp, %rdi

	vmovq	%xmm1, %r8
	vmovq	%xmm0, %rsp
	mov	%r9, 0*8(%r8)
	mov	%r10, 1*8(%r8)
	mov	%r11, 2*8(%r8)
	mov	%rbx, 3*8(%r8)
	mov	%rbp, 4*8(%r8)
	mov	%r12, 5*8(%r8)
	mov	%r13, 6*8(%r8)
	mov	%r14, 7*8(%r8)
	mov	%r15, 8*8(%r8)
	mov	%rdi, 9*8(%r8)

	pop	%r15
	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
.flint_mpn_mulhigh_10_end:
SIZE(flint_mpn_mulhigh_10, .flint_mpn_mulhigh_10_end)
.cfi_endproc
