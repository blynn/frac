#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "cf.h"

// sin 1 from Taylor series
static void *sin1_expansion(cf_t cf) {
  mpz_t p0, q0;  // p0/q0 = last convergent
  mpz_t p1, q1;  // p1/q1 = convergent
  mpz_t r0, s0;  // r0/s0 = lower bound for sin 1
  mpz_t r1, s1;  // r1/s1 = upper bound for sin 1
  mpq_t lim0, lim1;  // Related to error of lower and upper bounds:
                     //   error = 1 / (lim)
  mpz_t t0, t1, t2;  // Temporary variables.
  mpq_t tq0, tq1;

  mpz_init(r0); mpz_init(r1); mpz_init(s0); mpz_init(s1);
  mpz_init(p0); mpz_init(p1); mpz_init(q0); mpz_init(q1);
  mpz_init(t0); mpz_init(t1); mpz_init(t2);
  mpq_init(lim0); mpq_init(lim1);
  mpq_init(tq0); mpq_init(tq1);

  // Zero causes all sorts of problems.
  // Output it and forget about it.
  cf_put_int(cf, 0);

  // Initialize convergents.
  mpz_set_ui(p0, 1);
  mpz_set_ui(q1, 1);

  // sin 1 = 1 - 1/6 + ...
  mpz_set_ui(r0, 5);
  mpz_set_ui(s0, 6);
  mpz_set_ui(r1, 1);
  mpz_set_ui(s1, 1);
  mpq_set_ui(lim0, 1, 120);
  mpq_set_ui(lim1, 1, 6);

  int k = 1; // Number of Taylor series terms we've computed,
             // numbered from 0, unlike continued fraction terms.
  int i = 2; // Number of convergent we're computing.

  void increase_limit() {
    k++;
    mpz_mul_ui(r1, r0, 2 * k);
    mpz_mul_ui(r1, r1, 2 * k + 1);
    mpz_mul_ui(s1, s0, 2 * k);
    mpz_mul_ui(s1, s1, 2 * k + 1);
    mpz_add_ui(r1, r1, 1);
    k++;
    mpz_mul_ui(r0, r1, 2 * k);
    mpz_mul_ui(r0, r0, 2 * k + 1);
    mpz_mul_ui(s0, s1, 2 * k);
    mpz_mul_ui(s0, s0, 2 * k + 1);
    mpz_sub_ui(r0, r0, 1);
    mpz_set(mpq_denref(lim1), s0);
    mpz_mul_ui(mpq_denref(lim0), mpq_denref(lim1), 2 * k + 2);
    mpz_mul_ui(mpq_denref(lim0), mpq_denref(lim0), 2 * k + 3);
    gmp_printf("lower: %Zd / %Zd, %Qd\n", r0, s0, lim0);
    gmp_printf("upper: %Zd / %Zd, %Qd\n", r1, s1, lim1);
  }
  while (cf_wait(cf)) {
    for(;;) {
      printf("find x\n");
      // Find largest x such that (p1 x + p0)/(q1 x + q0)
      // <= r0/s0 for odd i, and >= r1/s1 for even i
      //  =>  x < (q0 r - p0 s)/(p1 s - q1 r)
      // for appropriate choice of r, s.
      mpz_ptr r, s;
      mpq_ptr lim;
      if (i & 1) {
	r = r0; s = s0; lim = lim0;
      } else {
	r = r1; s = s1; lim = lim1;
      }
      mpz_mul(t0, q0, r);
      mpz_mul(t1, p0, s);
      mpz_sub(t2, t0, t1);
      mpz_mul(t0, p1, s);
      mpz_mul(t1, q1, r);
      mpz_sub(t1, t0, t1);
      mpz_fdiv_q(t0, t2, t1);
      gmp_printf("soln: %Zd = %Zd/%Zd\n", t0, t2, t1);

      // Then x is the next convergent provided it is within limits.
      if (!mpz_sgn(t0)) {
	increase_limit();
      } else {
	// Check if t1/t2 is within the convergent bound
	mpz_mul(t1, p1, t0);
	mpz_add(t1, t1, p0);
	mpz_mul(t2, q1, t0);
	mpz_add(t2, t2, q0);

	mpq_set_num(tq0, t1);
	mpq_set_den(tq0, t2);
	mpq_set_num(tq1, r);
	mpq_set_den(tq1, s);
	gmp_printf("check: %Qd - %Qd\n", tq0, tq1);
	mpq_sub(tq0, tq0, tq1);
	if (mpq_sgn(tq0) < 0) {
	  mpq_neg(tq0, tq0);
	}
	// tq1 = 1 / 2 t2^2
	mpz_set_ui(mpq_numref(tq1), 1);
	mpz_mul(mpq_denref(tq1), t2, t2);
	mpz_mul_ui(mpq_denref(tq1), mpq_denref(tq1), 2);

	gmp_printf("check: %Qd + %Qd\n", tq0, lim);
	mpq_add(tq0, tq0, lim);
	gmp_printf("%Qd > %Qd?\n", tq0, tq1);
	if (mpq_cmp(tq0, tq1) >= 0) {
	  increase_limit();
	} else break;
      }
    };
    gmp_printf("term: %Zd\n", t0);
    cf_put(cf, t0);
    i++;
    mpz_set(p0, p1);
    mpz_set(q0, q1);
    mpz_set(p1, t1);
    mpz_set(q1, t2);
  }
  return NULL;
}

int main(int argc, char *argv[]) {
  cf_t x = cf_new_const(sin1_expansion);
  cf_t conv = cf_new_convergent(x);
  mpz_t z;
  mpz_init(z);
  for (int i = 0; i < 20; i++) {
    cf_get(z, conv); gmp_printf("%d: %Zd/", i, z);
    cf_get(z, conv); gmp_printf("%Zd\n", z);
  }
  return 0;
}
