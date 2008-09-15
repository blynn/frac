#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "cf.h"

// sin 1 from Taylor series
static void *sin1_expansion(cf_t cf) {
  mpz_t p0, q0;  // p0/q0 = last convergent
  mpz_t p1, q1;  // p1/q1 = convergent
  mpq_t r0;      // lower bound for sin 1
  mpq_t r1;      // upper bound for sin 1
  mpz_t t0, t1, t2;  // Temporary variables.
  mpq_t tq0;

  mpq_init(r0); mpq_init(r1);
  mpz_init(p0); mpz_init(p1); mpz_init(q0); mpz_init(q1);
  mpz_init(t0); mpz_init(t1); mpz_init(t2);
  mpq_init(tq0);

  // Zero causes all sorts of problems.
  // Output it and forget about it.
  cf_put_int(cf, 0);

  // Initialize convergents.
  mpz_set_ui(p0, 1);
  mpz_set_ui(q1, 1);

  // sin 1 = 1 - 1/6 + ...
  mpq_set_ui(r0, 5, 6);
  mpq_set_ui(r1, 1, 1);

  int k = 1; // Number of Taylor series terms we've computed,
             // numbered from 0, unlike continued fraction terms.
  int i = 2; // Number of convergent we're computing.

  void increase_limit() {
    k++;
    mpz_mul_ui(mpq_numref(r1), mpq_numref(r0), (2 * k) * (2 * k + 1));
    mpz_mul_ui(mpq_denref(r1), mpq_denref(r0), (2 * k) * (2 * k + 1));
    mpz_add_ui(mpq_numref(r1), mpq_numref(r1), 1);
    k++;
    mpz_mul_ui(mpq_numref(r0), mpq_numref(r1), (2 * k) * (2 * k + 1));
    mpz_mul_ui(mpq_denref(r0), mpq_denref(r1), (2 * k) * (2 * k + 1));
    mpz_sub_ui(mpq_numref(r0), mpq_numref(r0), 1);
  }
  while (cf_wait(cf)) {
    for(;;) {
      printf("find x\n");
      // Find largest x such that (p1 x + p0)/(q1 x + q0)
      // <= r0 for odd i, and >= r1 for even i
      //  =>  x < (q0 num(r) - p0 den(r))/(p1 den(r) - q1 num(r))
      // for the appropriate r.
      mpq_ptr r = i & 1 ? r0 : r1;
      mpz_mul(t0, q0, mpq_numref(r));
      mpz_mul(t1, p0, mpq_denref(r));
      mpz_sub(t2, t0, t1);
      mpz_mul(t0, p1, mpq_denref(r));
      mpz_mul(t1, q1, mpq_numref(r));
      mpz_sub(t1, t0, t1);
      mpz_fdiv_q(t0, t2, t1);

      // Then x is the next convergent provided substituting x + 1 overshoots
      // the other bound. (If not, we need better bounds to decide.)
      if (!mpz_sgn(t0)) {
	increase_limit();
      } else {
	mpz_add_ui(t1, t0, 1);
	mpz_set(mpq_numref(tq0), p0);
	mpz_addmul(mpq_numref(tq0), p1, t1);
	mpz_set(mpq_denref(tq0), q0);
	mpz_addmul(mpq_denref(tq0), q1, t1);
	if (i & 1) {
	  if (mpq_cmp(tq0, r1) > 0) break;
	  // Check (p1(x+1) + p0)/(q1(x+1) + q0) > r1
	} else {
	  if (mpq_cmp(tq0, r0) < 0) break;
	  // Check (p1(x+1) + p0)/(q1(x+1) + q0) < r0
	}
	increase_limit();
      }
    };
    cf_put(cf, t0);
    i++;
    mpz_addmul(p0, p1, t0);
    mpz_addmul(q0, q1, t0);
    mpz_swap(p1, p0);
    mpz_swap(q1, q0);
  }

  mpq_clear(r0); mpq_clear(r1);
  mpz_clear(p0); mpz_clear(p1); mpz_clear(q0); mpz_clear(q1);
  mpz_clear(t0); mpz_clear(t1); mpz_clear(t2);
  mpq_clear(tq0);
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
