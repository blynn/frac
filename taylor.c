#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "cf.h"

// sin 1 from Taylor series
static void *sin1_expansion(cf_t cf) {
  int k = 1;

  mpz_t r, s;
  mpz_init(r); mpz_init(s);

  mpz_set_ui(r, 1);
  mpz_set_ui(s, 1);
  int i;
  for (i = 1; i <= k ; i++ ) {
    mpz_mul_ui(s, s, (2 * i) * (2 * i + 1));
    mpz_mul_ui(r, r, (2 * i) * (2 * i + 1));
    if (i & 1) {
      mpz_sub_ui(r, r, 1);
    } else {
      mpz_add_ui(r, r, 1);
    }
  }

  mpz_t limit;
  mpz_init(limit);
  // Compute limit = (2k + 1)! / 2
  mpz_mul_ui(limit, s, k);
  mpz_mul_ui(limit, s, 2 * k + 1);

  mpz_t p0, p1;
  mpz_t q0, q1;
  mpz_init(p0); mpz_init(p1);
  mpz_init(q0); mpz_init(q1);
  mpz_set_ui(p1, 1);
  mpz_set_ui(q0, 1);

  mpz_t a, b;
  mpz_init(a); mpz_init(b);
  mpz_t t0, t1;
  mpz_init(t0); mpz_init(t1);

  mpz_set(a, r);
  mpz_set(b, s);

  mpz_t t2;
  mpz_init(t2);
  // Pump out start of expansion.
  i = 0;
  for (;;) {
    mpz_fdiv_qr(t0, t1, a, b);

    mpz_set(t2, q0);
    mpz_addmul(q0, q1, t0);
    mpz_mul(a, q0, q0);
    if (mpz_cmp(a, limit) >= 0) {
      mpz_set(q0, t2);
      break;
    }
    mpz_addmul(p0, p1, t0);
    mpz_swap(p0, p1);
    mpz_swap(q0, q1);
    i++;
    // gmp_printf("%Zd (%Zd) -> %Zd / %Zd\n", t0, t1, p1, q1);
    cf_put(cf, t0);
    // if (!mpz_sgn(t1)) die("sin(1) is rational!\n");
    mpz_set(a, b);
    mpz_set(b, t1);
  }

  void increase_k() {
    k++;
    mpz_mul_ui(limit, limit, 2 * k);
    mpz_mul_ui(limit, limit, 2 * k + 1);
    mpz_mul_ui(s, s, (2 * k) * (2 * k + 1));
    mpz_mul_ui(r, r, (2 * k) * (2 * k + 1));
    if (k & 1) {
      mpz_sub_ui(r, r, 1);
    } else {
      mpz_add_ui(r, r, 1);
    }
gmp_printf("r/s = %Zd/%Zd\n", r, s);
  }
  void increase_limit() {
    if (i & 1) {
      // Last output was for a lower bound.
      // We need an upper bound of sin 1.
      increase_k();
      if (k & 1) increase_k();
      printf("lower\n");
    } else {
      // Last output was for an upper bound.
      printf("upper\n");
      // We need a lower bound of sin 1.
      increase_k();
      if (!(k & 1)) increase_k();
    }
  }

  increase_limit();
  while(cf_wait(cf)) {
    for(;;) {
      printf("find x\n");
      // Find largest x such that (p1 x + p0)/(q1 x + q0) R r/s
      // where R is < for even i, and > for odd i
      // i.e. (p1 s - q1 r)x R q0 r - p0 s
      mpz_mul(t0, q0, r);
      mpz_mul(t1, p0, s);
      mpz_sub(a, t0, t1);
      mpz_mul(t0, p1, s);
      mpz_mul(t1, q1, r);
      mpz_sub(b, t0, t1);
      mpz_fdiv_q(t0, a, b);

      // Then x is the next convergent provided it is within limits.
      mpz_set(t2, q0);
      mpz_addmul(q0, q1, t0);
      mpz_mul(a, q0, q0);
      gmp_printf("%Zd > %Zd?\n", a, limit);
      if (mpz_cmp(a, limit) >= 0) {
	increase_limit();
	mpz_set(q0, t2);
      } else break;
    };
    gmp_printf("term: %Zd\n", t0);
    cf_put(cf, t0);
    i++;
    mpz_addmul(p0, p1, t0);
    mpz_swap(p0, p1);
    mpz_swap(q0, q1);
  }

  mpz_clear(t2);
  return NULL;
}

int main(int argc, char *argv[]) {
  cf_t x = cf_new_const(sin1_expansion);
  cf_t conv = cf_new_convergent(x);
  mpz_t z;
  mpz_init(z);
  for (int i = 0; i < 10; i++) {
    cf_get(z, conv); gmp_printf("%d: %Zd/", i, z);
    cf_get(z, conv); gmp_printf("%Zd\n", z);
  }
  return 0;
}
