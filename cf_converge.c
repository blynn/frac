#include <stdio.h>
#include <gmp.h>
#include "cf.h"

static void *convergent(cf_t cf) {
  cf_t input = cf_data(cf);
  mpz_t pold, p, pnew;
  mpz_t qold, q, qnew;
  mpz_init(pold); mpz_init(p); mpz_init(pnew);
  mpz_init(qold); mpz_init(q); mpz_init(qnew);
  mpz_set_ui(pold, 0); mpz_set_ui(p, 1);
  mpz_set_ui(qold, 1); mpz_set_ui(q, 0);
  mpz_t denom;
  mpz_init(denom);
  while(cf_wait(cf)) {
    cf_signal(input);
    cf_get(denom, input);
    mpz_mul(pnew, p, denom);
    mpz_add(pnew, pnew, pold);
    mpz_mul(qnew, q, denom);
    mpz_add(qnew, qnew, qold);

    cf_put(cf, pnew);
    cf_put(cf, qnew);
    mpz_set(pold, p); mpz_set(p, pnew);
    mpz_set(qold, q); mpz_set(q, qnew);
  }
  mpz_clear(denom);
  mpz_clear(pold); mpz_clear(p); mpz_clear(pnew);
  mpz_clear(qold); mpz_clear(q); mpz_clear(qnew);
  return NULL;
}
cf_t cf_new_convergent(cf_t a) {
  return cf_new(convergent, a);
}
