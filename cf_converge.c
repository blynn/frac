#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "cf.h"

struct mobius_data_s {
  cf_t input;
  mpz_t a, b, c, d;
};
typedef struct mobius_data_s *mobius_data_ptr;
static void *mobius(cf_t cf) {
  mobius_data_ptr md = cf_data(cf);
  cf_t input = md->input;

  mpz_t pold, p, pnew;
  mpz_t qold, q, qnew;
  mpz_init(pold); mpz_init(p); mpz_init(pnew);
  mpz_init(qold); mpz_init(q); mpz_init(qnew);
  mpz_set(pold, md->b); mpz_set(p, md->a);
  mpz_set(qold, md->d); mpz_set(q, md->c);
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

  mpz_clear(md->a);
  mpz_clear(md->b);
  mpz_clear(md->c);
  mpz_clear(md->d);
  free(md);
  return NULL;
}
cf_t cf_new_convergent_mobius(cf_t x, mpz_t a, mpz_t b, mpz_t c, mpz_t d) {
  mobius_data_ptr md = malloc(sizeof(*md));
  mpz_init(md->a);
  mpz_init(md->b);
  mpz_init(md->c);
  mpz_init(md->d);
  mpz_set(md->a, a);
  mpz_set(md->b, b);
  mpz_set(md->c, c);
  mpz_set(md->d, d);
  md->input = x;
  return cf_new(mobius, md);
}

cf_t cf_new_convergent(cf_t x) {
  mpz_t one, zero;
  mpz_init(one);
  mpz_init(zero);
  mpz_set_ui(one, 1);
  mpz_set_ui(zero, 0);
  return cf_new_convergent_mobius(x, one, zero, zero, one);
  mpz_clear(one);
  mpz_clear(zero);
}
