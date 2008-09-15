#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "cf.h"

// e = [2; 1, 2, 1, 1, 4, 1, ...]
static void *e_expansion(cf_t cf) {
  int even = 2;
  cf_put_int(cf, even);

  while(cf_wait(cf)) {
    cf_put_int(cf, 1);
    cf_put_int(cf, even);
    even += 2;
    cf_put_int(cf, 1);
  }
  return NULL;
}

cf_t cf_new_e() {
  return cf_new_const(e_expansion);
}

// 4/pi = 1 + 1/(3 + 4/(5 + 9/(7 + 16/(9 + ...))))
static void *pi_arctan_sequence(cf_t cf) {
  mpz_t num, denom;
  mpz_init(num);
  mpz_init(denom);

  mpz_set_ui(denom, 1);
  cf_put(cf, denom);
  mpz_set_ui(num, 1);
  cf_put(cf, num);
  while(cf_wait(cf)) {
    mpz_add_ui(denom, denom, 2);
    cf_put(cf, denom);
    mpz_add(num, num, denom);
    cf_put(cf, num);
  }

  mpz_clear(num);
  mpz_clear(denom);
  return NULL;
}

static void *regularized_pi(cf_t cf) {
  mpz_t a, b, c, d;
  mpz_init(a); mpz_init(b); mpz_init(c); mpz_init(d);
  mpz_set_ui(a, 0); mpz_set_ui(b, 4);
  mpz_set_ui(c, 1); mpz_set_ui(d, 0);
  cf_t nonregpi = cf_new(pi_arctan_sequence, NULL);
  cf_t conv = cf_new_nonregular_to_cf(nonregpi, a, b, c, d);
  mpz_t z;
  mpz_init(z);
  while(cf_wait(cf)) {
    cf_get(z, conv);
    cf_put(cf, z);
  }
  mpz_clear(z);
  cf_free(conv);
  cf_free(nonregpi);
  mpz_clear(a); mpz_clear(b); mpz_clear(c); mpz_clear(d);
  return NULL;
}

cf_t cf_new_pi() {
  return cf_new_const(regularized_pi);
}

// tan 1 = [1; 1, 1, 3, 1, 5, ...] 
static void *tan1_expansion(cf_t cf) {
  int odd = 1;

  while(cf_wait(cf)) {
    cf_put_int(cf, 1);
    cf_put_int(cf, odd);
    odd += 2;
  }

  return NULL;
}

cf_t cf_new_tan1() {
  return cf_new_const(tan1_expansion);
}

// exp(z) = 1 + z/(1 - z/(2 + z/(3 - z/(2 + z/(5 - z/(2 + z/ ...))))))
void *exp_expansion(cf_t cf) {
  mpz_ptr z = cf_data(cf);
  mpz_t minusz;
  int odd = 1;
  mpz_init(minusz);
  mpz_neg(minusz, z);
  cf_put_int(cf, 1);
  while(cf_wait(cf)) {
    cf_put(cf, z);
    cf_put_int(cf, odd);
    odd += 2;
    cf_put(cf, minusz);
    cf_put_int(cf, 2);
  }
  mpz_clear(minusz);
  mpz_clear(z);
  free(z);
  return NULL;
}

// tanh n = z/(1 + z^2/(3 + z^2/(5 + z^2/...)))
static void *gauss_tanh_expansion(cf_t cf) {
  mpz_ptr z = cf_data(cf);
  mpz_t z2;
  mpz_t odd;
  mpz_init(odd); mpz_init(z2);
  mpz_set_ui(odd, 1);
  mpz_set_ui(z2, 0);
  cf_put(cf, z2);
  cf_put(cf, z);
  mpz_mul(z2, z, z);
  while(cf_wait(cf)) {
    cf_put(cf, odd);
    mpz_add_ui(odd, odd, 2);
    cf_put(cf, z2);
  }
  mpz_clear(odd); mpz_clear(z2);
  mpz_clear(z);
  free(z);
  return NULL;
}

// tan n = z/(1 - z^2/(3 - z^2/(5 - z^2/...)))
static void *gauss_tan_expansion(cf_t cf) {
  mpz_ptr z = cf_data(cf);
  mpz_t z2;
  mpz_t odd;
  mpz_init(odd); mpz_init(z2);
  mpz_set_ui(odd, 1);
  mpz_set_ui(z2, 0);
  cf_put(cf, z2);
  cf_put(cf, z);
  mpz_mul(z2, z, z);
  mpz_neg(z2, z2);
  while(cf_wait(cf)) {
    cf_put(cf, odd);
    mpz_add_ui(odd, odd, 2);
    cf_put(cf, z2);
  }
  mpz_clear(odd); mpz_clear(z2);
  mpz_clear(z);
  free(z);
  return NULL;
}

cf_t cf_new_one_arg(void *(*fun)(cf_t), mpz_t z) {
  mpz_ptr p = malloc(sizeof(*p));
  mpz_init(p);
  mpz_set(p, z);
  return cf_new(fun, p);
}

struct funarg_s {
  void *(*fun)(cf_t);
  mpz_t arg;
};
typedef struct funarg_s *funarg_ptr;

static void *one_arg_nonreg(cf_t cf) {
  funarg_ptr p = cf_data(cf);
  mpz_t a, b, c, d;
  mpz_init(a); mpz_init(b); mpz_init(c); mpz_init(d);
  mpz_set_ui(a, 1); mpz_set_ui(b, 0);
  mpz_set_ui(c, 0); mpz_set_ui(d, 1);
  mpz_ptr copy = malloc(sizeof(*copy));
  mpz_init(copy);
  mpz_set(copy, p->arg);
  cf_t nonreg = cf_new(p->fun, copy);
  cf_t conv = cf_new_nonregular_to_cf(nonreg, a, b, c, d);
  mpz_t z;
  mpz_init(z);
  while(cf_wait(cf)) {
    cf_get(z, conv);
    cf_put(cf, z);
  }
  mpz_clear(z);
  cf_free(conv);
  cf_free(nonreg);
  mpz_clear(a); mpz_clear(b); mpz_clear(c); mpz_clear(d);

  mpz_clear(p->arg);
  free(p);
  return NULL;
}

cf_t cf_new_one_arg_nonreg(void *(*fun)(cf_t), mpz_t z) {
  funarg_ptr p = malloc(sizeof(*p));
  p->fun = fun;
  mpz_init(p->arg);
  mpz_set(p->arg, z);
  return cf_new(one_arg_nonreg, p);
}

cf_t cf_new_epow(mpz_t pow) {
  return cf_new_one_arg_nonreg(exp_expansion, pow);
}

cf_t cf_new_tanh(mpz_t z) {
  return cf_new_one_arg_nonreg(gauss_tanh_expansion, z);
}

cf_t cf_new_tan(mpz_t z) {
  return cf_new_one_arg_nonreg(gauss_tan_expansion, z);
}
