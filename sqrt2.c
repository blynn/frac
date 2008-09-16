#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "cf.h"

static void *sqrt2(cf_t cf) {
  cf_put_int(cf, 1);
  while(cf_wait(cf)) {
    cf_put_int(cf, 2);
  }
  return NULL;
}

cf_t cf_new_sqrt2() {
  return cf_new_const(sqrt2);
}

int main() {
  cf_t x, conv;
  x = cf_new_sqrt2();
  conv = cf_new_convergent(x);

  mpz_t p, q;
  mpz_init(p);
  mpz_init(q);
  for (int i = 0; i < 10; i++) {
    cf_get(p, conv);
    cf_get(q, conv);
    gmp_printf("p/q = %Zd/%Zd\n", p, q);
  }
  cf_free(conv);
  cf_free(x);

  x = cf_new_sqrt2();
  mpz_t a, b, c, d;
  mpz_init(a); mpz_init(b); mpz_init(c); mpz_init(d);
  mpz_set_si(a, -7);
  mpz_set_si(b, 5);
  mpz_set_si(c, -12);
  mpz_set_si(d, 13);
  cf_t dec = cf_new_mobius_to_decimal(x, a, b, c, d);
  cf_get(p, dec);
  putchar(cf_sign(dec) > 0 ? ' ' : '-');
  gmp_printf("%Zd.", p);
  for (int i = 0; i < 10; i++) {
    cf_get(p, dec);
    gmp_printf("%Zd", p);
  }
  putchar('\n');
  cf_free(x);
  mpz_clear(p);
  mpz_clear(q);
  return 0;
}
