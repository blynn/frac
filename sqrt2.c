#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "cf.h"

static void *sqrt2(cf_t cf) {
  mpz_t z;
  mpz_init(z);
  mpz_set_ui(z, 1);
  if (!cf_wait(cf)) goto finish;
  cf_put(cf, z);
  mpz_set_ui(z, 2);
  while(cf_wait(cf)) {
    cf_put(cf, z);
  }
finish:
  mpz_clear(z);
  return NULL;
}

cf_t cf_new_sqrt2() {
  return cf_new(sqrt2, NULL);
}

int main() {
  cf_t a, b;
  a = cf_new_sqrt2();
  b = cf_new_convergent(a);

  mpz_t p, q;
  mpz_init(p);
  mpz_init(q);
  for (int i = 0; i < 10; i++) {
    cf_signal(b);
    cf_get(p, b);
    cf_get(q, b);
    gmp_printf("p/q = %Zd/%Zd\n", p, q);
  }
  cf_free(a);
  mpz_clear(p);
  mpz_clear(q);
  return 0;
}
