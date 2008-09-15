#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

int main(int argc, char *argv[]) {
  int k = 5;
  if (argc > 1) k = atoi(argv[1]);

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
  gmp_printf("%Zd / %Zd\n", r, s);

  mpz_t limit;
  mpz_init(limit);
  // Compute limit = (2k + 1)! / 2
  mpz_mul_ui(limit, s, k);
  mpz_mul_ui(limit, s, 2 * k + 1);
  gmp_printf("%Zd\n", limit);

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
  for (;;) {
    mpz_fdiv_qr(t0, t1, a, b);

    mpz_addmul(p0, p1, t0);
    mpz_addmul(q0, q1, t0);
    mpz_mul(a, q0, q0);
    if (mpz_cmp(a, limit) >= 0) {
      printf("limit exceeded\n"); break;
    }
    mpz_swap(p0, p1);
    mpz_swap(q0, q1);
    gmp_printf("%Zd (%Zd) -> %Zd / %Zd\n", t0, t1, p1, q1);
    if (!mpz_sgn(t1)) break;
    mpz_set(a, b);
    mpz_set(b, t1);
  }

  mpz_clear(limit);
  mpz_clear(r); mpz_clear(s);
  mpz_clear(t0); mpz_clear(t1);
  mpz_clear(a); mpz_clear(b);
  return 0;
}
