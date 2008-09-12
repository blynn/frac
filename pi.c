// Compute (simple) continued fraction expansion for pi.
// Ben Lynn

#include <stdio.h>
#include <gmp.h>

int main() {
  mpz_t num, denom;

  mpz_t pold, p, pnew;
  mpz_t qold, q, qnew;
  mpz_t t0, t1;

  void recur() {
    mpz_mul(pnew, p, denom);
    mpz_addmul(pnew, pold, num);
    mpz_mul(qnew, q, denom);
    mpz_addmul(qnew, qold, num);

    mpz_gcd(t0, pnew, qnew);
    mpz_gcd(t1, p, q);
    mpz_gcd(t0, t0, t1);
    if (mpz_cmp_ui(t0, 1)) {
      mpz_divexact(pnew, pnew, t0);
      mpz_divexact(qnew, qnew, t0);
      mpz_divexact(p, p, t0);
      mpz_divexact(q, q, t0);
    }

    // If we just want convergents:
    // gmp_printf("p/q = %Zd/%Zd\n", pnew, qnew);
    mpz_fdiv_qr(pold, t0, p, q);
    // mpz_fdiv_qr(qold, t1, pnew, qnew);
    mpz_mul(qold, pold, qnew);

    if (mpz_cmp(qold, pnew) < 0) {
      mpz_add(qold, qold, qnew);
      //if (!mpz_cmp(pold, qold)) {
      if (mpz_cmp(qold, pnew) > 0) {  // It should never be 0 for pi.
	gmp_printf(" %Zd", pold);
	mpz_sub(t1, qold, pnew);
	mpz_sub(t1, qnew, t1);
	mpz_set(pold, q);
	mpz_set(qold, t0);
	mpz_set(p, qnew);
	mpz_set(q, t1);
	return;
      }
    }
    mpz_set(pold, p); mpz_set(p, pnew);
    mpz_set(qold, q); mpz_set(q, qnew);
  }

  mpz_init(num); mpz_init(denom);
  mpz_init(pold); mpz_init(p); mpz_init(pnew);
  mpz_init(qold); mpz_init(q); mpz_init(qnew);
  mpz_init(t0); mpz_init(t1);

  mpz_set_ui(pold, 4); mpz_set_ui(p, 0);
  mpz_set_ui(qold, 0); mpz_set_ui(q, 1);

  mpz_set_ui(num, 1); mpz_set_ui(denom, 1);
  recur();

  mpz_set_ui(num, 1); mpz_set_ui(denom, 3);
  recur();

  int i;
  for (i = 1; i < 5000; i++) {
    mpz_add_ui(denom, denom, 2);
    mpz_add_ui(num, num, (i << 1) + 1);
    recur();
  }
  printf("\n");

  return 0;
}
