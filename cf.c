// Demand channels. See squint paper by McIlroy.
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <gmp.h>

struct channel_s {
  void *data;
  struct channel_s *next;
};
typedef struct channel_s channel_t[1];
typedef struct channel_s *channel_ptr;

struct cf_s {
  // Each continued fraction is a separate thread.
  pthread_t thread;

  // Signal 'demand' to compute next term.
  // 'ack' acknowledges that the signal was received,
  // to prevent a new 'demand' signal while processing a previous one.
  pthread_cond_t demand, ack;
  pthread_mutex_t demand_mu, ack_mu;

  // Demand channel.
  pthread_cond_t read_cond;
  pthread_mutex_t chan_mu;
  channel_ptr chan, next;

  int quitflag;
  void *data;
};

typedef struct cf_s *cf_t;

int cf_wait(cf_t cf) {
  // Wait for 'demand' signal.
  pthread_cond_wait(&cf->demand, &cf->demand_mu);
  // Acknowledge it, allowing future 'demand' signals.
  pthread_mutex_lock(&cf->ack_mu);
  pthread_cond_signal(&cf->ack);
  pthread_mutex_unlock(&cf->ack_mu);
  if (cf->quitflag) {
    pthread_mutex_unlock(&cf->demand_mu);
    return 0;
  }
  return 1;
}

void cf_signal(cf_t cf) {
  pthread_mutex_lock(&cf->demand_mu);
  pthread_cond_signal(&cf->demand);
  pthread_mutex_unlock(&cf->demand_mu);
  pthread_cond_wait(&cf->ack, &cf->ack_mu);
}

void cf_free(cf_t cf) {
  cf->quitflag = 1;
  pthread_mutex_lock(&cf->demand_mu);
  pthread_cond_signal(&cf->demand);
  pthread_mutex_unlock(&cf->demand_mu);
  pthread_cond_wait(&cf->ack, &cf->ack_mu);
  pthread_join(cf->thread, NULL);
  pthread_mutex_destroy(&cf->demand_mu);
  pthread_mutex_destroy(&cf->ack_mu);
  pthread_cond_destroy(&cf->demand);
  pthread_cond_destroy(&cf->ack);
  free(cf);
}

void cf_put(cf_t cf, mpz_t z) {
  // TODO: Block or something if there's a lot on the queue.
  channel_ptr cnew = malloc(sizeof(*cnew));
  size_t count = (mpz_sizeinbase(z, 2) + 8 - 1) / 8;
  unsigned char *uc = malloc(count + 4);
  cnew->data = uc;
  uc[0] = count >> (8 * 3);
  uc[1] = (count >> (8 * 2)) & 255;
  uc[2] = (count >> 8) & 255;
  uc[3] = count & 255;
  mpz_export(uc + 4, NULL, 1, 1, 1, 0, z);
  cnew->next = NULL;
  pthread_mutex_lock(&cf->chan_mu);
  if (cf->chan) {
    cf->next->next = cnew;
  } else {
    // Channel is empty so send signal in case someone is waiting for it.
    cf->chan = cnew;
    pthread_cond_signal(&cf->read_cond);
  }
  cf->next = cnew;
  pthread_mutex_unlock(&cf->chan_mu);
}

void cf_get(mpz_t z, cf_t cf) {
  pthread_mutex_lock(&cf->chan_mu);
  // Wait for signal if channel is empty.
  if (!cf->chan) {
    pthread_cond_wait(&cf->read_cond, &cf->chan_mu);
  }
  channel_ptr c = cf->chan;
  cf->chan = c->next;
  unsigned char *uc = c->data;
  size_t count = uc[3]
               + (uc[2] << 8)
               + (uc[1] << (8 * 2))
               + (uc[0] << (8 * 3));
  mpz_import(z, count, 1, 1, 1, 0, uc + 4);
  free(c->data);
  free(c);
  pthread_mutex_unlock(&cf->chan_mu);
}

cf_t cf_new(void *(*func)(cf_t)) {
  cf_t cf = malloc(sizeof(*cf));
  pthread_attr_t attr;
  pthread_cond_init(&cf->demand, NULL);
  pthread_mutex_init(&cf->demand_mu, NULL);
  pthread_cond_init(&cf->ack, NULL);
  pthread_mutex_init(&cf->ack_mu, NULL);
  pthread_cond_init(&cf->read_cond, NULL);
  pthread_mutex_init(&cf->chan_mu, NULL);
  pthread_mutex_lock(&cf->demand_mu);
  pthread_mutex_lock(&cf->ack_mu);
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  cf->chan = NULL;
  cf->next = NULL;
  cf->quitflag = 0;

  pthread_create(&cf->thread, &attr, (void*(*)(void *)) func, cf);
  pthread_attr_destroy(&attr);
  return cf;
}

static void *sqrt2(cf_t cf) {
  mpz_t z;
  mpz_init(z);
  mpz_set_ui(z, 1);
  cf_put(cf, z);
  mpz_set_ui(z, 2);
  while(cf_wait(cf)) {
    cf_put(cf, z);
  }
  mpz_clear(z);
  return NULL;
}

cf_t cf_new_sqrt2() {
  return cf_new(sqrt2);
}

int main() {
  mpz_t z;
  mpz_init(z);
  cf_t a;
  a = cf_new_sqrt2();
  for (int i = 0; i < 10; i++) {
    cf_signal(a);
    cf_get(z, a);
    gmp_printf("a: %Zd\n", z);
  }
  cf_free(a);
  mpz_clear(z);
  return 0;
}
