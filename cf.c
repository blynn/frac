// Demand channels. See squint paper by McIlroy.
//
// TODO: Free unread channels on exit.
// TODO: Handle messy thread problems. What happens if a thread quits
// but then another tries to signal and read its channel?
// TODO: What if the continued fraction terminates?
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

  // Signal demand to compute next term.
  pthread_cond_t demand;

  // Demand channel.
  pthread_cond_t read_cond;
  pthread_mutex_t chan_mu;
  channel_ptr chan, next;

  int quitflag;
  void *data;
};

typedef struct cf_s *cf_t;

void *cf_data(cf_t cf) {
  return cf->data;
}

int cf_wait(cf_t cf) {
  pthread_mutex_lock(&cf->chan_mu);
  if (cf->chan) {
    // If there is still unread output, wait for 'demand' signal.
    pthread_cond_wait(&cf->demand, &cf->chan_mu);
  }
  pthread_mutex_unlock(&cf->chan_mu);
  if (cf->quitflag) {
    return 0;
  }
  return 1;
}

void cf_free(cf_t cf) {
  cf->quitflag = 1;
  pthread_cond_signal(&cf->demand);
  pthread_join(cf->thread, NULL);
  pthread_cond_destroy(&cf->demand);
  free(cf);
}

void cf_put(cf_t cf, mpz_t z) {
  // TODO: Block or something if there's a large backlog on the queue.
  channel_ptr cnew = malloc(sizeof(*cnew));
  if (!mpz_sgn(z)) {
    unsigned char *uc = malloc(4);
    uc[0] = uc[1] = uc[2] = uc[3] = 0;
    cnew->data=uc;
  } else {
    size_t count = (mpz_sizeinbase(z, 2) + 8 - 1) / 8;
    unsigned char *uc = malloc(count + 4);
    cnew->data = uc;
    uc[0] = count >> (8 * 3);
    uc[1] = (count >> (8 * 2)) & 255;
    uc[2] = (count >> 8) & 255;
    uc[3] = count & 255;
    mpz_export(uc + 4, NULL, 1, 1, 1, 0, z);
  }
  cnew->next = NULL;
  pthread_mutex_lock(&cf->chan_mu);
  if (cf->chan) {
    cf->next->next = cnew;
  } else {
    // Channel is empty. Now that we're populating it, send signal
    // in case someone is waiting for data.
    cf->chan = cnew;
    pthread_cond_signal(&cf->read_cond);
  }
  cf->next = cnew;
  pthread_mutex_unlock(&cf->chan_mu);
}

void cf_get(mpz_t z, cf_t cf) {
  pthread_mutex_lock(&cf->chan_mu);
  if (!cf->chan) {
    // If channel is empty, send demand signal and wait for read signal.
    pthread_cond_signal(&cf->demand);
    pthread_cond_wait(&cf->read_cond, &cf->chan_mu);
  }
  channel_ptr c = cf->chan;
  cf->chan = c->next;
  pthread_mutex_unlock(&cf->chan_mu);
  unsigned char *uc = c->data;
  size_t count = uc[3]
               + (uc[2] << 8)
               + (uc[1] << (8 * 2))
               + (uc[0] << (8 * 3));
  if (count) mpz_import(z, count, 1, 1, 1, 0, uc + 4);
  else mpz_set_ui(z, 0);
  free(c->data);
  free(c);
}

cf_t cf_new(void *(*func)(cf_t), void *data) {
  cf_t cf = malloc(sizeof(*cf));
  pthread_attr_t attr;
  pthread_cond_init(&cf->demand, NULL);
  pthread_cond_init(&cf->read_cond, NULL);
  pthread_mutex_init(&cf->chan_mu, NULL);
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  cf->chan = NULL;
  cf->next = NULL;
  cf->quitflag = 0;
  cf->data = data;
  pthread_create(&cf->thread, &attr, (void*(*)(void *)) func, cf);
  pthread_attr_destroy(&attr);
  return cf;
}
