#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include <sys/random.h>
#include <sys/time.h>

#include "params.h"
#include "meds.h"

#include <x86intrin.h>

double osfreq(void);

long long cpucycles(void)
{
  //unsigned long long result;
  //asm volatile(".byte 15;.byte 49;shlq $32,%%rdx;orq %%rdx,%%rax"
  //    : "=a" (result) ::  "%rdx");
  //return result;

  return __rdtsc();
}

int main(int argc, char *argv[])
{
  printf("name: %s\n", MEDS_name);

  long long keygen_time = 0;
  long long sign_time = 0;
  long long verify_time = 0;

  int rounds = 16;

  if (argc > 1)
    rounds = atoi(argv[1]);

  uint8_t key_seed[16] = {'d', 'e', 'a', 'd', 'b', 'e', 'e', 'f', 'd', 'e', 'a', 'd', 'b', 'e', 'e', 'f'};
  uint8_t sig_seed[16] = {'d', 'e', 'a', 'd', 'b', 'e', 'e', 'f', 'd', 'e', 'a', 'd', 'b', 'e', 'e', 'b'};

  char msg[16] = "TestTestTestTest";

#ifdef MEDS_MAX_PATH_LEN
  printf("max path len: %i\n", MEDS_MAX_PATH_LEN);
#endif

  uint8_t sk[SK_BYTES];
  uint8_t pk[MEDS_PK_BYTES];

  uint8_t sig[MEDS_SIG_BYTES];


  printf("m:    %i\n", MEDS_m);
  printf("n:    %i\n", MEDS_n);
  printf("q:    %i\n", MEDS_p);
  printf("k:    %i\n", MEDS_k);
  printf("s:    %i\n", MEDS_s);
  printf("sec:  %i\n", MEDS_sec_bytes);
  printf("t:    %i\n", MEDS_t);
  printf("w:    %i\n", MEDS_w);

  printf("pk:   %i bytes\n", MEDS_PK_BYTES);

  printf("sig:  %i bytes\n", MEDS_SIG_BYTES);


#ifdef HAVE_SEED_TREE
  printf("HAVE_SEED_TREE: True\n");
#else
  printf("HAVE_SEED_TREE: False\n");
#endif

  // struct timeval start_time, end_time;
  // gettimeofday(&start_time, NULL);


  // Run at least twice to get binary into caches...
  for (int round = 0; round < rounds; round++)
  {
    getrandom(key_seed, sizeof(key_seed), GRND_NONBLOCK) == -1 ?
      perror("getrandom") : "";

    getrandom(sig_seed, sizeof(sig_seed), GRND_NONBLOCK) == -1 ?
      perror("getrandom") : "";


    keygen_time = -cpucycles();
    keygen(key_seed, sizeof(key_seed), sk, sizeof(sk), pk, sizeof(pk));
    keygen_time += cpucycles();

    //for (int i = 0; i < sizeof(sk); i++)
    //  printf("%02x", sk[i]);
    //printf("\n");

    //for (int i = 0; i < sizeof(pk); i++)
    //  printf("%02x", pk[i]);
    //printf("\n");

    sign_time = -cpucycles();
    sign(sig_seed, sizeof(sig_seed), sk, sizeof(sk), msg, 4, sig, sizeof(sig));
    sign_time += cpucycles();

    //for (int i = 0; i < sizeof(sig); i++)
    //  printf("%02x", sig[i]);
    //printf("\n");

    verify_time = -cpucycles();
    verify(pk, sizeof(pk), msg, 4, sig, sizeof(sig));
    verify_time += cpucycles();


    //keygen_time /= rounds;
    //sign_time /= rounds;
    //verify_time /= rounds;

    double freq = osfreq() / 1000;

    printf("F: %f\n", freq);

    printf("%f (%llu cycles)  ", keygen_time / freq, keygen_time);
    printf("%f (%llu cycles)  ", sign_time / freq, sign_time);
    printf("%f (%llu cycles)  \n", verify_time / freq, verify_time);
    // print("sign:   {0:1.4f}".format(sig_time))
    // print("verify: {0:1.4f}".format(verify_time))
  }

  //gettimeofday(&end_time, NULL);

  //printf("seconds : %ld\nmicro seconds : %ld\n",
  //  end_time.tv_sec - start_time.tv_sec, end_time.tv_usec - start_time.tv_usec);

  return 0;
}

