#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include <sys/random.h>

#include "params.h"
#include "meds.h"

double osfreq(void);

long long cpucycles(void)
{
  unsigned long long result;
  asm volatile(".byte 15;.byte 49;shlq $32,%%rdx;orq %%rdx,%%rax"
      : "=a" (result) ::  "%rdx");
  return result;
}

int main(int argc, char *argv[])
{
  printf("paramter set: %s\n\n", MEDS_name);

  long long time = 0;
  long long keygen_time = 0xfffffffffffffff;
  long long sign_time = 0xfffffffffffffff;
  long long verify_time = 0xfffffffffffffff;

  int rounds = 1;

  if (argc > 1)
    rounds = atoi(argv[1]);

  uint8_t key_seed[16] = {'d', 'e', 'a', 'd', 'b', 'e', 'e', 'f', 'd', 'e', 'a', 'd', 'b', 'e', 'e', 'f'};
  uint8_t sig_seed[16] = {'d', 'e', 'a', 'd', 'b', 'e', 'e', 'f', 'd', 'e', 'a', 'd', 'b', 'e', 'e', 'e'};

  char msg[17] = "TestTestTestTest\0";

  //getrandom(key_seed, sizeof(key_seed), GRND_NONBLOCK) == -1 ?
  //  perror("getrandom") : "";

  //getrandom(sig_seed, sizeof(sig_seed), GRND_NONBLOCK) == -1 ?
  //  perror("getrandom") : "";

#ifdef MEDS_MAX_PATH_LEN
  printf("max path len: %i\n", MEDS_MAX_PATH_LEN);
#endif


  // Run at least twice to get binary into caches...
  for (int round = 0; round < rounds; round++)
  {
    uint8_t sk[SK_BYTES] = {0};
    uint8_t pk[MEDS_PK_BYTES] ={0};

    time = -cpucycles();
    keygen(key_seed, sizeof(key_seed), sk, sizeof(sk), pk, sizeof(pk));
    time += cpucycles();

    if (time < keygen_time) keygen_time = time;

    //for (int i = 0; i < sizeof(sk); i++)
    //  printf("%02x", sk[i]);
    //printf("\n");

    //for (int i = 0; i < sizeof(pk); i++)
    //  printf("%02x", pk[i]);
    //printf("\n");

    if (round == 0)
      printf("pk:  %li bytes\n", sizeof(pk));

    uint8_t sig[MEDS_SIG_BYTES] = {0};

    time = -cpucycles();
    sign(sig_seed, sizeof(sig_seed), sk, sizeof(sk), msg, 4, sig, sizeof(sig));
    time += cpucycles();

    if (time < sign_time) sign_time = time;

    //for (int i = 0; i < sizeof(sig); i++)
    //  printf("%02x", sig[i]);
    //printf("\n");

    if (round == 0)
      printf("sig: %li bytes\n", sizeof(sig));


    time = -cpucycles();
    char *ret = verify(pk, sizeof(pk), msg, 4, sig, sizeof(sig));
    time += cpucycles();

    if (time < verify_time) verify_time = time;

    if (round == 0)
    {
      printf("\n");
      if (ret)
        printf("success\n");
      else
        printf("!!! FAILED !!!\n");
      printf("\n");

      //if (rounds > 0)
      //{
      //  keygen_time = 0;
      //  sign_time = 0;
      //  verify_time = 0;
      //}
    }
  }

  double freq = osfreq();

  //keygen_time /= rounds;
  //sign_time /= rounds;
  //verify_time /= rounds;

  printf("Time:\n");
  printf("keygen: %f   (%llu cycles)\n", keygen_time / freq, keygen_time);
  printf("sign:   %f   (%llu cycles)\n", sign_time / freq, sign_time);
  printf("verify: %f   (%llu cycles)\n", verify_time / freq, verify_time);
  // print("sign:   {0:1.4f}".format(sig_time))
  // print("verify: {0:1.4f}".format(verify_time))

  return 0;
}

