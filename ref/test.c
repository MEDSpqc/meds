#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include "randombytes.h"

#include "params.h"
#include "api.h"
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

  //unsigned char entropy_input[48] = {0xB5, 0x11,  0};
  //unsigned char entropy_input[48] = {0x82, 0x9F, 0};
  unsigned char entropy_input[48] = {0};

  if (argc > 1)
  {
    uint64_t val = atol(argv[1]);

    for (int i = 0; i < 8; i++)
      entropy_input[i] = (val >> (i*8)) & 0xff;
  }

  printf("seed: ");
  for (int i = sizeof(entropy_input)-1; i >= 0; i--)
    printf("%X", entropy_input[i]);
  printf("\n");


  randombytes_init(entropy_input, NULL, 256);

  char msg[4] = "Test";

  printf("pk:  %i bytes\n", CRYPTO_PUBLICKEYBYTES);
  printf("sk:  %i bytes\n", CRYPTO_SECRETKEYBYTES);
  printf("sig: %i bytes\n", CRYPTO_BYTES);
  printf("\n");

  for (int round = 0; round < rounds; round++)
  {
    uint8_t sk[CRYPTO_SECRETKEYBYTES] = {0};
    uint8_t pk[CRYPTO_PUBLICKEYBYTES] = {0};

    time = -cpucycles();
    crypto_sign_keypair(pk, sk);
    time += cpucycles();

    if (time < keygen_time) keygen_time = time;

    uint8_t sig[CRYPTO_BYTES + sizeof(msg)] = {0};
    unsigned long long sig_len = sizeof(sig);

    time = -cpucycles();
    crypto_sign(sig, &sig_len, (const unsigned char *)msg, sizeof(msg), sk);
    time += cpucycles();

    if (time < sign_time) sign_time = time;

    unsigned char msg_out[4];
    unsigned long long msg_out_len = sizeof(msg_out);


    time = -cpucycles();
    int ret = crypto_sign_open(msg_out, &msg_out_len, sig, sizeof(sig), pk);
    time += cpucycles();

    if (time < verify_time) verify_time = time;

    if (ret == 0)
      printf("success\n");
    else
    {
      printf("!!! FAILED !!!\n");
      return -1;
    }
  }

  double freq = osfreq();

  printf("\n");
  printf("Time (min of %i runs):\n", rounds);
  printf("keygen: %f   (%llu cycles)\n", keygen_time / freq, keygen_time);
  printf("sign:   %f   (%llu cycles)\n", sign_time / freq, sign_time);
  printf("verify: %f   (%llu cycles)\n", verify_time / freq, verify_time);

  return 0;
}

