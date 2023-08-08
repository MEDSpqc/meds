from Crypto.Cipher import AES

DRBG_ctx_key = [0]*32
DRBG_ctx_v = [0]*16
DRBG_ctx_reseed_counter = 1


def randombytes_init(entropy_input, personalization_string, security_strength):
  global DRBG_ctx_key
  global DRBG_ctx_v
  global DRBG_ctx_reseed_counter

  assert len(entropy_input) == 48
  assert security_strength == 256

  seed_material = bytes(entropy_input)
    
  if personalization_string:
    for i in range(48):
      seed_material[i] ^= personalization_string[i]

  DRBG_ctx_key = [0]*32
  DRBG_ctx_v = [0]*16

  DRBG_ctx_key, DRBG_ctx_v = AES256_CTR_DRBG_Update(seed_material, DRBG_ctx_key, DRBG_ctx_v)

  DRBG_ctx_reseed_counter = 1

def randombytes(xlen):
  global DRBG_ctx_key
  global DRBG_ctx_v
  global DRBG_ctx_reseed_counter

  block = bytes([0]*16)
  i = 0

  ret = bytes()
  
  while xlen > 0:
    for j in range(15, -1, -1):
      if DRBG_ctx_v[j] == 0xff:
        DRBG_ctx_v[j] = 0x00
      else:
        DRBG_ctx_v[j] += 1
        break

    cipher = AES.new(bytes(DRBG_ctx_key), AES.MODE_ECB)
    block = cipher.encrypt(bytes(DRBG_ctx_v))


    if xlen > 15:
      ret = ret + block
      i += 16
      xlen -= 16
    else:
      ret = ret + block[:xlen]
      xlen = 0

  DRBG_ctx_key, DRBG_ctx_v = AES256_CTR_DRBG_Update(None, DRBG_ctx_key, DRBG_ctx_v)
  DRBG_ctx_reseed_counter += 1

  return ret


def AES256_CTR_DRBG_Update(provided_data, Key, V):
  temp = bytes()
 
  for i in range(3):
    for j in range(15, -1, -1):
      if V[j] == 0xff:
        V[j] = 0x00
      else:
        V[j] += 1
        break
    
    cipher = AES.new(bytes(Key), AES.MODE_ECB)
    temp = temp + cipher.encrypt(bytes(V))

  temp = bytearray(temp)

  if provided_data != None:
    for i in range(48):
      temp[i] ^= provided_data[i];

  return temp[:32], temp[32:48]

if __name__ == "__main__":

  seed = [6, 21, 80, 35, 77, 21, 140, 94, 201, 85, 149, 254, 4, 239, 122, 37, 118, 127, 46, 36, 204, 43, 196, 121, 208, 157, 134, 220, 154, 188, 253, 231, 5, 106, 140, 38, 111, 158, 249, 126, 208, 133, 65, 219, 210, 225, 255, 161]
  
  randombytes_init(seed, None, 256) #list(range(48)), None, 256)
  
  print("".join([f"{int(v)} " for v in randombytes(48)]))
  
