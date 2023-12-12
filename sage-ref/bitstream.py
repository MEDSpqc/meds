class BitStream:
  def __init__(self, data = None):
    if data:
      self.data = data
    else:
      self.data = bytearray()

    self.bit_pos = 0
    self.byte_pos = 0

  def finalize(self):
    if self.bit_pos > 0:
      self.bit_pos = 8

  def write(self, data, data_len):
    if (self.bit_pos + data_len < 8):
      if self.bit_pos == 0:
        self.data.append(0)

      self.data[self.byte_pos] |= data << self.bit_pos

      self.bit_pos += data_len

      if (self.bit_pos == 7):
        self.bit_pos = 0
        self.byte_pos += 1

      return 0

    if (self.bit_pos > 0):
      self.data[self.byte_pos] |= (data << self.bit_pos) & 0xFF

      data >>= 8 - self.bit_pos
      data_len -= 8 - self.bit_pos

      self.bit_pos = 0
      self.byte_pos += 1

    while (data_len >= 8):
      self.data.append(data & 0xFF)
      self.byte_pos += 1

      data >>= 8
      data_len -= 8

      self.bit_pos = data_len

    if (data_len > 0):
      self.data.append(data)

      self.bit_pos = data_len

    return 0

  def read(self, data_len):
    data = 0

    if (self.bit_pos + data_len < 8):
      data = (self.data[self.byte_pos] >> self.bit_pos) & ((1 << data_len)-1)

      self.bit_pos += data_len

      if (self.bit_pos == 7):
        self.bit_pos = 0
        self.byte_pos += 1

      return data

    off = 0

    if (self.bit_pos > 0):
      data = self.data[self.byte_pos] >> self.bit_pos

      off = 8-self.bit_pos
      data_len -= 8-self.bit_pos

      self.bit_pos = 0
      self.byte_pos += 1

    while (data_len >= 8):
      data |= self.data[self.byte_pos] << off
      off += 8
      self.byte_pos += 1
      data_len -= 8
      self.bit_pos = 0

    if (data_len > 0):
      data |= (self.data[self.byte_pos] & ((1 << data_len)-1)) << off

      self.bit_pos = data_len

    return data

if __name__ == "__main__":
  import random

  num = 10

  widths = [random.randint(1, 20) for _ in range(num)] 
  data = [random.randint(0, 1<<w - 1) for w in widths]

  print(data)
  print(widths)
  print()

  bs = BitStream()

  for v, w in zip(data, widths):
    bs.write(v, w)
    print(v, w, bin(v)[2:].zfill(w), bin(int(bs.data.hex(), 16))[2:].zfill(len(bs.data)*8))

  print()

  print(bs.data.hex())


  bs = BitStream(bs.data)

  ret = []

  acc = 0

  for v, w in zip(data, widths):
    val = bs.read(w)
    print("x " if val != v else '  ', val, v, w, acc)
    ret.append(val)

    acc += w

  print()
  print("OK!" if ret == data else "ERROR!")

