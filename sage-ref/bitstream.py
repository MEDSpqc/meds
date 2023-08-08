class BitStream:
  def __init__(self, data = None):
    if data:
      self.data = data
      self.byte_pos = 0
    else:
      self.data = bytearray()

    self.bit_pos = 0

  def finalize(self):
    if self.bit_pos > 0:
      self.bit_pos = 8

  def write(self, data, data_len):
    byte_pos = len(self.data) - 1

    if (self.bit_pos + data_len < 8):
      self.data[byte_pos] |= data << self.bit_pos

      self.bit_pos += data_len

      if (self.bit_pos == 7):
        self.bit_pos = 0
        self.byte_pos += 1

      return 0

    if (self.bit_pos > 0):
      self.data[byte_pos] |= (data << self.bit_pos) & 0xFF

      data >>= 8 - self.bit_pos
      data_len -= 8 - self.bit_pos

      self.bit_pos = 0

    while (data_len >= 8):
      self.data.append(data & 0xFF)

      data >>= 8
      data_len -= 8

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
      data |= (self.data[self.byte_pos] & ((1 <<  data_len)-1)) << off

      self.bit_pos = data_len

    return data

if __name__ == "__main__":

  data = [6001, 5101, 3590, 598, 3684, 626, 3614, 7813, 8083, 6243]

  bs = BitStream()

  for v in data:
    bs.write(v, 13)

  print(bs.data.hex())


  bs = BitStream(bs.data)

  for i in range(len(data)):
     print(bs.read(13))

