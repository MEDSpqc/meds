import math

class Node:
  def __init__(self, height, value = object, f = lambda a, b : (a, a), h=0, i=0):
    self.height = height

    self.left  = None
    self.right = None

    self.update(value, f, h, i)

  def update(self, value, f = lambda a : (a, a), h=0, i=0):
    self.value = value

    if self.height > 0:
      lv, rv = f(self.value, (1<<h) + i - 1)

      self.left  = Node(self.height-1, lv, f, h+1, i<<1)
      self.right = Node(self.height-1, rv, f, h+1, (i<<1) + 1)

  def __str__(self):
    if self.value:
      if type(self.value) == bytes:
        return str(self.value.hex())
    else:
      return "--"


  def __getitem__(self, key):
    assert key < 2**self.height

    if self.height == 0:
      return self.value

    if key < 2**(self.height-1):
      return self.left[key]
    else:
      return self.right[key & ((1 << (self.height-1)) - 1)]

  def __setitem__(self, key, value):
    assert key < 2**self.height

    self.value = None

    if self.height == 0:
      self.value = value

      return

    if key < 2**(self.height-1):
      self.left[key] = value
    else:
      self.right[key & ((1 << (self.height-1)) - 1)] = value


  def delete(self, key):
    self.value = None

    if self.height == 0:
      return True

    if key < 2**(self.height-1):
      return self.left.delete(key)
    else:
      return self.right.delete(key & ((1 << (self.height-1)) - 1))

  def path(self, leafs, lshift, level=0, offset=0):
    if not self.value:
      if self.height > 0:
        return self.left.path(leafs, lshift-1, level+1, offset<<1) + self.right.path(leafs, lshift-1, level+1, (offset<<1) + 1)
      else:
        return []
    else:
      if offset < (leafs / (1 << lshift)):
        return [self.value]
      else:
        return []

  def patch(self, path, f = lambda a, b : (a, a), h=0, i=0):
    try:
      if not self.value:
        if self.height > 0:
          path = self.left.patch(path, f, h+1, i<<1)
          path = self.right.patch(path, f, h+1, (i<<1) + 1)
      else:
        self.update(path[0], f, h, i)
        path = path[1:]
    except IndexError: pass

    return path


class SeedTree:
  def __init__(self, leafs, value = object, f = lambda a, b : (a, a)):
    self.height = math.ceil(math.log(leafs,2))
    self.leafs = leafs
    self.root = Node(self.height, value, f)

  def __getitem__(self, key):
    return self.root[key]

  def __setitem__(self, key, value):
    self.root[key] = value

  def delete(self, key):
    return self.root.delete(key)

  def path(self):
    return self.root.path(self.leafs, self.height)

  def patch(self, path, f = lambda a : (a, a)):
    self.root.patch(path, f)

  def tree(self):
    stack = [self.root]

    i = 0

    while stack[i].left:
      stack.append(stack[i].left)
      stack.append(stack[i].right)

      i = i + 1

    ret = ""

    for i in range(self.height+1):
      for j in range(1<<i):
        ret += "   "*((1<<(self.height-i))-1)
        try:
          ret += "  " + str(stack[0])[:2] + "  "
        except TypeError:
          ret += "  XX  "
        ret += "   "*((1<<(self.height-i))-1)

        stack = stack[1:]

      ret += "\n\n"

    return ret


if __name__ == "__main__":
  import random
  from Crypto.Hash import SHAKE256

  leafs = 14

  def hash_pair(value, _):
    shake = SHAKE256.new()
    shake.update(value)
    digest = shake.read(256>>3)

    return digest[:16], digest[16:]

  seed = random.randbytes(16)

  shake = SHAKE256.new()
  shake.update(seed)
  seed = shake.read(16)


  ## signer side
  tree = SeedTree(leafs, seed, hash_pair)

  print(tree.tree())

  for i in range(leafs):
    print(tree[i].hex())


  tree.delete(2)
  tree.delete(8)
  tree.delete(12)

  print()

  print(tree.tree())

  for i in range(leafs):
    print(tree[i].hex() if type(tree[i]) == bytes else "--")


  print()

  path = tree.path()

  print("\npath:\n")
  print("\n".join([v.hex() for v in path]))


  ## verifier side
  tree2 = SeedTree(leafs)

  tree2.delete(2)
  tree2.delete(8)
  tree2.delete(12)


  print(tree2.tree())

  tree2.patch(path, hash_pair)

  print(tree2.tree())

  for i in range(leafs):
    print(tree2[i].hex() if type(tree2[i]) == bytes else "--")

