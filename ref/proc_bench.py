#!/usr/bin/python

import sys
import parse

keygen_time = []
sign_time = []
verify_time = []

keygen_cyc = []
sign_cyc = []
verify_cyc = []

for line in sys.stdin:
  if line.startswith("Signature verification failed!"):
    sys.stderr.write("\n\n  Signature verification failed!\n\n\n")
    continue
  if ":" in line:
    if line.startswith("name:"):
      name = parse.parse("name: {name}", line.strip())['name']
    if line.startswith("pk:"):
      pk = parse.parse("pk:   {pk:d} bytes", line.strip())['pk']
    continue
    print(line.strip())
  else:
    parsed = parse.parse("{keygen_time:f} ({keygen_cyc:d} cycles)  {sign_time:f} ({sign_cyc:d} cycles)  {verify_time:f} ({verify_cyc:d} cycles)", line.strip())

    keygen_time.append(parsed['keygen_time'])
    sign_time.append(parsed['sign_time'])
    verify_time.append(parsed['verify_time'])

    keygen_cyc.append(parsed['keygen_cyc']/1000000)
    sign_cyc.append(parsed['sign_cyc']/1000000)
    verify_cyc.append(parsed['verify_cyc']/1000000)

keygen_time.sort()
sign_time.sort()
verify_time.sort()

keygen_cyc.sort()
sign_cyc.sort()
verify_cyc.sort()

print(f"{name} & {keygen_time[len(keygen_time)>>1]:f} & {keygen_cyc[len(keygen_cyc)>>1]}" + 
      f" & {sign_time[len(sign_time)>>1]:f} & {sign_cyc[len(sign_cyc)>>1]}" + 
      f" & {verify_time[len(verify_time)>>1]} & {verify_cyc[len(verify_cyc)>>1]} \\\\")

