#!/usr/bin/env python3
import sys
import re
import os
import subprocess

# fd = os.open('foo.sam', os.O_RDWR | os.O_CREAT)
# os.dup2(fd, 4)

nm = nu = ne = reads = 0

sam = sys.stderr

for line in sys.stdin:

    if line.startswith("@"):
        # os.write(4, line.encode())
        continue

    fields = line.strip().split("\t", 11)
    flag = int(fields[1])
    # Skip read if it fails any of the checks
    if int(fields[1]) & 0x4 or int(fields[1]) & 0x100 or int(fields[1]) & 0x800:
        continue
    # Read is not mapped
    if (flag & 0x900 != 0):
        nm += 1
        continue
    # Too many missmatches
    if not re.findall(r"XM:i:[012]", fields[11]):
        ne += 1
        continue
    #Reads are not uniqely mapped
    if not re.findall(r"NH:i:1", fields[11]):
        nu += 1
        continue
    XM = re.findall(r"XM:i:[012]", fields[11])
    NH = re.findall(r"NH:i:1", fields[11])
    # print(XM,NH)
    sign =  "-" if flag & 0x10 else "+"
    contacts_line=f"{fields[2]}\t{fields[3]}\t{int(fields[3]) + len(fields[9])}\t{fields[0]}\t1\t{sign}\n"
    sys.stdout.write(contacts_line)
    # os.write(4, line.encode())
    reads += 1

# Print statistics to standard error
sys.stderr.write(f'Total reads: {reads}\n')
sys.stderr.write(f'Unmapped reads: {nm}\n')
sys.stderr.write(f'Non-unique reads: {nu}\n')
sys.stderr.write(f'Reads with >2 mismatches: {ne}\n')

# os.close(fd)
# p = subprocess.Popen(['samtools', 'view', '-S', '-b', '-o', 'contacts.bam'], stdin=4)
# p = subprocess.Popen(['wc', '-l'], stdin=4)
# p.wait()