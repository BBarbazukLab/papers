import os
import sys

os.chdir(sys.argv[1])

old_files = [i for i in os.listdir(sys.argv[1]) if "Cast" in i]
new_files = [i.replace('-','_') for i in old_files]

for x, i in enumerate(new_files):
    os.rename(old_files[x], i)

