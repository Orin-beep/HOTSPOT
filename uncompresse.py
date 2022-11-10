import os
import sys

path = sys.argv[1]  # models/mdl1
fs = os.listdir(path)
for f in fs:
    os.system(f"tar -xvf {path}/{f} -C {path}")
    id = f[1:f.rfind('.tar')]
    part = path + '/part' + id
    for m in os.listdir(part):
        os.system(f'mv {part}/{m} models/')
    os.system(f'rm -r {part}')
    os.system(f'rm {path}/{f}')
os.system(f'rm -r {path}')


