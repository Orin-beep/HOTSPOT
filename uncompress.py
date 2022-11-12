import os
import sys


model_path = sys.argv[1]    # models/
tars = os.listdir(model_path)   # models1.tar, models2.tar
tars = [model_path+'/'+x for x in tars]

for tar in tars:
    os.system(f'tar -xvf {tar} -C {model_path}')
    os.system(f'rm {tar}')

folders = os.listdir(model_path)
folders = [model_path+'/'+x for x in folders]
for folder in folders:
    os.system(f'mv {folder}/* {model_path}')
    os.system(f'rm -r {folder}')



