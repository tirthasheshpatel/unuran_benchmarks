import os

files = [file for file in os.listdir() if not os.path.isdir(file)]
files.remove('helper.py')

for file in files:
    print(f"\n\nIn file {file}\n\n")
    os.system(f'cat {file} | grep \#define')
