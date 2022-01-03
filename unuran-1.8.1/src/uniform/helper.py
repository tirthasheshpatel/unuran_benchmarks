import os
import re

files = os.listdir()
files.remove('.deps')
files.remove('helper.py')

for file in files:
    os.system(f'cat {file} | grep \#define')
