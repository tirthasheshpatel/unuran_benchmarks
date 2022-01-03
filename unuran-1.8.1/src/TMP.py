import os
from sys import stderr

dirs = [i for i in os.listdir('.') if os.path.isdir(i)]
dirs.remove("tests")

stderr.write(f"{dirs}")

for dir in dirs:
    files = os.listdir(dir)
    for file in files:
        print(f"File {dir + '/' + file}: ", end="")
        os.system(f"cat '{dir + '/' + file}' | grep '#include'")
