
import subprocess
import sys

def install(package):
    subprocess.check_call([sys.executable, "-m", "pip", "install", package])

print("Installing packages...")

lst = ["numpy", "matplotlib", "pandas", "sklearn", "deepwalk", "node2vec", "networkx", "random"]

for i in lst:
    if not i in sys.modules.keys():
        install(i)
        print(i+" installed")
    else:
        print(i+" not installed")

