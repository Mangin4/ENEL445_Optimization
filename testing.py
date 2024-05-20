from os import environ
N_THREADS = '8'
environ['OMP_NUM_THREADS'] = N_THREADS
environ['OPENBLAS_NUM_THREADS'] = N_THREADS
environ['MKL_NUM_THREADS'] = N_THREADS
environ['VECLIB_MAXIMUM_THREADS'] = N_THREADS
environ['NUMEXPR_NUM_THREADS'] = N_THREADS
import numpy as np
from numpy.random import rand
import time

start = time.time()

n = 8000

data1 = rand(n, n)
data2 = rand(n, n)

result = data1.dot(data2)

duration = time.time()-start
print(duration)