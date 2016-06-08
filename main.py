import numpy as np
import sys
# sys.settrace
import fmodule
import matplotlib.pyplot as plt

# print fmodule
# print
# print fmodule.__doc__
# print
# print fmodule.solver.__doc__
# print
# print fmodule.solve.solve_run.__doc__

n=10
input=np.asfortranarray(np.zeros((n,n,n)))

result = fmodule.solver(input,n)

print 'Error: ',np.linalg.norm(np.asarray(input-np.asfortranarray(np.zeros((n,n,n)))).reshape(-1),
                     ord=1)

