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

n=128
input=np.asfortranarray(np.zeros((n,n,n)))
output=np.asfortranarray(np.zeros((n,n,n)))

# print input

fmodule.solver(input,n)

# print input

print 'Error: ',np.linalg.norm(np.asarray(input-output).reshape(-1),
                               ord=1)

