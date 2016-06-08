import numpy as np
import sys
sys.settrace
import fmodule
import matplotlib.pyplot as plt

print fmodule
print
print fmodule.__doc__
print
print fmodule.solver.__doc__
# print
# print fmodule.solve.solve_run.__doc__

n=10
input=np.asfortranarray(np.zeros((n,n,n)))
input=input*1

result = fmodule.solver(input)

# print input
# print result

