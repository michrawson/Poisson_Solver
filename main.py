import numpy as np
import math
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

n=8
l=10
error = np.zeros(l)
for j in range(l):
    n=n+8
    print 'n',n
    input=np.asfortranarray(np.zeros((n,n,n)))
    output=np.asfortranarray(np.zeros((n,n,n)))

    hx=10.0/n

    for i3 in range(1,n+1):
        x3 = hx*(i3-n/2)
        for i2 in range(1,n+1):
            x2 = hx*(i2-n/2)
            for i1 in range(1,n+1):
                x1 = hx*(i1-n/2)
                r2 = x1*x1+x2*x2+x3*x3
                r = np.sqrt(r2)
                input[i1-1,i2-1,i3-1] = 1.0/np.pi/np.sqrt(np.pi)*np.exp(-r2)
                if r==0:
                    output[i1-1,i2-1,i3-1] = 2.0/np.sqrt(np.pi)
                else:
                    output[i1-1,i2-1,i3-1] = math.erf(r)/r

    # print input

    fmodule.solver(input,n)

    # print input

    print 'Error: ',np.max(np.abs(np.asarray(input-output).reshape(-1)))

    error[j]=np.max(np.abs(np.asarray(input-output).reshape(-1)))

fig = plt.figure()

ax = fig.add_subplot(111)
ax.set_title('Error')

ax.set_xlabel('number points')
ax.set_ylabel('log(error)')

ax.axis([10, 100, -25, -1])

ax.plot(range(10,10*l+1,10), np.log(error))
plt.show()
