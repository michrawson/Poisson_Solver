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

def test1(n,hx):
    density=np.asfortranarray(np.zeros((n,n,n)))
    potential=np.asfortranarray(np.zeros((n,n,n)))

    for i3 in range(1,n+1):
        x3 = hx*(i3-n/2)
        for i2 in range(1,n+1):
            x2 = hx*(i2-n/2)
            for i1 in range(1,n+1):
                x1 = hx*(i1-n/2)

                r2 = x1*x1+x2*x2+x3*x3
                r = np.sqrt(r2)

                potential[i1-1,i2-1,i3-1] = np.exp(-r2) # True Charge

                density[i1-1,i2-1,i3-1] = np.exp(-r2)*(-6.0+4.0*r2) # Density
    return potential, density


def test2(n,hx):
    density=np.asfortranarray(np.zeros((n,n,n)))
    potential=np.asfortranarray(np.zeros((n,n,n)))

    for i3 in range(1,n+1):
        x3 = hx*(i3-n/2)
        for i2 in range(1,n+1):
            x2 = hx*(i2-n/2)
            for i1 in range(1,n+1):
                x1 = hx*(i1-n/2)
                r2 = x1*x1+x2*x2+x3*x3
                r = np.sqrt(r2)
                density[i1-1,i2-1,i3-1] = -4.0/np.sqrt(np.pi)*np.exp(-r2) # Density
                if r==0:
                    potential[i1-1,i2-1,i3-1] = 2.0/np.sqrt(np.pi) # True Charge
                else:
                    potential[i1-1,i2-1,i3-1] = math.erf(r)/r # True Charge
    return potential, density


n=32


rang = 10.0
hx=rang/n

potential, density = test2(n,hx)

potential_approx = fmodule.solve(hx,density,n=n)

print 'Error: ',np.max(np.abs(np.asarray(potential-potential_approx).reshape(-1)))

sys.exit()




# L[j]=n
# error[j]=np.max(np.abs(np.asarray(input-output).reshape(-1)))
#
# fig = plt.figure()
#
# ax = fig.add_subplot(111)
# ax.set_title('Error - Trap.')
#
# ax.set_xlabel('number points')
# ax.set_ylabel('log(error)')
#
# axes = plt.gca()
# axes.set_xlim([0,100])
# axes.set_ylim([-12,-2])
#
# log_error = np.log10(error)
#
# ax.plot(L, log_error)
# plt.show()
