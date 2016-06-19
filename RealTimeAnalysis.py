import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

x = np.array([25,50,75,100,125,150,175,200,225,250])
y = np.array([321,534,690,1012,1382,1976,3004,4242,7040,10243]) # rounded mean values
u = np.array([20,30,30,40,100,100,100,150,200,1000]) # rounded (nearest ten/hundred) standard deviations

def expon(x,A,B,C):
    return A*np.e**(B*x) + C

params = [4,0.0308,+300]
popt, pcov = curve_fit(expon,x,y,params, absolute_sigma = True)
xopt = np.linspace(0,275,1000)
yopt = expon(xopt,*popt)

dymin = (y - expon(x, *popt))/u #chi-squared value using optimal params
min_chisq = sum(dymin*dymin)
dof = len(x) - len(popt) #degrees of freedom

print 'Chi square: ', min_chisq
print 'Number of degrees of freedom: ',dof
print 'Chi sqaure per degree of freedom: ', min_chisq/dof
print ''


plt.plot(x,y,'b+',label='Measured values')
plt.plot(xopt,yopt,'g',label='Model fit')
plt.xlim(0,275)
#plt.ylim(1e0,1e4)
plt.xlabel("Chain Length",fontsize=20)
plt.ylabel("Time [s]",fontsize=20)
plt.yscale('log')
plt.legend()
plt.show()
