# program to model the pressure profile below the seafloor after an
# instantaneous change in sea level
# allows variable porosity and permeability
# written by Hugh Daigle, 8 September 2021
import numpy as np # import numpy
import matplotlib.pyplot as plt #import pyplot

dsf=800 #seafloor depth in m
rhow=1024 #water density in kg/m3
kbw = 2.25 #water bulk modulus in GPa
dz=20 #vertical resolution in m
maxiter=300000 #maximum number of time steps
nz=41 #%number of vertical steps
mu=0.001 #water viscosity

depth=np.arange(0,dz*(nz-1),dz)
P_init=list(map(lambda num : (num+dsf)*rhow*9.80665, depth)) #initial pressure in Pa
P_star_init=[0]*nz #initial overpresure in Pa (assumed to be zero)
phi=list(map(lambda num : 0.775*np.exp(-num/1251), depth)) #porosity from Kominz et al 2011 (https://doi.org/10.2110/jsr.2011.60)
k=list(map(lambda num : 10**((8.3887*num)-20.862), phi)) #permeability in m^2 from Daigle and Screaton 2015 (https://doi.org/10.1111/gfl.12090)
alpha=max(map(lambda num1, num2 : (num1*kbw*1e9)/(num2*mu), k, phi)) #Finite difference stability criterion
dt=np.floor(0.5*dz*dz/alpha) #set time step

dsf_new=650 #new seafloor depth in m
P_steady=list(map(lambda num : (num+dsf_new)*rhow*9.80665, depth)) #new steady-state pressure in Pa
P_star_old=list(map(lambda num1, num2 : num1-num2, P_init, P_steady))
P_star_new=P_star_old
P_star_old[0]=0

iteration=1
for i in range (1,maxiter):
    P_star_new[0]=0
    for j in range (1,nz-2):
        P_star_new[j]=((P_star_old[j])+((dt*kbw*1e9/mu)*(((1/dz)*(((0.5*((k[j+1]/phi[j+1])+(k[j]/phi[j])))*(((P_star_old[j+1])-(P_star_old[j]))/dz))-((0.5*((k[j-1]/phi[j-1])+(k[j]/phi[j])))*(((P_star_old[j])-(P_star_old[j-1]))/dz)))))))      
    P_star_new[-1]=P_star_new[-2]
    P_star_old=P_star_new
    if np.floor(dt*iteration/(365*24*3600))==5:
        P5=list(map(lambda num1, num2 : num1+num2, P_steady, P_star_new))
    if max(map(lambda num1, num2 : np.absolute(num1/num2), P_star_new, P_steady))<0.001:
        break #stop when pressure profile deviates by no more than 0.1% from the steady-state profile
    iteration=iteration+1
    print(iteration)
    
time=dt*iteration/(365*24*3600) #in years
print(time)
P_new=list(map(lambda num1, num2 : num1+num2, P_steady, P_star_new))
plt.figure(figsize=(12,6))
plt.subplot(131)
plt.plot(P_init,depth,color='black',label="Initial")
plt.plot(P5,depth,color='green',label='5 years')
plt.plot(P_new,depth,color='red',label="Final")
plt.plot(P_steady,depth,color='red',linestyle='dashed',label="Steady state")
plt.xlabel('Pore pressure (Pa)')
plt.ylabel('Depth (mbsf)')
plt.legend()
plt.gca().invert_yaxis()
plt.subplot(132)
plt.plot(k,depth,color='black')
plt.xscale('log')
plt.xlabel('Permeability (m^2)')
plt.gca().invert_yaxis()
plt.subplot(133)
plt.plot(phi,depth,color='black')
plt.xlabel('Porosity')
plt.gca().invert_yaxis()
plt.show()