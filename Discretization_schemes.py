import numpy as np
import matplotlib.pyplot as plt

k=[0.1,1,5]
ii=10000
dt=1000
f=0.0001

for n in range(len(k)):
	fig, ax = plt.subplots(figsize=(16, 10))
	plt.xlabel('velocity U-component', fontsize=16)
	plt.ylabel('velocity V-component', fontsize=16)
	u=np.zeros((ii))
	v=np.zeros((ii))
	u[0]=2
	for i in range(ii-1):					# Explicit Euler Method
		u[i+1]=f*v[i]*dt+u[i]
		v[i+1]=-f*u[i]*dt+f*k[n]*(1-u[i]**2)*v[i]*dt+v[i]
	plt.xlim(min(u)-1,max(u)+1)
	plt.ylim(min(v)-1,max(v)+1)
	plt.plot(u,v, label='Explicit Euler Method (forward method)')
	u=np.zeros((ii))
	v=np.zeros((ii))
	u[-1]=2
	for i in reversed(range(ii-1)):					# Implicit Euler Method - not correct now!!!
		u[i]=-f*v[i+1]*dt+u[i+1]		# Should be corrcted!!!
		v[i]=f*u[i+1]*dt-f*k[n]*(1-u[i+1]**2)*v[i+1]*dt+v[i+1]
	#plt.plot(u,v, label='Implicit Euler Method (backwards method)')
	u=np.zeros((ii))
	v=np.zeros((ii))
	u[0]=2
	for i in range(ii-1):					# Semi-implicit Euler Method
		u[i+1]=f*v[i]*dt+u[i]		#Rough estimation of the next term !!!!
		v[i+1]=-f*u[i]*dt+f*k[n]*(1-u[i]**2)*v[i]*dt+v[i]	#Rough estimation of the next term !!!!
		u[i+1]=(f*v[i]+f*v[i+1])*dt*0.5+u[i]
		v[i+1]=((-f*u[i]+f*k[n]*(1-u[i]**2)*v[i])+(-f*u[i+1]+f*k[n]*(1-u[i+1]**2)*v[i+1]))*dt*0.5+v[i]
	plt.plot(u,v, label='Semi-Implicit Euler Method (trapezoidal method)')
	u=np.zeros((ii))
	v=np.zeros((ii))
	u[0]=2; a=0.75
	for i in range(ii-1):					# General two-points Euler Method
		#### v[n+1]=v[n]-f*dt*u[n] - the solution 
		u[i+1]=f*v[i]*dt+u[i]	#Rough estimation of the next term !!!!
		v[i+1]=-f*u[i]*dt+f*k[n]*(1-u[i]**2)*v[i]*dt+v[i]	#Rough estimation of the next term !!!!
		u[i+1]=((f*v[i])*(1-a)+(f*v[i+1])*a)*dt+u[i]
		v[i+1]=((1-a)*(-f*u[i]+f*k[n]*(1-u[i]**2)*v[i])+a*(-f*u[i+1]+f*k[n]*(1-u[i+1]**2)*v[i+1]))*dt+v[i]
	plt.xlim(min(u)-1,max(u)+1)
	plt.ylim(min(v)-1,max(v)+1)
	plt.plot(u,v, label='General Two-point scheme (alpha=0.75)')
	ax.annotate('k=%s' %(k[n]), xy=(0,0),  xycoords='figure fraction',xytext=(910, 100), family='Courier New, monospace',textcoords='offset points',ha="left", va="bottom",fontsize=24, bbox=dict(facecolor='lightgrey', edgecolor='none', pad=5.0))
	plt.legend(loc=3)
	fig.savefig(str(n), dpi=300)


