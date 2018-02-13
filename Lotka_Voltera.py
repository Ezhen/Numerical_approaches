"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""" Written by Evgeny Ivanov. Numerical methods exam, 20/01/2017. The task about two concentrations """
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""



from pylab import *; import matplotlib.pyplot as plt;import sys,os,shutil

# Dimensions and the time step #
time=arange(0,2401,1.)
depth=arange(0,101,1)
h=depth[-1]-depth[0]; hh=h*h
k=0.01
folder = os.path.abspath("10")

# Declaration of concentration arrays #
c1=zeros((len(time),len(depth)))
c2=zeros((len(time),len(depth)))

# Initial conditions #
for j in range(len(depth)):
	z=depth[j]
	c1[0,j]=(1+cos(pi*z/h))+1/2.
	#c2[0,j]=(1+sin(pi*z/h))+1/2.
	c2[0,j]=2.0
c1[0,0]=c1[0,1]; c1[0,-1]=c1[0,-2]
c2[0,0]=c2[0,1]; c2[0,-1]=c2[0,-2]



a=0.1 #0.1/1/10

# A loop through time and depth
fig, ax = plt.subplots(figsize=(16, 10))
#ax.set_title('Wave propagation',family='Courier New, monospace',fontsize=20, y=0.88)
plt.xlim(-0.1,3.0)
plt.ylim(0,100)
plt.xlabel('Concentration [units]', fontsize=16)
plt.ylabel('Depth [m]', fontsize=16)
lines=[]
lines.append(ax.plot(c1[0],depth,   'r',label='C1',linewidth=3))
lines.append(ax.plot(c2[0],depth,  'g', label='C2',linewidth=3))
legend(loc=3)
ax.lines
file_name = os.path.abspath(folder+"/tmp"+str(0)+".png")
ax.annotate('time: %s' %(0), xy=(0,0),  xycoords='figure fraction',xytext=(910, 100), family='Courier New, monospace',textcoords='offset points',ha="left", va="bottom",fontsize=24, bbox=dict(facecolor='lightgrey', edgecolor='none', pad=5.0))
fig.savefig(file_name, dpi=100)
ax.lines.pop(0);ax.lines.pop(0)
for i in range(len(time)-1):
	for j in range(1,len(depth)-1):
		dz=depth[j]-depth[j-1]; zz=dz*dz
		dt=zz/k/1.99 #1.99 2.0
		p1=c1[i,j+1]+c1[i,j-1]-2*c1[i,j]
		p2=c2[i,j+1]+c2[i,j-1]-2*c2[i,j]
		# Solving the system of differential equations
		c1[i+1,j] = ((a*k/hh) *c1[i,j]*(1.-c2[i,j]) + k*p1/zz)*dt + c1[i,j]
		c2[i+1,j] = ((-1.*a*k/hh) *c2[i,j]*(1.-c1[i,j]) + k*p2/zz)*dt + c2[i,j]
		# Apply boundary conditions
		c1[i+1,0]=c1[i+1,1]; c1[i+1,-1]=c1[i+1,-2]
		c2[i+1,0]=c2[i+1,1]; c2[i+1,-1]=c2[i+1,-2]
	print i, sum(c1[i+1]), sum(c2[i+1])
	if i%20==0:
		lines=[]
		lines.append(ax.plot(c1[i],depth, 'r', label='C1',linewidth=3))
		lines.append(ax.plot(c2[i],depth, 'g', label='C2',linewidth=3))
		ax.lines
		file_name = os.path.abspath(folder+"/tmp"+str(i/20)+".png")
		ax.annotate('time: %s' %(i), xy=(0,0),  xycoords='figure fraction',xytext=(910, 100), family='Courier New, monospace',textcoords='offset points',ha="left", va="bottom",fontsize=24, bbox=dict(facecolor='lightgrey', edgecolor='none', pad=5.0))
		fig.savefig(file_name, dpi=50)
   		ax.lines.pop(0);ax.lines.pop(0)

file_name = os.path.abspath("10//tmp%d.png")
os.system("ffmpeg -framerate 20/1 -f image2 -y -i "+file_name+" -r 24 -bit_rate 1800 -vb 20M A01_error.mpeg")


# Plotting of results
plt.subplot(211)
a=plt.contourf(time,depth,c1.T)
cbar1=plt.colorbar(a)
cbar1.set_clim(-0.4, 3.1)
plt.ylabel('Depth [m]')
gca().invert_yaxis()
plt.title('C1')
plt.subplot(212)
b=plt.contourf(time,depth,c2.T)
cbar2=plt.colorbar(b)
cbar2.set_clim(-0.4,3.1)
plt.xlabel('Time [s]')
plt.ylabel('Depth [m]')
gca().invert_yaxis()
plt.title('C2')

plt.show()
