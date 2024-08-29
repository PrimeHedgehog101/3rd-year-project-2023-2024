from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
from matplotlib.collections import LineCollection
###############################################################################
#N=8
Positions = [(-1,0),(1,0),(1,2),(-1,2)] # (x,y)
Strengths = [(0,-5.92384392),(5.92384392,0),(0,5.92384392),(-5.92384392,0)] # (G_r,G_i)
Positions_p = [(1.,0),(2.00001,1),(3,0)] # Peturbed positions

# Positions = [(0,1),(0.9483,0.31104),(0.5843,-0.809),(-0.5843,-0.809),(-0.9483,0.31104)] # (x,y)
# Strengths = [(1,0),(1,0),(1,0),(1,0),(1,0)] # (G_r,G_i)

# Positions = [(0,-1),(0,1),(1,-2),(1,2)] # (x,y)
# Strengths = [(-1,0),(1,0),(-1,0),(1,0)] # (G_r,G_i)

Nstep = 9000 # Used in Runge-Cutter
h     = 0.005

x0_min = -10 # creates the min/max for the mesh grid used to plot streamlines
x0_max = -10

y0_min = -28
y0_max = 6

Num_of_blue_lines_squared = 30 #Number of streamlines (squared)

listx = np.linspace(x0_min,x0_max,Num_of_blue_lines_squared) #Creates mesh grid for streamlines
listy = np.linspace(y0_min,y0_max,Num_of_blue_lines_squared)
#listy = np.delete(listy,7)
listx = np.linspace(x0_min,x0_max,1)

listx2 = [1.1,0.9,1,1] # Manual initial positions for flow lines in fig1
listy2 = [0,0,0.1,-0.1]

t = np.linspace(0,30,20000)

#Manual_IC = 'True' # If true: uses manual positions, if false: use mesh grid
Manual_IC = 'False'

#colour_plot = 'True'
colour_plot = 'False'

N = len(Positions)

stream = 'on'
#stream = 'off'

stream_power = 1
stream_angle = np.pi/4

###############################################################################

a=[i//N * 2 for i in range(N**2)]
b=[(i%N) * 2 for i in range(N**2)]
c=list(np.array(a)+1)
d=list(np.array(b)+1)
g=list(np.zeros(N*N))
D=list(np.zeros(2*N))
J=list(range(0,2*N,2))
G=[i for j in Strengths for i in j]
G1=[G[2*i]+G[2*i+1]*1j for i in range(N)]
x0 = [i for j in Positions for i in j]
x0 = x0[:2*N]


# This function sets up the problem to the N vortices problem
def odes2(x,t):
    # Defining new variable to lessen code
    for i in range(0,N*N):
            g[i]=(x[a[i]]**2)+(x[b[i]]**2)+(x[c[i]]**2)+(x[d[i]]**2)-(2*x[a[i]]*x[b[i]])-(2*x[c[i]]*x[d[i]]) 
    
    # Reproducing the summation in terms of x,y
    for j in J:
        Sum=0
        Sum2=0
        for i in range(0,N):
            if g[int(i+(N*j/2))]!=0:
                p=(x[j+1]-x[2*i+1])*(1/(g[int(i+(N*j/2))]))*(-1/(2*np.pi))*G[2*i] #y
                q=-(x[j]-x[2*i])*(1/(g[int(i+(N*j/2))]))*(-1/(2*np.pi))*G[2*i+1] #x
                Sum += p+q 
            else:
                None
            D[j]=Sum 
            if stream == 'on':
                D[j]=D[j]+(stream_power*np.cos(stream_angle))
            if stream == 'off':
                None
        for i in range(0,N):
            if g[int(i+(N*j/2))]!=0:
                p=(x[j]-x[2*i])*(1/(g[int(i+(N*j/2))]))*(1/(2*np.pi))*G[2*i] #x
                q=(x[j+1]-x[2*i+1])*(1/(g[int(i+(N*j/2))]))*(1/(2*np.pi))*G[2*i+1] #y
                Sum2 += p+q 
            if stream == 'off':
                None
            D[j+1]=Sum2 
            if stream == 'on':
                D[j+1]=D[j+1]+(stream_power*np.sin(stream_angle))
            else:
                None
             
    return D  


# Plots the path of each of the N vortices 
def plot():
    x = odeint(odes2,x0,t)
    fig, ax = plt.subplots(figsize=(12,8),dpi=100) #(8,6)
    #plt.xlabel('x')
    #plt.ylabel('y')
    #plt.grid()
    plt.tick_params(axis='both', which='both', direction='in', length=3, width=0.7)
    plt.tick_params(axis='x', which='both', top=True)
    plt.tick_params(axis='y', which='both', right=True)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    #plt.axis([-5,350,-6.2,6.2])
    
    if colour_plot == 'True':
        for i in range(N):
            points = np.array([x[:,2*i],x[:,2*i+1]]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            dt = t[1]-t[0]
            dydt = np.gradient(x[:,2*i+1], dt)
            dxdt = np.gradient(x[:,2*i], dt)
            grad_magnitude = np.sqrt(dydt**2+dxdt**2)
            #norm = plt.Normalize(grad_magnitude.min(), grad_magnitude.max())
            norm = plt.Normalize(0, 0.1)
            lc = LineCollection(segments, cmap='plasma', norm=norm)
            lc.set_array(grad_magnitude)
            lc.set_linewidth(3)
            line = ax.add_collection(lc)
            plt.plot(x0[2*i],x0[2*i+1],marker='o',color='black')
            #cbar.ax.tick_params(labelsize=10)
        cbar=fig.colorbar(line, ax=ax)
        cbar.ax.tick_params(labelsize=18)
        
    if colour_plot == 'False':
        for i in range(N):
            plt.plot(x[:,2*i],x[:,2*i+1],color='black')
            plt.plot(x0[2*i],x0[2*i+1],marker='o',color='black')
            #xxx1 = np.linspace(0,1,50)
            #yyy1 = np.linspace(1,1,50)
            #yyy2 = np.linspace(0,0,50)
            
        
        
        #plt.plot(xxx1,yyy1,linestyle='--',color='black',alpha=1,linewidth=0.8,dashes=(10,6))
        #plt.plot(xxx1,yyy2,linestyle='--',color='black',alpha=1,linewidth=0.8,dashes=(10,6))
        #plt.plot(yyy1,xxx1,linestyle='--',color='black',alpha=1,linewidth=0.8,dashes=(10,6))
        #plt.plot(yyy2,xxx1,linestyle='--',color='black',alpha=1,linewidth=0.8,dashes=(10,6))
        
            #print(np.sqrt((x[:,2*i][1001]-x0[2*i])**2+(x[:,2*i+1][1001]-x0[2*i+1])**2))
    #print(np.gradient(x[:,1], dt))
    #fig.savefig("Vortex_h2.pdf", bbox_inches='tight')
    #fig.savefig("Vortex_h3.jpg", bbox_inches='tight')
    #print(max(grad_magnitude))
    
# Plots the x position against time. Shows periodicity of trajectory of plot()
def plot_x():
    x = odeint(odes2,x0,t)
    fig, ax = plt.subplots(figsize=(8,8),dpi=200) #(8,6)
    #plt.xlabel('x')
    #plt.ylabel('y')
    #plt.grid()
    plt.tick_params(axis='both', which='both', direction='in', length=3, width=0.7)
    plt.tick_params(axis='x', which='both', top=True)
    plt.tick_params(axis='y', which='both', right=True)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.axis([-1.5,2.5,-1.5,2.5])
    
    plt.plot(t,x[:,0],color='black')
    plt.plot(t,x[:,2],color='red')
    plt.plot(t,x[:,4],color='blue')
    #plt.plot(x0[2*i],x0[2*i+1],marker='o',color='black')
    
# plots the distance between z_1 and z_2 over time
def plot_dis():
    x = odeint(odes2,x0,t)
    fig, ax = plt.subplots(figsize=(8,8),dpi=200) #(8,6)
    #plt.xlabel('x')
    #plt.ylabel('y')
    #plt.grid()
    plt.tick_params(axis='both', which='both', direction='in', length=3, width=0.7)
    plt.tick_params(axis='x', which='both', top=True)
    plt.tick_params(axis='y', which='both', right=True)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.axis([-1.5,2.5,-1.5,2.5])
    
    plt.plot(np.sqrt((x[:,0]-x[:,2])**2+(x[:,1]-x[:,3])**2),color='black')
    #plt.plot(t,x[:,2],color='red')
    #plt.plot(t,x[:,4],color='blue')
    #plt.plot(x0[2*i],x0[2*i+1],marker='o',color='black')

def Ham(time):
    x = odeint(odes2,x0,t)
    H_v = 0
    for i in range(N):
        for j in range(N):    
            if i != j:
                H_v += (-1/(2*np.pi))*G1[j]*G1[i]*np.log(abs((x[:,2*j][time]+x[:,2*j+1][time]*1j)-(x[:,2*i][time]+x[:,2*i+1][time]*1j)))
                #print((x[:,k][time]+x[:,k+1][time]*1j)-(x[:,2*i][time]+x[:,2*i+1][time]*1j))
    return H_v
        
        
        


# Sets up the x,y velocity feilds for N vortices assuming vortices are stationary
def velfun(x,y,t): 
    xvelocity = 0
    yvelocity = 0
    if stream == 'off':
        for i in range(N): 
            gamma = x**2 + y**2 + x0[2*i]**2 + x0[2*i+1]**2 - 2*x*x0[2*i] - 2*y*x0[2*i+1]
            xvelocity += -1/(2*np.pi)*(1/gamma)*(G[2*i]*(y-x0[2*i+1])-G[2*i+1]*(x-x0[2*i]))
            yvelocity += 1/(2*np.pi)*(1/gamma)*(G[2*i]*(x-x0[2*i])+G[2*i+1]*(y-x0[2*i+1]))
    else:
        for i in range(N): 
            gamma = x**2 + y**2 + x0[2*i]**2 + x0[2*i+1]**2 - 2*x*x0[2*i] - 2*y*x0[2*i+1]
            xvelocity += -1/(2*np.pi)*(1/gamma)*(G[2*i]*(y-x0[2*i+1])-G[2*i+1]*(x-x0[2*i])) 
            yvelocity += 1/(2*np.pi)*(1/gamma)*(G[2*i]*(x-x0[2*i])+G[2*i+1]*(y-x0[2*i+1]))
        xvelocity = xvelocity + (stream_power*np.cos(stream_angle))
        yvelocity = yvelocity + (stream_power*np.sin(stream_angle))
    return xvelocity, yvelocity 


# Carrys out Runge-Cutter algorithm
def loop(x0,y0,h): 
    
    t = np.linspace(0,Nstep*h,Nstep) 
    x = np.zeros(Nstep) 
    y = np.zeros(Nstep) 
    x[0]=x0
    y[0]=y0    
    
# Runge-Kutter 2 method:
    
    for k in range(1,Nstep): 
            
        t1    = t[k-1] 
        x1   = x[k-1] 
        y1    = y[k-1] 
        ux,uy = velfun(x1,y1,t1) 
        xp    = x[k-1] + h*ux 
        yp    = y[k-1] + h*uy 
        
        uxp,uyp = velfun(xp,yp,t1+h) 
        
        uxa = 0.5*(ux + uxp) 
        uya = 0.5*(uy + uyp)  
        x[k] = x[k-1] + h*uxa 
        y[k] = y[k-1] + h*uya 
        
# Stops chaos like behavour by restricting computation near certain points;
    
        if -1.08<x[k]<-0.92 and -0.08<y[k]<0.08:
            
            # This ensures both arrays the same size
            x=np.trim_zeros(x,'b')
            y=np.trim_zeros(y,'b')
            break
        
        # if -1.05<x[k]<-0.95 and 1.95<y[k]<2.05:
            
        #     # This ensures both arrays the same size
        #     x=np.trim_zeros(x,'b')
        #     y=np.trim_zeros(y,'b')
        #     break
        
        
    for o,b in zip(x[10:],y[10:]):
        if abs(o-x0)<0.009 and abs(b-y0)<0.009:
            positionx = list(x).index(o)
            positiony = list(y).index(b)
            x = x[:positionx]
            y = y[:positiony]
            break

    return x,y


# Plots streamlines of stationary vortices
def fig1(): 
    fig1, ax = plt.subplots(figsize=(8,6),dpi=200)
    if Manual_IC == 'False':
        for i in listx: 
            for j in listy: 
                
                # Stops division by zero at (0,0)
                if not (i==1 and j==0):
                    #Calls the function which produces x,y arrays and plots them.
                    x,y=loop(i,j,h) 
                    plt.plot(x,y,'black',linewidth=1.2) 
            
        for i in range(N):
            plt.plot(x0[2*i],x0[2*i+1],marker='o',color='black')
      
    if Manual_IC == 'True':
        for i,j in zip(listx2,listy2): 
                
            # Stops division by zero at (0,0)
            if not (i==1 and j==0):
                #Calls the function which produces x,y arrays and plots them.
                    x,y=loop(i,j,h) 
                    plt.plot(x,y,'black',linewidth=1.2) 
            
        for i in range(N):
            plt.plot(x0[2*i],x0[2*i+1],marker='o',color='black')
        
    
    # Graphing design.
    xx = [-0.5843,0.5843]
    xx2 = np.linspace(0.5843,0.9483,10)
    xx3 = np.linspace(-0.9483,-0.5843,10)
    xx4 = np.linspace(-0.9483,0,10)
    xx5 = np.linspace(0,0.9483,10)
    yy1 = [-0.809,-0.809]
    yy2 = 3.077*xx2 - 2.606
    yy3 = -3.077*xx3 - 2.606
    yy4 = 0.7265*xx4 + 1
    yy5 = -0.7265*xx5 + 1
    #plt.plot(xx,yy1,linestyle='--',color='black',alpha=1,linewidth=0.8,dashes=(10,6))
    #plt.plot(xx2,yy2,linestyle='--',color='black',alpha=1,linewidth=0.8,dashes=(10,6))
    #plt.plot(xx3,yy3,linestyle='--',color='black',alpha=1,linewidth=0.8,dashes=(10,6))
    #plt.plot(xx4,yy4,linestyle='--',color='black',alpha=1,linewidth=0.8,dashes=(10,6))
    #plt.plot(xx5,yy5,linestyle='--',color='black',alpha=1,linewidth=0.8,dashes=(10,6))
    #plt.xlabel(r'$x$') 
    #plt.ylabel(r'$y$') 
    #ax.spines['top'].set_visible(False)
    #ax.spines['right'].set_visible(False)
    #plt.grid() 
    plt.axis([-10,10,-4,6]) 
    plt.tick_params(axis='both', which='both', direction='in', length=3, width=0.7)
    plt.tick_params(axis='x', which='both', top=True)
    plt.tick_params(axis='y', which='both', right=True)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.title('Streamlines of flow')
    fig1.savefig("Vortex_x2.pdf", bbox_inches='tight')
    

###############################################################################
# A=[]
# for k in range(0,N):
#     B=[]
#     for i in range(0,N):
#         if i==k:
#             B.append(0)
#         else:
#             B.append(1/((x0[2*k]-x0[2*i])+(x0[2*k+1]-x0[2*i+1])*1j))
#     A.append(B)
# A = np.array(A)
# eigv , eigV = LA.eig(A)

# eigV=eigV.T
# for val, vec in zip(eigv,eigV):
#     assert np.allclose(np.dot(A,vec),val*vec)

# def near(a,b,rtol=1e-5,atol=1e-8):
#     return np.abs(a-b)<(atol+rtol*np.abs(b))
 
# eigV_01 = eigV[near(eigv,0)]
# eigV_0=[]
# for i in eigV_01:
#     for j in i:
#         eigV_0.append(j)
# eigV_0R = np.around(eigV_0,6) 

# if np.size(eigV_0) == N:
    
#     eigV_0N = 1/(eigV_0[N-1]) * np.array(eigV_0)
#     eigV_0NR = np.around(eigV_0N,6) 
# else:
#     eigV_0N = "   There is no zero eigenvalue :(   "
#     eigV_0NR = "   There is no zero eigenvalue :(   "

#This solves stationary equilbria in a uniform stream, A is normal matrix and B 
#is uniform stream matrix with factor of -2pi i.
A=[]
for k in range(0,N):
    B=[]
    for i in range(0,N):
        if i==k:
            B.append(0)
        else:
            B.append(1/((x0[2*k]-x0[2*i])+(x0[2*k+1]-x0[2*i+1])*1j))
    A.append(B)
A = np.array(A)
B = 2*np.pi*1j*-stream_power*(np.cos(stream_angle)-1j*np.sin(stream_angle))*np.array([1,1,1,1])
solution = np.linalg.solve(A,B)

# A=[]
# for k in range(0,N):
#     B=[]
#     for i in range(0,N):
#         if i==k:
#             B.append(0)
#         else:
#             B.append(1/((x0[2*k]-x0[2*i])+(x0[2*k+1]-x0[2*i+1])*1j))
#     A.append(B)
# A = np.array(A)
# B = 2*np.pi*1j*np.array([-1j+2,-1+2,1j+2,1+2])
# solution = np.linalg.solve(A,B)
###############################################################################
AA=[]
for j in range(N):
    sum3=0
    BB=[]
    for i in range(N):
        if i==j:
            for k in range(N):
                if j != k:
                    w = (G[2*k]+G[2*k+1]*1j)/(x0[2*j]-x0[2*k]+(x0[2*j+1]-x0[2*k+1])*1j)**2
                    sum3 += w
                else:
                    None
            BB.append(sum3)
        else:
            q = -(G[2*i]+G[2*i+1]*1j)/(x0[2*j]-x0[2*i]+(x0[2*j+1]-x0[2*i+1])*1j)**2
            BB.append(q)
    AA.append(BB)

eigv2 , eigV2 = LA.eig(1j*np.array(AA))  
eigV2 = eigV2.T               
eigv2R = np.around(eigv2,10)     
eigV2R = np.around(eigV2,1) 
AA = 1j*np.array(AA)
###############################################################################

def SE(): # Stationary equilbruim, finds strengths s.t SE and plots motion. Allows pertubation.
    QQ=[]
    for i in range(N):
        QQ.append(np.real(eigV_0[i]))
        QQ.append(np.imag(eigV_0[i]))
    G=QQ
    x0p = [i for j in Positions_p for i in j]
    def odes2(x,t):
        
        for i in range(0,N*N):
                g[i]=(x[a[i]]**2)+(x[b[i]]**2)+(x[c[i]]**2)+(x[d[i]]**2)-(2*x[a[i]]*x[b[i]])-(2*x[c[i]]*x[d[i]])

        for j in J:
            Sum=0
            Sum2=0
            for i in range(0,N):
                if g[int(i+(N*j/2))]!=0:
                    p=(x[j+1]-x[2*i+1])*(1/(g[int(i+(N*j/2))]))*(-1/(2*np.pi))*G[2*i] #y
                    q=-(x[j]-x[2*i])*(1/(g[int(i+(N*j/2))]))*(-1/(2*np.pi))*G[2*i+1] #x
                    Sum += p+q
                else:
                    None
                D[j]=Sum
            for i in range(0,N):
                if g[int(i+(N*j/2))]!=0:
                    p=(x[j]-x[2*i])*(1/(g[int(i+(N*j/2))]))*(1/(2*np.pi))*G[2*i] #x
                    q=(x[j+1]-x[2*i+1])*(1/(g[int(i+(N*j/2))]))*(1/(2*np.pi))*G[2*i+1] #y
                    Sum2 += p+q
                else:
                    None
                D[j+1]=Sum2
                 
        return D  
    x = odeint(odes2,x0p,t) 
    
    fig = plt.figure(dpi=200)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid()
    #plt.axis([5.99,6.01,4.99,5.01])
    for i in range(N):
        plt.plot(x[:,2*i],x[:,2*i+1],color='black')
        if x0[2*i]-x0p[2*i]!=0 or x0[2*i+1]-x0p[2*i+1]!=0:
            plt.plot(x0[2*i],x0[2*i+1],marker='o',color='blue')
            plt.plot(x0p[2*i],x0p[2*i+1],marker='o',color='black')
        else:
            plt.plot(x0[2*i],x0[2*i+1],marker='o',color='black')
            


    



#if 0.95<x[k]<1.05 and -0.05<y[k]<0.05:
    
   # # This ensures both arrays the same size
   # x=np.trim_zeros(x,'b')
   # y=np.trim_zeros(y,'b')
    
    #break


    
        