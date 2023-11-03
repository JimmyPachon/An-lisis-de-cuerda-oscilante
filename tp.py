



import matplotlib.animation as animation
import math
import numpy as np
from matplotlib import pyplot as plt
from scipy.linalg import toeplitz
import pylab
from numpy import linalg as la


#c es la raíz del cociente entre la tensión y la densidad 
c=1.0
#L es la longitud de la cuerda en el reposo
L=1.0


#N es el número de partes móviles
N=50



#h es la distancia horizontal entre partes móviles para N masas
h=L/(N+1)


#--------------------- Item a ----------------------------------------

#Arma una matriz tridiagonal, B para el caso de extremo fijo en función del número de masas

def B(X):
    
  a = np.zeros(X)
  a[0] = 2
  a[1] = -1
  first_row1 = np.array(a)
  first_col1 = np.array(a)
  B = toeplitz(first_col1, first_row1)
  
  return B

def C(X):

   #C será la matriz tridiagonal para el caso de extremo libre.
   b = np.zeros(X)
   b[0] = 2
   b[1] = -1
   first_row2 = np.array(b)
   first_col2 = np.array(b)
   C = toeplitz(first_col2, first_row2)
   #Se modifica la ultima fila por la condición de extremo libre
   (C[X-1])[X-1]=1.0
   
   return C



    
#Función que da una lista de las frecuencias angulares dada una matriz y el número de partes móviles
def W(X,Y):
    # Halla los autovalores y autovectores de la matriz
    mu, v = la.eig(X)
    
    # w da una lista con las frecuencias angulares 
    w=np.zeros(Y)
    
    for i in range(0,Y):
       w[i]=(c*(Y+1)/L)*math.sqrt(mu[i])
   
    return w



#Función que da una lista de las longitudes de onda dada una matriz y el número de partes móviles

def Lam(X,Y):
    
    w=W(X,Y)
    #lamb da una lista de longitudes de onda
    lamb=np.zeros(Y)
    
    for i in range(0,Y):
     lamb[i]=((math.pi)*2*c)/w[i]
     
    return lamb


# --------------------------- Item b ----------------------------------
#EXTREMO CERRADO

#min_lamb es la lista con los valores mínimos de longitud de onda para cada N


min_lamb=np.zeros(N-1)

for i in range(2,N+1):
    
      min_lamb[i-2]=min(Lam(B(i),i))
      

# h_aux es una lista con las distancias entre partes móviles para cada N

h_aux = np.zeros(N-1)


for i in range(2,N+1):
    
    h_aux[i-2]= L/(i+1)   
    


#Gráficos pedidos

x= h_aux
y= min_lamb

z=np.zeros(N-1)

for i in range(0,N-1):
    
    z[i]= math.pi*h_aux[i]


plt.subplot(221)
plt.plot(x,y,'o',linewidth=1.0)
plt.plot(x,z,'g',linewidth=1.0)
plt.xlabel('h')
plt.ylabel('Lambda mínimo')
plt.title('Extremos')
plt.legend(['Cerrado','pi*h'])

    


#EXTREMO ABIERTO

#min_lamb2 es la lista con los valores mínimos de longitud de onda para cada N



min_lamb2=np.zeros(N-1)

for i in range(2,N+1):
    
      min_lamb2[i-2]=min(Lam(C(i),i))



#Gráficos pedidos

y2= min_lamb2

plt.subplot(223)
plt.plot(x,y2,'o',linewidth=1.0)
plt.plot(x,z,'g',linewidth=1.0)
plt.xlabel('h')
plt.ylabel('Lambda mínimo')
plt.legend(['Abierto','pi*h'])

# Acá se compara 2L/(N+1) con 2L/(N) en función de N

#N_fila es una lista con números que van de 1 hasta N

N_fila=np.zeros(N)

for i in range(1,N+1):
  
    N_fila[i-1]=i

#lambda_N es una lista con los valores lamba=2L/N

lambda_N=np.zeros(N)

for i in range(1,N+1):
    
    lambda_N[i-1]=2*L/(N_fila[i-1])
    
#lambda_N1 es una lista con los valores lamba=pi*L/(N+1)
    
lambda_N1=np.zeros(N)
    
for i in range(1,N+1):
    
    lambda_N1[i-1]=math.pi*L/(N_fila[i-1]+1)
    
    
plt.subplot(222)
plt.plot(N_fila,lambda_N,'b',linewidth=1.0)
plt.plot(N_fila,lambda_N1,'r',linewidth=1.0)
plt.xlabel('N')
plt.ylabel('Lambda')
plt.legend(['2L/N','pi*L/(N+1)'])





#---------------------------------- Item c------------------------------------------

#ln es la lista de longitudes de onda modales para n masas en el caso de extremos cerrados
l2=np.sort(Lam(B(2),2))
l3=np.sort(Lam(B(3),3))
l5=np.sort(Lam(B(5),5))
l7=np.sort(Lam(B(7),7))
l11=np.sort(Lam(B(11),11))

#Con este código se toman los cocientes para cada caso.
#for i in range(0,len(l11)):
#     print(l3[2]/l11[i])
     
     
#Gráfico para representar el aliasing espacial
   
x2= np.array(range(400))*0.01
k= np.zeros(len(x2))
g= np.zeros(len(x2))

for i in range(len(x2)):
    k[i]=np.sin((2*math.pi/l2[0])*x2[i])
    
for i in range(len(x2)):
    g[i]=np.sin((2*math.pi/l5[1])*x2[i])
    
plt.subplot(224)
plt.plot(x2,k,'b',linewidth=2.0)
plt.plot(x2,g,'r',linewidth=2.0)


#--------------------------------- Item d --------------------------------

c=1.0
N=5
n= np.zeros(N)
h=L/(N+1)


for i in range(1,N+1):
    n[i-1]=i
    


nh= np.zeros(N)

for i in range(0,N):
    nh[i]=h*n[i]
    
    
#Matriz identidad NxN multiplicada por 2
I2 = [] 
for i in range(N):
	fila = []
	for j in range(N):
		val = 2 if i==j else 0
		fila.append(val)
	I2.append(fila)

    

#Psi0 son las condiciones iniciales de desplazamiento


Psi0=np.zeros(N)

for i in range(0,N):
    
    if nh[i]<=(L/2) :
        
        Psi0[i]=0.2*nh[i]/L
    
    else:
        
        Psi0[i]=0.2*(1-(nh[i]/L))
        

        
    
        
t= 0.0
tmax=10.0
tau=0.05




r2 = (c**2)*(tau**2)/(h**2)

Psi1 = Psi0





Psi_aux=np.zeros((int(tmax//tau)+1,N))    

while (t<=tmax):
    
                Psi_aux[int(t//tau)]= Psi1
                Psi_new = (I2-(r2)*B(N)).dot(Psi1)-Psi0
                Psi0 = Psi1
                Psi1 = Psi_new
                t = t + tau

        
Psi_extremos=np.zeros((int(tmax//tau)+1,N+2))

for i in range(0,len(Psi_aux)):
    
    for j in range (1,N+1):
        
        Psi_extremos[i][j]=Psi_aux[i][j-1]
        
        
    Psi_extremos[i][0]=0   
    Psi_extremos[i][N+1]=0
    


# Animación 1 , sólo activar una a la vez
                
fig=plt.figure()
ax=fig.gca()
    
nh_extremos=np.zeros(len(nh)+2)

for i in range(0,len(nh)):
    nh_extremos[i+1]= nh[i]
    
nh_extremos[len(nh)+1]=1.0



def actualizar(i):
        
   ax.clear()
   
   ax.plot(nh_extremos,Psi_extremos[i])
   
        
   plt.title("Animación de las 5 cuentas")
   
   plt.ylim(-0.2,0.2)
   plt.xlim(-0.1,1.1)
   

ani=animation.FuncAnimation(fig,actualizar,range(int(tmax//tau)))
plt.show()

"""
#--------- Masa 3 en función del tiempo ----------------------------------

# Animación 2, sólo dejar activa una a la vez

fig2=plt.figure()
ax2=fig2.gca()
  
Psi_tercer_masa=np.zeros(int(tmax//tau))

for i in range(0,len(Psi_tercer_masa)):
    Psi_tercer_masa[i]=Psi_extremos[i][3]


t_lista=np.zeros(int(tmax//tau))

for i in range(0,len(t_lista)):
    t_lista[i]= i*tau
    
    
    
def actualizar2(i):
        
   ax2.clear()
   
   ax2.plot(t_lista[:i],Psi_tercer_masa[:i])
        
   plt.title("Animación de la cuenta 3")
     
   plt.ylim(-0.2,0.2)
   plt.xlim(-0.1,tmax)

ani2=animation.FuncAnimation(fig2,actualizar2,range(int(tmax//tau)))
plt.show()
"""
