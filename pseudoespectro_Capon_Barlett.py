import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
#%matplotlib inline

def trans(matriz):
    return np.conj(matriz).T

M=2            #número de elementos en el arreglo de antenas
D=1          #número de fuentes de las que se desea el AoA
sig2=0.1       #varianza del ruido
th1=60  #ángulo de la fuente transmisora
#th2=60
th1=th1*np.pi/180
#th2=th2*pi/180
a1=1
#a2=[1]
a=1
for i in range(2,M+1):
    a1=np.matrix([a1,np.exp(1j*i*np.pi*np.sin(th1))])
    #vector principal del arreglo
 #   a2=[a2 exp(1j*i*pi*sin(th2))]
A=trans(a1)
Rss=1
Rxx= A*Rss*trans(A)+sig2*np.eye(M)
theta= range(181)
th=np.ones(180)
P_Bartlett=np.ones(180)
P_Capon=np.ones(180)
P_MUSIC= np.ones(180)

for k in range(180):
    th[k]=theta[k]*np.pi/180
    a=1
    for jj in range(2,M+1):
        a =np.matrix([a,np.exp(1j*jj*np.pi*np.sin(th[k]))])
    P_Bartlett[k]=np.real(a*Rxx*trans(a))
    P_Capon[k]=np.real(np.linalg.inv(a*np.linalg.inv(Rxx)*trans(a)))
    Dia,V=np.linalg.eig(Rxx)
    index = np.argsort(abs(Dia[0]))
    #Dia=ix*np.eye(2)
    #V[1,0]=V[1,0]*-1
    EN=V[:,0]
    P_MUSIC[k]=np.real(np.linalg.inv(np.absolute(((a*EN*trans(EN)*trans(a))))))
plt.figure(1)
plt.subplot(3,1,1)
plt.plot(th*180/np.pi,(P_Bartlett/max(P_Bartlett)),'-b')
plt.grid(True)
plt.xlabel('Angle(?)')
plt.ylabel('|P(\theta)| (dB)')
plt.title('Pseudoespectro con Algoritmo de Bartlett')
plt.figure(2)
plt.subplot(3,1,2)
plt.plot(th*180/np.pi,(P_Capon/max(P_Capon)),'-g')
plt.grid(True)
plt.xlabel('Angle(?)')
plt.ylabel('|P(\theta)| (dB)')
plt.title('Pseudoespectro con Algoritmo de Capon')
plt.figure(3)
plt.subplot(3,1,3)
plt.plot(th*180/np.pi,(P_MUSIC/max(P_MUSIC)),'-r')
plt.grid(True)
plt.xlabel('Angle(?)')
plt.ylabel('|P(\theta)| (dB)')
plt.title('Pseudoespectro con Algoritmo de MUSIC')
plt.figure(4)
plt.plot(th*180/np.pi,(P_Bartlett/max(P_Bartlett)),th*180/np.pi,(P_Capon/max(P_Capon)),th*180/np.pi,(P_MUSIC/max(P_MUSIC)))
plt.legend(['Bartlett','Capon','MUSIC'])
plt.title('Pseudoespectro con Algoritmo de MUSIC, Bartlett & Capon')
plt.grid(True)
plt.xlabel('Angle(?)')
plt.ylabel('|P(\theta)| (dB)')