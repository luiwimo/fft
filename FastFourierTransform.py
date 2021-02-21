from math import *
from numpy import *
from matplotlib import *
import numpy as np
import matplotlib.pyplot as plt
import sys 
from numpy import cos, linspace, pi
from pylab import plot, show, title, xlabel, ylabel, subplot
from scipy import fft, arange
from scipy import signal




fm=100.0
T=1.0/fm
f=5
t=np.arange(0,1,T)

xe = t + np.random.normal(size=100)
#xe=cos(2*pi*f*t)
#xe=signal.square(2 * np.pi * 5 * t)

x=[0]*len(xe)
for ks in range(len(x)):
	x[ks]=xe[ks]

y=np.fft.fft(x)



def complete(x):                #Completa el vector

	N=len(x)                    #Obtiene el tamano de x
	m=round(np.log2(N))         #Calcula el log_2 de N
	f=2**m                      #Calcula 2^m
	g=2**(m+1)                  #Calcula 2^(m+1)

	if N<=f:                    #Si la longitud del vector x
           						#es menor que 2^m               
		m=m
		while N<f:

			x.append(0)         #Agrega un cero al final de x
			N=len(x)            #mientras N sea menor que 2^m

								#Si la longitud del vector x
	else:                       #es mayor que 2^m

		m=m+1
		while N<g:

			x.append(0)         #Agrega un cero al final de x
			N=len(x)            #mientras N sea menor que 2^(m+1)
	m=int(m)                    #Convierte m a entero

	return x,N,m                #Regresa los valores x,N,m


def BitRev(xc,N,m):             #Reordena los elementos del vector
	i=0      					#Inicializa i=0 
	xo=[0]*N                    #Inicializa xo de tamano N con ceros

	for i in range(N):          #Inicia ciclo for desde i=0 hasta N
		ind=bin(i)              #bin() convierte "i" a binario 
	                            #y devuelve una cadena con terminacion "b"
		num=ind.split("b")      #corta la cadena hasta encontrar "b"
		l=num[1]                #Obtiene el elemento 1 del arreglo num que
		                        #es donde se encuentra el numero binario
		n=len(l)                #Calcula el tamano del numero binario

		while n<m:              #mientras el tamano de la cadena del numero  
			l="0"+l             #binario mas pequeno sea menor que la del 
		   	n=len(l)            #numero mas grande se agregan ceros

		o=0						#Inicializa o=0 
		lo=[0]*m;               #Inicializa lo de tamano m con ceros
		for o in range(m):      #Inicia ciclo for desde o=0 hasta m 
			lo[o]=l[o]          #Se divide la cadena y se guarda en lista
		lo.reverse()            #Esta funcion invierte el orden de la lista
 
		o=0                     #Inicializa o=0 
		ls=""                   #Se declara una cadena ls
		for o in range(m):      #Inicia ciclo for desde o=0 hasta m 
			ls=ls+lo[o]         #Se guardan los elementos de la lista
 								#en una cadena

		r=int(ls,2)             #Conversion de binario a decimal entero

		k=0						#Inicializa k=0 
		while k<=r:             #El ciclo while reordena el vector original
			if k==r:            #recibido con los nuevos indices obtenidos
				xo[i]=xc[k]	    #del bit-reversal
			k=k+1

	return xo                   #Regresa vector reordenado

def wn(N):                      #Calcula los valores de Wn necesarios
	
	wnd=[0]*(N/2)               #Inicializa wnd de tamano N/2

	for i in range(N/2):        #Calcula los valores necesarios 
	                            #de la exponencial compleja
		wnd[i]=round(np.cos((2*pi*i)/N),4)-1j*round(np.sin((2*pi*i)/N),4)

	return wnd                  #Regresa los valores calculados

def fftm(xo,wnd,N):             #Aplica el algoritmo de Danielson-Lanczos
	pnt=2                       #Se define el numero de puntos inicial 
	z=[0]*N                     #Se inicializa z con ceros de tamano N

	while pnt<=N:               #El ciclo se repite conforme pnt se incrementa
							    #en multiplos de 2 hasta llegar a N
		rep=N/pnt               #Se asegura que Puntos * Repeticiones = N
		med=pnt/2               
		nx=0
		ind=0

		while rep>0:            #Se realiza la transformada de pnt puntos rep veces
			for j in range(2):

				lm=0;
				for k in range(med):

					if lm>=N:
						lm=0
					else:  
						lm=lm
		         #A continuacion se tienen tres casos:
		             #Cuando se resta el factor de giro
		             #Cuando se suma el factor de giro
		             #Cuando el factor de giro vale 1 

					for n in range(2):
						if n==1:
							if j==1:
								z[ind]=z[ind]-(xo[k+med+nx]*wnd[lm])
							else:
								z[ind]=z[ind]+(xo[k+med+nx]*wnd[lm])
						else:
							z[ind]=xo[k+nx]
					ind=ind+1

					lm=lm+(N/(2*med))				

			rep=rep-1
			nx=nx+pnt

		xo=z
		z=[0]*N
		pnt=2*pnt
		Nxo=len(xo)
	return xo,Nxo

def GraphEspect(xora,fm):
 n = len(xora) # longitud de la senal
 k = arange(n)
 T = n/fm
 frq = k/T # 2 lados del rango de frecuancia
 frq = frq[range(n/2)] # Un lado del rango de frecuencia
 
 Y=[0]*n
 i=0
 for i in range(n): 
     Y[i]= xora[i]/n # fft calcula la normalizacion

 Y = Y[0:n/2]
 plot(frq,Y,'r') # grafica el espectro de frecuencia
 xlabel('Frecuencia (Hz)')
 ylabel('|Y(f)|')



xc,N,m=complete(x)

xo=BitRev(xc,N,m)

wnd=wn(N)

xor,Nxor=fftm(xo,wnd,N)


xora=[0]*Nxor
ya=[0]*Nxor

for a in range(Nxor):
	xora[a]=abs(xor[a])


for a in range(len(y)):
	ya[a]=abs(y[a])
 

subplot(2,2,1)

plot(t,xe)
xlabel('Tiempo')
ylabel('Amplitud')
plt.grid()

subplot(2,2,2)
GraphEspect(xora,fm)
plt.grid()

subplot(2,1,2)
GraphEspect(ya,fm)
plt.grid()

plt.show()
sys.exit()
