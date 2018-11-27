# -*- coding: utf-8 -*-
"""
Created on Sun Apr 08 10:24:36 2018

@author: Xinwei
"""

from numpy import *

# initial state
def iniWaveFcn(hdim):    
    y = zeros(hdim+1) 
    y[0] = 1.0 
    return y
 
# basis of sub-Hilbert space   
def getBasis(N,m):

	 hdim = (N-abs(m))/2+1
    n1vec = []
    n0vec = []
    nm1vec = []
    if m>=0:        
        for ii in range(0,hdim):
            nm1vec.append(ii)
            n1vec.append(ii + m)
            n0vec.append(N - 2*ii - m)
    else:
        for ii in range(0,hdim):
            n1vec.append(ii)
            nm1vec.append(ii-m)
            n0vec.append(N-2*ii+m)
            
    return n1vec,n0vec,nm1vec
        
   
# Hamiltonian in each subspace	
def getHam(N,m):
    
    hdim = (N-abs(m))/2+1
    diag = zeros(hdim+1)
    eiag = zeros(hdim) 
    
	#diagonal elements	
    for ii in range(0, hdim):
        nm1 = nm1vec[ii]           		
		  n1 = n1vec[ii]
		  n0 = n0vec[ii]
		  diag[ii] = c/2.d0/N*((2*n0-1)*(n1+nm1)+(n1-nm1)**2)

	#off-diagonal	
    for ii in range(1, hdim):
		  nm1 = nm1vec[ii]
		  n1 = n1vec[ii]
		  n0 = n0vec[ii]
		  eiag[ii-1] = c/N*sqrt(float(n0+1)*(n0+2)*n1*nm1)	

    return diag,eiag

    

# probability of jump  
def decayPro(N,m,y,n1vec,n0vec,nm1vec,dt):
    
    hdim = (N-abs(m))/2+1
    
    tmp = map(abs, y)
    n0mean = n0vec*tmp**2
    n1mean = n1vec*tmp**2
    nm1mean = nm1vec*tmp**2
            
    dp1 = gamma*n1mean*dt
	 dp0 = gamma*n0mean*dt
	 dpm1 = gamma*nm1mean*dt
    return dp1,dp0,dpm1

# quantum jump
def jumpProcess(flag,y,N,m):
    
   tmp = y
	
	if(flap==1):
		m = m-1
		loss = n1vec		
	elif(flap==0):
		loss = n0vec		
	else:
		m = m+1
		loss = nm1vec

	
	ntot = ntot-1
	hdim = (ntot-abs(m))/2+1
	
	if(hdim<hdim0) then
	
		if(flap.ne.0) then
			tmp(1:hdim) = y(2:hdim0)
			tmp1 = loss
			loss(1:hdim) = tmp1(2:hdim0)
	
	Getbasis(N,m)		
	y = zeros(hdim)
   y = tmp*sqrt(loss)
   
   return N,m,y
    
# non-Hermitian evolution    
def nonHermiEov(psi,ham):
    # diagonalization
    evec = eig(ham)
    psi = evec*psi
    psi = exp(-1i*eval*dt-gamma/2.0)*psi
    psi = evec*psi
    return psi


# Monte carlo wave function method   
def mcwfEvo(ntot):

    
    tmax = 2*pi*5.d0
    dt = 1.d-3
    tstep = tmax/dt+1
    y = init()
    Getbasis(N,m)
    q = qmax
    ham = GetHam(hdim)
    rr = np.random.rand([tstep,3])
    N = ntot
    m = 0
    
    for tidx in range(tstep+1):
        		r1 = rr[tidx,1]
		      r0 = rr[tidx,2]
		      rm1 = rr[tidx,3]
            energy()  
            dp1,dp2,dp3 = decay_pro(n1mean,n0mean,nm1mean)
            
            if(r1>dp1.and.r0>dp0.and.rm1>dpm1):
                time = tidx*dt
                q = qmax*(1.d0-2.d0*time/tmax)
                GetHam(N,m)
                N,m,psi = NonHermiEov()
            elif(r1<dp1):
                N,m,psi = jump(1)	
            elif(r2<dp2):
                N,m,psi = jump(0)
            elif(r3<dp3):
                N,m,psi = jump(-1)
                
    		   norm = sum(abs(psi)**2)
		    psi = psi/sqrt(norm)
              
    psif = psi
    return psif
                                        