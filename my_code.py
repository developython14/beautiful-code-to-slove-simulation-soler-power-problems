import numpy as np
import math
from numpy import linalg
from numpy.lib.function_base import average 
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error


def rendement_electrique (Tc ,Tcrefrence=298, nrefrence = 0.17,betapv = 1/(270-25) ):
    return nrefrence*(1-betapv*(Tc-Tcrefrence))
    
def solve (Qnf): 
    Vvent = []
    with open(r'/content/vv1.txt') as f :
        x = f.readlines()
    for i in range (0,len(x)) :
        x[i] = x[i].replace(',' ,'.')
        Vvent.append(x[i][:len(x[i])-1])

  
    Vvent = [float(i) for i in Vvent ]
    Vvent = [float(i) for i in Vvent ]
    Vvent = Vvent[6:len(Vvent):5]
    Vvent[-1] = 0.08
    Vvent = np.array(Vvent)

    
    g = []
    with open(r'/content/g1.txt') as f :
       x = f.readlines()
    for i in range (0,len(x)) :
        x[i] = x[i].replace(',' ,'.')
        g.append(x[i][:len(x[i])-1])


    G= [float(i) for i in g ]
    G = np.array(G)

    Tamb = []
    with open(r'/content/tamb.txt') as f :
        x = f.readlines()
    for i in range (0,len(x)) :
        x[i] = x[i].replace(',' ,'.')
        Tamb.append(x[i][:len(x[i])-1])


    Tamb= [float(i)+273 for i in Tamb[:-1] ]
    Tamb = np.array(Tamb)
    l=0.79 
    Ltube=1.62 
    Dh=0.012 
    kc=0.036 
    kt=0.033 
    ec=0.0003 
    et=0.0005 
    kail=230 
    eail=0.01 
    cps=897 
    rops=2700    
    ktub=387 
    rotube=8960 
    ctub=385 
    etub=0.002 
    kiso=0.039 
    eiso=0.07 
    roiso=40 
    segma=5.67*10**-8 
    av=0.06 
    ac=0.85 
    at=0.5 
    tovi=0.84  
    bettapv=1/(270-25)
    Tcref=298 
    nuref=0.17
    Ntube=12 
    L=1.58 
    pi =math.pi
    sec = pi*(0.022**2)/4
    sec1 = pi*(0.012**2)/4
    V=Qnf/sec
    s=L*l 
    s2=pi*Dh*(1/4)*L*Ntube 
    s1=s-s2 
    sl=3.14*Dh*L*Ntube 
    s3=pi*Dh*((1/2)+(1/4))*L*Ntube 
    s4=(pi+2)*Dh*(1/2)*L*Ntube+s1 
    cc=677 
    rot=1200 
    ct=1250 
    cpail=897 
    roail=2700 
    cpiso=1500 
    beta=0.88 
    epsv=0.83   
    kv=1 
    roc=2702 
    rov=2530 
    cv=500 
    ev=0.0032 
    kps=230   
    eps=0.001 
    sail=eail*L*6 
    mtub=rotube*pi*(0.012/2)**2-(0.011/2)**2*Ltube*Ntube 
    mc=roc*s*ec 
    mv=rov*s*ev 
    mt=rot*s*et 
    miso=roiso*s*eiso 
    mail=roail*eail*L*0.04*6 
    mps=rops*(s4)*eps



    Tvold =[None]*167
    Tnfe= [25+273 ]*167
    Tcold=[None]*167
    Ttold=[None]*167
    Ttubeold=[None]*167
    Teauold=[None]*167
    Tpsold=[None]*167 
    Tailold=[None]*167
    Tisoold=[None]*167
    Tv=[None] *167
    Tc=[None]*167
    Tt=[None] *167
    Ttube=[None] *167
    Teau=[None] *167
    Tps=[None] *167
    Tail=[None] *167
    Tiso=[None] *167
    Tciel=[None] *167
    knf=[None] *167
    cpnf=[None] *167
    munf=[None] *167
    ronf=[None] *167
    pr= [None] *167
    Re= [None] *167
    Hr_vciel  =[None] *167
    Hcv_v =[None] *167
    Hcv_nftub=[None] *167
    Teau2 = [None] *167
    puissance = [None] *167
    ntheau = [None] *167
    Tvold[0]=320 
    Tcold[0]=310 
    Ttold[0]=300 
    Ttubeold[0]=299 
    Teauold[0]=298 
    Tpsold[0]=292 
    Tailold[0]=292 
    Tisoold[0]=292 
    mmax=96 
    t =np.arange(167)

    dt=300 
    tol=10**-5 
    err=tol+1 
    m=0 
    Tv[m]=Tvold[m] 
    Tc[m]=Tcold[m] 
    Tt[m]=Ttold[m] 
    Ttube[m]=Ttubeold[m] 
    Teau[m]=Teauold[m] 
    Tps[m]=Tpsold[m] 
    Tail[m]=Tailold[m] 
    Tiso[m]=Tisoold[m] 
    iter1=0 
    Hcd_vc=1/(ev/kv+ec/kc) #%conduction entre le verre et la cellule 
    Hcd_ct=1/(ec/kc+et/kt) #%conduction entre la cellule et tedlar
    Hcd_ttub=1/(et/kt+etub/ktub) #%conduction entre tedlar et le tube  
    Hcd_tps=1/(eps/kps+et/kt) #%conduction entre la plaque  et tedlar 
    Hcd_tubps=1/(etub/ktub+eps/kps) #%conduction entre le tube et la plaque 
    Hcd_psail=1/(eps/kps+eail/kail) #%conduction entre la plaque et l'aillette
    Hcd_psiso=1/(eps/kps+eiso/kiso) #%cnduction entre la plaque  et l'isolant
    Hcd_ailiso=1/(eail/kail+eiso/kiso)  #%conduction entre l'ailette et l'isolant
    tol=10**-5 
    err=tol+1
    while err>tol:
        iter1=iter1+1 
        m=0 
        Tciel[m]=0.0552*Tamb[m]**1.5 
        knf[m] =-7843*10**-9*Teauold[m]**2 +0.0062*Teauold[m] - 0.54 
        ronf[m] =-0.003 *Teauold[m]**2 + 1.505*Teauold[m] + 816.781 
        cpnf[m] =-0.0000643 *Teauold[m]**3 + 0.0552 *Teauold[m]**2 - 20.86 *Teauold[m] + 6719.637
        munf[m] = (2414*10**-8)*10**(247.8/(Teauold[m]-140)) 
        pr[m]=cpnf[m]*munf[m]/knf[m] 
        Re[m]=ronf[m]*V*Dh/munf[m] 
    #%les coefficient de transfert 
    #%1-RAYONNEMENT   
    #% Verre - voûte céleste
        Hr_vciel[m]=segma*epsv*(Tvold[m]**2+Tciel[m]**2)*(Tvold[m]+Tciel[m]) 
    #%2 CONDUCTION 
    #%1] le verre et l'aire    
        Hcv_v[m]=5.7+3.8*float(Vvent[m])  
    # %2 le nano-fluide et le tube
        if  Re[m]<=2300:
            Hcv_nftub[m]=4.36*knf[m]/Dh 
        elif  Re[m]>2300:
            Hcv_nftub[m]=0.023*(knf[m]/Dh)*(Re[m]**0.8)*(pr[m]**0.4) 
        Tv[m]=(Tciel[m]*Hr_vciel[m]+Tcold[m]*Hcd_vc+Tamb[m]*Hcv_v[m]+av*G[m])/(Hr_vciel[m]+Hcd_vc+Hcv_v[m]) 
        Tc[m]=(Tv[m]*Hcd_vc+Ttold[m]*Hcd_ct-nuref*(1+bettapv*Tcref)*G[m]*tovi*beta+tovi*ac*G[m]*beta)/(Hcd_vc+Hcd_ct-nuref*bettapv*G[m]*tovi*beta) 
        Tt[m]=(Tcold[m]*Hcd_ct*s+Tpsold[m]*Hcd_tps*s1+Ttubeold[m]*Hcd_ttub*s2+tovi*s*at*G[m]*(1-beta))/(Hcd_ct*s+Hcd_tps*s1+Hcd_ttub*s2) 
        Ttube[m] = (Ttold[m]*Hcd_ttub*s2 + Tpsold[m]*Hcd_tubps*s3 +Teauold[m]*Hcv_nftub[m]*sl)/(Hcd_tubps*s3 + Hcd_ttub*s2 + Hcv_nftub[m]*sl) 
        A =(Hcv_nftub[m]*sl/Ltube)/(ronf[m]*V*sec*cpnf[m]) 
        Teau[m]=Ttubeold[m]+Ttubeold[m]*(math.exp(-A*Ltube)-1)/(A*Ltube )-Tnfe[m]*(math.exp(-A*Ltube)-1)/(A*Ltube)
        Tps[m]=(Ttold[m]*Hcd_tps*s1+Ttubeold[m]*Hcd_tubps*s3+Tisoold[m]*Hcd_psiso*s4)/(Hcd_tps*s1+Hcd_tubps*s3+Hcd_psiso*s4) 
        Tail[m]=(Tpsold[m]*Hcd_psail*sail+Tisoold[m]*Hcd_ailiso*sail)/(Hcd_psail*sail+Hcd_ailiso*sail) 
        Tiso[m]=(Tpsold[m]*Hcd_psiso*s4+Tailold[m]*Hcd_ailiso*sail)/(Hcd_psiso*s4+Hcd_ailiso*sail) 
        Teau2[m]=Ttube[m]*(1-math.exp(-A*Ltube))+Tnfe[m]*math.exp(-A*Ltube) 
        puissance[m]=(ronf[m]*V*sec*cpnf[m])*(Teau2[m]-Tnfe[m]) 
        ntheau[m]=puissance[m]/(G[m]*s)
        erreur = [None]*8 
        erreur[1]=Tv[m]-Tvold[m] 
        erreur[2]=Tc[m]-Tcold[m] 
        erreur[3]=Tt[m]-Ttold[m] 
        erreur[4]=Ttube[m]-Ttubeold[m] 
        erreur[5]=Teau[m]-Teauold[m] 
        erreur[6]=Tps[m]-Tpsold[m] 
        erreur[7]=Tail[m]-Tailold[m] 
        erreur =list(map(abs,erreur[1:]))
        err=max(erreur)
        m=0
        Tvold[m]=Tv[m] 
        Tcold[m]=Tc[m] 
        Ttold[m]=Tt[m] 
        Ttubeold[m]=Ttube[m] 
        Teauold[m]=Teau[m] 
        Tpsold[m]=Tps[m] 
        Tailold[m]=Tail[m] 
        Tisoold[m]=Tiso[m] 

    m=0
    Tvold[m+1]=Tv[m] 
    Tcold[m+1]=Tc[m] 
    Ttold[m+1]=Tt[m] 
    Ttubeold[m+1]=Ttube[m] 
    Teauold[m+1]=Teau[m] 
    Tpsold[m+1]=Tps[m] 
    Tailold[m+1]=Tps[m]
    Tisoold[m+1]=Tiso[m] 
    tol=10**-5 
    for m in range(1,mmax):
        err=tol+1
        iter1 = 0
        while err>tol  :
            iter1 = iter1+1
            Tciel[m]=0.0552*Tamb[m]**1.5 
            knf[m] =-7843*10**-9*Teauold[m]**2 +0.0062*Teauold[m] - 0.54 
            ronf[m] =-0.003 *Teauold[m]**2 + 1.505*Teauold[m] + 816.781 
            cpnf[m] = -0.0000643 *Teauold[m]**3 + 0.0552 *Teauold[m]**2 - 20.86 *Teauold[m] + 6719.637 
            munf[m] = (2414*10**-8)*10**(247.8/(Teauold[m]-140)) 
            Re[m]=ronf[m]*V*Dh/munf[m] 
            pr[m]=cpnf[m]*munf[m]/knf[m] 
            Hr_vciel[m]=segma*epsv*(Tvold[m]**2+Tciel[m]**2)*(Tvold[m]+Tciel[m]) 
            Hcv_v[m]=5.7+3.8*float(Vvent[m])  
            if Re[m]<=2300 :
                Hcv_nftub[m]=4.36*knf[m]/Dh 
            elif  Re[m]>2300 :
                Hcv_nftub[m]=0.023*(knf[m]/Dh)*(Re[m]**0.8)*(pr[m]**(0.4)) 
            Tv[m]=(Tv[m-1]*mv*cv+dt*(Tciel[m]*Hr_vciel[m]*s+Tcold[m]*(Hcd_vc*s)+Tamb[m]*Hcv_v[m]*s+av*G[m]*s))/(mv*cv+dt*s*(Hr_vciel[m]+Hcd_vc+Hcv_v[m])) 
            Tc[m]=(Tc[m-1]*mc*cc+dt*(Tvold[m]*Hcd_vc*s+Ttold[m]*Hcd_ct*s-nuref*(1+bettapv*Tcref)*tovi*beta*G[m]*s+tovi*ac*G[m]*beta*s))/(mc*cc+dt*s*(Hcd_vc+Hcd_ct-nuref*bettapv*G[m]*tovi*beta)) 
            Tt[m]=(Tt[m-1]*mt*ct+dt*(Tcold[m]*Hcd_ct*s+Tpsold[m]*Hcd_tps*s1+Ttubeold[m]*Hcd_ttub*s2+tovi*s*at*G[m]*(1-beta)))/(mt*ct+dt*(Hcd_ct*s+Hcd_tps*s1+Hcd_ttub*s2)) 
            Ttube[m]=(Ttube[m-1]*mtub*ctub+dt*(Ttold[m]*Hcd_ttub*s2+Tpsold[m]*Hcd_tubps*s3+ Teauold[m]*Hcv_nftub[m]*sl))/(mtub*ctub+dt*(Hcd_tubps*s3+Hcd_ttub*s2+Hcv_nftub[m]*sl)) 
            A =(Hcv_nftub[m]*sl/L)/(ronf[m]*V*sec*cpnf[m]) 
            Teau[m]=Ttubeold[m]+Ttubeold[m]*(math.exp(-A*Ltube)-1)/(A*Ltube )-Tnfe[m]*(math.exp(-A*Ltube)-1)/(A*Ltube) 
            Tps[m]=(Tps[m-1]*cps*mps+dt*(Ttold[m]*Hcd_tps*s1+Ttubeold[m]*Hcd_tubps*s3+Tisoold[m]*Hcd_psiso*s4))/(mps*cps+dt*(Hcd_tps*s1+Hcd_tubps*s3+Hcd_psiso*s4)) 
            Tail[m]=(Tail[m-1]*cpail*mail+dt*(Tpsold[m]*Hcd_psail*sail+Tisoold[m]*Hcd_ailiso*sail))/(mail*cpail+dt*(Hcd_psail*sail+Hcd_ailiso*sail)) 
            Tiso[m]=(Tiso[m-1]*miso*cpiso+dt*(Tpsold[m]*Hcd_psiso*s4+Tailold[m]*Hcd_ailiso*sail))/(miso*cpiso+dt*(Hcd_psiso*s4+Hcd_ailiso*sail)) 
            Teau2[m]=Ttube[m]*(1-math.exp(-A*Ltube))+Tnfe[m]*math.exp(-A*Ltube) 
            puissance[m]=(ronf[m]*V*sec1*cpnf[m]*(Teau2[m]-Tnfe[m]))
            ntheau[m]=puissance[m]/(G[m]*s) 
            erreur = [None]*9
            erreur[1]=Tv[m]-Tvold[m] 
            erreur[2]=Tc[m]-Tcold[m] 
            erreur[3]=Tt[m]-Ttold[m] 
            erreur[4]=Ttube[m]-Ttubeold[m] 
            erreur[5]=Teau[m]-Teauold[m] 
            erreur[6]=Tps[m]-Tpsold[m] 
            erreur[7]=Tail[m]-Tailold[m] 
            erreur[8]=Tiso[m]-Tisoold[m] 
            erreur =list(map(abs,erreur[1:]))
            err=max(erreur)
            Tvold[m]=Tv[m] 
            Tcold[m]=Tc[m] 
            Ttold[m]=Tt[m] 
            Ttubeold[m]=Ttube[m] 
            Teauold[m]=Teau[m] 
            Tpsold[m]=Tps[m] 
            Tailold[m]=Tail[m] 
            Tisoold[m]=Tiso[m]
        Tvold[m+1]=Tv[m] 
        Tcold[m+1]=Tc[m] 
        Ttold[m+1]=Tt[m] 
        Ttubeold[m+1]=Ttube[m] 
        Teauold[m+1]=Teau[m] 
        Tpsold[m+1]=Tps[m] 
        Tailold[m+1]=Tail[m] 
        Tisoold[m+1]=Tiso[m]
    G_s = [i*s for i in G]
    red_t_t = np.array(ronf[:96])*V*sec*np.array(cpnf[:96])*(np.array(Teau2[:96])-np.array(Tnfe[:96]))
    rend_t = red_t_t/np.array(G_s[:-1])
    rendement_exg_th  = np.array(ronf[:96])*V*sec*np.array(cpnf[:96])*((np.array(Teau2[:96])-np.array(Tnfe[:96]))-np.array(Tamb[:96])*np.log(np.array(Teau2[:96])/np.array(Tnfe[:96])))
    rendement_exg_th = rendement_exg_th/(G[:96]*s*(1-np.array(Tamb)[:96]/5778))
    rend_e = [rendement_electrique(i) for i in Tc[:96]]
    rendement__exergetique_elec_th = np.array(rend_e)/(1-(np.array(Tamb[:-1])/5778))

    return Hcv_nftub

import numpy as np
import math
from numpy import linalg
from numpy.lib.function_base import average 
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error

plt.style.use('bmh')
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'Ubuntu'
plt.rcParams['font.monospace'] = 'Ubuntu Mono'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 10
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['figure.titlesize'] = 12

def rendement_electrique (Tc ,Tcrefrence=298, nrefrence = 0.17,betapv = 1/(270-25) ):
    return nrefrence*(1-betapv*(Tc-Tcrefrence))
    
def solve (Qnf): 
    Vvent = []
    with open(r'/content/vv1.txt') as f :
        x = f.readlines()
    for i in range (0,len(x)) :
        x[i] = x[i].replace(',' ,'.')
        Vvent.append(x[i][:len(x[i])-1])

  
    Vvent = [float(i) for i in Vvent ]
    Vvent = [float(i) for i in Vvent ]
    Vvent = Vvent[6:len(Vvent):5]
    Vvent[-1] = 0.08
    Vvent = np.array(Vvent)

    
    g = []
    with open(r'/content/g1.txt') as f :
       x = f.readlines()
    for i in range (0,len(x)) :
        x[i] = x[i].replace(',' ,'.')
        g.append(x[i][:len(x[i])-1])


    G= [float(i) for i in g ]
    G = np.array(G)

    Tamb = []
    with open(r'/content/tamb.txt') as f :
        x = f.readlines()
    for i in range (0,len(x)) :
        x[i] = x[i].replace(',' ,'.')
        Tamb.append(x[i][:len(x[i])-1])


    Tamb= [float(i)+273 for i in Tamb[:-1] ]
    Tamb = np.array(Tamb)
    l=0.79 
    Ltube=1.62 
    Dh=0.012 
    kc=0.036 
    kt=0.033 
    ec=0.0003 
    et=0.0005 
    kail=230 
    eail=0.01 
    cps=897 
    rops=2700    
    ktub=387 
    rotube=8960 
    ctub=385 
    etub=0.002 
    kiso=0.039 
    eiso=0.07 
    roiso=40 
    segma=5.67*10**-8 
    av=0.06 
    ac=0.85 
    at=0.5 
    tovi=0.84  
    bettapv=1/(270-25)
    Tcref=298 
    nuref=0.17
    Ntube=12 
    L=1.58 
    pi =math.pi
    sec = pi*(0.022**2)/4
    sec1 = pi*(0.012**2)/4
    V=Qnf/sec
    s=L*l 
    s2=pi*Dh*(1/4)*L*Ntube 
    s1=s-s2 
    sl=3.14*Dh*L*Ntube 
    s3=pi*Dh*((1/2)+(1/4))*L*Ntube 
    s4=(pi+2)*Dh*(1/2)*L*Ntube+s1 
    cc=677 
    rot=1200 
    ct=1250 
    cpail=897 
    roail=2700 
    cpiso=1500 
    beta=0.88 
    epsv=0.83   
    kv=1 
    roc=2702 
    rov=2530 
    cv=500 
    ev=0.0032 
    kps=230   
    eps=0.001 
    sail=eail*L*6 
    mtub=rotube*pi*(0.012/2)**2-(0.011/2)**2*Ltube*Ntube 
    mc=roc*s*ec 
    mv=rov*s*ev 
    mt=rot*s*et 
    miso=roiso*s*eiso 
    mail=roail*eail*L*0.04*6 
    mps=rops*(s4)*eps



    Tvold =[None]*167
    Tnfe= [25+273 ]*167
    Tcold=[None]*167
    Ttold=[None]*167
    Ttubeold=[None]*167
    Teauold=[None]*167
    Tpsold=[None]*167 
    Tailold=[None]*167
    Tisoold=[None]*167
    Tv=[None] *167
    Tc=[None]*167
    Tt=[None] *167
    Ttube=[None] *167
    Teau=[None] *167
    Tps=[None] *167
    Tail=[None] *167
    Tiso=[None] *167
    Tciel=[None] *167
    knf=[None] *167
    cpnf=[None] *167
    munf=[None] *167
    ronf=[None] *167
    #
    keau=[None] *167
    cpeau=[None] *167
    mueau=[None] *167
    roeau=[None] *167
    #
    pr= [None] *167
    Re= [None] *167
    Hr_vciel  =[None] *167
    Hcv_v =[None] *167
    Hcv_nftub=[None] *167
    Teau2 = [None] *167
    puissance = [None] *167
    ntheau = [None] *167
    Tvold[0]=320 
    Tcold[0]=310 
    Ttold[0]=300 
    Ttubeold[0]=299 
    Teauold[0]=298 
    Tpsold[0]=292 
    Tailold[0]=292 
    Tisoold[0]=292 
    mmax=96 
    t =np.arange(167)

    dt=300 
    tol=10**-5 
    err=tol+1 
    m=0 
    Tv[m]=Tvold[m] 
    Tc[m]=Tcold[m] 
    Tt[m]=Ttold[m] 
    Ttube[m]=Ttubeold[m] 
    Teau[m]=Teauold[m] 
    Tps[m]=Tpsold[m] 
    Tail[m]=Tailold[m] 
    Tiso[m]=Tisoold[m] 
    iter1=0 
    Hcd_vc=1/(ev/kv+ec/kc) #%conduction entre le verre et la cellule 
    Hcd_ct=1/(ec/kc+et/kt) #%conduction entre la cellule et tedlar
    Hcd_ttub=1/(et/kt+etub/ktub) #%conduction entre tedlar et le tube  
    Hcd_tps=1/(eps/kps+et/kt) #%conduction entre la plaque  et tedlar 
    Hcd_tubps=1/(etub/ktub+eps/kps) #%conduction entre le tube et la plaque 
    Hcd_psail=1/(eps/kps+eail/kail) #%conduction entre la plaque et l'aillette
    Hcd_psiso=1/(eps/kps+eiso/kiso) #%cnduction entre la plaque  et l'isolant
    Hcd_ailiso=1/(eail/kail+eiso/kiso)  #%conduction entre l'ailette et l'isolant
    tol=10**-5 
    err=tol+1
    while err>tol:
        iter1=iter1+1 
        m=0 
        Tciel[m]=0.0552*Tamb[m]**1.5 
        knf[m] =-7843*10**-9*Teauold[m]**2 +0.0062*Teauold[m] - 0.54 
        ronf[m] =-0.003 *Teauold[m]**2 + 1.505*Teauold[m] + 816.781 
        cpnf[m] =-0.0000643 *Teauold[m]**3 + 0.0552 *Teauold[m]**2 - 20.86 *Teauold[m] + 6719.637
        munf[m] = (2414*10**-8)*10**(247.8/(Teauold[m]-140)) 
        pr[m]=cpnf[m]*munf[m]/knf[m] 
        Re[m]=ronf[m]*V*Dh/munf[m] 
    #%les coefficient de transfert 
    #%1-RAYONNEMENT   
    #% Verre - voûte céleste
        Hr_vciel[m]=segma*epsv*(Tvold[m]**2+Tciel[m]**2)*(Tvold[m]+Tciel[m]) 
    #%2 CONDUCTION 
    #%1] le verre et l'aire    
        Hcv_v[m]=5.7+3.8*float(Vvent[m])  
    # %2 le nano-fluide et le tube
        if  Re[m]<=2300:
            Hcv_nftub[m]=4.36*knf[m]/Dh 
        elif  Re[m]>2300:
            Hcv_nftub[m]=0.023*(knf[m]/Dh)*(Re[m]**0.8)*(pr[m]**0.4) 
        Tv[m]=(Tciel[m]*Hr_vciel[m]+Tcold[m]*Hcd_vc+Tamb[m]*Hcv_v[m]+av*G[m])/(Hr_vciel[m]+Hcd_vc+Hcv_v[m]) 
        Tc[m]=(Tv[m]*Hcd_vc+Ttold[m]*Hcd_ct-nuref*(1+bettapv*Tcref)*G[m]*tovi*beta+tovi*ac*G[m]*beta)/(Hcd_vc+Hcd_ct-nuref*bettapv*G[m]*tovi*beta) 
        Tt[m]=(Tcold[m]*Hcd_ct*s+Tpsold[m]*Hcd_tps*s1+Ttubeold[m]*Hcd_ttub*s2+tovi*s*at*G[m]*(1-beta))/(Hcd_ct*s+Hcd_tps*s1+Hcd_ttub*s2) 
        Ttube[m] = (Ttold[m]*Hcd_ttub*s2 + Tpsold[m]*Hcd_tubps*s3 +Teauold[m]*Hcv_nftub[m]*sl)/(Hcd_tubps*s3 + Hcd_ttub*s2 + Hcv_nftub[m]*sl) 
        A =(Hcv_nftub[m]*sl/Ltube)/(ronf[m]*V*sec*cpnf[m]) 
        Teau[m]=Ttubeold[m]+Ttubeold[m]*(math.exp(-A*Ltube)-1)/(A*Ltube )-Tnfe[m]*(math.exp(-A*Ltube)-1)/(A*Ltube)
        Tps[m]=(Ttold[m]*Hcd_tps*s1+Ttubeold[m]*Hcd_tubps*s3+Tisoold[m]*Hcd_psiso*s4)/(Hcd_tps*s1+Hcd_tubps*s3+Hcd_psiso*s4) 
        Tail[m]=(Tpsold[m]*Hcd_psail*sail+Tisoold[m]*Hcd_ailiso*sail)/(Hcd_psail*sail+Hcd_ailiso*sail) 
        Tiso[m]=(Tpsold[m]*Hcd_psiso*s4+Tailold[m]*Hcd_ailiso*sail)/(Hcd_psiso*s4+Hcd_ailiso*sail) 
        Teau2[m]=Ttube[m]*(1-math.exp(-A*Ltube))+Tnfe[m]*math.exp(-A*Ltube) 
        puissance[m]=(ronf[m]*V*sec*cpnf[m])*(Teau2[m]-Tnfe[m]) 
        ntheau[m]=puissance[m]/(G[m]*s)
        erreur = [None]*8 
        erreur[1]=Tv[m]-Tvold[m] 
        erreur[2]=Tc[m]-Tcold[m] 
        erreur[3]=Tt[m]-Ttold[m] 
        erreur[4]=Ttube[m]-Ttubeold[m] 
        erreur[5]=Teau[m]-Teauold[m] 
        erreur[6]=Tps[m]-Tpsold[m] 
        erreur[7]=Tail[m]-Tailold[m] 
        erreur =list(map(abs,erreur[1:]))
        err=max(erreur)
        m=0
        Tvold[m]=Tv[m] 
        Tcold[m]=Tc[m] 
        Ttold[m]=Tt[m] 
        Ttubeold[m]=Ttube[m] 
        Teauold[m]=Teau[m] 
        Tpsold[m]=Tps[m] 
        Tailold[m]=Tail[m] 
        Tisoold[m]=Tiso[m] 

    m=0
    Tvold[m+1]=Tv[m] 
    Tcold[m+1]=Tc[m] 
    Ttold[m+1]=Tt[m] 
    Ttubeold[m+1]=Ttube[m] 
    Teauold[m+1]=Teau[m] 
    Tpsold[m+1]=Tps[m] 
    Tailold[m+1]=Tps[m]
    Tisoold[m+1]=Tiso[m] 
    tol=10**-5 
    for m in range(1,mmax):
        err=tol+1
        iter1 = 0
        while err>tol  :
            iter1 = iter1+1
            Tciel[m]=0.0552*Tamb[m]**1.5 
            knf[m] =-7843*10**-9*Teauold[m]**2 +0.0062*Teauold[m] - 0.54 
            ronf[m] =-0.003 *Teauold[m]**2 + 1.505*Teauold[m] + 816.781 
            cpnf[m] = -0.0000643 *Teauold[m]**3 + 0.0552 *Teauold[m]**2 - 20.86 *Teauold[m] + 6719.637 
            munf[m] = (2414*10**-8)*10**(247.8/(Teauold[m]-140)) 
            Re[m]=ronf[m]*V*Dh/munf[m] 
            pr[m]=cpnf[m]*munf[m]/knf[m] 
            Hr_vciel[m]=segma*epsv*(Tvold[m]**2+Tciel[m]**2)*(Tvold[m]+Tciel[m]) 
            Hcv_v[m]=5.7+3.8*float(Vvent[m])  
            if Re[m]<=2300 :
                Hcv_nftub[m]=4.36*knf[m]/Dh 
            elif  Re[m]>2300 :
                Hcv_nftub[m]=0.023*(knf[m]/Dh)*(Re[m]**0.8)*(pr[m]**(0.4)) 
            Tv[m]=(Tv[m-1]*mv*cv+dt*(Tciel[m]*Hr_vciel[m]*s+Tcold[m]*(Hcd_vc*s)+Tamb[m]*Hcv_v[m]*s+av*G[m]*s))/(mv*cv+dt*s*(Hr_vciel[m]+Hcd_vc+Hcv_v[m])) 
            Tc[m]=(Tc[m-1]*mc*cc+dt*(Tvold[m]*Hcd_vc*s+Ttold[m]*Hcd_ct*s-nuref*(1+bettapv*Tcref)*tovi*beta*G[m]*s+tovi*ac*G[m]*beta*s))/(mc*cc+dt*s*(Hcd_vc+Hcd_ct-nuref*bettapv*G[m]*tovi*beta)) 
            Tt[m]=(Tt[m-1]*mt*ct+dt*(Tcold[m]*Hcd_ct*s+Tpsold[m]*Hcd_tps*s1+Ttubeold[m]*Hcd_ttub*s2+tovi*s*at*G[m]*(1-beta)))/(mt*ct+dt*(Hcd_ct*s+Hcd_tps*s1+Hcd_ttub*s2)) 
            Ttube[m]=(Ttube[m-1]*mtub*ctub+dt*(Ttold[m]*Hcd_ttub*s2+Tpsold[m]*Hcd_tubps*s3+ Teauold[m]*Hcv_nftub[m]*sl))/(mtub*ctub+dt*(Hcd_tubps*s3+Hcd_ttub*s2+Hcv_nftub[m]*sl)) 
            A =(Hcv_nftub[m]*sl/L)/(ronf[m]*V*sec*cpnf[m]) 
            Teau[m]=Ttubeold[m]+Ttubeold[m]*(math.exp(-A*Ltube)-1)/(A*Ltube )-Tnfe[m]*(math.exp(-A*Ltube)-1)/(A*Ltube) 
            Tps[m]=(Tps[m-1]*cps*mps+dt*(Ttold[m]*Hcd_tps*s1+Ttubeold[m]*Hcd_tubps*s3+Tisoold[m]*Hcd_psiso*s4))/(mps*cps+dt*(Hcd_tps*s1+Hcd_tubps*s3+Hcd_psiso*s4)) 
            Tail[m]=(Tail[m-1]*cpail*mail+dt*(Tpsold[m]*Hcd_psail*sail+Tisoold[m]*Hcd_ailiso*sail))/(mail*cpail+dt*(Hcd_psail*sail+Hcd_ailiso*sail)) 
            Tiso[m]=(Tiso[m-1]*miso*cpiso+dt*(Tpsold[m]*Hcd_psiso*s4+Tailold[m]*Hcd_ailiso*sail))/(miso*cpiso+dt*(Hcd_psiso*s4+Hcd_ailiso*sail)) 
            Teau2[m]=Ttube[m]*(1-math.exp(-A*Ltube))+Tnfe[m]*math.exp(-A*Ltube) 
            puissance[m]=(ronf[m]*V*sec1*cpnf[m]*(Teau2[m]-Tnfe[m]))
            ntheau[m]=puissance[m]/(G[m]*s) 
            erreur = [None]*9
            erreur[1]=Tv[m]-Tvold[m] 
            erreur[2]=Tc[m]-Tcold[m] 
            erreur[3]=Tt[m]-Ttold[m] 
            erreur[4]=Ttube[m]-Ttubeold[m] 
            erreur[5]=Teau[m]-Teauold[m] 
            erreur[6]=Tps[m]-Tpsold[m] 
            erreur[7]=Tail[m]-Tailold[m] 
            erreur[8]=Tiso[m]-Tisoold[m] 
            erreur =list(map(abs,erreur[1:]))
            err=max(erreur)
            Tvold[m]=Tv[m] 
            Tcold[m]=Tc[m] 
            Ttold[m]=Tt[m] 
            Ttubeold[m]=Ttube[m] 
            Teauold[m]=Teau[m] 
            Tpsold[m]=Tps[m] 
            Tailold[m]=Tail[m] 
            Tisoold[m]=Tiso[m]
        Tvold[m+1]=Tv[m] 
        Tcold[m+1]=Tc[m] 
        Ttold[m+1]=Tt[m] 
        Ttubeold[m+1]=Ttube[m] 
        Teauold[m+1]=Teau[m] 
        Tpsold[m+1]=Tps[m] 
        Tailold[m+1]=Tail[m] 
        Tisoold[m+1]=Tiso[m]
    G_s = [i*s for i in G]
    red_t_t = np.array(ronf[:96])*V*sec*np.array(cpnf[:96])*(np.array(Teau2[:96])-np.array(Tnfe[:96]))
    rend_t = red_t_t/np.array(G_s[:-1])
    rendement_exg_th  = np.array(ronf[:96])*V*sec*np.array(cpnf[:96])*((np.array(Teau2[:96])-np.array(Tnfe[:96]))-np.array(Tamb[:96])*np.log(np.array(Teau2[:96])/np.array(Tnfe[:96])))
    rendement_exg_th = abs(rendement_exg_th/(G[:96]*s*(1-np.array(Tamb)[:96]/5778)))
    rend_e = [rendement_electrique(i) for i in Tc[:96]]
    rendement__exergetique_elec_th = np.array(rend_e)/(1-(np.array(Tamb[:-1])/5778))

    return  Hcv_nftub
def solve_al (Qnf):

    Vvent = []
    with open(r'/content/vv1.txt') as f :
        x = f.readlines()
    for i in range (0,len(x)) :
        x[i] = x[i].replace(',' ,'.')
        Vvent.append(x[i][:len(x[i])-1])

  
    Vvent = [float(i) for i in Vvent ]
    Vvent = [float(i) for i in Vvent ]
    Vvent = Vvent[6:len(Vvent):5]
    Vvent[-1] = 0.08
    Vvent = np.array(Vvent)

    
    g = []
    with open(r'/content/g1.txt') as f :
       x = f.readlines()
    for i in range (0,len(x)) :
        x[i] = x[i].replace(',' ,'.')
        g.append(x[i][:len(x[i])-1])


    G= [float(i) for i in g ]
    G = np.array(G)

    Tamb = []
    with open(r'/content/tamb.txt') as f :
        x = f.readlines()
    for i in range (0,len(x)) :
        x[i] = x[i].replace(',' ,'.')
        Tamb.append(x[i][:len(x[i])-1])


    Tamb= [float(i)+273 for i in Tamb[:-1] ]
    Tamb = np.array(Tamb)
    l=0.79 
    Ltube=1.62 
    Dh=0.012 
    kc=0.036 
    kt=0.033 
    ec=0.0003 
    et=0.0005 
    kail=230 
    eail=0.01 
    cps=897 
    rops=2700    
    ktub=387 
    rotube=8960 
    ctub=385 
    etub=0.002 
    kiso=0.039 
    eiso=0.07 
    roiso=40 
    segma=5.67*10**-8 
    av=0.06 
    ac=0.85 
    at=0.5 
    tovi=0.84  
    bettapv=1/(270-25)
    Tcref=298 
    nuref=0.17
    Ntube=12 
    L=1.58 
    pi =math.pi
    sec = pi*(0.022**2)/4
    sec1 = pi*(0.012**2)/4
    V=Qnf/sec
    s=L*l 
    s2=pi*Dh*(1/4)*L*Ntube 
    s1=s-s2 
    sl=3.14*Dh*L*Ntube 
    s3=pi*Dh*((1/2)+(1/4))*L*Ntube 
    s4=(pi+2)*Dh*(1/2)*L*Ntube+s1 
    cc=677 
    rot=1200 
    ct=1250 
    cpail=897 
    roail=2700 
    cpiso=1500 
    beta=0.88 
    epsv=0.83   
    kv=1 
    roc=2702 
    rov=2530 
    cv=500 
    ev=0.0032 
    kps=230   
    eps=0.001 
    sail=eail*L*6 
    mtub=rotube*pi*(0.012/2)**2-(0.011/2)**2*Ltube*Ntube 
    mc=roc*s*ec 
    mv=rov*s*ev 
    mt=rot*s*et 
    miso=roiso*s*eiso 
    mail=roail*eail*L*0.04*6 
    mps=rops*(s4)*eps



    Tvold =[None]*167
    Tnfe= [25+273 ]*167
    Tcold=[None]*167
    Ttold=[None]*167
    Ttubeold=[None]*167
    Teauold=[None]*167
    Tpsold=[None]*167 
    Tailold=[None]*167
    Tisoold=[None]*167
    Tv=[None] *167
    Tc=[None]*167
    Tt=[None] *167
    Ttube=[None] *167
    Teau=[None] *167
    Tps=[None] *167
    Tail=[None] *167
    Tiso=[None] *167
    Tciel=[None] *167
    knf=[None] *167
    cpnf=[None] *167
    munf=[None] *167
    ronf=[None] *167
    #
    keau=[None] *167
    cpeau=[None] *167
    mueau=[None] *167
    roeau=[None] *167
    #
    pr= [None] *167
    Re= [None] *167
    Hr_vciel  =[None] *167
    Hcv_v =[None] *167
    Hcv_nftub=[None] *167
    Teau2 = [None] *167
    puissance = [None] *167
    ntheau = [None] *167
    Tvold[0]=320 
    Tcold[0]=310 
    Ttold[0]=300 
    Ttubeold[0]=299 
    Teauold[0]=298 
    Tpsold[0]=292 
    Tailold[0]=292 
    Tisoold[0]=292 
    mmax=96 
    t =np.arange(167)

    dt=300 
    tol=10**-5 
    err=tol+1 
    m=0 
    Tv[m]=Tvold[m] 
    Tc[m]=Tcold[m] 
    Tt[m]=Ttold[m] 
    Ttube[m]=Ttubeold[m] 
    Teau[m]=Teauold[m] 
    Tps[m]=Tpsold[m] 
    Tail[m]=Tailold[m] 
    Tiso[m]=Tisoold[m] 
    iter1=0 
    Hcd_vc=1/(ev/kv+ec/kc) #%conduction entre le verre et la cellule 
    Hcd_ct=1/(ec/kc+et/kt) #%conduction entre la cellule et tedlar
    Hcd_ttub=1/(et/kt+etub/ktub) #%conduction entre tedlar et le tube  
    Hcd_tps=1/(eps/kps+et/kt) #%conduction entre la plaque  et tedlar 
    Hcd_tubps=1/(etub/ktub+eps/kps) #%conduction entre le tube et la plaque 
    Hcd_psail=1/(eps/kps+eail/kail) #%conduction entre la plaque et l'aillette
    Hcd_psiso=1/(eps/kps+eiso/kiso) #%cnduction entre la plaque  et l'isolant
    Hcd_ailiso=1/(eail/kail+eiso/kiso)  #%conduction entre l'ailette et l'isolant
    tol=10**-5 
    err=tol+1
    while err>tol:
        iter1=iter1+1 
        m=0 
        Tciel[m]=0.0552*Tamb[m]**1.5 
        knf[m] =-7843*10**-9*Teauold[m]**2 +0.0062*Teauold[m] - 0.54 
        ronf[m] =-0.003 *Teauold[m]**2 + 1.505*Teauold[m] + 816.781 
        cpnf[m] =-0.0000643 *Teauold[m]**3 + 0.0552 *Teauold[m]**2 - 20.86 *Teauold[m] + 6719.637
        munf[m] = (2414*10**-8)*10**(247.8/(Teauold[m]-140)) 
        pr[m]=cpnf[m]*munf[m]/knf[m] 
        Re[m]=ronf[m]*V*Dh/munf[m] 
    #%les coefficient de transfert 
    #%1-RAYONNEMENT   
    #% Verre - voûte céleste
        Hr_vciel[m]=segma*epsv*(Tvold[m]**2+Tciel[m]**2)*(Tvold[m]+Tciel[m]) 
    #%2 CONDUCTION 
    #%1] le verre et l'aire    
        Hcv_v[m]=5.7+3.8*float(Vvent[m])  
    # %2 le nano-fluide et le tube
        if  Re[m]<=2300:
            Hcv_nftub[m]=4.36*knf[m]/Dh 
        elif  Re[m]>2300:
            Hcv_nftub[m]=0.023*(knf[m]/Dh)*(Re[m]**0.8)*(pr[m]**0.4) 
        Tv[m]=(Tciel[m]*Hr_vciel[m]+Tcold[m]*Hcd_vc+Tamb[m]*Hcv_v[m]+av*G[m])/(Hr_vciel[m]+Hcd_vc+Hcv_v[m]) 
        Tc[m]=(Tv[m]*Hcd_vc+Ttold[m]*Hcd_ct-nuref*(1+bettapv*Tcref)*G[m]*tovi*beta+tovi*ac*G[m]*beta)/(Hcd_vc+Hcd_ct-nuref*bettapv*G[m]*tovi*beta) 
        Tt[m]=(Tcold[m]*Hcd_ct*s+Tpsold[m]*Hcd_tps*s1+Ttubeold[m]*Hcd_ttub*s2+tovi*s*at*G[m]*(1-beta))/(Hcd_ct*s+Hcd_tps*s1+Hcd_ttub*s2) 
        Ttube[m] = (Ttold[m]*Hcd_ttub*s2 + Tpsold[m]*Hcd_tubps*s3 +Teauold[m]*Hcv_nftub[m]*sl)/(Hcd_tubps*s3 + Hcd_ttub*s2 + Hcv_nftub[m]*sl) 
        A =(Hcv_nftub[m]*sl/Ltube)/(ronf[m]*V*sec*cpnf[m]) 
        Teau[m]=Ttubeold[m]+Ttubeold[m]*(math.exp(-A*Ltube)-1)/(A*Ltube )-Tnfe[m]*(math.exp(-A*Ltube)-1)/(A*Ltube)
        Tps[m]=(Ttold[m]*Hcd_tps*s1+Ttubeold[m]*Hcd_tubps*s3+Tisoold[m]*Hcd_psiso*s4)/(Hcd_tps*s1+Hcd_tubps*s3+Hcd_psiso*s4) 
        Tail[m]=(Tpsold[m]*Hcd_psail*sail+Tisoold[m]*Hcd_ailiso*sail)/(Hcd_psail*sail+Hcd_ailiso*sail) 
        Tiso[m]=(Tpsold[m]*Hcd_psiso*s4+Tailold[m]*Hcd_ailiso*sail)/(Hcd_psiso*s4+Hcd_ailiso*sail) 
        Teau2[m]=Ttube[m]*(1-math.exp(-A*Ltube))+Tnfe[m]*math.exp(-A*Ltube) 
        puissance[m]=(ronf[m]*V*sec*cpnf[m])*(Teau2[m]-Tnfe[m]) 
        ntheau[m]=puissance[m]/(G[m]*s)
        erreur = [None]*8 
        erreur[1]=Tv[m]-Tvold[m] 
        erreur[2]=Tc[m]-Tcold[m] 
        erreur[3]=Tt[m]-Ttold[m] 
        erreur[4]=Ttube[m]-Ttubeold[m] 
        erreur[5]=Teau[m]-Teauold[m] 
        erreur[6]=Tps[m]-Tpsold[m] 
        erreur[7]=Tail[m]-Tailold[m] 
        erreur =list(map(abs,erreur[1:]))
        err=max(erreur)
        m=0
        Tvold[m]=Tv[m] 
        Tcold[m]=Tc[m] 
        Ttold[m]=Tt[m] 
        Ttubeold[m]=Ttube[m] 
        Teauold[m]=Teau[m] 
        Tpsold[m]=Tps[m] 
        Tailold[m]=Tail[m] 
        Tisoold[m]=Tiso[m] 

    m=0
    Tvold[m+1]=Tv[m] 
    Tcold[m+1]=Tc[m] 
    Ttold[m+1]=Tt[m] 
    Ttubeold[m+1]=Ttube[m] 
    Teauold[m+1]=Teau[m] 
    Tpsold[m+1]=Tps[m] 
    Tailold[m+1]=Tps[m]
    Tisoold[m+1]=Tiso[m] 
    tol=10**-5 
    knp = 40
    ronp = 3970
    cpnp = 765
    dnp=10e-8
    deau=0.434*10e-9
    phi= 0.1
    for m in range(1,mmax):
        err=tol+1
        iter1 = 0
        while err>tol  :
            iter1 = iter1+1
            Tciel[m]=0.0552*Tamb[m]**1.5 
            #proprites for the pure water 
            keau[m] =-7843*10**-9*Teauold[m]**2 +0.0062*Teauold[m] - 0.54 
            roeau[m] =-0.003 *Teauold[m]**2 + 1.505*Teauold[m] + 816.781 
            cpeau[m] = -0.0000643 *Teauold[m]**3 + 0.0552 *Teauold[m]**2 - 20.86 *Teauold[m] + 6719.637 
            mueau[m] = (2414*10**-8)*10**(247.8/(Teauold[m]-140))
            # correlation for the nanofluids :
            knf[m] =(knp+2*keau[m]-2*phi*(keau[m]-knp))/(knp+2*keau[m]+phi*(keau[m]-knp))*keau[m]
            ronf[m] =(1-phi)*roeau[m] +phi*ronp
            cpnf[m] = (1-phi)*cpeau[m] +phi*cpnp
            munf[m] = mueau[m]*(1/(1-34.87*(dnp/deau)**-0.3*phi**1.03))
            ######
            Re[m]=ronf[m]*V*Dh/munf[m] 
            pr[m]=cpnf[m]*munf[m]/knf[m] 
            Hr_vciel[m]=segma*epsv*(Tvold[m]**2+Tciel[m]**2)*(Tvold[m]+Tciel[m]) 
            Hcv_v[m]=5.7+3.8*float(Vvent[m])  
            if Re[m]<=2300 :
                Hcv_nftub[m]=4.36*knf[m]/Dh 
            elif  Re[m]>2300 :
                Hcv_nftub[m]=0.023*(knf[m]/Dh)*(Re[m]**0.8)*(pr[m]**(0.4)) 
            Tv[m]=(Tv[m-1]*mv*cv+dt*(Tciel[m]*Hr_vciel[m]*s+Tcold[m]*(Hcd_vc*s)+Tamb[m]*Hcv_v[m]*s+av*G[m]*s))/(mv*cv+dt*s*(Hr_vciel[m]+Hcd_vc+Hcv_v[m])) 
            Tc[m]=(Tc[m-1]*mc*cc+dt*(Tvold[m]*Hcd_vc*s+Ttold[m]*Hcd_ct*s-nuref*(1+bettapv*Tcref)*tovi*beta*G[m]*s+tovi*ac*G[m]*beta*s))/(mc*cc+dt*s*(Hcd_vc+Hcd_ct-nuref*bettapv*G[m]*tovi*beta)) 
            Tt[m]=(Tt[m-1]*mt*ct+dt*(Tcold[m]*Hcd_ct*s+Tpsold[m]*Hcd_tps*s1+Ttubeold[m]*Hcd_ttub*s2+tovi*s*at*G[m]*(1-beta)))/(mt*ct+dt*(Hcd_ct*s+Hcd_tps*s1+Hcd_ttub*s2)) 
            Ttube[m]=(Ttube[m-1]*mtub*ctub+dt*(Ttold[m]*Hcd_ttub*s2+Tpsold[m]*Hcd_tubps*s3+ Teauold[m]*Hcv_nftub[m]*sl))/(mtub*ctub+dt*(Hcd_tubps*s3+Hcd_ttub*s2+Hcv_nftub[m]*sl)) 
            A =(Hcv_nftub[m]*sl/L)/(ronf[m]*V*sec*cpnf[m]) 
            Teau[m]=Ttubeold[m]+Ttubeold[m]*(math.exp(-A*Ltube)-1)/(A*Ltube )-Tnfe[m]*(math.exp(-A*Ltube)-1)/(A*Ltube) 
            Tps[m]=(Tps[m-1]*cps*mps+dt*(Ttold[m]*Hcd_tps*s1+Ttubeold[m]*Hcd_tubps*s3+Tisoold[m]*Hcd_psiso*s4))/(mps*cps+dt*(Hcd_tps*s1+Hcd_tubps*s3+Hcd_psiso*s4)) 
            Tail[m]=(Tail[m-1]*cpail*mail+dt*(Tpsold[m]*Hcd_psail*sail+Tisoold[m]*Hcd_ailiso*sail))/(mail*cpail+dt*(Hcd_psail*sail+Hcd_ailiso*sail)) 
            Tiso[m]=(Tiso[m-1]*miso*cpiso+dt*(Tpsold[m]*Hcd_psiso*s4+Tailold[m]*Hcd_ailiso*sail))/(miso*cpiso+dt*(Hcd_psiso*s4+Hcd_ailiso*sail)) 
            Teau2[m]=Ttube[m]*(1-math.exp(-A*Ltube))+Tnfe[m]*math.exp(-A*Ltube) 
            puissance[m]=(ronf[m]*V*sec1*cpnf[m]*(Teau2[m]-Tnfe[m]))
            ntheau[m]=puissance[m]/(G[m]*s) 
            erreur = [None]*9
            erreur[1]=Tv[m]-Tvold[m] 
            erreur[2]=Tc[m]-Tcold[m] 
            erreur[3]=Tt[m]-Ttold[m] 
            erreur[4]=Ttube[m]-Ttubeold[m] 
            erreur[5]=Teau[m]-Teauold[m] 
            erreur[6]=Tps[m]-Tpsold[m] 
            erreur[7]=Tail[m]-Tailold[m] 
            erreur[8]=Tiso[m]-Tisoold[m] 
            erreur =list(map(abs,erreur[1:]))
            err=max(erreur)
            Tvold[m]=Tv[m] 
            Tcold[m]=Tc[m] 
            Ttold[m]=Tt[m] 
            Ttubeold[m]=Ttube[m] 
            Teauold[m]=Teau[m] 
            Tpsold[m]=Tps[m] 
            Tailold[m]=Tail[m] 
            Tisoold[m]=Tiso[m]
        Tvold[m+1]=Tv[m] 
        Tcold[m+1]=Tc[m] 
        Ttold[m+1]=Tt[m] 
        Ttubeold[m+1]=Ttube[m] 
        Teauold[m+1]=Teau[m] 
        Tpsold[m+1]=Tps[m] 
        Tailold[m+1]=Tail[m] 
        Tisoold[m+1]=Tiso[m]
    G_s = [i*s for i in G]
    red_t_t = np.array(ronf[:96])*V*sec*np.array(cpnf[:96])*(np.array(Teau2[:96])-np.array(Tnfe[:96]))
    rend_t = red_t_t/np.array(G_s[:-1])
    rendement_exg_th  = np.array(ronf[:96])*V*sec*np.array(cpnf[:96])*((np.array(Teau2[:96])-np.array(Tnfe[:96]))-np.array(Tamb[:96])*np.log(np.array(Teau2[:96])/np.array(Tnfe[:96])))
    rendement_exg_th = abs(rendement_exg_th/(G[:96]*s*(1-np.array(Tamb)[:96]/5778)))
    rend_e = [rendement_electrique(i) for i in Tc[:96]]
    rendement__exergetique_elec_th = np.array(rend_e)/(1-(np.array(Tamb[:-1])/5778))

    return Hcv_nftub

def solve_zn (Qnf):

    Vvent = []
    with open(r'/content/vv1.txt') as f :
        x = f.readlines()
    for i in range (0,len(x)) :
        x[i] = x[i].replace(',' ,'.')
        Vvent.append(x[i][:len(x[i])-1])

  
    Vvent = [float(i) for i in Vvent ]
    Vvent = [float(i) for i in Vvent ]
    Vvent = Vvent[6:len(Vvent):5]
    Vvent[-1] = 0.08
    Vvent = np.array(Vvent)

    
    g = []
    with open(r'/content/g1.txt') as f :
       x = f.readlines()
    for i in range (0,len(x)) :
        x[i] = x[i].replace(',' ,'.')
        g.append(x[i][:len(x[i])-1])


    G= [float(i) for i in g ]
    G = np.array(G)

    Tamb = []
    with open(r'/content/tamb.txt') as f :
        x = f.readlines()
    for i in range (0,len(x)) :
        x[i] = x[i].replace(',' ,'.')
        Tamb.append(x[i][:len(x[i])-1])


    Tamb= [float(i)+273 for i in Tamb[:-1] ]
    Tamb = np.array(Tamb)
    l=0.79 
    Ltube=1.62 
    Dh=0.012 
    kc=0.036 
    kt=0.033 
    ec=0.0003 
    et=0.0005 
    kail=230 
    eail=0.01 
    cps=897 
    rops=2700    
    ktub=387 
    rotube=8960 
    ctub=385 
    etub=0.002 
    kiso=0.039 
    eiso=0.07 
    roiso=40 
    segma=5.67*10**-8 
    av=0.06 
    ac=0.85 
    at=0.5 
    tovi=0.84  
    bettapv=1/(270-25)
    Tcref=298 
    nuref=0.17
    Ntube=12 
    L=1.58 
    pi =math.pi
    sec = pi*(0.022**2)/4
    sec1 = pi*(0.012**2)/4
    V=Qnf/sec
    s=L*l 
    s2=pi*Dh*(1/4)*L*Ntube 
    s1=s-s2 
    sl=3.14*Dh*L*Ntube 
    s3=pi*Dh*((1/2)+(1/4))*L*Ntube 
    s4=(pi+2)*Dh*(1/2)*L*Ntube+s1 
    cc=677 
    rot=1200 
    ct=1250 
    cpail=897 
    roail=2700 
    cpiso=1500 
    beta=0.88 
    epsv=0.83   
    kv=1 
    roc=2702 
    rov=2530 
    cv=500 
    ev=0.0032 
    kps=230   
    eps=0.001 
    sail=eail*L*6 
    mtub=rotube*pi*(0.012/2)**2-(0.011/2)**2*Ltube*Ntube 
    mc=roc*s*ec 
    mv=rov*s*ev 
    mt=rot*s*et 
    miso=roiso*s*eiso 
    mail=roail*eail*L*0.04*6 
    mps=rops*(s4)*eps



    Tvold =[None]*167
    Tnfe= [25+273 ]*167
    Tcold=[None]*167
    Ttold=[None]*167
    Ttubeold=[None]*167
    Teauold=[None]*167
    Tpsold=[None]*167 
    Tailold=[None]*167
    Tisoold=[None]*167
    Tv=[None] *167
    Tc=[None]*167
    Tt=[None] *167
    Ttube=[None] *167
    Teau=[None] *167
    Tps=[None] *167
    Tail=[None] *167
    Tiso=[None] *167
    Tciel=[None] *167
    knf=[None] *167
    cpnf=[None] *167
    munf=[None] *167
    ronf=[None] *167
    #
    keau=[None] *167
    cpeau=[None] *167
    mueau=[None] *167
    roeau=[None] *167
    #
    pr= [None] *167
    Re= [None] *167
    Hr_vciel  =[None] *167
    Hcv_v =[None] *167
    Hcv_nftub=[None] *167
    Teau2 = [None] *167
    puissance = [None] *167
    ntheau = [None] *167
    Tvold[0]=320 
    Tcold[0]=310 
    Ttold[0]=300 
    Ttubeold[0]=299 
    Teauold[0]=298 
    Tpsold[0]=292 
    Tailold[0]=292 
    Tisoold[0]=292 
    mmax=96 
    t =np.arange(167)

    dt=300 
    tol=10**-5 
    err=tol+1 
    m=0 
    Tv[m]=Tvold[m] 
    Tc[m]=Tcold[m] 
    Tt[m]=Ttold[m] 
    Ttube[m]=Ttubeold[m] 
    Teau[m]=Teauold[m] 
    Tps[m]=Tpsold[m] 
    Tail[m]=Tailold[m] 
    Tiso[m]=Tisoold[m] 
    iter1=0 
    Hcd_vc=1/(ev/kv+ec/kc) #%conduction entre le verre et la cellule 
    Hcd_ct=1/(ec/kc+et/kt) #%conduction entre la cellule et tedlar
    Hcd_ttub=1/(et/kt+etub/ktub) #%conduction entre tedlar et le tube  
    Hcd_tps=1/(eps/kps+et/kt) #%conduction entre la plaque  et tedlar 
    Hcd_tubps=1/(etub/ktub+eps/kps) #%conduction entre le tube et la plaque 
    Hcd_psail=1/(eps/kps+eail/kail) #%conduction entre la plaque et l'aillette
    Hcd_psiso=1/(eps/kps+eiso/kiso) #%cnduction entre la plaque  et l'isolant
    Hcd_ailiso=1/(eail/kail+eiso/kiso)  #%conduction entre l'ailette et l'isolant
    tol=10**-5 
    err=tol+1
    while err>tol:
        iter1=iter1+1 
        m=0 
        Tciel[m]=0.0552*Tamb[m]**1.5 
        knf[m] =-7843*10**-9*Teauold[m]**2 +0.0062*Teauold[m] - 0.54 
        ronf[m] =-0.003 *Teauold[m]**2 + 1.505*Teauold[m] + 816.781 
        cpnf[m] =-0.0000643 *Teauold[m]**3 + 0.0552 *Teauold[m]**2 - 20.86 *Teauold[m] + 6719.637
        munf[m] = (2414*10**-8)*10**(247.8/(Teauold[m]-140)) 
        pr[m]=cpnf[m]*munf[m]/knf[m] 
        Re[m]=ronf[m]*V*Dh/munf[m] 
    #%les coefficient de transfert 
    #%1-RAYONNEMENT   
    #% Verre - voûte céleste
        Hr_vciel[m]=segma*epsv*(Tvold[m]**2+Tciel[m]**2)*(Tvold[m]+Tciel[m]) 
    #%2 CONDUCTION 
    #%1] le verre et l'aire    
        Hcv_v[m]=5.7+3.8*float(Vvent[m])  
    # %2 le nano-fluide et le tube
        if  Re[m]<=2300:
            Hcv_nftub[m]=4.36*knf[m]/Dh 
        elif  Re[m]>2300:
            Hcv_nftub[m]=0.023*(knf[m]/Dh)*(Re[m]**0.8)*(pr[m]**0.4) 
        Tv[m]=(Tciel[m]*Hr_vciel[m]+Tcold[m]*Hcd_vc+Tamb[m]*Hcv_v[m]+av*G[m])/(Hr_vciel[m]+Hcd_vc+Hcv_v[m]) 
        Tc[m]=(Tv[m]*Hcd_vc+Ttold[m]*Hcd_ct-nuref*(1+bettapv*Tcref)*G[m]*tovi*beta+tovi*ac*G[m]*beta)/(Hcd_vc+Hcd_ct-nuref*bettapv*G[m]*tovi*beta) 
        Tt[m]=(Tcold[m]*Hcd_ct*s+Tpsold[m]*Hcd_tps*s1+Ttubeold[m]*Hcd_ttub*s2+tovi*s*at*G[m]*(1-beta))/(Hcd_ct*s+Hcd_tps*s1+Hcd_ttub*s2) 
        Ttube[m] = (Ttold[m]*Hcd_ttub*s2 + Tpsold[m]*Hcd_tubps*s3 +Teauold[m]*Hcv_nftub[m]*sl)/(Hcd_tubps*s3 + Hcd_ttub*s2 + Hcv_nftub[m]*sl) 
        A =(Hcv_nftub[m]*sl/Ltube)/(ronf[m]*V*sec*cpnf[m]) 
        Teau[m]=Ttubeold[m]+Ttubeold[m]*(math.exp(-A*Ltube)-1)/(A*Ltube )-Tnfe[m]*(math.exp(-A*Ltube)-1)/(A*Ltube)
        Tps[m]=(Ttold[m]*Hcd_tps*s1+Ttubeold[m]*Hcd_tubps*s3+Tisoold[m]*Hcd_psiso*s4)/(Hcd_tps*s1+Hcd_tubps*s3+Hcd_psiso*s4) 
        Tail[m]=(Tpsold[m]*Hcd_psail*sail+Tisoold[m]*Hcd_ailiso*sail)/(Hcd_psail*sail+Hcd_ailiso*sail) 
        Tiso[m]=(Tpsold[m]*Hcd_psiso*s4+Tailold[m]*Hcd_ailiso*sail)/(Hcd_psiso*s4+Hcd_ailiso*sail) 
        Teau2[m]=Ttube[m]*(1-math.exp(-A*Ltube))+Tnfe[m]*math.exp(-A*Ltube) 
        puissance[m]=(ronf[m]*V*sec*cpnf[m])*(Teau2[m]-Tnfe[m]) 
        ntheau[m]=puissance[m]/(G[m]*s)
        erreur = [None]*8 
        erreur[1]=Tv[m]-Tvold[m] 
        erreur[2]=Tc[m]-Tcold[m] 
        erreur[3]=Tt[m]-Ttold[m] 
        erreur[4]=Ttube[m]-Ttubeold[m] 
        erreur[5]=Teau[m]-Teauold[m] 
        erreur[6]=Tps[m]-Tpsold[m] 
        erreur[7]=Tail[m]-Tailold[m] 
        erreur =list(map(abs,erreur[1:]))
        err=max(erreur)
        m=0
        Tvold[m]=Tv[m] 
        Tcold[m]=Tc[m] 
        Ttold[m]=Tt[m] 
        Ttubeold[m]=Ttube[m] 
        Teauold[m]=Teau[m] 
        Tpsold[m]=Tps[m] 
        Tailold[m]=Tail[m] 
        Tisoold[m]=Tiso[m] 

    m=0
    Tvold[m+1]=Tv[m] 
    Tcold[m+1]=Tc[m] 
    Ttold[m+1]=Tt[m] 
    Ttubeold[m+1]=Ttube[m] 
    Teauold[m+1]=Teau[m] 
    Tpsold[m+1]=Tps[m] 
    Tailold[m+1]=Tps[m]
    Tisoold[m+1]=Tiso[m] 
    tol=10**-5 
    knp = 13
    ronp = 5600
    cpnp = 495
    dnp=50e-9
    deau=0.434*10e-9
    phi= 0.1
    for m in range(1,mmax):
        err=tol+1
        iter1 = 0
        while err>tol  :
            iter1 = iter1+1
            Tciel[m]=0.0552*Tamb[m]**1.5 
            #proprites for the pure water 
            keau[m] =-7843*10**-9*Teauold[m]**2 +0.0062*Teauold[m] - 0.54 
            roeau[m] =-0.003 *Teauold[m]**2 + 1.505*Teauold[m] + 816.781 
            cpeau[m] = -0.0000643 *Teauold[m]**3 + 0.0552 *Teauold[m]**2 - 20.86 *Teauold[m] + 6719.637 
            mueau[m] = (2414*10**-8)*10**(247.8/(Teauold[m]-140))
            # correlation for the nanofluids :
            knf[m] =(knp+2*keau[m]-2*phi*(keau[m]-knp))/(knp+2*keau[m]+phi*(keau[m]-knp))*keau[m]
            ronf[m] =(1-phi)*roeau[m] +phi*ronp
            cpnf[m] = (1-phi)*cpeau[m] +phi*cpnp
            munf[m] = mueau[m]*(1/(1-34.87*(dnp/deau)**-0.3*phi**1.03))
            ######
            Re[m]=ronf[m]*V*Dh/munf[m] 
            pr[m]=cpnf[m]*munf[m]/knf[m] 
            Hr_vciel[m]=segma*epsv*(Tvold[m]**2+Tciel[m]**2)*(Tvold[m]+Tciel[m]) 
            Hcv_v[m]=5.7+3.8*float(Vvent[m])  
            if Re[m]<=2300 :
                Hcv_nftub[m]=4.36*knf[m]/Dh 
            elif  Re[m]>2300 :
                Hcv_nftub[m]=0.023*(knf[m]/Dh)*(Re[m]**0.8)*(pr[m]**(0.4)) 
            Tv[m]=(Tv[m-1]*mv*cv+dt*(Tciel[m]*Hr_vciel[m]*s+Tcold[m]*(Hcd_vc*s)+Tamb[m]*Hcv_v[m]*s+av*G[m]*s))/(mv*cv+dt*s*(Hr_vciel[m]+Hcd_vc+Hcv_v[m])) 
            Tc[m]=(Tc[m-1]*mc*cc+dt*(Tvold[m]*Hcd_vc*s+Ttold[m]*Hcd_ct*s-nuref*(1+bettapv*Tcref)*tovi*beta*G[m]*s+tovi*ac*G[m]*beta*s))/(mc*cc+dt*s*(Hcd_vc+Hcd_ct-nuref*bettapv*G[m]*tovi*beta)) 
            Tt[m]=(Tt[m-1]*mt*ct+dt*(Tcold[m]*Hcd_ct*s+Tpsold[m]*Hcd_tps*s1+Ttubeold[m]*Hcd_ttub*s2+tovi*s*at*G[m]*(1-beta)))/(mt*ct+dt*(Hcd_ct*s+Hcd_tps*s1+Hcd_ttub*s2)) 
            Ttube[m]=(Ttube[m-1]*mtub*ctub+dt*(Ttold[m]*Hcd_ttub*s2+Tpsold[m]*Hcd_tubps*s3+ Teauold[m]*Hcv_nftub[m]*sl))/(mtub*ctub+dt*(Hcd_tubps*s3+Hcd_ttub*s2+Hcv_nftub[m]*sl)) 
            A =(Hcv_nftub[m]*sl/L)/(ronf[m]*V*sec*cpnf[m]) 
            Teau[m]=Ttubeold[m]+Ttubeold[m]*(math.exp(-A*Ltube)-1)/(A*Ltube )-Tnfe[m]*(math.exp(-A*Ltube)-1)/(A*Ltube) 
            Tps[m]=(Tps[m-1]*cps*mps+dt*(Ttold[m]*Hcd_tps*s1+Ttubeold[m]*Hcd_tubps*s3+Tisoold[m]*Hcd_psiso*s4))/(mps*cps+dt*(Hcd_tps*s1+Hcd_tubps*s3+Hcd_psiso*s4)) 
            Tail[m]=(Tail[m-1]*cpail*mail+dt*(Tpsold[m]*Hcd_psail*sail+Tisoold[m]*Hcd_ailiso*sail))/(mail*cpail+dt*(Hcd_psail*sail+Hcd_ailiso*sail)) 
            Tiso[m]=(Tiso[m-1]*miso*cpiso+dt*(Tpsold[m]*Hcd_psiso*s4+Tailold[m]*Hcd_ailiso*sail))/(miso*cpiso+dt*(Hcd_psiso*s4+Hcd_ailiso*sail)) 
            Teau2[m]=Ttube[m]*(1-math.exp(-A*Ltube))+Tnfe[m]*math.exp(-A*Ltube) 
            puissance[m]=(ronf[m]*V*sec1*cpnf[m]*(Teau2[m]-Tnfe[m]))
            ntheau[m]=puissance[m]/(G[m]*s) 
            erreur = [None]*9
            erreur[1]=Tv[m]-Tvold[m] 
            erreur[2]=Tc[m]-Tcold[m] 
            erreur[3]=Tt[m]-Ttold[m] 
            erreur[4]=Ttube[m]-Ttubeold[m] 
            erreur[5]=Teau[m]-Teauold[m] 
            erreur[6]=Tps[m]-Tpsold[m] 
            erreur[7]=Tail[m]-Tailold[m] 
            erreur[8]=Tiso[m]-Tisoold[m] 
            erreur =list(map(abs,erreur[1:]))
            err=max(erreur)
            Tvold[m]=Tv[m] 
            Tcold[m]=Tc[m] 
            Ttold[m]=Tt[m] 
            Ttubeold[m]=Ttube[m] 
            Teauold[m]=Teau[m] 
            Tpsold[m]=Tps[m] 
            Tailold[m]=Tail[m] 
            Tisoold[m]=Tiso[m]
        Tvold[m+1]=Tv[m] 
        Tcold[m+1]=Tc[m] 
        Ttold[m+1]=Tt[m] 
        Ttubeold[m+1]=Ttube[m] 
        Teauold[m+1]=Teau[m] 
        Tpsold[m+1]=Tps[m] 
        Tailold[m+1]=Tail[m] 
        Tisoold[m+1]=Tiso[m]
    G_s = [i*s for i in G]
    red_t_t = np.array(ronf[:96])*V*sec*np.array(cpnf[:96])*(np.array(Teau2[:96])-np.array(Tnfe[:96]))
    rend_t = red_t_t/np.array(G_s[:-1])
    rendement_exg_th  = np.array(ronf[:96])*V*sec*np.array(cpnf[:96])*((np.array(Teau2[:96])-np.array(Tnfe[:96]))-np.array(Tamb[:96])*np.log(np.array(Teau2[:96])/np.array(Tnfe[:96])))
    rendement_exg_th = abs(rendement_exg_th/(G[:96]*s*(1-np.array(Tamb)[:96]/5778)))
    rend_e = [rendement_electrique(i) for i in Tc[:96]]
    rendement__exergetique_elec_th = np.array(rend_e)/(1-(np.array(Tamb[:-1])/5778))

    return Hcv_nftub

def solve_cu (Qnf):

    Vvent = []
    with open(r'/content/vv1.txt') as f :
        x = f.readlines()
    for i in range (0,len(x)) :
        x[i] = x[i].replace(',' ,'.')
        Vvent.append(x[i][:len(x[i])-1])

  
    Vvent = [float(i) for i in Vvent ]
    Vvent = [float(i) for i in Vvent ]
    Vvent = Vvent[6:len(Vvent):5]
    Vvent[-1] = 0.08
    Vvent = np.array(Vvent)

    
    g = []
    with open(r'/content/g1.txt') as f :
       x = f.readlines()
    for i in range (0,len(x)) :
        x[i] = x[i].replace(',' ,'.')
        g.append(x[i][:len(x[i])-1])


    G= [float(i) for i in g ]
    G = np.array(G)

    Tamb = []
    with open(r'/content/tamb.txt') as f :
        x = f.readlines()
    for i in range (0,len(x)) :
        x[i] = x[i].replace(',' ,'.')
        Tamb.append(x[i][:len(x[i])-1])


    Tamb= [float(i)+273 for i in Tamb[:-1] ]
    Tamb = np.array(Tamb)
    l=0.79 
    Ltube=1.62 
    Dh=0.012 
    kc=0.036 
    kt=0.033 
    ec=0.0003 
    et=0.0005 
    kail=230 
    eail=0.01 
    cps=897 
    rops=2700    
    ktub=387 
    rotube=8960 
    ctub=385 
    etub=0.002 
    kiso=0.039 
    eiso=0.07 
    roiso=40 
    segma=5.67*10**-8 
    av=0.06 
    ac=0.85 
    at=0.5 
    tovi=0.84  
    bettapv=1/(270-25)
    Tcref=298 
    nuref=0.17
    Ntube=12 
    L=1.58 
    pi =math.pi
    sec = pi*(0.022**2)/4
    sec1 = pi*(0.012**2)/4
    V=Qnf/sec
    s=L*l 
    s2=pi*Dh*(1/4)*L*Ntube 
    s1=s-s2 
    sl=3.14*Dh*L*Ntube 
    s3=pi*Dh*((1/2)+(1/4))*L*Ntube 
    s4=(pi+2)*Dh*(1/2)*L*Ntube+s1 
    cc=677 
    rot=1200 
    ct=1250 
    cpail=897 
    roail=2700 
    cpiso=1500 
    beta=0.88 
    epsv=0.83   
    kv=1 
    roc=2702 
    rov=2530 
    cv=500 
    ev=0.0032 
    kps=230   
    eps=0.001 
    sail=eail*L*6 
    mtub=rotube*pi*(0.012/2)**2-(0.011/2)**2*Ltube*Ntube 
    mc=roc*s*ec 
    mv=rov*s*ev 
    mt=rot*s*et 
    miso=roiso*s*eiso 
    mail=roail*eail*L*0.04*6 
    mps=rops*(s4)*eps



    Tvold =[None]*167
    Tnfe= [25+273 ]*167
    Tcold=[None]*167
    Ttold=[None]*167
    Ttubeold=[None]*167
    Teauold=[None]*167
    Tpsold=[None]*167 
    Tailold=[None]*167
    Tisoold=[None]*167
    Tv=[None] *167
    Tc=[None]*167
    Tt=[None] *167
    Ttube=[None] *167
    Teau=[None] *167
    Tps=[None] *167
    Tail=[None] *167
    Tiso=[None] *167
    Tciel=[None] *167
    knf=[None] *167
    cpnf=[None] *167
    munf=[None] *167
    ronf=[None] *167
    #
    keau=[None] *167
    cpeau=[None] *167
    mueau=[None] *167
    roeau=[None] *167
    #
    pr= [None] *167
    Re= [None] *167
    Hr_vciel  =[None] *167
    Hcv_v =[None] *167
    Hcv_nftub=[None] *167
    Teau2 = [None] *167
    puissance = [None] *167
    ntheau = [None] *167
    Tvold[0]=320 
    Tcold[0]=310 
    Ttold[0]=300 
    Ttubeold[0]=299 
    Teauold[0]=298 
    Tpsold[0]=292 
    Tailold[0]=292 
    Tisoold[0]=292 
    mmax=96 
    t =np.arange(167)

    dt=300 
    tol=10**-5 
    err=tol+1 
    m=0 
    Tv[m]=Tvold[m] 
    Tc[m]=Tcold[m] 
    Tt[m]=Ttold[m] 
    Ttube[m]=Ttubeold[m] 
    Teau[m]=Teauold[m] 
    Tps[m]=Tpsold[m] 
    Tail[m]=Tailold[m] 
    Tiso[m]=Tisoold[m] 
    iter1=0 
    Hcd_vc=1/(ev/kv+ec/kc) #%conduction entre le verre et la cellule 
    Hcd_ct=1/(ec/kc+et/kt) #%conduction entre la cellule et tedlar
    Hcd_ttub=1/(et/kt+etub/ktub) #%conduction entre tedlar et le tube  
    Hcd_tps=1/(eps/kps+et/kt) #%conduction entre la plaque  et tedlar 
    Hcd_tubps=1/(etub/ktub+eps/kps) #%conduction entre le tube et la plaque 
    Hcd_psail=1/(eps/kps+eail/kail) #%conduction entre la plaque et l'aillette
    Hcd_psiso=1/(eps/kps+eiso/kiso) #%cnduction entre la plaque  et l'isolant
    Hcd_ailiso=1/(eail/kail+eiso/kiso)  #%conduction entre l'ailette et l'isolant
    tol=10**-5 
    err=tol+1
    while err>tol:
        iter1=iter1+1 
        m=0 
        Tciel[m]=0.0552*Tamb[m]**1.5 
        knf[m] =-7843*10**-9*Teauold[m]**2 +0.0062*Teauold[m] - 0.54 
        ronf[m] =-0.003 *Teauold[m]**2 + 1.505*Teauold[m] + 816.781 
        cpnf[m] =-0.0000643 *Teauold[m]**3 + 0.0552 *Teauold[m]**2 - 20.86 *Teauold[m] + 6719.637
        munf[m] = (2414*10**-8)*10**(247.8/(Teauold[m]-140)) 
        pr[m]=cpnf[m]*munf[m]/knf[m] 
        Re[m]=ronf[m]*V*Dh/munf[m] 
    #%les coefficient de transfert 
    #%1-RAYONNEMENT   
    #% Verre - voûte céleste
        Hr_vciel[m]=segma*epsv*(Tvold[m]**2+Tciel[m]**2)*(Tvold[m]+Tciel[m]) 
    #%2 CONDUCTION 
    #%1] le verre et l'aire    
        Hcv_v[m]=5.7+3.8*float(Vvent[m])  
    # %2 le nano-fluide et le tube
        if  Re[m]<=2300:
            Hcv_nftub[m]=4.36*knf[m]/Dh 
        elif  Re[m]>2300:
            Hcv_nftub[m]=0.023*(knf[m]/Dh)*(Re[m]**0.8)*(pr[m]**0.4) 
        Tv[m]=(Tciel[m]*Hr_vciel[m]+Tcold[m]*Hcd_vc+Tamb[m]*Hcv_v[m]+av*G[m])/(Hr_vciel[m]+Hcd_vc+Hcv_v[m]) 
        Tc[m]=(Tv[m]*Hcd_vc+Ttold[m]*Hcd_ct-nuref*(1+bettapv*Tcref)*G[m]*tovi*beta+tovi*ac*G[m]*beta)/(Hcd_vc+Hcd_ct-nuref*bettapv*G[m]*tovi*beta) 
        Tt[m]=(Tcold[m]*Hcd_ct*s+Tpsold[m]*Hcd_tps*s1+Ttubeold[m]*Hcd_ttub*s2+tovi*s*at*G[m]*(1-beta))/(Hcd_ct*s+Hcd_tps*s1+Hcd_ttub*s2) 
        Ttube[m] = (Ttold[m]*Hcd_ttub*s2 + Tpsold[m]*Hcd_tubps*s3 +Teauold[m]*Hcv_nftub[m]*sl)/(Hcd_tubps*s3 + Hcd_ttub*s2 + Hcv_nftub[m]*sl) 
        A =(Hcv_nftub[m]*sl/Ltube)/(ronf[m]*V*sec*cpnf[m]) 
        Teau[m]=Ttubeold[m]+Ttubeold[m]*(math.exp(-A*Ltube)-1)/(A*Ltube )-Tnfe[m]*(math.exp(-A*Ltube)-1)/(A*Ltube)
        Tps[m]=(Ttold[m]*Hcd_tps*s1+Ttubeold[m]*Hcd_tubps*s3+Tisoold[m]*Hcd_psiso*s4)/(Hcd_tps*s1+Hcd_tubps*s3+Hcd_psiso*s4) 
        Tail[m]=(Tpsold[m]*Hcd_psail*sail+Tisoold[m]*Hcd_ailiso*sail)/(Hcd_psail*sail+Hcd_ailiso*sail) 
        Tiso[m]=(Tpsold[m]*Hcd_psiso*s4+Tailold[m]*Hcd_ailiso*sail)/(Hcd_psiso*s4+Hcd_ailiso*sail) 
        Teau2[m]=Ttube[m]*(1-math.exp(-A*Ltube))+Tnfe[m]*math.exp(-A*Ltube) 
        puissance[m]=(ronf[m]*V*sec*cpnf[m])*(Teau2[m]-Tnfe[m]) 
        ntheau[m]=puissance[m]/(G[m]*s)
        erreur = [None]*8 
        erreur[1]=Tv[m]-Tvold[m] 
        erreur[2]=Tc[m]-Tcold[m] 
        erreur[3]=Tt[m]-Ttold[m] 
        erreur[4]=Ttube[m]-Ttubeold[m] 
        erreur[5]=Teau[m]-Teauold[m] 
        erreur[6]=Tps[m]-Tpsold[m] 
        erreur[7]=Tail[m]-Tailold[m] 
        erreur =list(map(abs,erreur[1:]))
        err=max(erreur)
        m=0
        Tvold[m]=Tv[m] 
        Tcold[m]=Tc[m] 
        Ttold[m]=Tt[m] 
        Ttubeold[m]=Ttube[m] 
        Teauold[m]=Teau[m] 
        Tpsold[m]=Tps[m] 
        Tailold[m]=Tail[m] 
        Tisoold[m]=Tiso[m] 

    m=0
    Tvold[m+1]=Tv[m] 
    Tcold[m+1]=Tc[m] 
    Ttold[m+1]=Tt[m] 
    Ttubeold[m+1]=Ttube[m] 
    Teauold[m+1]=Teau[m] 
    Tpsold[m+1]=Tps[m] 
    Tailold[m+1]=Tps[m]
    Tisoold[m+1]=Tiso[m] 
    tol=10**-5 
    knp = 32.9
    ronp = 6310
    cpnp = 551
    dnp=50e-9
    deau=0.434*10e-9
    phi= 0.1
    for m in range(1,mmax):
        err=tol+1
        iter1 = 0
        while err>tol  :
            iter1 = iter1+1
            Tciel[m]=0.0552*Tamb[m]**1.5 
            #proprites for the pure water 
            keau[m] =-7843*10**-9*Teauold[m]**2 +0.0062*Teauold[m] - 0.54 
            roeau[m] =-0.003 *Teauold[m]**2 + 1.505*Teauold[m] + 816.781 
            cpeau[m] = -0.0000643 *Teauold[m]**3 + 0.0552 *Teauold[m]**2 - 20.86 *Teauold[m] + 6719.637 
            mueau[m] = (2414*10**-8)*10**(247.8/(Teauold[m]-140))
            # correlation for the nanofluids :
            knf[m] =(knp+2*keau[m]-2*phi*(keau[m]-knp))/(knp+2*keau[m]+phi*(keau[m]-knp))*keau[m]
            ronf[m] =(1-phi)*roeau[m] +phi*ronp
            cpnf[m] = (1-phi)*cpeau[m] +phi*cpnp
            munf[m] = mueau[m]*(1/(1-34.87*(dnp/deau)**-0.3*phi**1.03))
            ######
            Re[m]=ronf[m]*V*Dh/munf[m] 
            pr[m]=cpnf[m]*munf[m]/knf[m] 
            Hr_vciel[m]=segma*epsv*(Tvold[m]**2+Tciel[m]**2)*(Tvold[m]+Tciel[m]) 
            Hcv_v[m]=5.7+3.8*float(Vvent[m])  
            if Re[m]<=2300 :
                Hcv_nftub[m]=4.36*knf[m]/Dh 
            elif  Re[m]>2300 :
                Hcv_nftub[m]=0.023*(knf[m]/Dh)*(Re[m]**0.8)*(pr[m]**(0.4)) 
            Tv[m]=(Tv[m-1]*mv*cv+dt*(Tciel[m]*Hr_vciel[m]*s+Tcold[m]*(Hcd_vc*s)+Tamb[m]*Hcv_v[m]*s+av*G[m]*s))/(mv*cv+dt*s*(Hr_vciel[m]+Hcd_vc+Hcv_v[m])) 
            Tc[m]=(Tc[m-1]*mc*cc+dt*(Tvold[m]*Hcd_vc*s+Ttold[m]*Hcd_ct*s-nuref*(1+bettapv*Tcref)*tovi*beta*G[m]*s+tovi*ac*G[m]*beta*s))/(mc*cc+dt*s*(Hcd_vc+Hcd_ct-nuref*bettapv*G[m]*tovi*beta)) 
            Tt[m]=(Tt[m-1]*mt*ct+dt*(Tcold[m]*Hcd_ct*s+Tpsold[m]*Hcd_tps*s1+Ttubeold[m]*Hcd_ttub*s2+tovi*s*at*G[m]*(1-beta)))/(mt*ct+dt*(Hcd_ct*s+Hcd_tps*s1+Hcd_ttub*s2)) 
            Ttube[m]=(Ttube[m-1]*mtub*ctub+dt*(Ttold[m]*Hcd_ttub*s2+Tpsold[m]*Hcd_tubps*s3+ Teauold[m]*Hcv_nftub[m]*sl))/(mtub*ctub+dt*(Hcd_tubps*s3+Hcd_ttub*s2+Hcv_nftub[m]*sl)) 
            A =(Hcv_nftub[m]*sl/L)/(ronf[m]*V*sec*cpnf[m]) 
            Teau[m]=Ttubeold[m]+Ttubeold[m]*(math.exp(-A*Ltube)-1)/(A*Ltube )-Tnfe[m]*(math.exp(-A*Ltube)-1)/(A*Ltube) 
            Tps[m]=(Tps[m-1]*cps*mps+dt*(Ttold[m]*Hcd_tps*s1+Ttubeold[m]*Hcd_tubps*s3+Tisoold[m]*Hcd_psiso*s4))/(mps*cps+dt*(Hcd_tps*s1+Hcd_tubps*s3+Hcd_psiso*s4)) 
            Tail[m]=(Tail[m-1]*cpail*mail+dt*(Tpsold[m]*Hcd_psail*sail+Tisoold[m]*Hcd_ailiso*sail))/(mail*cpail+dt*(Hcd_psail*sail+Hcd_ailiso*sail)) 
            Tiso[m]=(Tiso[m-1]*miso*cpiso+dt*(Tpsold[m]*Hcd_psiso*s4+Tailold[m]*Hcd_ailiso*sail))/(miso*cpiso+dt*(Hcd_psiso*s4+Hcd_ailiso*sail)) 
            Teau2[m]=Ttube[m]*(1-math.exp(-A*Ltube))+Tnfe[m]*math.exp(-A*Ltube) 
            puissance[m]=(ronf[m]*V*sec1*cpnf[m]*(Teau2[m]-Tnfe[m]))
            ntheau[m]=puissance[m]/(G[m]*s) 
            erreur = [None]*9
            erreur[1]=Tv[m]-Tvold[m] 
            erreur[2]=Tc[m]-Tcold[m] 
            erreur[3]=Tt[m]-Ttold[m] 
            erreur[4]=Ttube[m]-Ttubeold[m] 
            erreur[5]=Teau[m]-Teauold[m] 
            erreur[6]=Tps[m]-Tpsold[m] 
            erreur[7]=Tail[m]-Tailold[m] 
            erreur[8]=Tiso[m]-Tisoold[m] 
            erreur =list(map(abs,erreur[1:]))
            err=max(erreur)
            Tvold[m]=Tv[m] 
            Tcold[m]=Tc[m] 
            Ttold[m]=Tt[m] 
            Ttubeold[m]=Ttube[m] 
            Teauold[m]=Teau[m] 
            Tpsold[m]=Tps[m] 
            Tailold[m]=Tail[m] 
            Tisoold[m]=Tiso[m]
        Tvold[m+1]=Tv[m] 
        Tcold[m+1]=Tc[m] 
        Ttold[m+1]=Tt[m] 
        Ttubeold[m+1]=Ttube[m] 
        Teauold[m+1]=Teau[m] 
        Tpsold[m+1]=Tps[m] 
        Tailold[m+1]=Tail[m] 
        Tisoold[m+1]=Tiso[m]
    G_s = [i*s for i in G]
    red_t_t = np.array(ronf[:96])*V*sec*np.array(cpnf[:96])*(np.array(Teau2[:96])-np.array(Tnfe[:96]))
    rend_t = red_t_t/np.array(G_s[:-1])
    rendement_exg_th  = np.array(ronf[:96])*V*sec*np.array(cpnf[:96])*((np.array(Teau2[:96])-np.array(Tnfe[:96]))-np.array(Tamb[:96])*np.log(np.array(Teau2[:96])/np.array(Tnfe[:96])))
    rendement_exg_th = abs(rendement_exg_th/(G[:96]*s*(1-np.array(Tamb)[:96]/5778)))
    rend_e = [rendement_electrique(i) for i in Tc[:96]]
    rendement__exergetique_elec_th = np.array(rend_e)/(1-(np.array(Tamb[:-1])/5778))

    return Hcv_nftub


tc_nfal = solve_al(5e-5)
