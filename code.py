import numpy as np
import math
from numpy import linalg
from numpy.lib.function_base import average 


def masse_verre(Ro_verre,s,ev):
    return  Ro_verre*s*ev       

def masse_cellule(Ro_cellule,s,ev):
    return  Ro_cellule*s*ev       

def c_verre(Tv):
    return 500
def c_cellule(Tc):
    return 677

def h_r_ciel_verre(emissivite_v,Tv,Tambiante):
    Tciel = 0.0552 * Tambiante**1.5
    return 5.67e-8 *emissivite_v*(Tv+Tciel)*(Tv**2+Tciel**2)


def h_conduction_cellule_verre(kv,ev,kc,ec):
    return 1/(ec/kc+ev/kv)

def h_conection_air_verre(Vvent):
    return 5.7+3.86*Vvent

def rendement_electrique (Tc ,Tcrefrence=298, nrefrence = 0.12,betapv = 0.0045 ):
    return nrefrence*(1-betapv*(Tc-Tcrefrence))

 
def h_conduction_cellule_teddler(kc,ec,kd,ed):
    return 1/(ec/kc+ed/kd)

 
def h_conduction_ailete_ps(el,ke,eps,kps):
    return 1/(el/ke+eps/kps)

def h_conduction_aillete_isolant(ke,el,kisolant,eisolant) : 
    return  1/(el/ke+eisolant/kisolant)


def h_conduction_tube_tedler(kted,eted,ktube,etube) : 
    return 1/(eted/kted+etube/ktube )

def h_conduction_tube_ps(kps,eps,ktube,etube) : 
    return 1/(etube/ktube + eps/kps)


def h_conduction_ps_tedler(kted,eted,kps,eps) : 
     return 1/(eted/kted + eps/kps)

def h_connduction_tube_isolant(kisolant,eisolant,ktube,etube):
    return 1 / (eisolant/kisolant + etube/ktube)

def h_connduction_ps_isolant(kisolant,eisolant,kps,eps):
    return 1/(eisolant/kisolant+eps/kps)

def h_convection_tube_nano(D,V,knanof,cpnano, ronanof,visnanof):
    re = V * ronanof * D / visnanof
    pr = cpnano * visnanof /knanof
    if re <=2300 :
        return 4.86*knanof/D
    elif re> 2300 :
        return 0.023 *knanof/D* pr **0.4 * re**0.8 
class  slower_probleme:
    def __init__(self,av,emissivite_v,G,Tnf,kv,kc,Tps,ktube,etube,kted,eted,Tisolante,Tailiette,ev,el,ke,eps,kps,ec,kisolant,eisolant,Vvent,Tc,Ro_verre,Ro_cellule,Tambiante,Tv,Ttube, tovi,beta,ac,at,V,D,visnanof,l,L,Ntube):
        self.L = L
        self.l = l
        self.s = self.l * self.L
        self.av = av
        self.emissivite_v = emissivite_v
        self.G = G
        self.el = el 
        self.ke = ke 
        self.eps = eps 
        self.kps = kps 
        self.D = 0.012
        self.kc = kc
        self.ev = ev
        self.kv = kv
        self.ec = ec 
        self.kisolant = kisolant
        self.eisolant = eisolant
        self.kted = kted
        self.eted =eted
        self.ktube = ktube
        self.etube = etube
        self.Vvent = Vvent
        self.tovi = 0.84  
        self.beta = beta
        self.ro_tedler=1200
        self.ro_tube=8960
        self.ro_plaque=2700
        self.ro_aillete=2700
        self.ro_isol=20
        self.c_tedler=1250
        self.c_tube=385
        self.c_plaque=897
        self.c_aillete=897
        self.c_isol=1500
        self.k_verre= 1
        self.k_cellue=0.036
        self.k_tedler=0.033
        self.k_tube=380
        self.k_plaque=230
        self.k_aillete=230
        self.k_isol = 0.040

        self.ac = ac
        self.at = at
        self.V = 0.01
        self.Tnf = Tnf
        self.kanonf = -7843e-9*self.Tnf**2 +0.0062*self.Tnf - 0.54
        self.ronanof =-0.003 *self.Tnf**2 + 1.505*self.Tnf + 816.781
        self.cpnano = -0.0000643 * self.Tnf**3 + 0.0552 *self.Tnf**2 - 20.86 *self.Tnf + 6719.637 
        self.visnanof = (2414e-8)*10**(247.8/(self.Tnf-140))
        self.Tv = Tv
        self.Tc = Tc
        self.Tt = Tt
        self.Tps = Tps
        self.Tambiante = Tambiante
        self.Ttube = Ttube
        self.Ltube = 1.62
        self.Tailiette =Tailiette
        self.Ntube=Ntube
        self.Ro_cellule=Ro_cellule
        self.Ro_verre=Ro_verre
        self.sl = math.pi*self.Ltube*self.D*12
        self.s2 = math.pi* self.D/3 * self.L*self.Ntube
        self.s1 = self.s - self.s2
        self.s3 = self.sl - self.s2
        self.s4 = self.s1 +self.s3 +self.D*L*12 
        self.Tps = Tps
        self.Tisolante = Tisolante 
        self.masse_verre = masse_verre (self.Ro_verre,self.s,self.ev)
        self.masse_cellule = masse_cellule (self.Ro_cellule,self.s,self.ec)
        self.h_r_ciel_verre = h_r_ciel_verre(self.emissivite_v,self.Tv,self.Tambiante)
        self.h_conduction_cellule_verre = h_conduction_cellule_verre(kv,ev,kc,ec)
        self.c_verre = c_verre(self.Tv)
        self.c_cellule = c_cellule(self.Tc)
        self.h_conection_air_verre = h_conection_air_verre(Vvent)
        self.rendement_electrique = rendement_electrique(self.Tc ,Tcrefrence=298, nrefrence = 0.12,betapv = 0.0045 )
        self.h_conduction_cellule_teddler = h_conduction_cellule_teddler(kc,ec,kted,eted)
        self.h_conduction_aillete_isolant = h_conduction_aillete_isolant(ke,el,kisolant,eisolant)
        self.h_conduction_tube_tedler = h_conduction_tube_tedler(kted,eted,ktube,etube)
        self.h_conduction_tube_ps =h_conduction_tube_ps(kps,eps,ktube,etube)
        self.h_connduction_tube_isolant = h_connduction_tube_isolant(kisolant,eisolant,ktube,etube)
        self.h_convection_tube_nano = h_convection_tube_nano(self.D,self.V,self.kanonf,self.cpnano, self.ronanof,self.visnanof)
        self.Tciel = 0.0552 * Tambiante**1.5
        self.sec= math.pi /4* self.D**2 
        self.h_conduction_ps_tedler= h_conduction_ps_tedler(kted,eted,kps,eps)
        self.h_conduction_ailete_ps = h_conduction_ailete_ps(el,ke,eps,kps)
        self.h_connduction_ps_isolant = h_connduction_ps_isolant(kisolant,eisolant,kps,eps)
    def construction_of_matrix_A(self,dt):
        a=np.zeros(64).reshape(8,8) 
        b = np.zeros(8)
        c = [self.c_verre,self.c_cellule, self.c_tedler,self.c_tube,self.c_plaque,self.c_aillete,self.c_isol]
        ros = [  self.Ro_verre,self.Ro_cellule, self.ro_tedler,self.ro_tube,self.ro_plaque,self.ro_aillete,self.ro_isol]
        eppiseur = [self.ev,self.ec,self.eted,self.etube,self.eps,self.el,self.eisolant]
        p=np.zeros(len(ros))
        for i in range(len(ros)) : 
            p[i] = dt/(c[i] * ros[i]*self.sl*eppiseur[i])
        p[3] = dt/(c[3] * ros[3]*self.sl*eppiseur[3])
        p[5] = dt / (c[i] * ros[i]*self.L*6*eppiseur[5]*0.04)
        a[0,0] =1+(self.h_r_ciel_verre * self.s + self.h_conduction_cellule_verre*self.s + self.h_conection_air_verre*self.s)*p[0]
        a[0,1] = -self.h_conduction_cellule_verre*self.s*p[0]
        a[1,1]= 1+p[1]*self.s*(self.h_conduction_cellule_verre*self.beta + self.h_conduction_cellule_teddler*self.beta -self.tovi*0.12*self.beta*self.G*0.0045 )
        a[1,0] = -self.s * p[1] *self.h_conduction_cellule_verre*self.beta 
        a[1,2] = -self.s*p[1]*self.h_conduction_cellule_teddler*self.beta
        a[2,1] = -self.h_conduction_cellule_teddler * self.s*p[2] * self.beta
        a[2,2] = 1 + p[2]*(self.h_conduction_cellule_teddler*self.s*self.beta +self.h_conduction_tube_tedler*self.s2 + self.h_conduction_ps_tedler*self.s1)
        a[2,3] = -(self.h_conduction_tube_tedler * p[2] *self.s2)
        a[2,5] = -self.h_conduction_ps_tedler * self.s1* p[2]
        a[3,2] = -self.h_conduction_tube_tedler * self.s2 * p[3] 
        a[3,3] = 1 + p[3]*(self.h_conduction_tube_ps*self.s3 +self.h_conduction_tube_tedler*self.s2 + self.h_convection_tube_nano*self.sl)
        a[3,4] = -self.h_convection_tube_nano*self.sl*p[3]
        a[3,5] = -self.h_conduction_tube_ps * self.s3 *p[3]
        A = (self.h_convection_tube_nano * self.sl /self.Ltube) / (self.ronanof*self.V * self.sec * self.cpnano)
        a[4,3] = -1-((math.exp(-A  * self.Ltube)-1)/(A* self.Ltube ))
        a[5,2] = -(self.h_conduction_ps_tedler*self.s1 * p[4] )
        a[4,4] = 1
        a[5,3] = -(self.h_conduction_tube_ps * self.s3 *p[4])
        a[5,6] = -(self.h_conduction_ailete_ps * self.el * self.L *6*p[4])
        a[5,7] = -(self.h_connduction_ps_isolant* self.s4 * p[4]  )
        a[5,5] = 1 + p[4]*(self.h_conduction_ps_tedler * self.s1 + self.h_conduction_tube_ps*self.s3 +self.h_conduction_ailete_ps * self.el *self.L*6 + self.h_connduction_ps_isolant * self.s4  )
        a[6,5] = -(self.h_conduction_ailete_ps*self.el * p[5] *6* self.L)
        a[6,6] = 1 +(self.h_conduction_ailete_ps  + self.h_conduction_aillete_isolant ) * p[5] *self.el * self.L*6
        a[6,7] = -(self.h_conduction_aillete_isolant * p[5] * self.el * self.L*6)
        a[7,6] = -(self.h_conduction_aillete_isolant * self.el*self.L*6) *p[6]
        a[7,5] = -(self.h_connduction_ps_isolant *self.s4) *p[6]
        a[7,7] = 1 + p[6]*(self.h_connduction_ps_isolant * self.s4 + self.h_conduction_aillete_isolant * self.el * self.L*6 )
        b[0] = self.Tv +p[0]*(self.h_r_ciel_verre*self.s*self.Tciel + self.av *self.s*self.G +self.h_conection_air_verre*self.s* self.Tambiante )
        b[1] = self.Tc +self.s*p[1]*(self.tovi*self.ac*self.G*self.beta - self.tovi*0.12* self.beta *self.G*(1+0.0045*298))
        b[2] = self.Tt + self.tovi*self.at*self.G*self.s *(1-self.beta) * p[2] 
        b[3] = self.Ttube 
        b[4] = -(298*math.exp(-A * self.Ltube)-1)/(A*self.Ltube)
        b[5] = self.Tps
        b[6] = self.Tailiette
        b[7] = self.Tisolante
        return a,b
    #def solove():

    #print("number of iteration :"  ,)
    #print("soluion is : " ),
    def slove(self,t_in ,t_f,pas ):
        z =int((t_f-t_in)/pas)
        solution = np.empty([z,8])
        solution[0] = [self.Tv ,self.Tc ,self.Tt,self.Ttube ,self.Tnf,self.Tps ,self.Tailiette,self.Tisolante ]
        for i in range (t_in+1,z,1):
            a,b = self.construction_of_matrix_A(pas)
            solve = crammer(a,b)
            self.Tv = solve[0]
            self.Tc = solve[1]
            self.Tt = solve[2]
            self.Ttube = solve[3]
            self.Tnf =solve[4]
            self.Tps = solve[5]
            self.Tailiette= solve[6]
            self.Tisolante = solve[7]
            self.kanonf = -7843e-9*self.Tnf**2 +0.0062*self.Tnf - 0.54
            self.ronanof =-0.003 *self.Tnf**2 + 1.505*self.Tnf + 816.781
            self.cpnano = -0.0000643 * self.Tnf**3 + 0.0552 *self.Tnf**2 - 20.86 *self.Tnf + 6719.637 
            self.visnanof = (2414e-8)*10**(247.8/(self.Tnf-140))
            self.h_r_ciel_verre = h_r_ciel_verre(self.emissivite_v,self.Tv,self.Tambiante)
            self.h_convection_tube_nano = h_convection_tube_nano(self.D,self.V,self.kanonf,self.cpnano, self.ronanof,self.visnanof)
            solution[i] = solve
        return solution  
   
def gauss_siedel(a,b,tol):
    x=np.ones(8)
    # methode de gausse siedel  
    error = tol+1
    k =0

    while error > tol  :
          m,n = a.shape
          xn = [1,1,1,1,1,1,1,1]
          for i in range(m):
              s = 0
              for j in range(n):
                  s = s+a[i,j]*x[j]
              s = s - a[i,i]*x[i]
              xn[i] = (b[i]-s ) / a[i,i]
          xn=np.array(xn)
          error = max(abs(xn-x))
          x = xn
          k = k+1
    #print("x :" , x)
    #print("k : " ,k)
    #print("x :" , x)
    #print('error' ,error)
    #print(np.dot(a,x)-b)
    return x 


def crammer(a,b) : 
    #ref = np.array([a,a,a])
    #det_ref = np.linalg.det(a)
    #x = np.array([5,1,4])
    #for i in range(len(a)): 
        #dm = ref[i]
        #dm[:,i] = b
        #x[i] = np.linalg.det(dm)/det_ref
    return np.linalg.solve(a,b)
 


l=0.79;
Ltube=1.62;
kc=0.036;
kted=0.033;
ec=0.0003;
et=0.0005;
cps=897;
rops=2700;
ktube=380;
rotube=8960;
ctube=385;
etube=2e-3;
kiso=0.039;
eiso=0.07
roiso=40
av=0.06
ac=0.85
at=0.8
tovi=0.84
Tcref=298
nuref=0.12
Ntube=12
L=1.58
beta=0.88
kv=1
Ro_cellule=2703
Ro_verre=2530
cv=200
ev=0.0032
kps=230
eps=0.001
G = 300
Tnf = 298
Tps = 298
Tailiette =298
Tc =298
Tambiante =298
Tv = 298
Tisolante = 298
Tt = 298
emissivite_v = 0.83
Vvent = 0.5
eted = 0.0005
D = 0.012
el = 0.001 
ke = 230 
kisolant  = 0.039
eisolant  = 0.07
Ttube = 291.65
V = 0.001
kd = 0.033
ed =0.0005
probleme = slower_probleme(av,emissivite_v,G,Tnf,kv,kc,Tps,ktube,etube,kted,eted,Tisolante,Tailiette,ev,el,ke,eps,kps,ec,kisolant,eisolant,Vvent,Tc,Ro_verre,Ro_cellule,Tambiante,Tv,Ttube, Tt,tovi,beta,ac,at,V,D,l,L,Ntube)
re = probleme.slove(0,36000,300)
print(re[0])
print(re[119])
# '''Data for plotting
#t = np.arange(0.0, 100, 1)
#s = re[:,0]
#s1 = re[:,1]
#s2 = re[:,2]
#s3 = re[:,3]
#s4 = re[:,4]
#s5 = re[:,5]
#s6 = re[:,6]
#s7 = re[:,7]
#fig, ax = plt.subplots(figsize=(16,9))
#ax.plot(t, s, label = 'tv')
#ax.plot(t,s1,label = 'tc')
#ax.plot(t,s2,label = 'tt')
#ax.plot(t,s3,label = 'ttube')
#ax.plot(t,s4,label = 'tn')
#ax.plot(t,s5,label = 'tp')
#ax.plot(t,s6,label = 'tai')
#ax.plot(t,s7,label = 'tiso')
#ax.set(xlabel='time (s)', ylabel='temerature (k)')
#ax.grid()
#plt.legend()
#plt.show()


def gauss_siedel1(a,b,tol):
    x=np.ones(3)
    # methode de gausse siedel  
    error = tol+1
    k =0
    xn=[4,5,8]
    while error > tol  :
          m,n = a.shape
          for i in range(m):
              s = 0
              for j in range(n):
                  s = s+a[i,j]*x[j]
              s = s - a[i,i]*x[i]
              xn[i] = (b[i]-s ) / a[i,i]
          xn=np.array(xn)
          error = max(abs(xn-x))
          x = xn
          k = k+1
    #print("x :" , x)
    #print("k : " ,k)
    #print("x :" , x)
    #print('error' ,error)
    #print(np.dot(a,x)-b)
    return x 


#schma explicite 