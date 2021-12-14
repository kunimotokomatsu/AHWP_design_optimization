import numpy as np
#import pylab as py
import sys
from scipy.optimize import curve_fit
import scipy.constants as constants
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c=constants.c #m/s
pi=np.pi
e=np.e
radeg = 180./pi
ep0 = constants.epsilon_0 #s4 A2/m3 kg
u0=constants.mu_0 #H/m=V/(A/s)/m
im=complex(0.,1.)
arcmin=1./60.
h=constants.h #m^2 kg s-1
kb=constants.k #m^2 kg s^-2 K^-1
Tcmb=2.7255 #K
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#refrerence Thomas Essinger-Hileman, arXiv:1301.6160v3 [physics.optics] 24 Mar 2014
#n1,no,ne=refractive index of air, ordinary, extraordinay
#theta1=incident angle [rad] (in x-z plane)
#chi=extraordinary angle [rad] (from x-axis)
#pol='p' or 's', 'p'=x, 's'=y
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def rms(x):
    x=np.array(x)
    return np.sqrt(np.sum(x**2)/float(len(x)))

def get_Trans(air,air_err,sample,sample_err):
    t=sample/air
    err=np.sqrt((sample_err/air)**2.+(sample*air_err/(np.array(air)**2.))**2.)
    return t,err

def read_txt2f(filename):
    import fileinput
    arr1 = []
    arr2 = []
    file = open(filename,'r')
    for i in file.readlines():
        i=i.split()
        arr1.append(float(i[0]))
        arr2.append(float(i[1]))
    file.close()
    return np.array(arr1),np.array(arr2)

def read_txt3f(filename):
    import fileinput
    arr1 = []
    arr2 = []
    arr3 = []
    file = open(filename,'r')
    for i in file.readlines():
        i=i.split()
        arr1.append(float(i[0]))
        arr2.append(float(i[1]))
        arr3.append(float(i[2]))
    file.close()
    return np.array(arr1),np.array(arr2),np.array(arr3)

def read_txt4f(filename):
    import fileinput
    arr1 = []
    arr2 = []
    arr3 = []
    arr4 = []
    file = open(filename,'r')
    for i in file.readlines():
        i=i.split()
        arr1.append(float(i[0]))
        arr2.append(float(i[1]))
        arr3.append(float(i[2]))
        arr4.append(float(i[3]))
    file.close()
    return np.array(arr1),np.array(arr2),np.array(arr3),np.array(arr4)

#def fit_func(t,a0,a1,a2,a3,a4):
#    return a0+a1*np.cos(2.*t+2.*a2)+a3*np.cos(4.*t+4.*a4)

def fit_func(theta,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16):
    return a0+a1*np.cos(theta+a2)+a3*np.cos(2.*theta+2.*a4)+a5*np.cos(3.*theta+3.*a6)+a7*np.cos(4.*theta+4.*a8)+a9*np.cos(5.*theta+5.*a10)+a11*np.cos(6.*theta+6.*a12)+a13*np.cos(7.*theta+7.*a14)+a15*np.cos(8.*theta+8.*a16)

#def fit_func(theta,a0,a1,a2,a3,a4,a5,a6,a7,a8):
#    return a0+a1*np.cos(2.*theta+2.*a2)+a3*np.cos(4.*theta+4.*a4)+a5*np.cos(6.*theta+6.*a6)+a7*np.cos(8.*theta+8.*a8)

#def fit_func(theta,a0,a1,a2,a3,a4,a5,a6,a7,a8):
#    return a0+a1*np.cos(2.*theta+2.*a2)+a3*np.cos(4.*theta+4.*a4)+a5*np.cos(6.*theta+6.*a6)+a7*np.cos(8.*theta+8.*a8)

def expi(x):return e**complex(0,x)


def ygrid():
    m = np.zeros((4,4))
    m[0,0]=m[1,1]=1.
    m[0,1]=m[1,0]=-1.
    return 0.5*m

def xgrid():
    m = np.zeros((4,4))
    m[0,0]=m[0,1]=m[1,0]=m[1,1]=1.
    return 0.5*m

def Jones_to_Mueller(J):
    j11=J[0][0]
    j12=J[0][1]
    j21=J[1][0]
    j22=J[1][1]
    m = np.zeros((4,4),'complex')
    m[0][0]=0.5*(j11*np.conjugate(j11)+j21*np.conjugate(j21)+j12*np.conjugate(j12)+j22*np.conjugate(j22))
    m[0][1]=-0.5*(j22*np.conjugate(j22)-j21*np.conjugate(j21)+j12*np.conjugate(j12)-j11*np.conjugate(j11))
    m[0][2]=0.5*(j21*np.conjugate(j22)+j22*np.conjugate(j21)+j11*np.conjugate(j12)+j12*np.conjugate(j11))
    m[0][3]=-0.5*im*(j21*np.conjugate(j22)-j22*np.conjugate(j21)+j11*np.conjugate(j12)-j12*np.conjugate(j11))
    
    m[1][0]=-0.5*(j22*np.conjugate(j22)+j21*np.conjugate(j21)-j12*np.conjugate(j12)-j11*np.conjugate(j11))
    m[1][1]=0.5*(j22*np.conjugate(j22)-j21*np.conjugate(j21)-j12*np.conjugate(j12)+j11*np.conjugate(j11))
    m[1][2]=-0.5*(j21*np.conjugate(j22)+j22*np.conjugate(j21)-j11*np.conjugate(j12)-j12*np.conjugate(j11))
    m[1][3]=0.5*im*(j21*np.conjugate(j22)-j22*np.conjugate(j21)-j11*np.conjugate(j12)+j21*np.conjugate(j11))
    
    m[2][0]=0.5*(j12*np.conjugate(j22)+j22*np.conjugate(j12)+j11*np.conjugate(j21)+j21*np.conjugate(j11))
    m[2][1]=-0.5*(j12*np.conjugate(j22)+j22*np.conjugate(j12)-j11*np.conjugate(j21)-j21*np.conjugate(j11))
    m[2][2]=0.5*(j11*np.conjugate(j22)+j22*np.conjugate(j11)+j12*np.conjugate(j21)+j21*np.conjugate(j12))
    m[2][3]=-0.5*im*(j11*np.conjugate(j22)-j22*np.conjugate(j11)-j12*np.conjugate(j21)+j21*np.conjugate(j12))
    
    m[3][0]=0.5*im*(j12*np.conjugate(j22)-j22*np.conjugate(j12)+j11*np.conjugate(j21)-j21*np.conjugate(j11))
    m[3][1]=-0.5*im*(j12*np.conjugate(j22)-j22*np.conjugate(j12)-j11*np.conjugate(j21)+j21*np.conjugate(j11))
    m[3][2]=0.5*im*(j11*np.conjugate(j22)-j22*np.conjugate(j11)+j12*np.conjugate(j21)-j21*np.conjugate(j12))
    m[3][3]=0.5*(j11*np.conjugate(j22)+j22*np.conjugate(j11)-j12*np.conjugate(j21)-j21*np.conjugate(j12))
    
    return m

def n_eff(no,ne,chi,n1,theta1):
    return ne*np.sqrt(1.+(1./(np.real(ne)**2.)-1./(np.real(no)**2.))*(np.real(n1)*np.sin(theta1)*np.cos(chi))**2.)

#def n_eff(no,ne,chi,n1,theta1):
#    return ne

def sin_refraction(n,n1,theta1):
    n=np.real(n)
    n1=np.real(n1)
    return n1*np.sin(theta1)/n

def cos_refraction(n,n1,theta1):
    return np.sqrt(1.-sin_refraction(n,n1,theta1)**2)

def transfer_matrix_1plate(freq,no,ne,thickness,chi,n1,theta1):
    neff=n_eff(no,ne,chi,n1,theta1)
    sin_o=sin_refraction(no,n1,theta1)
    sin_e=sin_refraction(neff,n1,theta1)
    cos_o=cos_refraction(no,n1,theta1)
    cos_e=cos_refraction(neff,n1,theta1)
    sin_c=np.sin(chi)
    cos_c=np.cos(chi)
    cos_eo=cos_e*cos_o+sin_e*sin_o
    
    cos_de=np.cos(2.*pi*freq*neff*thickness*cos_e/c)
    sin_de=im*np.sin(2.*pi*freq*neff*thickness*cos_e/c)
    cos_do=np.cos(2.*pi*freq*no*thickness*cos_o/c)
    sin_do=im*np.sin(2.*pi*freq*no*thickness*cos_o/c)
    
    a1=(sin_c**2)*sin_e*sin_o+cos_e*cos_o
    a2=(sin_c**2)*cos_e*cos_eo+(cos_c**2)*(cos_o**3)
    a3=a1*no**2+cos_c**2*sin_o*sin_e*ne**2
    a4=a1*no**2-sin_c**2*sin_o*sin_e*ne**2
    
    t=np.zeros((4,4),'complex')
    t[0][0]=(sin_c**2*a3*cos_do+cos_c**2*a4*cos_de)/(no**2*a1)
    t[0][1]=(sin_c**2*cos_e*cos_o*cos_eo*ne**2*no*sin_do+cos_c**2*cos_o**2*neff*sin_de*a4)/(no**2*ne**2*a2)
    t[0][2]=(cos_c*sin_c*a4*(cos_de-cos_do))/(no**2*a1)
    t[0][3]=cos_c*sin_c*(a4*neff*sin_de-cos_o**2*no*ne**2*sin_do)/(no**2*ne**2*a2)
    
    t[1][0]=(neff*sin_c**2*a3*sin_do+cos_c**2*cos_o**2*no*ne**2*sin_de)/(neff*no*cos_o*a1)
    t[1][1]=(sin_c**2*cos_e*cos_eo*cos_do+cos_c**2*cos_o**3*cos_de)/a2
    t[1][2]=cos_c*sin_c*(cos_o**2*ne**2*no*sin_de-a4*neff*sin_do)/(neff*no*a1*cos_o)
    t[1][3]=cos_c*sin_c*cos_o*(cos_de-cos_do)/a2
    
    t[2][0]=cos_c*sin_c*a3*(cos_de-cos_do)/(no**2*a1)
    t[2][1]=cos_c*sin_c*(neff*cos_o**2*a3*sin_de-cos_e*cos_o*cos_eo*ne**2*no*sin_do)/(ne**2*no**2*a2)
    t[2][2]=(cos_c**2*a4*cos_do+sin_c**2*a3*cos_de)/(no**2*a1)
    t[2][3]=(cos_c**2*cos_o**2*ne**2*no*sin_do+sin_c**2*a3*neff*sin_de)/(ne**2*no**2*a2)
    
    t[3][0]=cos_c*sin_c*(cos_e*cos_eo*ne**2*no*sin_de-neff*cos_o*a3*sin_do)/(neff*no*a1)
    t[3][1]=cos_c*sin_c*cos_e*cos_o**2*cos_eo*(cos_de-cos_do)/a2
    t[3][2]=(neff*cos_c**2*cos_o*a4*sin_do+sin_c**2*cos_e*cos_eo*no*ne**2*sin_de)/(neff*no*a1)
    t[3][3]=(sin_c**2*cos_e*cos_eo*cos_de+cos_c**2*cos_o**3*cos_do)/a2
    
    return t

def transfer_matrix_multilayer(freq,no_arr,ne_arr,thickness_arr,chi_arr,n1,theta1):
    gtmm=transfer_matrix_1plate(freq,no_arr[0],ne_arr[0],thickness_arr[0],chi_arr[0],n1,theta1)
    for i in range(1,len(thickness_arr)):
        gtmm=np.dot(transfer_matrix_1plate(freq,no_arr[i],ne_arr[i],thickness_arr[i],chi_arr[i],n1,theta1),gtmm)
    return gtmm

def Jones_matrix_multilayer(freq,no_arr,ne_arr,thickness_arr,chi_arr,n1,theta1):
    t=transfer_matrix_multilayer(freq,no_arr,ne_arr,thickness_arr,chi_arr,n1,theta1)
    n1=np.real(n1)
    a0=(t[0][0]*np.cos(theta1)+t[0][1]*n1)/np.cos(theta1)
    a1=(t[0][2]+t[0][3]*n1*np.cos(theta1))/np.cos(theta1)
    a2=(t[1][0]*np.cos(theta1)+t[1][1]*n1)/n1
    a3=(t[1][2]+t[1][3]*n1*np.cos(theta1))/n1
    a4=t[2][0]+t[2][1]*n1*np.cos(theta1)
    a5=t[2][2]+t[2][3]*n1*np.cos(theta1)
    a6=(t[3][0]*np.cos(theta1)+t[3][1]*n1)/(n1*np.cos(theta1))
    a7=(t[3][2]+t[3][3]*n1*np.cos(theta1))/(n1*np.cos(theta1))
    a8=1./((a0+a2)*(a5+a7)-(a1+a3)*(a4+a6))
    
    #print freq,chi_arr*radeg,a0,a1,a2,a3,a4,a5,a6,a7,a8
    
    j11=2.*a8*(a5+a7)
    j21=2.*a8*(-a4-a6)
    j12=2.*a8*(-a1-a3)
    j22=2.*a8*(a0+a2)
    return np.array([[j11,j12],[j21,j22]])

def Mueller_matrix_multilayer(freq,no_arr,ne_arr,thickness_arr,chi_arr,n1,theta1):
    j=Jones_matrix_multilayer(freq,no_arr,ne_arr,thickness_arr,chi_arr,n1,theta1)
    return Jones_to_Mueller(j)
    
def IVA_multilayer(freq,hwp_angle,no_arr,ne_arr,thickness_arr,chi_arr,n1,theta1,a_in,p_in,pol):
    if pol=='p':
        grid=xgrid()
    if pol=='s':
        a_in=a_in+90./radeg
        grid=ygrid()
    IVA=[]
    Sin = np.array([1.,p_in*np.cos(2.*a_in),p_in*np.sin(2.*a_in),0.])
    for a in hwp_angle:
        mueller = Mueller_matrix_multilayer(freq,no_arr,ne_arr,thickness_arr,np.array(chi_arr)+a,n1,theta1)
        mueller = np.dot(grid,mueller)
        Sout = np.dot(mueller,Sin)
        IVA.append(np.real(Sout[0]))
    return np.array(IVA)

def insert_gaps(gap_arr,no_arr,ne_arr,thickness_arr,angle_arr,n_gap):
    thickness_arr_gap=[thickness_arr[0]]
    angle_arr_gap=[angle_arr[0]]
    no_arr_gap=[no_arr[0]]
    ne_arr_gap=[ne_arr[0]]
    for i in range(1,len(thickness_arr)):
        thickness_arr_gap=np.r_[thickness_arr_gap,[gap_arr[i-1],thickness_arr[i]]]
        angle_arr_gap=np.r_[angle_arr_gap,[0.,angle_arr[i]]]
        no_arr_gap=np.r_[no_arr_gap,[n_gap,no_arr[i]]]
        ne_arr_gap=np.r_[ne_arr_gap,[n_gap,ne_arr[i]]]
    return [no_arr_gap,ne_arr_gap,thickness_arr_gap,angle_arr_gap]

def add_MFA_layers(MFAthickness,nofMFAlayer,no_arr,ne_arr,thickness_arr,chi_arr,n1):
    none=(no_arr[0]+ne_arr[0])*0.5
    print(none)
    no_MFA=n1+np.arange(1,nofMFAlayer+1,1)*(none-n1)/float(nofMFAlayer+1)
    ne_MFA=n1+np.arange(1,nofMFAlayer+1,1)*(none-n1)/float(nofMFAlayer+1)
    thickness_MFA=np.ones(nofMFAlayer)*MFAthickness/float(nofMFAlayer)
    chi_MFA=np.zeros(nofMFAlayer)
    
    no_arr=np.r_[no_MFA,no_arr,no_MFA[::-1]]
    ne_arr=np.r_[ne_MFA,ne_arr,ne_MFA[::-1]]
    chi_arr=np.r_[chi_MFA,chi_arr,chi_MFA[::-1]]
    thickness_arr=np.r_[thickness_MFA,thickness_arr,thickness_MFA[::-1]]
    return [no_arr,ne_arr,thickness_arr,chi_arr]

def poleff_phase_multilayer(freqin,no_arr,ne_arr,thickness_arr,chi_arr,n1,theta1,a_in,p_in,pol):
    amp0_pre=[]
    amp4_pre=[]
    phase_pre=[]
    hwp_angle=np.arange(0,361,10)/radeg
    for freq in freqin:
        IVA=IVA_multilayer(freq,hwp_angle,no_arr,ne_arr,thickness_arr,chi_arr,n1,theta1,a_in,p_in,pol)
        IVA_err=np.ones(len(IVA))*np.mean(IVA)*0.01
        paramater=np.zeros(17)
        paramater[0]=np.mean(IVA)
        paramater[3]=0.1*(np.max(IVA)-np.min(IVA))
        paramater[7]=0.4*(np.max(IVA)-np.min(IVA))
        fit_popt,fit_cov=curve_fit(fit_func,hwp_angle,IVA,p0=paramater)
        amp0_pre.append(np.abs(fit_popt[0]))
        amp4_pre.append(np.abs(fit_popt[7]))
        if fit_popt[7]<0.:phase_pre.append(fit_popt[8]+0.25*pi)
        else:phase_pre.append(fit_popt[8])
    #sys.stdout.write('\r%3d'%(100*len(poleff_pre)/len(freqin))+'%')
    #sys.stdout.flush()
    #sys.stdout.write('\n')
    return np.array(amp4_pre),np.array(phase_pre),np.array(amp0_pre)

def amp2_phase2_multilayer(freqin,no_arr,ne_arr,thickness_arr,chi_arr,n1,theta1,a_in,p_in,pol):
    amp0_pre=[]
    amp2_pre=[]
    phase2_pre=[]
    hwp_angle=np.arange(0,361,10)/radeg
    for freq in freqin:
        IVA=IVA_multilayer(freq,hwp_angle,no_arr,ne_arr,thickness_arr,chi_arr,n1,theta1,a_in,p_in,pol)
        IVA_err=np.ones(len(IVA))*np.mean(IVA)*0.01
        paramater=np.zeros(17)
        paramater[0]=np.mean(IVA)
        paramater[3]=0.1*(np.max(IVA)-np.min(IVA))
        paramater[7]=0.4*(np.max(IVA)-np.min(IVA))
        fit_popt,fit_cov=curve_fit(fit_func,hwp_angle,IVA,p0=paramater)
        amp0_pre.append(np.abs(fit_popt[0]))
        amp2_pre.append(np.abs(fit_popt[3]))
        if fit_popt[3]<0.:phase2_pre.append(fit_popt[4]+0.5*pi)
        else:phase2_pre.append(fit_popt[4])
        sys.stdout.write('\r%3d'%(100*len(amp2_pre)/len(freqin))+'%')
        sys.stdout.flush()
    sys.stdout.write('\n')
    return np.array(amp2_pre),np.array(phase2_pre),np.array(amp0_pre)

def amp_phase_multilayer(freqin,no_arr,ne_arr,thickness_arr,chi_arr,n1,theta1,a_in,p_in,pol):
    popt=[]
    cov=[]
    hwp_angle=np.arange(0,361,10)/radeg
    for freq in freqin:
        IVA=IVA_multilayer(freq,hwp_angle,no_arr,ne_arr,thickness_arr,chi_arr,n1,theta1,a_in,p_in,pol)
        IVA_err=np.ones(len(IVA))*np.mean(IVA)*0.01
        paramater=np.zeros(17)
        paramater[0]=np.mean(IVA)
        paramater[3]=0.1*(np.max(IVA)-np.min(IVA))
        paramater[7]=0.4*(np.max(IVA)-np.min(IVA))
        fit_popt,fit_cov=curve_fit(fit_func,hwp_angle,IVA,sigma=IVA_err,p0=paramater,absolute_sigma=True)
        popt.append(fit_popt)
        cov.append(fit_cov)
        #sys.stdout.write('\r%3d'%(100*len(popt)/len(freqin))+'%')
        #sys.stdout.flush()
    #sys.stdout.write('\n')
    return popt,cov

def transfer_matrix_1plate_normal(freq,no,ne,thickness,chi,n1):
    sin_c=np.sin(chi)
    cos_c=np.cos(chi)
    
    eu=np.sqrt(ep0/u0)
    
    d_e=2.*pi*freq*ne*thickness/c
    d_o=2.*pi*freq*no*thickness/c
    
    cos_de=np.cos(d_e)
    sin_de=im*np.sin(d_e)
    cos_do=np.cos(d_o)
    sin_do=im*np.sin(d_o)
    
    t=np.zeros((4,4),'complex')
    t[0][0]=sin_c**2*cos_do+cos_c**2*cos_de
    t[0][1]=(sin_c**2*ne*sin_do+cos_c**2*no*sin_de)/(no*ne*eu)
    t[0][2]=cos_c*sin_c*(cos_de-cos_do)
    t[0][3]=cos_c*sin_c*(no*sin_de-ne*sin_do)/(no*ne*eu)
    
    t[1][0]=(sin_c**2*no*sin_do+cos_c**2*ne*sin_de)*eu
    t[1][1]=sin_c**2*cos_do+cos_c**2*cos_de
    t[1][2]=cos_c*sin_c*eu*(ne*sin_de-no*sin_do)
    t[1][3]=cos_c*sin_c*(cos_de-cos_do)
    
    t[2][0]=cos_c*sin_c*(cos_de-cos_do)
    t[2][1]=cos_c*sin_c*(no*sin_de-ne*sin_do)/(ne*no*eu)
    t[2][2]=cos_c**2*cos_do+sin_c**2*cos_de
    t[2][3]=(cos_c**2*ne*sin_do+sin_c**2*no*sin_de)/(ne*no*eu)
    
    t[3][0]=cos_c*sin_c*eu*(ne*sin_de-no*sin_do)
    t[3][1]=cos_c*sin_c*(cos_de-cos_do)
    t[3][2]=(cos_c**2*no*sin_do+sin_c**2*ne*sin_de)*eu
    t[3][3]=sin_c**2*cos_de+cos_c**2*cos_do
    
    return t

def transfer_matrix_multilayer_normal(freq,no_arr,ne_arr,thickness_arr,chi_arr,n1):
    tm=transfer_matrix_1plate_normal(freq,no_arr[0],ne_arr[0],thickness_arr[0],chi_arr[0],n1)
    for i in range(1,len(thickness_arr)):
        tm=np.dot(transfer_matrix_1plate_normal(freq,no_arr[i],ne_arr[i],thickness_arr[i],chi_arr[i],n1),tm)
    return tm

def Jones_matrix_multilayer_normal(freq,no_arr,ne_arr,thickness_arr,chi_arr,n1):
    t=transfer_matrix_multilayer_normal(freq,no_arr,ne_arr,thickness_arr,chi_arr,n1)
    n1=np.real(n1)
    eu=np.sqrt(ep0/u0)
    a0=t[0][0]+t[0][1]*n1*eu
    a1=t[0][2]+t[0][3]*n1*eu
    a2=t[1][0]/(n1*eu)+t[1][1]
    a3=t[1][2]/(n1*eu)+t[1][3]
    a4=t[2][0]+t[2][1]*n1*eu
    a5=t[2][2]+t[2][3]*n1*eu
    a6=t[3][0]/(n1*eu)+t[3][1]
    a7=t[3][2]/(n1*eu)+t[3][3]
    a8=1./((a0+a2)*(a5+a7)-(a1+a3)*(a4+a6))
    
    #print freq,chi_arr*radeg,a0,a1,a2,a3,a4,a5,a6,a7,a8
    
    j11=2.*a8*(a5+a7)
    j21=2.*a8*(-a4-a6)
    j12=2.*a8*(-a1-a3)
    j22=2.*a8*(a0+a2)
    return np.array([[j11,j12],[j21,j22]])

def Jones_to_0f2f4f(freq,jones_arr,ain):
    amp4=[]
    phase4=[]
    amp2=[]
    phase2=[]
    amp0=[]
    cos_ain=np.cos(2*ain)
    sin_ain=np.sin(2*ain)
    for i in xrange(len(freq)):
        j=np.array(jones_arr[i])
        j_inv=np.conjugate(j)
        coeff_cos4=(cos_ain*(j[1][1]*j_inv[1][1]-j[0][0]*j_inv[1][1]-j[1][0]*j_inv[1][0]-j[0][1]*j_inv[1][0]-j[1][0]*j_inv[0][1]-j[0][1]*j_inv[0][1]-j[1][1]*j_inv[0][0]+j[0][0]*j_inv[0][0])+sin_ain*(j[1][0]*j_inv[1][1]+j[0][1]*j_inv[1][1]+j[1][1]*j_inv[1][0]-j[0][0]*j_inv[1][0]+j[1][1]*j_inv[0][1]-j[0][0]*j_inv[0][1]-j[1][0]*j_inv[0][0]-j[0][1]*j_inv[0][0]))/8
        
        coeff_sin4=(cos_ain*(-j[1][0]*j_inv[1][1]-j[0][1]*j_inv[1][1]-j[1][1]*j_inv[1][0]+j[0][0]*j_inv[1][0]-j[1][1]*j_inv[0][1]+j[0][0]*j_inv[0][1]+j[1][0]*j_inv[0][0]+j[0][1]*j_inv[0][0])+sin_ain*(j[1][1]*j_inv[1][1]-j[0][0]*j_inv[1][1]-j[1][0]*j_inv[1][0]-j[0][1]*j_inv[1][0]-j[1][0]*j_inv[0][1]-j[0][1]*j_inv[0][1]-j[1][1]*j_inv[0][0]+j[0][0]*j_inv[0][0]))/8
        
        coeff_cos2=((-j[1][1]*j_inv[1][1]+j[1][0]*j_inv[1][0]-j[0][1]*j_inv[0][1]+j[0][0]*j_inv[0][0])+cos_ain*(-j[1][1]*j_inv[1][1]-j[1][0]*j_inv[1][0]+j[0][1]*j_inv[0][1]+j[0][0]*j_inv[0][0])+sin_ain*(-j[0][1]*j_inv[1][1]-j[0][0]*j_inv[1][0]-j[1][1]*j_inv[0][1]-j[1][0]*j_inv[0][0]))/4
        
        coeff_sin2=((j[1][0]*j_inv[1][1]+j[1][1]*j_inv[1][0]+j[0][0]*j_inv[0][1]+j[0][1]*j_inv[0][0])+cos_ain*(j[0][1]*j_inv[1][1]+j[0][0]*j_inv[1][0]+j[1][1]*j_inv[0][1]+j[1][0]*j_inv[0][0])+sin_ain*(-j[1][1]*j_inv[1][1]-j[1][0]*j_inv[1][0]+j[0][1]*j_inv[0][1]+j[0][0]*j_inv[0][0]))/4
        
        const=(2*(j[1][1]*j_inv[1][1]+j[1][0]*j_inv[1][0]+j[0][1]*j_inv[0][1]+j[0][0]*j_inv[0][0])+cos_ain*(j[1][1]*j_inv[1][1]+j[0][0]*j_inv[1][1]-j[1][0]*j_inv[1][0]+j[0][1]*j_inv[1][0]+j[1][0]*j_inv[0][1]-j[0][1]*j_inv[0][1]+j[1][1]*j_inv[0][0]+j[0][0]*j_inv[0][0])+sin_ain*(-j[1][0]*j_inv[1][1]+j[0][1]*j_inv[1][1]-j[1][1]*j_inv[1][0]-j[0][0]*j_inv[1][0]+j[1][1]*j_inv[0][1]+j[0][0]*j_inv[0][1]-j[1][0]*j_inv[0][0]+j[0][1]*j_inv[0][0]))/8
        a4=np.sqrt(np.abs(coeff_cos4)**2+np.abs(coeff_sin4)**2)
        a2=np.sqrt(np.abs(coeff_cos2)**2+np.abs(coeff_sin2)**2)
        amp4.append(np.real(a4))
        if coeff_cos4<0:phase4.append(pi/4+np.arctan(np.real(coeff_sin4)/np.real(coeff_cos4))/4)
        else:phase4.append(np.arctan(np.real(coeff_sin4)/np.real(coeff_cos4))/4)
        amp2.append(np.real(a2))
        if coeff_cos2<0:phase2.append(pi/2+np.arctan(np.real(coeff_sin2)/np.real(coeff_cos2))/2)
        else:phase2.append(np.arctan(np.real(coeff_sin2)/np.real(coeff_cos2))/2)
        amp0.append(np.real(const))
    return [np.array(amp4),np.array(phase4),np.array(amp2),np.array(phase2),np.array(amp0)]



def fastcal_2f4f_normal(freq,no_arr,ne_arr,thickness_arr,chi_arr,ain,n1):
    amp4=[]
    phase4=[]
    amp2=[]
    phase2=[]
    amp0=[]
    cos_ain=np.cos(2*ain)
    sin_ain=np.sin(2*ain)
    for f in freq:
        j=Jones_matrix_multilayer_normal(f,no_arr,ne_arr,thickness_arr,np.array(chi_arr)-ain,n1)
        j_inv=np.conjugate(j)
        coeff_cos4=(cos_ain*(j[1][1]*j_inv[1][1]-j[0][0]*j_inv[1][1]-j[1][0]*j_inv[1][0]-j[0][1]*j_inv[1][0]-j[1][0]*j_inv[0][1]-j[0][1]*j_inv[0][1]-j[1][1]*j_inv[0][0]+j[0][0]*j_inv[0][0])+sin_ain*(j[1][0]*j_inv[1][1]+j[0][1]*j_inv[1][1]+j[1][1]*j_inv[1][0]-j[0][0]*j_inv[1][0]+j[1][1]*j_inv[0][1]-j[0][0]*j_inv[0][1]-j[1][0]*j_inv[0][0]-j[0][1]*j_inv[0][0]))/8.
        
        coeff_sin4=(cos_ain*(-j[1][0]*j_inv[1][1]-j[0][1]*j_inv[1][1]-j[1][1]*j_inv[1][0]+j[0][0]*j_inv[1][0]-j[1][1]*j_inv[0][1]+j[0][0]*j_inv[0][1]+j[1][0]*j_inv[0][0]+j[0][1]*j_inv[0][0])+sin_ain*(j[1][1]*j_inv[1][1]-j[0][0]*j_inv[1][1]-j[1][0]*j_inv[1][0]-j[0][1]*j_inv[1][0]-j[1][0]*j_inv[0][1]-j[0][1]*j_inv[0][1]-j[1][1]*j_inv[0][0]+j[0][0]*j_inv[0][0]))/8.
        
        coeff_cos2=((-j[1][1]*j_inv[1][1]+j[1][0]*j_inv[1][0]-j[0][1]*j_inv[0][1]+j[0][0]*j_inv[0][0])+cos_ain*(-j[1][1]*j_inv[1][1]-j[1][0]*j_inv[1][0]+j[0][1]*j_inv[0][1]+j[0][0]*j_inv[0][0])+sin_ain*(-j[0][1]*j_inv[1][1]-j[0][0]*j_inv[1][0]-j[1][1]*j_inv[0][1]-j[1][0]*j_inv[0][0]))/4.
        
        coeff_sin2=((j[1][0]*j_inv[1][1]+j[1][1]*j_inv[1][0]+j[0][0]*j_inv[0][1]+j[0][1]*j_inv[0][0])+cos_ain*(j[0][1]*j_inv[1][1]+j[0][0]*j_inv[1][0]+j[1][1]*j_inv[0][1]+j[1][0]*j_inv[0][0])+sin_ain*(-j[1][1]*j_inv[1][1]-j[1][0]*j_inv[1][0]+j[0][1]*j_inv[0][1]+j[0][0]*j_inv[0][0]))/4.
        
        const=(2*(j[1][1]*j_inv[1][1]+j[1][0]*j_inv[1][0]+j[0][1]*j_inv[0][1]+j[0][0]*j_inv[0][0])+cos_ain*(j[1][1]*j_inv[1][1]+j[0][0]*j_inv[1][1]-j[1][0]*j_inv[1][0]+j[0][1]*j_inv[1][0]+j[1][0]*j_inv[0][1]-j[0][1]*j_inv[0][1]+j[1][1]*j_inv[0][0]+j[0][0]*j_inv[0][0])+sin_ain*(-j[1][0]*j_inv[1][1]+j[0][1]*j_inv[1][1]-j[1][1]*j_inv[1][0]-j[0][0]*j_inv[1][0]+j[1][1]*j_inv[0][1]+j[0][0]*j_inv[0][1]-j[1][0]*j_inv[0][0]+j[0][1]*j_inv[0][0]))/8.
        a4=np.sqrt(np.abs(coeff_cos4)**2+np.abs(coeff_sin4)**2)
        a2=np.sqrt(np.abs(coeff_cos2)**2+np.abs(coeff_sin2)**2)
        amp4.append(np.real(a4))
        if coeff_cos4<0:phase4.append(pi/4+np.arctan(np.real(coeff_sin4)/np.real(coeff_cos4))/4)
        else:phase4.append(np.arctan(np.real(coeff_sin4)/np.real(coeff_cos4))/4)
        amp2.append(np.real(a2))
        if coeff_cos2<0:phase2.append(pi/2+np.arctan(np.real(coeff_sin2)/np.real(coeff_cos2))/2)
        else:phase2.append(np.arctan(np.real(coeff_sin2)/np.real(coeff_cos2))/2)
        amp0.append(np.real(const))
    return np.array(amp4),np.array(phase4),np.array(amp2),np.array(phase2),np.array(amp0)

def Mueller_matrix_multilayer_normal(freq,no_arr,ne_arr,thickness_arr,chi_arr,ain,n1):
    j=Jones_matrix_multilayer_normal(freq,no_arr,ne_arr,thickness_arr,np.array(chi_arr)-ain,n1)
    return Jones_to_Mueller(j)
'''
def phase_adjust_nf(p,n):#rad
    p=np.array(p)
    while(p[0]<0.):p[0]=p[0]+2.*pi/n
    while(p[0]>pi/n):p[0]=p[0]-2.*pi/n
    for i in range(0,len(p)):
        while((p[i]-p[0])<-pi/n):p[i]=p[i]+2.*pi/n
        while((p[i]-p[0])>0.9*pi/n):p[i]=p[i]-2.*pi/n
    return np.array(p)
'''
def phase_adjust_nf(p,n):#rad
    p=np.array(p)
    priod=2.*pi/float(n)
    for i in range(1,len(p)):
        while(p[i]-p[i-1]>0.8*priod):p[i]=p[i]-priod
        while(p[i]-p[i-1]<-0.8*priod):p[i]=p[i]+priod
    return np.array(p)

def phase_adjust_nf_w_amp_minus(p,n):#rad
    p=np.array(p)
    priod=pi/float(n)
    for i in range(1,len(p)):
        while(p[i]-p[i-1]>priod):
            #print p[i]*radeg, p[i-1]*radeg, (p[i]-p[i-1])*radeg
            p[i]=p[i]-priod
        while(p[i]-p[i-1]<-priod):p[i]=p[i]+priod
    return np.array(p)

#https://wiki.kek.jp/pages/viewpage.action?pageId=108104267
def litebird_lft_old_band(freq,data):
    band_width=[[34.e9,46.e9],[42.5e9,57.5e9],[53.1e9,66.9e9],[60.2e9,75.8e9],[69.e9,87.e9],[78.8e9,99.2e9],[88.5e9,111.5e9],[101.2e9,136.9e9],[119.e9,161.e9],[141.1e9,190.9e9],[165.8e9,224.3e9],[199.8e9,270.3e9]]
    band_data=[[],[],[],[],[],[],[],[],[],[],[],[]]
    band_freq=[[],[],[],[],[],[],[],[],[],[],[],[]]
    band_center=[40,50,60,69,78,89,100,119,140,166,195,235]
    freq=np.array(freq)
    data=np.array(data)
    for i in range(0,len(freq)):
        if freq[i] >=band_width[0][0] and freq[i] <=band_width[0][1]:
            band_data[0].append(data[i])
            band_freq[0].append(freq[i])
        
        if freq[i] >=band_width[1][0] and freq[i] <=band_width[1][1]:
            band_data[1].append(data[i])
            band_freq[1].append(freq[i])
        
        if freq[i] >=band_width[2][0] and freq[i] <=band_width[2][1]:
            band_data[2].append(data[i])
            band_freq[2].append(freq[i])
        
        if freq[i] >=band_width[3][0] and freq[i] <=band_width[3][1]:
            band_data[3].append(data[i])
            band_freq[3].append(freq[i])
        
        if freq[i] >=band_width[4][0] and freq[i] <=band_width[4][1]:
            band_data[4].append(data[i])
            band_freq[4].append(freq[i])
        
        if freq[i] >=band_width[5][0] and freq[i] <=band_width[5][1]:
            band_data[5].append(data[i])
            band_freq[5].append(freq[i])
        
        if freq[i] >=band_width[6][0] and freq[i] <=band_width[6][1]:
            band_data[6].append(data[i])
            band_freq[6].append(freq[i])
        
        if freq[i] >=band_width[7][0] and freq[i] <=band_width[7][1]:
            band_data[7].append(data[i])
            band_freq[7].append(freq[i])
        
        if freq[i] >=band_width[8][0] and freq[i] <=band_width[8][1]:
            band_data[8].append(data[i])
            band_freq[8].append(freq[i])
        
        if freq[i] >=band_width[9][0] and freq[i] <=band_width[9][1]:
            band_data[9].append(data[i])
            band_freq[9].append(freq[i])
        
        if freq[i] >=band_width[10][0] and freq[i] <=band_width[10][1]:
            band_data[10].append(data[i])
            band_freq[10].append(freq[i])
        
        if freq[i] >=band_width[11][0] and freq[i] <=band_width[11][1]:
            band_data[11].append(data[i])
            band_freq[11].append(freq[i])

    return np.array(band_center),np.array(band_freq),np.array(band_data),np.array(band_width)

def litebird_lft_band(freq,data):
    band_center=[40.,50.,60.,69.,78.,89.,100.,119.,140.]
    band_center=np.array(band_center)*1.e9
    band_frac=[0.3,0.3,0.23,0.23,0.23,0.23,0.23,0.3,0.3]
    band_frac=np.array(band_frac)*0.5
    band_width=[]
    band_freq=[]
    band_data=[]
    for i in range(len(band_center)):
        band_width.append([band_center[i]*(1.-band_frac[i]),band_center[i]*(1.+band_frac[i])])
        band_freq.append([])
        band_data.append([])
    freq=np.array(freq)
    data=np.array(data)
    for i in range(len(freq)):
        for j in range(len(band_center)):
            if freq[i] >=band_width[j][0] and freq[i] <=band_width[j][1]:
                band_data[j].append(data[i])
                band_freq[j].append(freq[i])
    return np.array(band_center),np.array(band_freq),np.array(band_data),np.array(band_width)

def litebird_band(freq,data):
    band_center=[40.,50.,60.,69.,78.,89.,100.,119.,140.,166.,195.,235.,280.,337.,402.]
    band_center=np.array(band_center)*1.e9
    band_frac=[0.3,0.3,0.23,0.23,0.23,0.23,0.23,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.23]
    band_frac=np.array(band_frac)*0.5
    band_width=[]
    band_freq=[]
    band_data=[]
    for i in range(len(band_center)):
        band_width.append([band_center[i]*(1.-band_frac[i]),band_center[i]*(1.+band_frac[i])])
        band_freq.append([])
        band_data.append([])
    freq=np.array(freq)
    data=np.array(data)
    for i in range(len(freq)):
        for j in range(len(band_center)):
            if freq[i] >=band_width[j][0] and freq[i] <=band_width[j][1]:
                band_data[j].append(data[i])
                band_freq[j].append(freq[i])            
    return np.array(band_center),np.array(band_freq),np.array(band_data),np.array(band_width)


def DLPF(time,data,err,freq_cut,rotation_time):
    
    n=int(np.max(time)/rotation_time)
    
    num_cut=np.where(np.array(time)<=n*rotation_time)[0][-1]
    
    data_fft=np.copy(data[:num_cut])
    err_fft=np.copy(err[:num_cut])
    time_fft=np.copy(time[:num_cut])
    
    fft=np.fft.fft(data_fft)
    fft_freq=np.fft.fftfreq(len(data_fft),time[1]-time[0])
    
    fft_filter=np.copy(fft*2.)
    fft_filter[(fft_freq==0)]=0.5*fft_filter[(fft_freq==0)]
    fft_filter[(fft_freq>freq_cut)]=0.
    fft_filter[(fft_freq<0.)]=0.
    
    ifft=np.real(np.fft.ifft(fft_filter))
    
    return [time_fft,ifft,err_fft]


def unit():
    u=np.zeros((4,4),'complex')
    u[0,0] = complex(1., 0.)
    u[1,1] = complex(1., 0.)
    u[2,2] = complex(1., 0.)
    u[3,3] = complex(1., 0.)
    return u

def Mueller_rot(theta):
    mr = np.zeros((4,4),'complex')
    mr[0,0]=1.
    mr[1,1]= np.cos(2.*theta)
    mr[1,2]= -np.sin(2.*theta)
    mr[2,1]= np.sin(2.*theta)
    mr[2,2]= np.cos(2.*theta)
    mr[3,3]=1.
    return mr

def Mueller_shwp_wo_refl(freq,no,ne,thickness):
    d=2.*pi*freq*thickness*np.abs(ne-no)/c
    mr = np.zeros((4,4),'complex')
    mr[0,0]=1.
    mr[1,1]=1.
    mr[2,2]=np.cos(d)
    mr[2,3]=-np.sin(d)
    mr[3,2]=np.sin(d)
    mr[3,3]=np.cos(d)
    return mr
        
def fit_4f(rho,a0,a4,phi4):
    return a0+a4*np.cos(4*(rho+phi4))

def poleff_phase_wo_refl(freqin,no_arr,ne_arr,thickness_arr,angle_arr,a_in,p_in):
    amp4=[]
    phase=[]
    amp0=[]
    for freq in freqin:
        m=unit()
        for i in range(0,len(thickness_arr)):
            m=np.dot(Mueller_rot(angle_arr[i]),m)
            m=np.dot(Mueller_shwp_wo_refl(freq,no_arr[i],ne_arr[i],thickness_arr[i]),m)
            m=np.dot(Mueller_rot(-angle_arr[i]),m)
        
        hwp_angle=np.arange(0,361,10)/radeg
        IVA=[]
        Sin = np.array([1.,p_in*np.cos(2.*a_in),p_in*np.sin(2.*a_in),0.])
        for a in hwp_angle:
            m_rot=unit()
            m_rot = np.dot(Mueller_rot(a),m_rot)
            m_rot = np.dot(m,m_rot)
            m_rot = np.dot(Mueller_rot(-a),m_rot)
            m_rot = np.dot(xgrid(),m_rot)
            Sout = np.dot(m_rot,Sin)
            IVA.append(np.abs(Sout[0]))
        init=[1.,1.,0]
        popt,cov=curve_fit(fit_4f,np.array(hwp_angle),np.array(IVA),p0=init,maxfev=10000)
        #print popt[0]
        amp0.append(np.abs(popt[0]))
        amp4.append(np.abs(popt[1]))
        if popt[1]>0:phase.append(popt[2])
        else:phase.append(popt[2]+0.25*pi)
    return np.array(amp4),np.array(phase),np.array(amp0)


def shwp_wo_refl_w_angle(freq,no,ne,thickness,angle):
    d=2.*pi*freq*thickness*np.abs(ne-no)/c
    mr = np.zeros((4,4))
    mr[0,0]=1.
    mr[1,1]=np.cos(d)*np.sin(2.*angle)**2+np.cos(2.*angle)**2
    mr[1,2]=(np.cos(d)-1.)*np.cos(2.*angle)*np.sin(2.*angle)
    mr[1,3]=-np.sin(d)*np.sin(2.*angle)
    mr[2,1]=(np.cos(d)-1.)*np.cos(2.*angle)*np.sin(2.*angle)
    mr[2,2]=np.cos(d)*np.cos(2.*angle)**2+np.sin(2.*angle)**2
    mr[2,3]=-np.sin(d)*np.cos(2.*angle)
    mr[3,1]=np.sin(d)*np.sin(2.*angle)
    mr[3,2]=np.sin(d)*np.cos(2.*angle)
    mr[3,3]=np.cos(d)
    return mr

def fastcal_4f_normal_wo_refl(freq,no_arr,ne_arr,thickness_arr,angle_arr,a_in,p_in):
    amp0=[]
    amp4=[]
    phase4=[]
    Sin = np.array([1.,p_in*np.cos(2.*a_in),p_in*np.sin(2.*a_in),0.])
    for i in range(0,len(freq)):
        m=unit()
        for j in range(0,len(angle_arr)):
            m=np.dot(shwp_wo_refl_w_angle(freq[i],no_arr[j],ne_arr[j],thickness_arr[j],angle_arr[j]),m)
        coeff_cos4=0.25*(Sin[1]*(m[1][1]-m[2][2])+Sin[2]*(m[2][1]+m[1][2]))
        coeff_sin4=0.25*(Sin[1]*(m[2][1]+m[1][2])+Sin[2]*(m[2][2]-m[1][1]))
        const=0.25*(2.*Sin[0]*m[0][0]+Sin[1]*(m[1][1]+m[2][2])+Sin[2]*(m[1][2]-m[2][1]))
        a4=np.sqrt(np.abs(coeff_cos4)**2+np.abs(coeff_sin4)**2)
        amp4.append(np.real(a4))
        amp0.append(np.real(const))
        if coeff_cos4<0:phase4.append(-(pi/4.+np.arctan(coeff_sin4/coeff_cos4)/4.))
        else:phase4.append(-np.arctan(coeff_sin4/coeff_cos4)/4.)
    return np.array(amp4),np.array(phase4),np.array(amp0)

def fastcal_4f_normal_wo_refl_2(freq,no_arr,ne_arr,thickness_arr,angle_arr,a_in,p_in):
    amp0=[]
    amp4=[]
    phase4=[]
    for i in range(0,len(freq)):
       m=unit()
       for j in range(0,len(angle_arr)):
           m=np.dot(shwp_wo_refl_w_angle(freq[i],no_arr[j],ne_arr[j],thickness_arr[j],angle_arr[j]),m)
       a4=0.25*np.sqrt((m[1][1]-m[2][2])**2+(m[1][2]+m[2][1])**2)
       p4=0.25*np.arctan((m[1][2]+m[2][1])/(m[1][1]-m[2][2]))
       amp4.append(np.real(a4))
       if a4<0:phase4.append(pi/4+p4)
       else:phase4.append(p4)
       amp0.append(np.NaN)
    return np.array(amp4),np.array(phase4),np.array(amp0)

def band_ave_poleff_tophat(poleff,phase):
    poleff_sum=poleff[0]
    phase_sum=phase[0]
    for i in range(1,len(poleff)):
        coeff_cos4=poleff_sum*np.cos(phase_sum*4.)+poleff[i]*np.cos(phase[i]*4.)
        coeff_sin4=-poleff_sum*np.sin(phase_sum*4.)-poleff[i]*np.sin(phase[i]*4.)
        poleff_sum=np.sqrt(coeff_cos4**2+coeff_sin4**2)
        if coeff_cos4<0:phase_sum= -(pi/4.+np.arctan(coeff_sin4/coeff_cos4)/4.)
        else:phase_sum= -np.arctan(coeff_sin4/coeff_cos4)/4.
    return np.real(poleff_sum/len(poleff)),np.real(phase_sum),np.real(np.max(phase)-np.min(phase))

def band_ave_poleff(poleff,phase, weight_func):
    poleff_sum=weight_func[0]*poleff[0]
    phase_sum=phase[0]
    for i in range(1,len(poleff)):
        coeff_cos4=poleff_sum*np.cos(phase_sum*4.)+weight_func[i]*poleff[i]*np.cos(phase[i]*4.)
        coeff_sin4=-poleff_sum*np.sin(phase_sum*4.)-weight_func[i]*poleff[i]*np.sin(phase[i]*4.)
        poleff_sum=np.sqrt(coeff_cos4**2+coeff_sin4**2)
        if coeff_cos4<0:phase_sum= -(pi/4.+np.arctan(coeff_sin4/coeff_cos4)/4.)
        else:phase_sum= -np.arctan(coeff_sin4/coeff_cos4)/4.
    return np.real(poleff_sum),np.real(phase_sum),np.real(np.max(phase)-np.min(phase))

def sum_cos4f(amp,phase, weight_func):
    amp_sum=weight_func[0]*amp[0]
    phase_sum=phase[0]
    for i in range(1,len(amp)):
        coeff_cos4=amp_sum*np.cos(phase_sum*4.)+weight_func[i]*amp[i]*np.cos(phase[i]*4.)
        coeff_sin4=-amp_sum*np.sin(phase_sum*4.)-weight_func[i]*amp[i]*np.sin(phase[i]*4.)
        amp_sum=np.sqrt(coeff_cos4**2+coeff_sin4**2)
        if coeff_cos4<0:phase_sum= -(pi/4.+np.arctan(coeff_sin4/coeff_cos4)/4.)
        else:phase_sum= -np.arctan(coeff_sin4/coeff_cos4)/4.
    return np.real(amp_sum),np.real(phase_sum),np.real(np.max(phase)-np.min(phase))
    
def sum_cos2f(amp,phase, weight_func):
    amp_sum=weight_func[0]*amp[0]
    phase_sum=phase[0]
    for i in range(1,len(amp)):
        coeff_cos2=amp_sum*np.cos(phase_sum*2.)+weight_func[i]*amp[i]*np.cos(phase[i]*2.)
        coeff_sin2=-amp_sum*np.sin(phase_sum*2.)-weight_func[i]*amp[i]*np.sin(phase[i]*2.)
        amp_sum=np.sqrt(coeff_cos2**2+coeff_sin2**2)
        if coeff_cos2<0:phase_sum= -(pi/2.+np.arctan(coeff_sin2/coeff_cos2)/2.)
        else:phase_sum= -np.arctan(coeff_sin2/coeff_cos2)/2.
    return np.real(amp_sum),np.real(phase_sum),np.real(np.max(phase)-np.min(phase))

def sum_trig(amp, phase):
    count=0
    if np.max(amp)>0.:
        for i in range(1,len(amp)):
            if amp[i]>0:
                if count==0:
                    amp_sum=amp[i]
                    phase_sum=phase[i]
                    count=1
                else:
                    coeff_cos=amp_sum*np.cos(phase_sum)+amp[i]*np.cos(phase[i])
                    coeff_sin=amp_sum*np.sin(phase_sum)+amp[i]*np.sin(phase[i])
                    amp_sum=np.sqrt(coeff_cos**2+coeff_sin**2)
                    if coeff_cos<0.:phase_sum=np.arctan(coeff_sin/coeff_cos)-pi
                    else:phase_sum=np.arctan(coeff_sin/coeff_cos)
    else:pass
    return np.real(amp_sum),np.real(phase_sum)


def fastcal_2f4f_normal_2(freq,no_arr,ne_arr,thickness_arr,chi_arr,ain,n1):
    amp4=[]
    phase4=[]
    amp2=[]
    phase2=[]
    amp0=[]
    for f in freq:
        m=Mueller_matrix_multilayer_normal(f,no_arr,ne_arr,thickness_arr,chi_arr,ain,n1)
        a4=0.25*np.sqrt((m[1][1]-m[2][2])**2+(m[1][2]+m[2][1])**2)
        a2_p=0.5*np.sqrt(m[0][1]**2+m[0][2]**2)
        a2_i=0.5*np.sqrt(m[1][0]**2+m[2][0]**2)
        p2_p=0.5*np.arctan(m[0][2]/m[0][1])
        p2_i=0.5*np.arctan(m[2][0]/m[1][0])
        p4=0.25*np.arctan((m[1][2]+m[2][1])/(m[1][1]-m[2][2]))
        a2p2=sum_cos2f([a2_i,a2_p],[p2_i,p2_p],[1,1])
        a2=a2p2[0]
        p2=a2p2[1]
        amp4.append(np.real(a4))
        if a4<0:phase4.append(pi/4+p4)
        else:phase4.append(p4)
        amp2.append(np.real(a2))
        if a2<0:phase2.append(pi/2+p2)
        else:phase2.append(p2)
        amp0.append(np.NaN)
    return np.array(amp4),np.array(phase4),np.array(amp2),np.array(phase2),np.array(amp0)
'''
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def unit():
    u=np.zeros((4,4),'complex')
    u[0,0] = complex(1., 0.)
    u[1,1] = complex(1., 0.)
    u[2,2] = complex(1., 0.)
    u[3,3] = complex(1., 0.)
    return u
    
def Mueller_rot(theta):
    mr = np.zeros((4,4),'complex')
    mr[0,0]=1.
    mr[1,1]= np.cos(2.*theta)
    mr[1,2]= -np.sin(2.*theta)
    mr[2,1]= np.sin(2.*theta)
    mr[2,2]= np.cos(2.*theta)
    mr[3,3]=1.
    return mr


def pauli_matrix(n):
    pm=np.zeros((2,2),'complex')
    if n==0:
    pm[0][0]=pm[1][1]=1.

    if n==1:
    pm[0][0]=1.
    pm[1][1]=-1.

    if n==2:
    pm[0][1]=pm[1][0]=1.

    if n==3:
    pm[0][1]=-im
    pm[1][0]=im
    return pm


def delta(no,ne,freq,thickness,chi,n1,theta1):
    k0=2.*pi*freq/c
    neff=n_eff(no,ne,chi,n1,theta1)
    cos_no=cos_refraction(no,n1,theta1)
    cos_ne=cos_refraction(neff,n1,theta1)
    delta_o=k0*no*thickness*cos_no
    delta_e=k0*neff*thickness*cos_ne
    d=np.zeros((4,4),'complex')
    d[0][0]=expi(delta_o)
    d[1][1]=expi(delta_e)
    d[2][2]=expi(-delta_o)
    d[3][3]=expi(-delta_e)
    return d

def rot_xy(chi):
    r_xy=np.zeros((3,3),'complex')
    r_xy[0][0]=np.cos(chi)
    r_xy[0][1]=-np.sin(chi)
    r_xy[1][0]=np.sin(chi)
    r_xy[1][1]=np.cos(chi)
    r_xy[2][2]=1.
    return r_xy

def Dielectric_matrix(no,ne,chi):
    ne=np.real(ne)
    no=np.real(no)
    dm=np.zeros((3,3),'complex')
    dm[0][0]=ep0*ne**2
    dm[1][1]=ep0*no**2
    dm[2][2]=ep0*no**2

    dm=np.dot(rot_xy(chi),dm)
    dm=np.dot(dm,rot_xy(-chi))
    return dm

def DHtoEH(no,ne,chi):
    dm_inv=np.linalg.inv(Dielectric_matrix(no,ne,chi))
    de=np.zeros((4,4),'complex')
    de[0][0]=dm_inv[0][0]
    de[0][2]=dm_inv[1][0]
    de[1][1]=1.
    de[2][0]=dm_inv[0][1]
    de[2][2]=dm_inv[1][1]
    de[3][3]=1.
    return de

def EHtoDH(no,ne,chi):
    return np.linalg.inv(DHtoEH(no,ne,chi))

def DDtoDH(no,ne,chi,n1,theta1):
    neff=n_eff(no,ne,chi,n1,theta1)
    cos_no=cos_refraction(no,n1,theta1)
    cos_ne=cos_refraction(neff,n1,theta1)
    sin_no=sin_refraction(no,n1,theta1)
    sin_ne=sin_refraction(neff,n1,theta1)

    D_no=np.array([-np.sin(chi)*cos_no,np.cos(chi)*cos_no,np.sin(chi)*sin_no])
    D_no=D_no/np.sqrt(np.sum(D_no**2))
    #print np.sqrt(np.sum(D_no**2))

    D_ne=np.array([-np.cos(chi)*cos_no*cos_ne,-np.sin(chi)*(sin_no*sin_ne+cos_no*cos_ne),np.cos(chi)*cos_no*sin_ne])
    D_ne=D_ne/np.sqrt(np.sum(D_ne**2))
    #print np.sqrt(np.sum(D_ne**2))

    H_no=np.array([-np.cos(chi)*cos_no**2,-np.sin(chi),np.cos(chi)*cos_no*sin_no])
    H_no=H_no/np.sqrt(np.sum(H_no**2))
    #print np.sqrt(np.sum(H_no**2))

    H_ne=np.array([(cos_no*cos_ne+sin_no*sin_ne)*cos_ne*np.sin(chi),-cos_no*np.cos(chi),-(cos_no*cos_ne+sin_no*sin_ne)*sin_ne*np.sin(chi)])
    H_ne=H_ne/np.sqrt(np.sum(H_ne**2))
    #print np.sqrt(np.sum(H_ne**2))

    dh=np.zeros((4,4),'complex')
    dh[0][0]=D_no[0]
    dh[0][1]=D_ne[0]
    dh[0][2]=D_no[0]
    dh[0][3]=D_ne[0]
    dh[1][0]=H_no[1]*c/np.real(no)
    dh[1][1]=H_ne[1]*c/np.real(neff)
    dh[1][2]=-H_no[1]*c/np.real(no)
    dh[1][3]=-H_ne[1]*c/np.real(neff)
    dh[2][0]=D_no[1]
    dh[2][1]=D_ne[1]
    dh[2][2]=D_no[1]
    dh[2][3]=D_ne[1]
    dh[3][0]=-H_no[0]*c/np.real(no)
    dh[3][1]=-H_ne[0]*c/np.real(neff)
    dh[3][2]=H_no[0]*c/np.real(no)
    dh[3][3]=H_ne[0]*c/np.real(neff)
    return dh

def DHtoDD(no,ne,chi,n1,theta1):
    return np.linalg.inv(DDtoDH(no,ne,chi,n1,theta1))

def transfer_matrix_1plate(freq,no,ne,thickness,chi,n1,theta1):
    g=DHtoEH(no,ne,chi)
    g=np.dot(g,DDtoDH(no,ne,chi,n1,theta1))
    g=np.dot(g,delta(no,ne,freq,thickness,chi,n1,theta1))
    g=np.dot(g,DHtoDD(no,ne,chi,n1,theta1))
    g=np.dot(g,EHtoDH(no,ne,chi))
    return g

def Mueller_matrix_multilayer(freq,no_arr,ne_arr,thickness_arr,chi_arr,n1,theta1):
    m=np.zeros((4,4),'complex')
    jones=Jones_matrix_multilayer(freq,no_arr,ne_arr,thickness_arr,chi_arr,n1,theta1)
    jones_d=np.conjugate(jones.T)
    for i in range(0,4):
        for j in range(0,4):
            k=np.dot(pauli_matrix(i),jones)
            k=np.dot(k,pauli_matrix(j))
            k=np.dot(k,jones_d)
            m[i][j]=0.5*np.trace(k)
    return m


'''




