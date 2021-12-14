import numpy as np
import matplotlib.pyplot as plt
import itertools
import matplotlib.cm as cm
import sys
import os
import glob
cmap = cm.get_cmap("jet")
plt.rcParams.update({'font.size': 12})
from lib_KK import litebird_lft_band #(freq,data)
from lib_KK import band_ave_poleff_tophat #(poleff,phase)
from lib_KK import fastcal_4f_normal_wo_refl #(freq,no_arr,ne_arr,thickness_arr,angle_arr,a_in,p_in)
from lib_KK import rms
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c=299792458. #m/s
pi=np.pi
e=np.e
radeg = (180./pi)
ep0 = 8.8542e-12 #s4 A2/m3 kg
u0=4.*pi*1e-7 #H/m=V/(A/s)/m
im=complex(0.,1.)
arcmin=1./60.
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
num=100

#sym or asym
patterns=['asym','sym']

layers=[3,5,7,9]
no=3.047
ne=3.361
#no=complex(no,0.5*no*1.e-4)
#ne=complex(ne,0.5*ne*1.e-4)
a_in = 45. #deg
p_in = 1. # 0to1

f_range=[[4.e9,191.e9],[34.e9,161.e9,],[84.e9,111.e9]]
f_step=1.e9

#thickness=c/((f_f+f_i)*np.abs(no-ne))
#thickness=4.896e-3
thickness=4.9e-3

output_dir='output'
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if len(glob.glob(output_dir))==0:os.mkdir(output_dir)

for f_r in f_range:
    for pattern in patterns:
        for layer in layers:
            random_seed=123456789
            np.random.seed(random_seed)
            random_angle_arr=np.random.rand(int((layer-1)*0.5),num)*180.
            
            f_i=f_r[0]
            f_f=f_r[1]
            freq=np.arange(f_i,f_f,f_step)
            no_arr = np.ones(layer)*no
            ne_arr = np.ones(layer)*ne
            min_best=0
            angle_best=np.zeros(layer)
            angle_arr=np.zeros(layer)
            ave_best=np.zeros(layer)
            thickness_arr=np.ones(layer)*thickness
            thickness_best=np.zeros(layer)
            count=0
            ave_best=0

            ave_arr=[]
            for i in range(0,num):
                #while(mean<0.995):
                for j in range(0,int((layer-1)/2)):
                    random_angle=random_angle_arr[j][i]
                    angle_arr[j]=random_angle
                    if pattern=='asym':
                        angle_arr[-j-1]=-random_angle #anti-symmetric
                    elif pattern=='sym':
                        angle_arr[-j-1]=random_angle #symmetric
                    else:print('invalid pattern name')
                angle_arr=np.round(angle_arr,3)
                m=fastcal_4f_normal_wo_refl(freq,no_arr,ne_arr,thickness_arr,angle_arr/radeg,a_in/radeg,p_in)
                modeff_band=2.*m[0]
                phase_band=m[1]
                ave=band_ave_poleff_tophat(modeff_band,phase_band)
                ave_current=ave[0]
                ave_arr.append(ave[0])
                if ave_current>ave_best:
                    ave_best=ave_current
                    angle_best=angle_arr
                    thickness_best=thickness_arr
                angle_arr=np.zeros(layer)
                count=count+1
                sys.stdout.write(str('  %i'%count)+' '+str(np.round(angle_best,1))+'\r')
                #sys.stdout.write(str('  %i'%count)+'/'+str('%i'%num)+' '+str('%.3f'%(ave_best[0]))+' '+str(angle_best)+'\r')
                sys.stdout.flush()
            print(str('  %i'%count)+' '+str(np.round(np.sum(ave_best),2))+' '+str(angle_best))

            freq_pre=np.arange(1.e9,200.e9,0.1e9)
            a_in=0.
            pre=fastcal_4f_normal_wo_refl(freq_pre,no_arr,ne_arr,thickness_best,angle_best/radeg,a_in/radeg,p_in)
            #pre_single=fastcal_4f_normal_wo_refl(freq,no_arr,ne_arr,thickness_best,np.array([157.339, 47.09, 0.,-47.09,-157.339])/radeg,a_in/radeg,p_in)

            #pre_single=fastcal_4f_normal_wo_refl(freq,no_arr,ne_arr,thickness_best,np.array([10., 0.,-10.])/radeg,a_in/radeg,p_in)
            pre_single=fastcal_4f_normal_wo_refl(freq_pre,np.array([no]),np.array([ne]),np.array([thickness]),np.array([0]),a_in/radeg,p_in)
            #pre_af=fastcal_4f_normal_wo_refl(freq,no_arr,ne_arr,thickness_minuit,angle_minuit/radeg,a_in/radeg,p_in)


            if layer==5:
                fig, ax = plt.subplots(figsize=(8.5,7.5))
                mappable=plt.scatter(random_angle_arr[0],random_angle_arr[1],c=ave_arr,s=3,cmap=cmap)
                plt.scatter(angle_best[0],angle_best[1],c="k",s=20, marker='*',label=r'Best ($\chi_{1}^{\rm Best}$, $\chi_{2}^{\rm Best}$)')
                plt.scatter(180-angle_best[0],180-angle_best[1],c="k",s=10, marker='o',label=r'Opposite ($180-\chi_{1}^{\rm Best}$, $180-\chi_{2}^{\rm Best}$)')
                plt.xlabel(r'$\chi_{1}$ [deg.]', fontsize=16)
                plt.ylabel(r'$\chi_{2}$ [deg.]', fontsize=16)
                plt.grid()
                plt.ylim(0,180)
                plt.xlim(0,180)
                fig.colorbar(mappable, ax=ax).set_label(r'Band averaged polarization efficiency $2A_{4}$')
                ax.legend(bbox_to_anchor=(0., 1.02, 1.05, 1.02), loc=3, ncol=2, mode="expand", borderaxespad=0)
                plt.tight_layout()
                plt.savefig(output_dir+'/'+pattern+'_'+str(layer)+'_run'+str(num)+'_'+str(int(f_i*1.e-9))+'_'+str(int(f_f*1.e-9))+'_angle_set_w_opposite.png')
                #plt.show()
                plt.close()

            fig, ax = plt.subplots(nrows=2,ncols=1,figsize=(8.5,8.5))
            plt.axes(ax[0])
            plt.ylim(0.0,1.05)
            plt.xlim(0,200)
            plt.plot([f_i*1.e-9,f_i*1.e-9],[0.,1.2],'m-')
            plt.plot([f_f*1.e-9,f_f*1.e-9],[0.,1.2],'m-')
            plt.plot([f_i*1.e-9,f_f*1.e-9],[0.9,0.9],'m-')
            plt.plot(freq_pre*1.e-9,2.*pre_single[0],'k--',label='Single',linewidth=0.5)
            plt.plot(freq_pre*1.e-9,2.*pre[0],'r-',label='angle '+str(angle_best)+' '+str(thickness_best[0]),linewidth=0.5)
            plt.ylabel('Mod. eff.', fontsize=16)
            plt.legend(fontsize=10)
            plt.xlabel('Frequency GHz', fontsize=16)
            plt.axes(ax[1])
            plt.xlim(0,200)
            plt.plot(freq_pre*1.e-9,pre_single[1]*radeg,'k--',label='Sangle',linewidth=0.5)
            plt.plot(freq_pre*1.e-9,pre[1]*radeg,'r-',label=str(layer)+' layer',linewidth=0.5)
            plt.ylabel('$\phi_{4}$', fontsize=16)
            plt.legend(fontsize=10)
            plt.xlabel('Frequency GHz', fontsize=16)
            plt.tight_layout()
            plt.savefig(output_dir+'/opt_'+pattern+'_'+str(layer)+'_run'+str(num)+'_'+str(int(f_i*1.e-9))+'_'+str(int(f_f*1.e-9))+'_odd.png')
            plt.close()
            print('END')

            np.savez(output_dir+'/opt_'+pattern+'_'+str(layer)+'_run'+str(num)+'_'+str(int(f_i*1.e-9))+'_'+str(int(f_f*1.e-9))+'_odd', \
            layer=layer, num=num, no=no, ne=ne, a_in=a_in, p_in=p_in, \
            freq_opt=freq, poleff=2.*pre[0], phi4=pre[1]*radeg, ave_arr=ave_arr, \
            random_seed=random_seed, random_angle_arr=random_angle_arr, angle_best=angle_best, ave_best=ave_best,thickness=thickness)

#plt.show()
sys.exit()