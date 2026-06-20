import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.animation as animation
import os, time
try:
    from numba import njit
    HAS_NUMBA=True
except ImportError:
    def njit(**kw): return lambda f: f
    HAS_NUMBA=False

GAMMA=5/3; ETA=0.5; JET_VEL=20.0; SIGMA_Z=0.1; SIGMA_PHI=0.8; A_MAG=0.8 #parametros #A-MAG es el radio de la region magnetizada
NR=50; NZ=70; R_MAX=8.0; Z_MAX=30.0; CFL=0.4; T_END=30.0; DT_OUT=999
dr=R_MAX/NR; dz=Z_MAX/NZ
R=(np.arange(NR)+0.5)*dr #rejilla
Z=(np.arange(NZ)+0.5)*dz
Rf=np.arange(NR+1)*dr
iRHO,iVR,iVZ,iVF,iBR,iBZ,iBF,iPRS=range(8); NV=8 #variable primaria: Q=[rho,vR,vz,vphi,BR,Bz,Bphi,prs] numerada para acceso rapido
P_AMB=1.0/GAMMA; Bz0=np.sqrt(2.0*SIGMA_Z*P_AMB) 
Bm=np.sqrt(2.0*SIGMA_PHI*P_AMB/(A_MAG**2*(0.5-2.0*np.log(A_MAG))))

def profile(r, t=0.0):
    # PLUTO usa cosh(R^4) para t<0.1 con campo toroidal, luego cosh(R^8)
    if SIGMA_PHI > 1e-9 and t < 0.1:
        return 1.0/np.cosh(np.minimum(r**4, 500.0))
    return 1.0/np.cosh(np.minimum(r**8, 500.0))

def jet_prim(r):
    bphi=np.where(r<A_MAG,-Bm*r/A_MAG,-Bm*A_MAG/np.maximum(r,1e-12)) #def campo magnetico
    p=P_AMB+Bm**2*(1.0-np.minimum(r**2/A_MAG**2,1.0)) #condicion de balance de presion
    Q=np.zeros((NV,len(r)))
    Q[iRHO]=ETA; Q[iVZ]=JET_VEL; Q[iBZ]=Bz0; Q[iBF]=bphi; Q[iPRS]=p
    return Q

def prim_to_cons(Q): #primitivas a conservativas
    U=Q.copy(); rho=Q[iRHO]
    U[iVR]=rho*Q[iVR]; U[iVZ]=rho*Q[iVZ]; U[iVF]=rho*Q[iVF]
    v2=Q[iVR]**2+Q[iVZ]**2+Q[iVF]**2
    B2=Q[iBR]**2+Q[iBZ]**2+Q[iBF]**2
    U[iPRS]=Q[iPRS]/(GAMMA-1.0)+0.5*rho*v2+0.5*B2
    return U

def cons_to_prim(U): #conservativas a primitivas
    Q=U.copy(); rho=np.maximum(U[iRHO],1e-10)
    Q[iVR]=U[iVR]/rho; Q[iVZ]=U[iVZ]/rho; Q[iVF]=U[iVF]/rho
    v2=Q[iVR]**2+Q[iVZ]**2+Q[iVF]**2
    B2=U[iBR]**2+U[iBZ]**2+U[iBF]**2
    p=(GAMMA-1.0)*(U[iPRS]-0.5*rho*v2-0.5*B2)
    Q[iPRS]=np.maximum(p,1e-10)
    return Q

def ptot(Q): #presion total (gas + magnetica)
    return Q[iPRS]+0.5*(Q[iBR]**2+Q[iBZ]**2+Q[iBF]**2)

def energy(Q): #energia total por unidad de volumen
    v2=Q[iVR]**2+Q[iVZ]**2+Q[iVF]**2
    B2=Q[iBR]**2+Q[iBZ]**2+Q[iBF]**2
    return Q[iPRS]/(GAMMA-1.0)+0.5*Q[iRHO]*v2+0.5*B2

def fast_speed(Q,axis): #velocidad de la onda magnetosonica rapida en direccion axis (0=R, 1=z)
    rho=np.maximum(Q[iRHO],1e-10); p=np.maximum(Q[iPRS],1e-10)
    bN=Q[iBR] if axis==0 else Q[iBZ]
    B2=Q[iBR]**2+Q[iBZ]**2+Q[iBF]**2
    cs2=GAMMA*p/rho; vA2=B2/rho; bN2=bN**2/rho
    disc=np.maximum((cs2+vA2)**2-4.0*cs2*bN2,0.0)
    cf=np.sqrt(0.5*(cs2+vA2+np.sqrt(disc)))
    v=Q[iVR] if axis==0 else Q[iVZ]
    return cf+np.abs(v)

def flux_R(Q):
    rho=Q[iRHO]; vR,vZ,vF=Q[iVR],Q[iVZ],Q[iVF]; bR,bZ,bF=Q[iBR],Q[iBZ],Q[iBF]
    pt=ptot(Q); vdB=vR*bR+vZ*bZ+vF*bF; e=energy(Q)
    F=np.empty_like(Q)
    F[iRHO]=rho*vR; F[iVR]=rho*vR**2+pt-bR**2; F[iVZ]=rho*vR*vZ-bR*bZ
    F[iVF]=rho*vR*vF-bR*bF; F[iBR]=0.0; F[iBZ]=vR*bZ-vZ*bR
    F[iBF]=vR*bF-vF*bR; F[iPRS]=(e+pt)*vR-bR*vdB
    return F

def flux_Z(Q):
    rho=Q[iRHO]; vR,vZ,vF=Q[iVR],Q[iVZ],Q[iVF]; bR,bZ,bF=Q[iBR],Q[iBZ],Q[iBF]
    pt=ptot(Q); vdB=vR*bR+vZ*bZ+vF*bF; e=energy(Q)
    G=np.empty_like(Q)
    G[iRHO]=rho*vZ; G[iVR]=rho*vZ*vR-bZ*bR; G[iVZ]=rho*vZ**2+pt-bZ**2
    G[iVF]=rho*vZ*vF-bZ*bF; G[iBR]=vZ*bR-vR*bZ; G[iBZ]=0.0
    G[iBF]=vZ*bF-vF*bZ; G[iPRS]=(e+pt)*vZ-bZ*vdB
    return G

def upwind_step_R(Q,dt):
    
    QL=Q[:,:-1,:]; QR=Q[:,1:,:]
    UL=prim_to_cons(QL); UR=prim_to_cons(QR)
    S=np.maximum(fast_speed(QL,axis=0),fast_speed(QR,axis=0))
    F_iface=0.5*(flux_R(QL)+flux_R(QR))-0.5*S[np.newaxis]*(UR-UL) #flujo de Rusanov 
    Ri=Rf[1:-1][np.newaxis,:,np.newaxis] #factor geometrico R en el flujo (solo para las interfaces internas)
    FR_full=np.zeros((NV,NR+1,NZ))
    FR_full[:,1:-1,:]=F_iface*Ri #flujo en las interfaces internas
    FR_full[:,0,:]=0.0 #no hay flujo en R=0 por simetria
    FR_full[:,-1,:]=flux_R(Q[:,-1,:])*Rf[-1] #flujo en la frontera externa (outflow)
    Rc=R[np.newaxis,:,np.newaxis] #geometrico
    return -(dt/(Rc*dr))*(FR_full[:,1:,:]-FR_full[:,:-1,:])

def upwind_step_Z(Q,dt):
    # Lax-Friedrichs (Rusanov) en z — primer orden, estable
    QL=Q[:,:,:-1]; QR=Q[:,:,1:]
    UL=prim_to_cons(QL); UR=prim_to_cons(QR)
    S=np.maximum(fast_speed(QL,axis=1),fast_speed(QR,axis=1))
    G_iface=0.5*(flux_Z(QL)+flux_Z(QR))-0.5*S[np.newaxis]*(UR-UL)
    G_full=np.zeros((NV,NR,NZ+1))
    G_full[:,:,1:-1]=G_iface
    G_full[:,:,0]=flux_Z(Q[:,:,0])
    G_full[:,:,-1]=flux_Z(Q[:,:,-1])
    return -(dt/dz)*(G_full[:,:,1:]-G_full[:,:,:-1])

def source_terms(Q,dt): #terminos geometricos cilindricos
    Rc=np.maximum(R,1e-10)[np.newaxis,:,np.newaxis]
    rho=Q[iRHO]; vR,vF=Q[iVR],Q[iVF]; bR,bF=Q[iBR],Q[iBF]
    S=np.zeros_like(Q)
    S[iVR]=dt*(rho*vF**2-bF**2)/Rc
    S[iVF]=-dt*(rho*vR*vF-bR*bF)/Rc
    S[iBF]=dt*(vR*bF-vF*bR)/Rc
    return S

def apply_bc(Q, t=0.0): #boquilla del jet
    Qj   = jet_prim(R)
    prof = profile(R, t)
    # BC z=0: inyeccion del jet mezclada con estado actual
    for v in range(NV):
        Q[v, :, 0] = Q[v, :, 0]*(1-prof) + Qj[v]*prof
    # R=0: simetria axial
    Q[:, 0, :]   = Q[:, 1, :]
    Q[iVR, 0, :] = 0.0
    Q[iBR, 0, :] = 0.0
    # outflow en R_max y z_max
    Q[:, -1, :] = Q[:, -2, :]
    Q[:, :, -1] = Q[:, :, -2]
    return Q

def init_state():
    Q=np.zeros((NV,NR,NZ))
    Q[iRHO]=1.0; Q[iBZ]=Bz0; Q[iPRS]=P_AMB
    return apply_bc(Q, 0.0)

def compute_dt(Q):
    cfR=fast_speed(Q,axis=0).max(); cfZ=fast_speed(Q,axis=1).max()
    return CFL*min(dr/max(cfR,1e-10),dz/max(cfZ,1e-10))

def L(Q,dt):
    return upwind_step_R(Q,dt)+upwind_step_Z(Q,dt)+source_terms(Q,dt)

def rk2_step(Q, dt, t=0.0):
    U=prim_to_cons(Q)
    U1=U+L(Q,dt); Q1=apply_bc(cons_to_prim(U1), t)
    U2=0.5*(U+U1+L(Q1,dt))
    return apply_bc(cons_to_prim(U2), t)

def error_L2(Q_py,Q_ref,var=iRHO):
    return np.sqrt(np.mean((Q_py[var]-Q_ref)**2))

def plot_all_vars(Q,t,save=False,outdir="."):
    var_names=["rho","vR","vz","vphi","BR","Bz","Bphi","prs"]
    R_full=np.concatenate([-R[::-1],R])
    fig,axes=plt.subplots(3,3,figsize=(15,12)); axes=axes.flatten()
    for idx,name in enumerate(var_names):
        field=Q[idx]; f_full=np.concatenate([field[::-1,:],field],axis=0)
        if idx in (iRHO,iPRS):
            cmap_v=plt.colormaps["inferno"]; pos=f_full[f_full>0]
            norm=mcolors.LogNorm(vmin=pos.min() if pos.size else 1e-4,vmax=f_full.max())
        else:
            cmap_v=plt.colormaps["RdBu_r"]; am=np.percentile(np.abs(f_full),99) or 1.0
            norm=mcolors.Normalize(vmin=-am,vmax=am)
        im=axes[idx].pcolormesh(R_full,Z,f_full.T,cmap=cmap_v,norm=norm,shading="auto")
        plt.colorbar(im,ax=axes[idx],label=name)
        axes[idx].set_title(name); axes[idx].set_xlabel("R"); axes[idx].set_ylabel("z")
    axes[-1].set_visible(False)
    fig.suptitle("Python (upwind+RK2)  t="+str(round(t,2)),fontsize=13)
    plt.tight_layout()
    if save:
        fname=os.path.join(outdir,"py_all_t"+str(round(t,2))+".png")
        plt.savefig(fname,dpi=120,bbox_inches="tight"); plt.close(); print("  -> "+fname)
    else:
        plt.show(); plt.close()



def run(outdir="output_mhd", save_times=(1.0, 4.0, 8.0)):
    os.makedirs(outdir,exist_ok=True)
    Q=init_state(); t=0.0; step=0; snaps=[]; t0=time.time()
    targets=sorted(save_times); ti=0
    print("Grid "+str(NR)+"x"+str(NZ)+"  T_END="+str(targets[-1])+"  CFL="+str(CFL))
    print("Guardando en t="+str(targets))
    while ti<len(targets):
        t_stop=targets[ti]
        while t<t_stop-1e-10:
            dt=min(compute_dt(Q),t_stop-t)
            Q=rk2_step(Q,dt,t); t+=dt; step+=1
            if step%500==0:
                print(str(step)+"  t="+str(round(t,3))+"  dt="+str(round(dt,6))+"  "+str(round(time.time()-t0,1))+"s")
        print("  -> snapshot t="+str(round(t,3)))
        snaps.append((t,Q.copy()))
        plot_all_vars(Q,t,save=True,outdir=outdir)
        ti+=1
    print("Fin: "+str(step)+" pasos en "+str(round(time.time()-t0,1))+"s")
    return snaps


