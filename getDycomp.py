from lib_QMout import *
import numpy 
import sys 

def leer_geom(data):
    nat=int(data[0].split()[0])
    geoms=np.zeros((nat,3))
    symbols=[]

    for j,line in enumerate(data[2:2+nat]):
        parts=line.split()
        symbols.append(parts[0])
        geoms[j,:]=list(map(float,parts[1:4]))

    return nat,symbols,geoms

ic0="./ICOND_00000/QM.out"
ic0dat=read_qmout(ic0)
total_stts,stts=getNstates(ic0dat)
niconds=50

rvals=np.zeros(niconds)
dyson_el=np.zeros((niconds,total_stts,total_stts),dtype=complex)
for i in range(1,niconds,1):
    QMoutfile = f"./ICOND_{i:05d}/QM.out"
    QMinfile = f"./ICOND_{i:05d}/QM.in"

    outdat=read_qmout(QMoutfile)
    indat=read_qmout(QMinfile)

    total_stts,stts=getNstates(outdat)
    dyson=read_Dyson(outdat,total_stts)

    for j in range(total_stts):
        row=dyson[j]
        for k in range(total_stts):
            re=row[2*k]
            im=row[2*k+1]
            dyson_el[i][j][k]=re+1j*im

    nat,symbols,geoms=leer_geom(indat)
    for j in range(nat):
        for k in range(j+1,nat):
            dx,dy,dz=abs(geoms[k,0]-geoms[j,0]),abs(geoms[k,1]-geoms[j,1]),abs(geoms[k,2]-geoms[j,2])
            rvals[i]=np.sqrt(dx*dx+dy*dy+dz*dz)

order=np.argsort(rvals)
r_sorted=rvals[order]
dyson_sorted=dyson_el[order,:,:]
for p in range(total_stts):
    for q in range(total_stts):
        fname=f"./Dyson_Comp/dyson_{p+1:02d}_{q+1:02d}.dat"
        with open(fname,"w") as f:
            for idx in range(niconds):
                if r_sorted[idx]==0.0:
                    continue
                R=r_sorted[idx]
                D=dyson_sorted[idx,p,q]
                Dnorm=np.sqrt(D.real*D.real+D.imag*D.imag)
                f.write(f"{R:14.8f} {Dnorm:20.10f}\n")
                # f.write(f"{R:14.8f} {D:20.10f}\n")
