from lib_QMout import *
import matplotlib.pyplot as plt
import os

plot_c={
    0: "red",
    1: "blue",
    2: "green",
    3: "yellow",
    4: "orange",
    5: "gray"
}

to_eV=27.4112

dirs=os.listdir("./")
iconds=[]

for i in dirs:
    if "ICOND_" in i:
        iconds.append(i)

iconds.sort()

icond0=iconds[0]+"/QM.out"
data=read_qmout(icond0)
ttl_stts,stts=getNstates(data)

icond0=iconds[0]+"/QM.in"
data=read_qmout(icond0)
nat=int(data[0])
nIC=int((nat*(nat-1))/2)
print(nIC)

pec_e=np.zeros((len(iconds),ttl_stts))
geoms=np.zeros((len(iconds),nIC))
points=np.zeros((len(iconds),nIC+ttl_stts))

for dirI in range(len(iconds)):
    if os.path.isdir(iconds[dirI]):
        file_path=iconds[dirI]+"/QM.out"
        if os.path.isfile(file_path):
            data=read_qmout(file_path)
            e_mat=list(getE(data,ttl_stts))
#
            e_mat=list(e_mat)
            for i in range(ttl_stts):
                e_list=list(e_mat[i])
                for j in range(ttl_stts*2):
                    if e_list[j] != 0.0:
                        pec_e[dirI][i]=e_list[j]

        file_path=iconds[dirI]+"/QM.in"
        if os.path.isfile(file_path):
            data=read_qmout(file_path)
            xyz=getxyz(data,nat)

            cntIC=0               
            for k in range(nat-1):
                at1=xyz[k]
                for m in range(k+1,nat):
                    at2=xyz[m]
                    val=0.0
                    for cord in range(3):
                        val+=(float(at1[cord])-float(at2[cord]))**2.
                    geoms[dirI][cntIC]=np.sqrt(val)
                    cntIC+=1

    points[dirI][:nIC]=geoms[dirI][:]
    points[dirI][nIC:]=pec_e[dirI][:]

np.savetxt("PECs_fromICONDs.dat",points,fmt="%.10f")

fig,ax=plt.subplots(figsize=(10,10))
counter=int(nIC)
minval=min(points[:,counter])

for i in range(len(points[:,counter])):
    minvalit=min(points[i,:])
    if minvalit < minval or minvalit == minval:
        minval = minvalit
        index=i

plt.vlines(points[index,0],-1,10,lw=2)
for mult in range(len(stts)):
    nmstt=int(stts[mult])*(mult+1)
    chck_rep=True
    for j in range(counter,counter+nmstt):
        if mult > 4:
            c=5
        
        if j > counter+int(stts[mult]):
            chck_rep=False

        if chck_rep and j == counter:
           # plt.plot(points[:,0],(points[:,j]-minval)*to_eV,lw=2,color=plot_c[mult],label="Mult= "+str(int((2*mult/2)+1)))
            plt.scatter(points[:,0],(points[:,j]-minval)*to_eV,lw=2,color=plot_c[mult],label="Mult= "+str(int((2*mult/2)+1)))
        elif chck_rep:
            #plt.plot(points[:,0],(points[:,j]-minval)*to_eV,lw=2,color=plot_c[mult])
            plt.scatter(points[:,0],(points[:,j]-minval)*to_eV,lw=2,color=plot_c[mult])

    counter+=nmstt
                    
ax.legend(fontsize=15,loc="upper right",frameon=False)
ax.set_xlim(min(points[:,0])-0.01,max(points[:,0])+0.01)
ax.set_ylim(-0.5,10)
fig.tight_layout()
plt.show()

