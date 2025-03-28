from lib_QMout import *
import matplotlib.pyplot as plt

plot_c={
    0: "red",
    1: "blue",
    2: "green",
    3: "yellow",
    4: "orange",
    5: "gray"
}

to_eV=27.4112

file_path="./QM.out"

data=read_qmout(file_path)
ttl_stts,stts=getNstates(data)
e_mat=list(getE(data,ttl_stts))
e_vals=get_SOC(e_mat,ttl_stts)

soc_eigval,soc_eigvec=np.linalg.eig(e_vals)
#soc_eigval=(soc_eigval-min(soc_eigval))*to_eV

np.savetxt("./Evals_states.dat",soc_eigval,fmt="%.8f")

#e_vals=(e_vals-min(e_vals))*to_eV
thrshld=float(input("Threshold: "))
fig,ax = plt.subplots(figsize=(5,10))
counter=0
nst=0
for i in range(len(stts)):
    nmstates=(i+1)*int(stts[i])
    chck_rep=True
    for j in range(counter,counter+nmstates):
        if i > 4:
            c=5

        if j != counter:
            diff=abs(soc_eigval[j]-soc_eigval[j-1])
            if diff < thrshld and chck_rep:
                nst+=1
            else:
                nst=0
        else:
            nst=0

        x1=nst*(1.0+0.1)
        x2=nst*(1.0+0.1)+1.0
            
        if j == counter+int(stts[i]):
            chck_rep=False

        if j == counter:
            plt.hlines(soc_eigval[j],x1,x2,lw=2,color=plot_c[i],label="Mult= "+str(int((2*i/2)+1)))
        else:
            plt.hlines(soc_eigval[j],x1,x2,lw=2,color=plot_c[i])

    counter+=nmstates

ax.legend(fontsize=15,frameon=False)
fig.tight_layout()
plt.show()
