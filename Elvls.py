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

e_vals=[]
e_mat=list(e_mat)
for i in range(ttl_stts):
    e_list=list(e_mat[i])
    for j in range(ttl_stts*2):
        if e_list[j] != 0.0:
            e_vals.append(e_list[j])

e_vals=np.array(e_vals)
e_vals=(e_vals-min(e_vals))*to_eV
np.savetxt("./Evals_states.dat",e_vals,fmt="%.8f")

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
            diff=abs(e_vals[j]-e_vals[j-1])
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
            plt.hlines(e_vals[j],x1,x2,lw=2,color=plot_c[i],label="Mult= "+str(int((2*i/2)+1)))
        else:
            plt.hlines(e_vals[j],x1,x2,lw=2,color=plot_c[i])

    counter+=nmstates

ax.legend(fontsize=15,frameon=False)
fig.tight_layout()
plt.show()
