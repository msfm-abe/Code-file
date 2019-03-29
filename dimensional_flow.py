import numpy as np,matplotlib.pyplot as plt
import os

length , epsilon = 35.0, 0.18 #実験系依存のパラメータ
diffusivity = 1.5

mu , sigma = 0, 0.1 #ガウスノイズのパラメータ

dimensionless_time = np.linspace(0.01,10,10000)
time = np.array([])
exit_flow = np.array([])

#-------------------------------------------------------------------------------

time = dimensionless_time*(epsilon*np.power(length,2))/diffusivity

for t in dimensionless_time:
    flow = 0
    for n in range(100):
        flow += np.power(-1,n)*(2*n+1)*np.exp(-np.power(n+0.5,2)*np.power(np.pi,2)*t)
    flow *= np.pi
    flow += np.random.normal(mu,sigma)
    exit_flow = np.append( exit_flow,flow*diffusivity/(epsilon*np.power(length,2)) )

#-------------------------------------------------------------------------------

folder = "./"+f"mu={mu},sigma={sigma}/"
file = folder+f"artificial_data_De={diffusivity}_length={length}_epsilon={epsilon}_mu={mu}_sigma={sigma}.dat"
string = folder+f"artificial_data_De={diffusivity}_length={length}_epsilon={epsilon}_mu={mu}_sigma={sigma}"

if not os.path.exists(folder):
    os.makedirs(folder)
# file = "./"+string+".dat"
with open(file,"w") as fileobj:
    for i in range(len(exit_flow)):
        fileobj.write( str(round(time[i],1)) )
        fileobj.write("\t")
        fileobj.write( str(round(exit_flow[i],7)) )
        fileobj.write("\n")

plt.plot(time,exit_flow)
# plt.legend(loc="upper right")
plt.grid(True)
plt.xlabel("Time [ms]")
plt.ylabel("Exit Flow [1/ms]")
plt.title("Time-Exit Flow")
#plt.savefig(string+".png")
plt.show()
