import numpy as np, matplotlib.pyplot as plt ,os

length, epsilon = 35.0 , 1.0 #実験系依存のパラメータ
diffusivity, k_a, k_d  = 1.5, 20., 5.

mu, sigma = 0 , 0.1 #ガウスノイズのパラメータ

dimensionless_time = np.linspace(0.01,10,10000)

dimensional_time = np.array([])
dimensional_exit_flow = np.array([])

##################################

def p_n(n):
    return (n+0.5)*np.pi

def r_plus(p_n, k_a, k_d):
    return ( -(p_n**2 + k_a + k_d)+np.sqrt((p_n**2 + k_a + k_d)**2 - 4*(p_n**2)*k_d) )/2.0

def r_minus(p_n, k_a, k_d):
    return ( -(p_n**2 + k_a + k_d)-np.sqrt((p_n**2 + k_a + k_d)**2 - 4*(p_n**2)*k_d) )/2.0

def A_n(p_n, r_plus, r_minus, k_a):
    return (r_plus + (p_n**2) + k_a)/(r_plus - r_minus)

def Flow(k_a, k_d, time):
    r_flow = np.zeros(len(time))
    for n in range(100):
        pn = p_n(n); r_p = r_plus(p_n=pn, k_a=k_a, k_d=k_d); r_m = r_minus(p_n=pn, k_a=k_a, k_d=k_d)
        A = A_n(p_n=pn, r_plus=r_p, r_minus=r_m, k_a=k_a)
        r_flow += np.power(-1.0,n)*(2.0*n+1.0)*( A*np.exp(r_m*time)+(1.0-A)*np.exp(r_p*time) )
    r_flow *= np.pi

    return r_flow

##################################

dimensional_time = dimensionless_time*(epsilon*np.power(length,2))/diffusivity

dimensionless_exit_flow = np.array([])
dimensionless_exit_flow = Flow(k_a=k_a, k_d=k_d, time=dimensionless_time)
for i in range(len(dimensionless_exit_flow)):
    dimensionless_exit_flow[i] += np.random.normal(mu, sigma)
dimensional_exit_flow = dimensionless_exit_flow*diffusivity/(epsilon*np.power(length,2))

folder = f"./Rmodel_mu={mu},sigma={sigma}"
file =  folder+f"/artificial_data_Rmodel_De={diffusivity}_ka={k_a}_kd={k_d}_length={length}_epsilon={epsilon}_mu={mu}_sigma={sigma}.dat"
string = folder+f"/artificial_data_Rmodel_De={diffusivity}_ka={k_a}_kd={k_d}_length={length}_epsilon={epsilon}_mu={mu}_sigma={sigma}.png"

if not os.path.exists(folder):
    os.makedirs(folder)

with open(file,"w") as fileobj:
    for i in range(len(dimensional_exit_flow)):
        fileobj.write( str(round(dimensional_time[i],1)) )
        fileobj.write("\t")
        fileobj.write( str(round(dimensional_exit_flow[i],7)) )
        fileobj.write("\n")

plt.plot(dimensional_time, dimensional_exit_flow)
plt.xlabel("Time[ms]")
plt.ylabel("Exit Flow[1/ms]")
plt.savefig(string)
plt.show()
