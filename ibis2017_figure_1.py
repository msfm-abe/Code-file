import numpy as np
import matplotlib.pyplot as plt

true_diffusivity, true_ka, true_kd = 1.5, 20.0, 5.
length, epsilon = 35.0 , 1.0
mu, sigma = 0 , 0.1

estimated_diffusivity, estimated_ka, estimated_kd = 1.8, 20., 4.


################################################################################

def reader(mu=mu, sigma=sigma, epsilon=epsilon, length=length, diffusivity=true_diffusivity, ka=true_ka, kd=true_kd):
    #folder = f"./mu={mu},sigma={sigma}/"
    #file = folder+f"artificial_data_De={diffusivity}_length={length}_epsilon={epsilon}_mu={mu}_sigma={sigma}.dat"

    folder = f"./Rmodel_mu={mu},sigma={sigma}/"
    file = folder+f"artificial_data_Rmodel_De={diffusivity}_ka={ka}_kd={kd}_length={length}_epsilon={epsilon}_mu={mu}_sigma={sigma}.dat"

    dimensional_time = np.array([])
    dimensional_flow = np.array([])

    for line in open(file,"r"):
        data = line.split("\t")
        dimensional_time = np.append(dimensional_time,float(data[0]))
        dimensional_flow = np.append(dimensional_flow,float(data[1]))

    return dimensional_time, dimensional_flow

################################################################################
#############################calculates flow curve

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

##############################

dimensionless_time = np.linspace(0.01,10,10000)
time = np.array([])
exit_flow = np.array([])
estimated_flow = np.array([])

time = dimensionless_time*(epsilon*np.power(length,2))/true_diffusivity

# for t in dimensionless_time:
#     flow = 0
#     for n in range(100):
#         flow += np.power(-1,n)*(2*n+1)*np.exp(-np.power(n+0.5,2)*np.power(np.pi,2)*t)
#     flow *= np.pi
#     exit_flow = np.append( exit_flow,flow*true_diffusivity/(epsilon*np.power(length,2)) )

dimensional_time, dimensional_flow = reader()

exit_flow = Flow(k_a=true_ka, k_d=true_kd, time=dimensionless_time)*true_diffusivity/(epsilon*np.power(length,2))

estimated_flow = Flow(k_a=estimated_ka, k_d=estimated_kd, time=dimensionless_time)*true_diffusivity/(epsilon*np.power(length,2))

# k_a = -0.05
# estimated_d_e = 1.54
# for t in dimensionless_time:
#     flow = 0
#     for n in range(100):
#         flow += np.exp(-k_a*t)*np.power(-1,n)*(2*n+1)*np.exp(-np.power(n+0.5,2)*np.power(np.pi,2)*t)
#     flow *= np.pi
#     estimated_flow = np.append( estimated_flow,flow*estimated_d_e/(epsilon*np.power(length,2)) )

#plt.plot(time, dimensional_flow, color = "0.7")
plt.plot(time, dimensional_flow, color = "0.7", label = "Experiment")
plt.plot(time, exit_flow, color = "r", label="True : $D_e=1.50, k_a=20.0, k_d=5.0$")
plt.plot(time , estimated_flow , color = "k", label=f"Estimated : $D_e={estimated_diffusivity}, k_a={estimated_ka}, k_d={estimated_kd}$")
plt.xlabel("Time[ms]")
plt.ylabel("Exit Flow[1/ms]")
plt.legend(loc = "upper right")
plt.grid(True)

#place = f"C:/Users/anbaigashi/Desktop/学会.発表/IBIS2017/資料/ibis2017_artificial_data_De={true_diffusivity}_length={length}_epsilon={epsilon}_mu={mu}_sigma={sigma}.png"
#plt.savefig(place)
plt.savefig("teireikai11_16_2017_fig4.png")
plt.show()
