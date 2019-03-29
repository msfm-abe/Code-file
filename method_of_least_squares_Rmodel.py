import numpy as np,matplotlib.pyplot as plt
import time

#start = time.time()
#e_time = np.array([])

mu , sigma = 0 , 0.005 #ノイズパラメータ
epsilon , length = 1.0 , 35.0 #実験系パラメータ
diffusivity = 1.5 # ターゲットパラメータ

d_e = np.round(np.linspace(1.0,3.0,3),2) #パラメータ探索範囲
k_a = np.round(np.linspace(0,1,3),2)
k_d = np.round(np.linspace(0,1,3),2)


#以下，関数定義
###############################################################################

def reader(mu , sigma , epsilon , length , diffusivity):
    folder = f"./mu={mu},sigma={sigma}/"
    file = folder+f"artificial_data_De={diffusivity}_length={length}_epsilon={epsilon}_mu={mu}_sigma={sigma}.dat"

    dimensional_time = np.array([])
    dimensional_flow = np.array([])
    #count = 0
    for line in open(file,"r"):
        data = line.split("\t")
        dimensional_time = np.append(dimensional_time,float(data[0]))
        dimensional_flow = np.append(dimensional_flow,float(data[1]))
        #print(dimensional_time,dimensional_flow)
        # count += 1
        # if count==1000:
        #     break

    return dimensional_time , dimensional_flow# , count

def dimensionless(dimensional_time , dimensional_flow  , diffusivity=diffusivity, epsilon=epsilon , length=length):
    dimensionless_time = np.array([])
    dimensionless_flow = np.array([])
    dimensionless_time = ( diffusivity/(epsilon*np.power(length,2)) )*dimensional_time
    dimensionless_flow = ( (epsilon*np.power(length,2))/diffusivity )*dimensional_flow

    return dimensionless_time , dimensionless_flow


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


def method_of_least_squares(dimensionless_flow, sdc):
    error = 0.0
    for artificial_data , standard_diffusion_curve in zip(dimensionless_flow,sdc):
        error += np.power(artificial_data - standard_diffusion_curve,2)

    return error

def file_writer(d_e, k_a, k_d, squared_error):
    #file = f"./mu={mu},sigma={sigma}/squared_error_De={diffusivity}_length={length}_epsilon={epsilon}_mu={mu}_sigma={sigma}.dat"
    file = "sample.dat"
    with open(file,"w") as fileobj:
        for i,de in enumerate(d_e):
            for j,a in enumerate(k_a):
                for k,d in enumerate(k_d):
                    fileobj.write(str(de)); fileobj.write("\t")
                    fileobj.write(str(a));fileobj.write("\t")
                    fileobj.write(str(d));fileobj.write("\t")
                    fileobj.write(str(squared_error[i,j,k]))
                    fileobj.write("\n")

    return None

def graph_plot(d_e , squared_error):
    folder = f"./mu={mu},sigma={sigma}/"
    file = folder+f"squared_error_De={diffusivity}_length={length}_epsilon={epsilon}_mu={mu}_sigma={sigma}.png"
    plt.plot(d_e,squared_error)
    plt.xlabel("Diffusivity[$mm^2/ms$]")
    plt.ylabel("Squared Error")
    plt.grid(True)
    plt.savefig(file)
    plt.show()

    return None

###############################################################################


dimensional_time , dimensional_flow = reader(mu=mu,sigma=sigma,epsilon=epsilon,length=length,diffusivity=diffusivity)

squared_error = np.array([])

for d in d_e:
    for ka in k_a:
        for kd in k_d:
            dimensionless_time , dimensionless_flow = dimensionless(diffusivity=d,dimensional_time=dimensional_time,dimensional_flow=dimensional_flow)
            #sdc = standard_diffusion_curve(dimensionless_time=dimensionless_time , k_a = ka)
            sdc = Flow(k_a=ka, k_d=kd, time=dimensionless_time)
            error = method_of_least_squares(dimensionless_flow=dimensionless_flow , sdc=sdc)
            squared_error = np.append(squared_error,error)

            # step = time.time()
            # elapsed_time = step-start
            # e_time = np.append(e_time,elapsed_time)

squared_error = squared_error.reshape(len(d_e), len(k_a), len(k_d))

file_writer(d_e=d_e, k_a=k_a, k_d=k_d, squared_error=squared_error)

#graph_plot(d_e=d_e , squared_error=squared_error)

# end = time.time()
# elapsed_time = end - start
# e_time = np.append(e_time,elapsed_time)
# print(f"elapsed time : {elapsed_time}[sec]")

#plt.plot(e_time)
#plt.show()

# for i in range(count):
#     print("time["+str(i)+"]=",dimensional_time[i],"flow["+str(i)+"]=",dimensional_flow[i])
#     print("d_time["+str(i)+"]=",dimensionless_time[i],"d_flow["+str(i)+"]=",dimensionless_flow[i])
#     print("d_time["+str(i)+"]=",dimensionless_time[i],"sdc["+str(i)+"]=",sdc[i])
#     print("\n")

#print(squared_error)

#graph_plot(dimensionless_time=dimensionless_time , dimensionless_flow=dimensionless_flow , sdc=sdc)
