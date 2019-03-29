import numpy as np,matplotlib.pyplot as plt

mu , sigma = 0 , 0.005 #ノイズパラメータ
epsilon , length = 1.0 , 35.0 #実験系パラメータ
diffusivity = 1.5 # ターゲットパラメータ

num_sample = 2 # サンプル数
num_time = 1000 # 時間の分割数
num_de = 20 # deの分割数

d_e = np.round(np.linspace(1.0,3.0,num_de),2) #パラメータ探索範囲

artificial_dimensionless_time = np.linspace(0.01,10,num_time)

#以下，関数定義
###############################################################################

def artificial_data(mu=mu , sigma=sigma , epsilon=epsilon , length=length , diffusivity=diffusivity , artificial_dimensionless_time=artificial_dimensionless_time):
    dimensional_time = np.array([])
    dimensional_flow = np.array([])
    dimensional_time = artificial_dimensionless_time*(epsilon*np.power(length,2))/diffusivity

    for time in artificial_dimensionless_time:
        flow = 0
        for n in range(100):
            flow += np.power(-1,n)*(2*n+1)*np.exp(-np.power(n+0.5,2)*np.power(np.pi,2)*time)
        flow *= np.pi
        flow += np.random.normal(mu,sigma)
        dimensional_flow = np.append( dimensional_flow,flow*diffusivity/(epsilon*np.power(length,2)) )

    return dimensional_time , dimensional_flow

def dimensionless(dimensional_time , dimensional_flow  , diffusivity=diffusivity, epsilon=epsilon , length=length):
    dimensionless_time = np.array([])
    dimensionless_flow = np.array([])
    dimensionless_time = ( diffusivity/(epsilon*np.power(length,2)) )*dimensional_time
    dimensionless_flow = ( (epsilon*np.power(length,2))/diffusivity )*dimensional_flow

    return dimensionless_time , dimensionless_flow

def standard_diffusion_curve(dimensionless_time):
    sdc = np.array([])
    for time in dimensionless_time:
        exit_flow = 0.0
        for n in range(100):
            exit_flow += np.pi*np.power(-1,n)*(2*n+1)*np.exp(-np.power(n+0.5,2)*np.power(np.pi,2)*time)
        sdc = np.append(sdc,exit_flow)

    return sdc

def method_of_least_squares(dimensionless_flow, sdc):
    error = 0.0
    for artificial_data , standard_diffusion_curve in zip(dimensionless_flow,sdc):
        error += np.power(artificial_data - standard_diffusion_curve,2)

    return error

###############################################################################


#dimensional_time , dimensional_flow = reader(mu=mu,sigma=sigma,epsilon=epsilon,length=length,diffusivity=diffusivity)

squared_error = np.array([])
min_de = np.array([])

for i in range(num_sample):
    dimensional_time , dimensional_flow = artificial_data(mu=mu , sigma=sigma , epsilon=epsilon , length=length , diffusivity=diffusivity , artificial_dimensionless_time=artificial_dimensionless_time)

    for d in d_e:
        dimensionless_time , dimensionless_flow = dimensionless(diffusivity=d,dimensional_time=dimensional_time,dimensional_flow=dimensional_flow)
        sdc = standard_diffusion_curve(dimensionless_time=dimensionless_time)
        error = method_of_least_squares(dimensionless_flow=dimensionless_flow , sdc=sdc)
        squared_error = np.append(squared_error,error)


squared_error = squared_error.reshape(num_sample,num_de)

#file1 = f"squared_error_sample_array_De={diffusivity}_length={length}_epsilon={epsilon}_mu={mu}_sigma={sigma}.dat"
s_file1 = "squared_error_sample_array_De={0}_length={1}_epsilon={2}_mu={3}_sigma={4}.dat"
file1 = s_file1.format(diffusivity,length,epsilon,mu,sigma)
np.savetxt(file1,squared_error)

#print(np.load("sample1.npy"))

for i in range(num_sample):
    temp = squared_error[i].min()
    pos = np.where(temp==squared_error[i])
    min_de = np.append(min_de,d_e[pos])

#file2 = f"squared_error_sample_hist_De={diffusivity}_length={length}_epsilon={epsilon}_mu={mu}_sigma={sigma}.dat"
s_file2 = "squared_error_sample_hist_De={0}_length={1}_epsilon={2}_mu={3}_sigma={4}.dat"
file2 = s_file2.format(diffusivity,length,epsilon,mu,sigma)
np.savetxt(file2,min_de)

#print(np.load("sample2.npy"))

#a = np.loadtxt("sample2.dat")

#plt.hist(a,range(1,3))
#plt.show()

# print("\n")
# print(d_e)
# print(min_de)
