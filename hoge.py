import numpy as np , matplotlib.pyplot as plt

def reader_(mu, sigma, epsilon, length, diffusivity):
#    folder = f"./mu={mu},sigma={sigma}/"
#    file = folder+f"artificial_data_De={diffusivity}_length={length}_epsilon={epsilon}_mu={mu}_sigma={sigma}.dat"

    s_folder = "./mu={0},sigma={1}/"
    folder = s_folder.format(mu, sigma)
    s_file = folder + "artificial_data_De={0}_length={1}_epsilon={2}_mu={3}_sigma={4}.dat"
#    s_file = "./artificial_data_De={0}_length={1}_epsilon={2}_mu={3}_sigma={4}.dat"
    file = s_file.format(diffusivity, length, epsilon, mu, sigma)

    dimensional_time = np.array([])
    dimensional_flow = np.array([])
    #count = 0
    for line in open(file,"r"):
        data = line.split("\t")
        dimensional_time = np.append(dimensional_time,float(data[0]))
        dimensional_flow = np.append(dimensional_flow,float(data[1]))

    return dimensional_time , dimensional_flow


def reader(mu, sigma, epsilon, length, diffusivity, Ka, Kd):
#    folder = f"./mu={mu},sigma={sigma}/"
#    file = folder+f"artificial_data_De={diffusivity}_length={length}_epsilon={epsilon}_mu={mu}_sigma={sigma}.dat"

    s_folder = "./Rmodel_mu={0},sigma={1}/"
    folder = s_folder.format(mu, sigma)
    s_file = folder + "artificial_data_Rmodel_De={0}_ka={1}_kd={2}_length={3}_epsilon={4}_mu={5}_sigma={6}.dat"
#    s_file = "./artificial_data_De={0}_length={1}_epsilon={2}_mu={3}_sigma={4}.dat"
    file = s_file.format(diffusivity, Ka, Kd, length, epsilon, mu, sigma)

    dimensional_time = np.array([])
    dimensional_flow = np.array([])
    #count = 0
    for line in open(file,"r"):
        data = line.split("\t")
        dimensional_time = np.append(dimensional_time,float(data[0]))
        dimensional_flow = np.append(dimensional_flow,float(data[1]))

    return dimensional_time , dimensional_flow

c , d = reader_(mu=0,sigma=0.005,epsilon=1.0,length=35.0,diffusivity=1.5)

a , b = reader(mu=0,sigma=0.005,epsilon=1.0,length=35.0,diffusivity=1.5,Ka=20.0,Kd=5.0)

plt.plot(a,b)
plt.plot(c,d)
plt.xlabel("Time[ms]")
plt.ylabel("Exit Flow[1/ms]")
#plt.savefig("ibis2017autum_sokuteirei1.png")
plt.show()
