import numpy as np
import matplotlib.pyplot as plot

mu, sigma = 0, 0.005
length, epsilon = 35.0, 1.0
De = 1.5
ka, kd = 20., 5.

################################################################################

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

# def reader_irmodel(mu, sigma, De, length, epsilon):
#     folder = f"./mu={mu},sigma={sigma}/"
# #    file = folder+f"artificial_data_De={diffusivity}_length={length}_epsilon={epsilon}_mu={mu}_sigma={sigma}.dat"
#     file = folder+f"squared_error_Rmodel_De={De}_length={length}_epsilon={epsilon}_mu={mu}_sigma={sigma}.dat"
#
#     diffusivity = np.array([]); adsorption_rate = np.array([]); desorption_rate = np.array([]); squared_error = np.array([])
#     for line in open(file,"r"):
#         data = line.split("\t")
#         diffusivity = np.append(diffusivity, float(data[0]))
#         adsorption_rate = np.append(adsorption_rate, float(data[1]))
#         desorption_rate = np.append(desorption_rate, float(data[2]))
#         squared_error = np.append(squared_error, float(data[3]))
#
#     return diffusivity, adsorption_rate, desorption_rate, squared_error
#
# def reader_rmodel(mu, sigma, De, ka, kd, length, epsilon):
#     file = f"./Rmodel_mu={mu},sigma={sigma}/squared_error_Rmodel_De={De}_ka={ka}_kd={kd}_length={length}_epsilon={epsilon}_mu={mu}_sigma={sigma}.dat"
#
#     diffusivity = np.array([]); adsorption_rate = np.array([]); desorption_rate = np.array([]); squared_error = np.array([])
#     for line in open(file,"r"):
#         data = line.split("\t")
#         diffusivity = np.append(diffusivity, float(data[0]))
#         adsorption_rate = np.append(adsorption_rate, float(data[1]))
#         desorption_rate = np.append(desorption_rate, float(data[2]))
#         squared_error = np.append(squared_error, float(data[3]))
#
#     return diffusivity, adsorption_rate, desorption_rate, squared_error

################################################################################

a, b = reader(mu=0,sigma=0.005,epsilon=1.0,length=35.0,diffusivity=1.5,Ka=20.0,Kd=5.0)
plt.plot(a,b)
plt.show()

# diffusivity, adsorption_rate, desorption_rate, squared_error = reader_irmodel(mu=mu, sigma=sigma, De=De, length=length, epsilon=epsilon)
# print(diffusivity)
# print("\n")
# print(squared_error)
