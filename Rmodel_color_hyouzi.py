import numpy as np
import matplotlib.pyplot as plt

mu, sigma = 0, 0.1
length, epsilon = 35.0, 1.0
De = 1.5
ka, kd = 20., 5.

################################################################################

def reader_irmodel(mu, sigma, De, length, epsilon):
    folder = f"./mu={mu},sigma={sigma}/"
#    file = folder+f"artificial_data_De={diffusivity}_length={length}_epsilon={epsilon}_mu={mu}_sigma={sigma}.dat"
    file = folder+f"squared_error_Rmodel_De={De}_length={length}_epsilon={epsilon}_mu={mu}_sigma={sigma}.dat"

    diffusivity = np.array([]); adsorption_rate = np.array([]); desorption_rate = np.array([]); squared_error = np.array([])
    for line in open(file,"r"):
        data = line.split("\t")
        diffusivity = np.append(diffusivity, float(data[0]))
        adsorption_rate = np.append(adsorption_rate, float(data[1]))
        desorption_rate = np.append(desorption_rate, float(data[2]))
        squared_error = np.append(squared_error, float(data[3]))

    return diffusivity, adsorption_rate, desorption_rate, squared_error

def reader_rmodel(mu, sigma, De, ka, kd, length, epsilon):
    file = f"./Rmodel_mu={mu},sigma={sigma}/squared_error_Rmodel_De={De}_ka={ka}_kd={kd}_length={length}_epsilon={epsilon}_mu={mu}_sigma={sigma}.dat"

    diffusivity = np.array([]); adsorption_rate = np.array([]); desorption_rate = np.array([]); squared_error = np.array([])
    for line in open(file,"r"):
        data = line.split("\t")
        diffusivity = np.append(diffusivity, float(data[0]))
        adsorption_rate = np.append(adsorption_rate, float(data[1]))
        desorption_rate = np.append(desorption_rate, float(data[2]))
        squared_error = np.append(squared_error, float(data[3]))

    return diffusivity, adsorption_rate, desorption_rate, squared_error

def squared_error_graph(squared_error, d_e):
    x = np.arange(0,1.05,0.05)
    y = np.arange(0,1.05,0.05)
    X, Y = np.meshgrid(x, y)

    squared_error = squared_error.reshape(21,21)

    plt.pcolor(X, Y, squared_error)
    plt.colorbar()
    d_e = round(d_e,2)
    plt.ylabel("Adsorption rate $k_{a}$")
    plt.xlabel("Desorption rate $k_{d}$")
    plt.title(f"Squared Error (Diffusivity={d_e})")
#    plt.savefig(f"hoge_De={d_e}.png")
    plt.savefig(f"./mu={mu},sigma={sigma}/squared_error_Rmodel_De={d_e}_length={length}_epsilon={epsilon}_mu={mu}_sigma={sigma}.png")
    plt.clf()

def squared_error_graph_rmodel(squared_error, d_e):
    x = np.arange(0,31,1)
    y = np.arange(0,31,1)
    X, Y = np.meshgrid(x, y)

    squared_error = squared_error.reshape(31,31)

    plt.pcolor(X, Y, squared_error)
    plt.colorbar()
    d_e = round(d_e,2)
    plt.ylabel("Adsorption rate $k_{a}$")
    plt.xlabel("Desorption rate $k_{d}$")
    plt.title(f"Squared Error (Diffusivity={d_e})")
#    plt.savefig(f"hoge_De={d_e}.png")
    plt.savefig(f"./Rmodel_mu={mu},sigma={sigma}/squared_error_Rmodel_De={d_e}_ka={ka}_kd={kd}_length={length}_epsilon={epsilon}_mu={mu}_sigma={sigma}.dat.png")
    plt.clf()

################################################################################

#diffusivity, adsorption_rate, desorption_rate, squared_error = reader_irmodel(mu=mu, sigma=sigma, De=De, length=length, epsilon=epsilon)

diffusivity, adsorption_rate, desorption_rate, squared_error = reader_rmodel(mu=mu, sigma=sigma, De=De, ka=ka, kd=kd, length=length, epsilon=epsilon)

s = np.array([])
de = 1.3
for i in range(1,51):
    #for j in range(441*(i-1),441*i):
    for j in range(961*(i-1),961*i):
        s = np.append(s, squared_error[j])

    #s = s.reshape(21,21)
    s = s.reshape(31,31)
    de = round(de,2)
#    squared_error_graph(squared_error=s, d_e=de)
    squared_error_graph_rmodel(squared_error=s, d_e=de)
    s = np.array([])
    de += 0.01
    if de>1.54 and de<1.56:
        de += 0.01
