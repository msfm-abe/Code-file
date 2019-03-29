import numpy as np

true_diffusivity = 1.5
length , epsilon = 35.0 , 1.0
mu , sigma = 0 , 0.05

diffusivity = np.array([])
ka = np.array([])
se = np.array([])

################################################################################

def reader(mu=mu , sigma=sigma , epsilon=epsilon , length=length , diffusivity=true_diffusivity):
    folder = f"./mu={mu},sigma={sigma}/"
    file = folder+f"squared_error_IRmodel_De={true_diffusivity}_length={length}_epsilon={epsilon}_mu={mu}_sigma={sigma}.dat"

    diffusivity = np.array([])
    ka = np.array([])
    se = np.array([])

    for line in open(file,"r"):
        data = line.split("\t")
        diffusivity = np.append(diffusivity,float(data[0]))
        ka = np.append(ka,float(data[1]))
        se = np.append(se, float(data[2]))

    return diffusivity, ka, se

################################################################################

diffusivity, ka, se = reader()

file = "hogehoge.dat"
count = 1
with open(file,"w") as fileobj:
    for i in range(len(diffusivity)):
        fileobj.write(str(diffusivity[i])); fileobj.write("\t")
        fileobj.write(str(ka[i])); fileobj.write("\t")
        fileobj.write(str(se[i])); fileobj.write("\n")
        count += 1
        if count==21:
            fileobj.write("\n")
            count = 1

# print(diffusivity)
# print("\n")
# print(ka)
# print("\n")
# print(se)
