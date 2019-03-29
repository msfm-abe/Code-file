import numpy as np,matplotlib.pyplot as plt

dimensionless_time = np.linspace(0.01,5.0,10000)
dimensionless_exit_flow = np.array([])
for t in dimensionless_time:
    flow = 0
    for n in range(100):
        flow += np.power(-1,n)*(2*n+1)*np.exp(-np.power(n+0.5,2)*np.power(np.pi,2)*t)
    flow *= np.pi
    dimensionless_exit_flow = np.append(dimensionless_exit_flow,flow)

#plt.yscale("log")
plt.plot(dimensionless_time,dimensionless_exit_flow)
plt.grid(True)
plt.xlabel("Dimensionless Time")
plt.ylabel("Dimensionless Exit Flow")
plt.title("Standard Diffusion Curve")
#plt.savefig("standard_diffusion_curve.eps")
plt.show()
