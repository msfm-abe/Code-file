#!/usr/bin/env python
# coding: utf-8


import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
from IPython.core.pylabtools import figsize
import math


# 読み込みファイル，計算に使用するデータ数，モンテカルロステップ数`T`，モンテカルロサンプル数`sample`を指定．
# 
# `epsilon` と `length` は今後の計算に必要．
data_file = 'TAP_artificial_data_De=1.5_ka=20_kd=5_sigma=5e-6_1-10^4ms.dat'
data_points = 10
#T = 819200
T = 5
sample = 1
#model = np.array([0, 0.5, 1.])
model = 0.1*np.arange(0, 11, dtype=float)
#model = 0.01*np.arange(60, 71, dtype=float)

true_diffusivity = 1.5
true_k_a = 20.0
true_k_d = 5.0
true_standard_deviation = 5e-6
epsilon = 1.0
length = 35.0


def p_n(n):
    return (n+0.5)*np.pi

def r_plus(p_n, k_a, k_d):
    return ( -(p_n**2 + k_a + k_d)+np.sqrt((p_n**2 + k_a + k_d)**2 - 4*(p_n**2)*k_d) )/2.0

def r_minus(p_n, k_a, k_d):
    return ( -(p_n**2 + k_a + k_d)-np.sqrt((p_n**2 + k_a + k_d)**2 - 4*(p_n**2)*k_d) )/2.0

def A_n(p_n, r_plus, r_minus, k_a):
    return (r_plus + (p_n**2) + k_a)/(r_plus - r_minus)

# def Flow(diffusivity, k_a, k_d, epsilon, length, time):
#     dimensionless_time = time * diffusivity/(epsilon*length**2)
    
#     f_flow = np.zeros(time.shape[0])
#     pn = p_n(0); r_p = r_plus(p_n=pn, k_a=k_a, k_d=k_d); r_m = r_minus(p_n=pn, k_a=k_a, k_d=k_d)
#     A = A_n(p_n=pn, r_plus=r_p, r_minus=r_m, k_a=k_a)
#     f_flow += np.power(-1.0, 0)*(2.0*0+1.0)*( A*np.exp(r_m*dimensionless_time)+(1.0-A)*np.exp(r_p*dimensionless_time) )
    
#     l_flow = np.zeros(time.shape[0])
#     pn = p_n(1); r_p = r_plus(p_n=pn, k_a=k_a, k_d=k_d); r_m = r_minus(p_n=pn, k_a=k_a, k_d=k_d)
#     A = A_n(p_n=pn, r_plus=r_p, r_minus=r_m, k_a=k_a)
#     l_flow = f_flow + np.power(-1.0, 1)*(2.0*1+1.0)*( A*np.exp(r_m*dimensionless_time)+(1.0-A)*np.exp(r_p*dimensionless_time) )
    
#     N = np.ones(time.shape[0], dtype=int)

#     for i in range(time.shape[0]):
#         while np.abs((l_flow[i]-f_flow[i])/f_flow[i]) > 10**-5:
#             N[i] += 1
#             f_flow[i] = l_flow[i]
#             pn = p_n(N[i]); r_p = r_plus(p_n=pn, k_a=k_a, k_d=k_d); r_m = r_minus(p_n=pn, k_a=k_a, k_d=k_d)
#             A = A_n(p_n=pn, r_plus=r_p, r_minus=r_m, k_a=k_a)
#             l_flow[i] += np.power(-1.0, N[i])*(2.0*N[i]+1.0)*( A*np.exp(r_m*dimensionless_time[i])+(1.0-A)*np.exp(r_p*dimensionless_time[i]) )
    
##     while np.abs((l_flow[0]-f_flow[0])/f_flow[0]) > 10**-5:
##         for i in range(time.shape[0]):
##             if np.abs((l_flow[i]-f_flow[i])/f_flow[i]) > 1.*10**-5:
##                 N[i] += 1
##                 f_flow[i] = l_flow[i]
##                 pn = p_n(N[i]); r_p = r_plus(p_n=pn, k_a=k_a, k_d=k_d); r_m = r_minus(p_n=pn, k_a=k_a, k_d=k_d)
##                 A = A_n(p_n=pn, r_plus=r_p, r_minus=r_m, k_a=k_a)
##                 l_flow[i] += np.power(-1.0, N[i])*(2.0*N[i]+1.0)*( A*np.exp(r_m*dimensionless_time[i])+(1.0-A)*np.exp(r_p*dimensionless_time[i]) )
    
# #    print(N)
#     l_flow *= np.pi
#     l_flow *= diffusivity/(epsilon*np.power(length, 2))

#     return l_flow

def Flow(diffusivity, k_a, k_d, epsilon, length, time):
    dimensionless_time = time * diffusivity/(epsilon*length**2)
    r_flow = np.zeros(len(time))
    for n in range(100):
        pn = p_n(n); r_p = r_plus(p_n=pn, k_a=k_a, k_d=k_d); r_m = r_minus(p_n=pn, k_a=k_a, k_d=k_d)
        A = A_n(p_n=pn, r_plus=r_p, r_minus=r_m, k_a=k_a)
        r_flow += np.power(-1.0,n)*(2.0*n+1.0)*( A*np.exp(r_m*dimensionless_time)+(1.0-A)*np.exp(r_p*dimensionless_time) )
    r_flow *= np.pi
    r_flow *= diffusivity/(epsilon*np.power(length, 2))

    return r_flow


dimensional_time = np.zeros(data_points)
dimensional_exit_flow = np.zeros(data_points)
i = 0
j = 0
for line in open(data_file, "r"):
    if i%(10000/data_points)==0:
        data = line.split("\t")
        dimensional_time[j] = data[0]
        dimensional_exit_flow[j] = data[1]
        j += 1
    else:
        pass
    i += 1


def energy(model, diffusivity_0, diffusivity_1, k_a, k_d, std,
           epsilon=epsilon, length=length, 
           time=dimensional_time, artificial_flow=dimensional_exit_flow):

    Energy = np.square(artificial_flow - ((1.-model)*Flow(diffusivity=diffusivity_0,
                                                     k_a=0., k_d=0.,
                                                     epsilon=epsilon, length=length, 
                                                     time=time) + \
                                          model*Flow(diffusivity=diffusivity_1,
                                                     k_a=k_a, k_d=k_d, 
                                                     epsilon=epsilon, length=length, 
                                                     time=time))
                      ).sum()
    
    Energy *= 1./(2.*std**2)
    Energy += (time.shape[0]/2.) * np.log(2.*np.pi*std**2) 
    
    return Energy
        

def r(model_cand, model_current,
       diffusivity_0_cand, diffusivity_1_cand, k_a_cand, k_d_cand, std_cand,
       diffusivity_0_current, diffusivity_1_current, k_a_current, k_d_current, std_current, 
       epsilon=epsilon, length=length, 
       time=dimensional_time, artificial_flow=dimensional_exit_flow):
    
    E_cand = energy(model=model_cand,
                    diffusivity_0=diffusivity_0_cand,
                    diffusivity_1=diffusivity_1_cand,
                    k_a=k_a_cand, k_d=k_d_cand,
                    std=std_cand, 
                    epsilon=epsilon, length=length, 
                    time=time, artificial_flow=artificial_flow)
    
    E_current = energy(model=model_current, 
                       diffusivity_0=diffusivity_0_current,
                       diffusivity_1=diffusivity_1_current,
                       k_a=k_a_current, k_d=k_d_current,
                       std=std_current,
                       epsilon=epsilon, length=length, 
                       time=time, artificial_flow=artificial_flow)
    
    
    exponent = -E_cand + E_current
    
    return np.exp(exponent), E_cand, E_current


def r_prime(model1, diffusivity1_0, diffusivity1_1, k_a1, k_d1, std1, 
            model2, diffusivity2_0, diffusivity2_1, k_a2, k_d2, std2, 
            epsilon=epsilon, length=length, time=dimensional_time, artificial_flow=dimensional_exit_flow):
    
    E1 = energy(model=model1, diffusivity_0=diffusivity2_0, diffusivity_1=diffusivity2_1,
                k_a=k_a2, k_d=k_d2, std=std2,
                epsilon=epsilon, length=length,
                time=dimensional_time, artificial_flow=dimensional_exit_flow)
    
    E2 = energy(model=model2, diffusivity_0=diffusivity1_0, diffusivity_1=diffusivity1_1,
                k_a=k_a1, k_d=k_d1, std=std1,
                epsilon=epsilon, length=length, 
                time=dimensional_time, artificial_flow=dimensional_exit_flow)
        
    E3 = energy(model=model1, diffusivity_0=diffusivity1_0, diffusivity_1=diffusivity1_1,
                k_a=k_a1, k_d=k_d1, std=std1,
                epsilon=epsilon, length=length, 
                time=dimensional_time, artificial_flow=dimensional_exit_flow)
        
    E4 = energy(model=model2, diffusivity_0=diffusivity2_0, diffusivity_1=diffusivity2_1,
                k_a=k_a2, k_d=k_d2, std=std2,
                epsilon=epsilon, length=length,
                time=dimensional_time, artificial_flow=dimensional_exit_flow)
    
    E = -E1-E2+E3+E4
    
    return np.exp(E)


def dlog_p(model, diffusivity_0, diffusivity_1, k_a, k_d, 
           std, time=dimensional_time, artificial_flow=dimensional_exit_flow, 
           epsilon=epsilon, length=length):
    
    dlogp = -1./(std**2)             *(
                 (Flow(diffusivity=diffusivity_0, k_a=0, k_d=0, time=time, epsilon=epsilon, length=length) \
                  - Flow(diffusivity=diffusivity_1, k_a=k_a, k_d=k_d, time=time, epsilon=epsilon, length=length)) \
                 *(artificial_flow - (1-model)*Flow(diffusivity=diffusivity_0,
                                                    k_a=0, k_d=0,
                                                    time=time, epsilon=epsilon, length=length) \
                   - model*Flow(diffusivity=diffusivity_1, 
                                k_a=k_a, k_d=k_d, 
                                time=time, epsilon=epsilon, length=length)\
                   )
              ).sum()
           
    return dlogp


def exchange_monte_calro(monte_calro_step, model, epsilon=epsilon, length=length, 
                         time=dimensional_time, artificial_flow=dimensional_exit_flow):

    diffusivity_0 = np.zeros([model.shape[0], monte_calro_step])
    diffusivity_1 = np.zeros([model.shape[0], monte_calro_step])
    k_a = np.zeros([model.shape[0], monte_calro_step])
    k_d = np.zeros([model.shape[0], monte_calro_step])
    sigma = np.zeros([model.shape[0], monte_calro_step])

    E = np.zeros([model.shape[0], monte_calro_step])
    dlogp = np.zeros([model.shape[0], monte_calro_step])
    
    pick = np.zeros([model.shape[0], 5])
    update = np.zeros([model.shape[0], 5])
    
    exchange_rate = np.zeros(model.shape[0] - 1)
    
    """initialize parameters"""
    
    for i in range(model.shape[0]):
        diffusivity_0[i, 0] = np.random.uniform(0, 50)
        diffusivity_1[i, 0] = np.random.uniform(0, 50)
        k_a[i, 0] = np.random.uniform(0, 50)
        k_d[i, 0] = np.random.uniform(0, 50)
        sigma[i, 0] = np.random.uniform(1e-6, 5e-5)
        E[i, 0] = energy(model[i], diffusivity_0[i, 0], diffusivity_1[i, 0],
                         k_a[i, 0], k_d[i, 0], sigma[i, 0])
        dlogp[i, 0] = dlog_p(model[i], diffusivity_0[i, 0], diffusivity_1[i, 0], k_a[i, 0], k_d[i, 0],
                             sigma[i, 0])
        
    t = 1

    while t < monte_calro_step:
        
        """update parameters"""
        
        for m in range(model.shape[0]):
            diffusivity_0[m, t] = diffusivity_0[m, t-1]; diffusivity_1[m, t] = diffusivity_1[m, t-1]; 
            k_a[m, t] = k_a[m, t-1]; k_d[m, t] = k_d[m, t-1]; sigma[m, t] = sigma[m, t-1]
            
            choice = np.random.randint(0, 5)
            if choice == 0:
                diffusivity_0[m, t] += stats.uniform.rvs(loc=-1e0, scale=2e0)
                if diffusivity_0[m, t] < 0:
                    diffusivity_0[m, t] = diffusivity_0[m, t-1]
            elif choice == 1:
                diffusivity_1[m, t] += stats.uniform.rvs(loc=-1e0, scale=2e0)
                if diffusivity_1[m, t] < 0:
                    diffusivity_1[m, t] = diffusivity_1[m, t-1]
            elif choice == 2:
                k_a[m, t] += stats.uniform.rvs(loc=-1e0, scale=2e0)
                if k_a[m, t] < 0:
                    k_a[m, t] = k_a[m, t-1]
            elif choice == 3:
                k_d[m, t] += stats.uniform.rvs(loc=-1e0, scale=2e0)
                if k_d[m, t] < 0:
                    k_d[m, t] = k_d[m, t-1]
            elif choice == 4:
                sigma[m, t] += stats.uniform.rvs(loc=-1e-5, scale=2e-5)
                if sigma[m, t] < 0:
                    sigma[m, t] = sigma[m, t-1]

            pick[m, choice] += 1

            transition_prob, e_cand, e_current = r(model_cand=model[m], model_current=model[m],
                                                   diffusivity_0_cand=diffusivity_0[m, t], 
                                                   diffusivity_1_cand=diffusivity_1[m, t],
                                                   k_a_cand=k_a[m, t], k_d_cand=k_d[m, t],
                                                   std_cand=sigma[m, t],
                                                   diffusivity_0_current=diffusivity_0[m, t-1],
                                                   diffusivity_1_current=diffusivity_1[m, t-1],
                                                   k_a_current=k_a[m, t-1], k_d_current=k_d[m, t-1],
                                                   std_current=sigma[m, t-1]) 
            
            
            if np.array([1.0, transition_prob]).min() > stats.uniform.rvs():
                update[m, choice] += 1
                E[m, t] = e_cand
                if model[m] == 0:
                    diff_1 = stats.uniform.rvs(loc=0., scale=50.)
                    adr = stats.uniform.rvs(loc=0., scale=50.); disr = stats.uniform.rvs(loc=0., scale=50.)
                    dlogp[m, t] = dlog_p(model[m], diffusivity_0[m, t], diff_1, adr, disr, sigma[m, t])
                    diffusivity_1[m, t] = diff_1; k_a[m, t] = adr; k_d[m, t] = disr
                elif model[m] == 1:
                    diff_0 = stats.uniform.rvs(loc=0., scale=50.)
                    dlogp[m, t] = dlog_p(model[m], diff_0, diffusivity_1[m, t], k_a[m, t], k_d[m, t], sigma[m, t])
                    diffusivity_0[m, t] = diff_0
                else:
                    dlogp[m, t] = dlog_p(model[m], diffusivity_0[m, t], diffusivity_1[m, t], k_a[m, t], k_d[m, t],
                                     sigma[m, t])
            else:
                diffusivity_0[m, t] = diffusivity_0[m, t-1]
                diffusivity_1[m, t] = diffusivity_1[m, t-1]
                k_a[m, t] = k_a[m, t-1]
                k_d[m, t] = k_d[m, t-1]
                sigma[m, t] = sigma[m, t-1]
                E[m, t] = E[m, t-1]
                dlogp[m, t] = dlogp[m, t-1]
                
        """exchange monte calro"""
        
#        x = np.arange(1, model.shape[0])
#        y = np.random.permutation(x) 
        
#        for m in y:
#            exchange_prob = r_prime(model1=model[m-1], diffusivity1_0=diffusivity_0[m-1, t], 
#                                    diffusivity1_1=diffusivity_1[m-1, t],
#                                    k_a1=k_a[m-1, t], k_d1=k_d[m-1, t],
#                                    std1=sigma[m-1, t], 
#                                    model2=model[m], diffusivity2_0=diffusivity_0[m, t],
#                                    diffusivity2_1=diffusivity_1[m, t],
#                                    k_a2=k_a[m, t], k_d2=k_d[m, t],
#                                    std2=sigma[m, t])
            
#            if np.array([1.0, exchange_prob]).min() > stats.uniform.rvs():
#                diffusivity_0[m-1, t], diffusivity_0[m, t] = diffusivity_0[m, t], diffusivity_0[m-1, t]
#                diffusivity_1[m-1, t], diffusivity_1[m, t] = diffusivity_1[m, t], diffusivity_1[m-1, t]
#                k_a[m-1, t], k_a[m, t] = k_a[m, t], k_a[m-1, t]
#                k_d[m-1, t], k_d[m, t] = k_d[m, t], k_d[m-1, t]
#                sigma[m-1, t], sigma[m, t] = sigma[m, t], sigma[m-1, t]
#                E[m-1, t], E[m, t] = E[m, t], E[m-1, t]
#                dlogp[m-1, t], dlogp[m, t] = dlogp[m, t], dlogp[m-1, t]

#                exchange_rate[m-1] += 1

        t += 1
                
    return diffusivity_0, diffusivity_1, k_a, k_d, sigma, E, dlogp, pick, update, exchange_rate


def file_write(w_file, model, diffusivity_0, diffusivity_1, k_a, k_d, sigma, E, dlogp, pick, update, exchange):
    
    with open(w_file, "w") as fileobj:
        
        fileobj.write(str('#model parametes M\n'))
        for i in range(model.shape[0]):
            fileobj.write('#M={}'.format((model[i]))); fileobj.write(str('\t'))
        fileobj.write(str('\n'))
        fileobj.write(str("#acceptance rate\n"))
        for m in range(model.shape[0]):
            for i in range(5):
                fileobj.write(str(pick[m, i])); fileobj.write(str("\t"))
            fileobj.write(str("\n"))
        fileobj.write(str("\n"))
        fileobj.write(str('#update \n'))
        for m in range(model.shape[0]):
            for i in range(5):
                fileobj.write(str(update[m, i])); fileobj.write(str("\t"))
            fileobj.write(str("\n"))
        fileobj.write(str("\n"))

        fileobj.write(str("#exchange rate\n"))
        for m in range(model.shape[0]-1):
            fileobj.write(str(exchange[m])); fileobj.write(str("\t"))
        fileobj.write(str("\n")); fileobj.write(str("\n"))

        fileobj.write(str('#De^M=0, De^M=1, ka, kd, sigma, energy, log bayes f\n'))
        for m in range(model.shape[0]):
            fileobj.write(str("#m={}\n".format(model[m])))
            for i in range(T):
                fileobj.write(str(diffusivity_0[m, i])); fileobj.write(str("\t"))
                fileobj.write(str(diffusivity_1[m, i])); fileobj.write(str("\t"))
                fileobj.write(str(k_a[m, i])); fileobj.write(str("\t"))
                fileobj.write(str(k_d[m, i])); fileobj.write(str("\t"))
                fileobj.write(str(sigma[m, i])); fileobj.write(str("\t"))
                fileobj.write(str(E[m, i])); fileobj.write(str("\t"))
                fileobj.write(str(dlogp[m, i])); fileobj.write(str("\t"))
                fileobj.write(str("\n"))
            fileobj.write(str("\n"))
            
    return 


for i in range(sample):
    w_file = './gakkai/EM_TAP-sigma_koukannashi_datapoints={}_T={}mcs_model={}-{}_replica={}_koushinhaba-m1p1m1e-5p1e-5_{}_0.dat'.format(data_points, T, model[0], model[-1], model.shape[0], data_file)
    
    diffusivity_0, diffusivity_1, k_a, k_d, sigma, E, dlogp, pick, update, exchange    = exchange_monte_calro(monte_calro_step=T, model=model)
    
    file_write(w_file, model, diffusivity_0, diffusivity_1, k_a, k_d, sigma, E, dlogp, pick, update, exchange)
