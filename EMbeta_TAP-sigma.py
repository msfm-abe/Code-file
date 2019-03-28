#!/usr/bin/env python
# coding: utf-8


import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
from IPython.core.pylabtools import figsize
import math


# 読み込みファイル，計算に使用するデータ数，モンテカルロステップ数`T`，逆温度$\beta$(配列),モンテカルロサンプル数`sample`を指定．
# 
# `epsilon` と `length` は今後の計算に必要．
data_file = 'TAP_artificial_data_De=1.5_ka=20_kd=5_sigma=5e-6_1-10^4ms.dat'
data_points = 10
#T = 204800
T = 102400
#T = 5
sample = 1

temp = np.array([3., 4., 6., 7., 8., 9.])
#temp = np.array([512, 1024])
beta = 1/temp

#beta = np.array([1, 0.5, 0.25, 0.2, 0.])
#temp = 1/beta
#temp[-1] = 1e6

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


def r(beta, model_cand, model_current,
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
    
    return np.exp(beta*exponent), E_cand, E_current


def r_prime(beta_l, beta_lplus1, model_l, diffusivity0_l, diffusivity1_l, k_a_l, k_d_l, std_l, 
            model_lplus1, diffusivity0_lplus1, diffusivity1_lplus1, k_a_lplus1, k_d_lplus1, std_lplus1, 
            epsilon=epsilon, length=length, time=dimensional_time, artificial_flow=dimensional_exit_flow):
    
    E1 = energy(model=model_l, diffusivity_0=diffusivity0_lplus1, diffusivity_1=diffusivity1_lplus1,
                k_a=k_a_lplus1, k_d=k_d_lplus1, std=std_lplus1,
                epsilon=epsilon, length=length,
                time=dimensional_time, artificial_flow=dimensional_exit_flow)
    
    E2 = energy(model=model_lplus1, diffusivity_0=diffusivity0_l, diffusivity_1=diffusivity1_l,
                k_a=k_a_l, k_d=k_d_l, std=std_l,
                epsilon=epsilon, length=length, 
                time=dimensional_time, artificial_flow=dimensional_exit_flow)
        
    E3 = energy(model=model_l, diffusivity_0=diffusivity0_l, diffusivity_1=diffusivity1_l,
                k_a=k_a_l, k_d=k_d_l, std=std_l,
                epsilon=epsilon, length=length, 
                time=dimensional_time, artificial_flow=dimensional_exit_flow)
        
    E4 = energy(model=model_lplus1, diffusivity_0=diffusivity0_lplus1, diffusivity_1=diffusivity1_lplus1,
                k_a=k_a_lplus1, k_d=k_d_lplus1, std=std_lplus1,
                epsilon=epsilon, length=length,
                time=dimensional_time, artificial_flow=dimensional_exit_flow)
    
    E = -beta_l*E1-beta_lplus1*E2+beta_l*E3+beta_lplus1*E4
    
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


def exchange_monte_calro(monte_calro_step, beta, model, epsilon=epsilon, length=length, 
                         time=dimensional_time, artificial_flow=dimensional_exit_flow):
    
    diffusivity_0 = np.zeros([2*beta.shape[0], monte_calro_step])
    diffusivity_1 = np.zeros([2*beta.shape[0], monte_calro_step])
    k_a = np.zeros([2*beta.shape[0], monte_calro_step])
    k_d = np.zeros([2*beta.shape[0], monte_calro_step])
    sigma = np.zeros([2*beta.shape[0], monte_calro_step])
    
    Energy = np.zeros([2*beta.shape[0], monte_calro_step])
    dlogp = np.zeros([2*beta.shape[0], monte_calro_step])
    
    pick = np.zeros([2*beta.shape[0], 5]); update = np.zeros([2*beta.shape[0], 5])
    exchange_count = np.zeros(2*beta.shape[0]-1)
    
    """initialize parameters"""
    
    for i in range(2*beta.shape[0]):
        diffusivity_0[i, 0] = np.random.uniform(0, 50)
        diffusivity_1[i, 0] = np.random.uniform(0, 50)
        k_a[i, 0] = np.random.uniform(0, 50)
        k_d[i, 0] = np.random.uniform(0, 50)
        sigma[i, 0] = np.random.uniform(1e-6, 5e-5)

        if i < beta.shape[0]:
            model = 0
        else:
            model = 1
        Energy[i, 0] = energy(model, diffusivity_0[i, 0], diffusivity_1[i, 0], k_a[i, 0], k_d[i, 0], sigma[i, 0])
        dlogp[i, 0] = dlog_p(model, diffusivity_0[i, 0], diffusivity_1[i, 0], k_a[i, 0], k_d[i, 0], sigma[i, 0])

    t = 1

    while t < monte_calro_step:
        """update parameters"""
        for b in range(2*beta.shape[0]):
            diffusivity_0[b, t] = diffusivity_0[b, t-1]; diffusivity_1[b, t] = diffusivity_1[b, t-1]
            k_a[b, t] = k_a[b, t-1]; k_d[b, t] = k_d[b, t-1]; sigma[b, t] = sigma[b, t-1]
            Energy[b, t] = Energy[b, t-1]

            choice = np.random.randint(0, 5); pick[b, choice] += 1
            if choice == 0:
                diffusivity_0[b, t] += stats.uniform.rvs(loc=-1e0, scale=2e0)
                if diffusivity_0[b, t] < 0:
                    diffusivity_0[b, t] = diffusivity_0[b, t-1]
            elif choice == 1:
                diffusivity_1[b, t] += stats.uniform.rvs(loc=-1e0, scale=2e0)
                if diffusivity_1[b, t] < 0:
                    diffusivity_1[b, t] = diffusivity_1[b, t-1]
            elif choice == 2:
                k_a[b, t] += stats.uniform.rvs(loc=-1e0, scale=2e0)
                if k_a[b, t] < 0:
                    k_a[b, t] = k_a[b, t-1]
            elif choice == 3:
                k_d[b, t] += stats.uniform.rvs(loc=-1e0, scale=2e0)
                if k_d[b, t] < 0:
                    k_d[b, t] = k_d[b, t-1]
            elif choice == 4:
                sigma[b, t] += stats.uniform.rvs(loc=-1e-5, scale=2e-5)
                if sigma[b, t] < 0:
                    sigma[b, t] = sigma[b, t-1]
            
            if b+1 <= beta.shape[0]:
                model = 0
                inv_temp = beta[b]
            else:
                model = 1
                inv_temp = beta[2*beta.shape[0]-1-b]
            
#             print('model=',model, 'beta=', inv_temp)
            transition_prob, e_cand, e_current = r(inv_temp, model, model,
                                                   diffusivity_0[b, t], diffusivity_1[b, t],
                                                   k_a[b, t], k_d[b, t], sigma[b, t],
                                                   diffusivity_0[b, t-1], diffusivity_1[b, t-1],
                                                   k_a[b, t-1], k_d[b, t-1], sigma[b, t-1])
            
            if np.array([1.0, transition_prob]).min() > stats.uniform.rvs():
                Energy[b, t] = e_cand
                update[b, choice] += 1
                if model == 0:
                    diff_1 = stats.uniform.rvs(loc=0., scale=50.)
                    adr = stats.uniform.rvs(loc=0., scale=50.); disr = stats.uniform.rvs(loc=0., scale=50.)
                    dlogp[b, t] = dlog_p(model, diffusivity_0[b, t], diff_1, adr, disr, sigma[b, t])
                    diffusivity_1[b, t] = diff_1; k_a[b, t] = adr; k_d[b, t] = disr
                elif model == 1:
                    diff_0 = stats.uniform.rvs(loc=0., scale=50.)
                    dlogp[b, t] = dlog_p(model, diff_0, diffusivity_1[b, t], k_a[b, t], k_d[b, t], sigma[b, t])
                    diffusivity_0[b, t] = diff_0
            else:
                diffusivity_0[b, t] = diffusivity_0[b, t-1]; diffusivity_1[b, t] = diffusivity_1[b, t-1]
                k_a[b, t] = k_a[b, t-1]; k_d[b, t] = k_d[b, t-1]; sigma[b, t] = sigma[b, t-1]
                Energy[b, t] = Energy[b, t-1]; dlogp[b, t] = dlogp[b, t-1]
        
        """exchange parameters"""
        y = np.random.permutation(2*beta.shape[0]-1)
        
        for b in y:
            if b < beta.shape[0]-1:
                beta_l, beta_lplus1 = beta[b], beta[b+1]
                model_l, model_lplus1 = 0., 0.
            elif b == beta.shape[0]-1:
                beta_l, beta_lplus1 = beta[b], beta[b]
                model_l, model_lplus1 = 0., 1.0
            elif b > beta.shape[0]-1:
                beta_l, beta_lplus1 = beta[2*beta.shape[0]-1-b], beta[2*beta.shape[0]-1-b-1]
                model_l, model_lplus1 = 1.0, 1.0
            
            exchange_prob = r_prime(beta_l, beta_lplus1,
                                    model_l, diffusivity_0[b, t], diffusivity_1[b, t],
                                    k_a[b, t], k_d[b, t], sigma[b, t], 
                                    model_lplus1, diffusivity_0[b+1, t], diffusivity_1[b+1, t],
                                    k_a[b+1, t], k_d[b+1, t], sigma[b+1, t])
            
            if np.array([1.0, exchange_prob]).min() > stats.uniform.rvs():
                diffusivity_0[b, t], diffusivity_0[b+1, t] = diffusivity_0[b+1, t], diffusivity_0[b, t]
                diffusivity_1[b, t], diffusivity_1[b+1, t] = diffusivity_1[b+1, t], diffusivity_1[b, t]
                k_a[b, t], k_a[b+1, t] = k_a[b+1, t], k_a[b, t]; k_d[b, t], k_d[b+1, t] = k_d[b+1, t], k_d[b, t]
                sigma[b, t], sigma[b+1, t] = sigma[b+1, t], sigma[b, t]
                Energy[b, t], Energy[b+1, t] = Energy[b+1, t], Energy[b, t]
                dlogp[b, t], dlogp[b+1, t] = dlogp[b+1, t], dlogp[b, t]
                exchange_count[b] += 1
 
        t += 1
                
    return diffusivity_0, diffusivity_1, k_a, k_d, sigma, Energy, dlogp, pick, update, exchange_count


def file_write(w_file, beta, model, diffusivity_0, diffusivity_1, k_a, k_d, sigma, E, dlogp, pick, update, exchange):
    
    with open(w_file, "w") as fileobj:
                
        fileobj.write(str('#model parametes M, temperature \n'))
        for i in range(model.shape[0]):
            for j in range(beta.shape[0]):
                if i==0 and j==0:
                    fileobj.write('#(M, Temp)=({}, {})'.format(model[i], temp[j])); fileobj.write(str('\t'))
                elif i==0:
                    fileobj.write('({}, {})'.format(model[i], temp[j])); fileobj.write(str('\t'))
                elif i==1:
                    fileobj.write('({}, {})'.format(model[i], temp[beta.shape[0]-1-j])); fileobj.write(str('\t'))
        fileobj.write(str('\n \n'))

        fileobj.write(str("#acceptance rate\n"))
        for m in range(2*beta.shape[0]):
            for i in range(5):
                fileobj.write(str(pick[m, i])); fileobj.write(str("\t"))
            fileobj.write(str("\n"))
        fileobj.write(str("\n"))
        fileobj.write(str('#update \n'))
        for m in range(2*beta.shape[0]):
            for i in range(5):
                fileobj.write(str(update[m, i])); fileobj.write(str("\t"))
            fileobj.write(str("\n"))
        fileobj.write(str("\n"))

        fileobj.write(str("#exchange rate\n"))
        for m in range(2*beta.shape[0]-1):
            fileobj.write(str(exchange[m])); fileobj.write(str("\t"))
        fileobj.write(str("\n")); fileobj.write(str("\n"))

        fileobj.write(str('#De^M=0, De^M=1, ka, kd, sigma, energy, log bayes f\n'))
        for b in range(2*beta.shape[0]):
            if b < beta.shape[0]:
                fileobj.write(str("#M={}, beta={}\n".format(0, beta[b])))
            else:
                fileobj.write(str("#M={}, beta={}\n".format(1, beta[2*beta.shape[0]-1-b])))
            for i in range(T):
                fileobj.write(str(diffusivity_0[b, i])); fileobj.write(str("\t"))
                fileobj.write(str(diffusivity_1[b, i])); fileobj.write(str("\t"))
                fileobj.write(str(k_a[b, i])); fileobj.write(str("\t"))
                fileobj.write(str(k_d[b, i])); fileobj.write(str("\t"))
                fileobj.write(str(sigma[b, i])); fileobj.write(str("\t"))
                fileobj.write(str(E[b, i])); fileobj.write(str("\t"))
                fileobj.write(str(dlogp[b, i])); fileobj.write(str("\t"))
                fileobj.write(str("\n"))
            fileobj.write(str("\n"))
            
    return


for i in range(sample):
    w_file = './gakkai/EMbeta_TAP-sigma_datapoints={}_T={}mcs_Temp={}-{}_replica={}_koushinhaba-m1p1m1e-5p1e-5_{}_4.dat'.format(data_points, T, temp[0], temp[-1], temp.shape[0], data_file)
    
    diffusivity_0, diffusivity_1, k_a, k_d, sigma, E, dlogp, pick, update, exchange = exchange_monte_calro(T, beta, [0, 1])
    
    file_write(w_file, beta, np.array([0, 1]), diffusivity_0, diffusivity_1, k_a, k_d, sigma, E, dlogp, pick, update, exchange)
