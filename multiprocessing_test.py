from multiprocessing import Pool
import numpy as np

# def nijyou(x):
#     print(np.power(x,2))
#     return None
#
# if __name__ == "__main__":
#     p = Pool(4)
#     p.map(nijyou, range(1000))

def kakezan(a , b):
    return a*b

def wrapper_kakezan(args):
    return kakezan(*args)

if __name__ == "__main__":
    tutumimono = [[i, 3] for i in range(100000)]
    p = Pool(processes=2)
    print( p.map(wrapper_kakezan, tutumimono) )
