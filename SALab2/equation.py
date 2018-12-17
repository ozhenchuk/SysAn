import numpy as np
from copy import deepcopy

eps = 1e-5

def Conjugate_Grad(A, b, x0, iter=-1, flag_it=True, flag_c=False): ###b, x0 - vector-rows
    n = len(x0[0])
    print('A= ', A)
    print('b= ', b)
    print('x0= ', x0)
    print('n= ', n)
    x = x0.T
    r = np.dot(A,x) - b.T
    print('r= ', r)
    print('r.T= ', r.T)
    p = -r
    print('r.T*r= ', r.T*r)
    print('(r.T,r)= ', np.dot(r.T,r))
    r_k_norm = float(np.dot(r.T,r))
    print('r_k_norm= ', r_k_norm)
    for i in range(2*n):
        Ap = np.dot(A,p)
        alpha = r_k_norm / float(np.dot(p.T,Ap))
        x += alpha*p
        r += alpha*Ap
        r_kplus1_norm = float(np.dot(r.T,r))
        beta = r_kplus1_norm/r_k_norm
        r_k_norm = r_kplus1_norm
        res = np.linalg.norm(np.dot(A,x)-b.T)
        ###if r_kplus1_norm < eps:
        if res < eps:
            ##print('Itr:', i)
            ##print(i, ' ', res, ' IRRTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT')
            break
        p = beta * p - r
    ret = x.T
    print('Al-b= ', res)
    if iter == -1:
        iter = 100
    if(res > eps and iter > 0 and flag_it): ##iter > 0 and 
        ##print('IRRTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT')
        ret = Conjugate_Grad(A, b, ret, iter-1)
    if flag_c:
        for i in range(len(b[0])):
            mat_line1 = []
            mat_line1.append(A[i])
            print(b[0][i], ' ', np.dot(mat_line1,ret.T), ' ')
    return ret
