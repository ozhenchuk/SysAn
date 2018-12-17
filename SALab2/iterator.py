import numpy as np
from copy import deepcopy

from polinoms import Chebyshev, Legender, Lagger, Hermit
from equation import Conjugate_Grad

class Iterator:
    def __init__(self, N, n, m, filename, polinom_type, p, lambda_flag):
        self.eps = 1e-10
        self.p = p
        self.p_type = polinom_type
        self.p_flag = True
        if self.p_type in [2,3]:
            self.p_flag = False
        self.lambda_flag = lambda_flag
        self.infile = filename
        file_in = open(filename, 'r')
        ##fout_inp = open('fout_inp.txt','w+')
        self.N = N
        ##fout_inp.write(str(self.N)+" "+str(range(self.N)))
        self.n = []
        for i in range(3):
            self.n.append(n[i])
        self.x = [0.]*3
        for i in range(3):
            self.x[i]=[0.]*self.n[i]
            for j in range(n[i]):
                self.x[i][j] = [0.]*self.N
        for i in range(self.N):
            line = file_in.readline().split(" ")
            count = 0
            for j in range(3):
                for k in range(self.n[j]):
                    self.x[j][k][i] = float(line[count])
                    ##fout_inp.write(str(self.x[j][k][i]) + " ")
                    count += 1 
            ##fout_inp.write(str(i)+"\n")
        self.min_x = np.zeros((3, max(self.n)))
        self.max_x = np.zeros((3, max(self.n)))
        for i in range(3):
            for j in range(self.n[i]):
                self.max_x[i][j] = max(self.x[i][j])
                self.min_x[i][j] = min(self.x[i][j])
        self.m = m
        self.y = np.zeros((self.m, self.N), "f") 
        for i in range(N):
            line = file_in.readline().split()
            for j in range(m):
                self.y[j][i] = float(line[j])
                ##fout_inp.write(str(self.y[j][i]) + " ")
            ##fout_inp.write("\n")
        file_in.close() 
        self.y_cnt = deepcopy(self.y)
        self.min_y = np.zeros(self.m)
        self.max_y = np.zeros(self.m) 
        for i in range(self.m):
            self.max_y[i] = max(self.y[i])
            self.min_y[i] = min(self.y[i]) 
        if polinom_type == 0:
            self.polinom = Chebyshev
        elif polinom_type == 1:
            self.polinom = Legender
        elif polinom_type == 2:
            self.polinom = Lagger
        elif polinom_type == 3:
            self.polinom = Hermit
        ##fout_inp.close()

    def find_Bq(self, index = 0):
        bq = deepcopy(self.y[index])
        return np.matrix(bq)

    def find_Psi3(self, p, index = 0):
        F = np.zeros((self.N, len(self.x[0])*(p[0] + 1) + len(self.x[1])*(p[1] + 1) + len(self.x[2])*(p[2] + 1)))
        for i in range(3):
            for j in range(self.n[i]):
                for k in range(p[i] + 1):
                    for l in range(self.N):
                        tmp_index = j*(p[i] + 1)
                        for r in range(i):
                            tmp_index += (p[r] + 1)*len(self.x[r])
                        F[l][k + tmp_index] = self.polinom(k, self.x[i][j][l])
                        ##F[l][tmp_index] = 1.
        bq = self.find_Bq(index) 
        A = np.dot(F.T, F)
        ##print('DIM A = ', A)
        b = (np.dot(F.T, bq.T)).T
        ##print('DIM b = ', b)
        x01 = np.ones((1, len(self.x[0])*(p[0] + 1) + len(self.x[1])*(p[1] + 1) + len(self.x[2])*(p[2] + 1)))
        ##print('DIM x01 = ', x01)
        lambda_coef = Conjugate_Grad(A, b, x01, 10, self.p_flag)
        ##print('DIM lambda_coef = ', lambda_coef)
        return [F, lambda_coef]
    def find_Psi1(self, x, p, index = 0):
        F = np.zeros((self.N, len(x)*(p + 1)))
        for j in range(len(x)):
            for k in range(p + 1):
                for l in range(self.N):
                    tmp_index = j*(p + 1)
                    F[l][k + tmp_index] = self.polinom(k, x[j][l])
                    ##F[l][tmp_index] = 1.
        bq = (self.find_Bq(index))/3 
        A = np.dot(F.T, F)
        ##print('DIM A = ', A)
        b = (np.dot(F.T, bq.T)).T
        ##print('DIM b = ', b)
        x01 = np.ones((1, len(x)*(p + 1)))
        ##print('DIM x01 = ', x01)
        lambda_coef = Conjugate_Grad(A, b, x01, 10, self.p_flag)
        ##print('DIM lambda_coef = ', lambda_coef)
        return [F, lambda_coef]
 
    def find_Ef(self, index = 0):
        print('Lambda_flag ',self.lambda_flag)
        if (self.lambda_flag == 0):
            set_a = [0.]*3
            m = 0
            set_lambda = [0.]*3
            for i in range(3):
                m += self.n[i]
            for i in range(3):
                calc = np.zeros((self.N, self.n[i]))
                tmp = self.find_Psi1(self.x[i], self.p[i], index)
                set_lambda[i] = deepcopy(tmp[1][0])
                for j in range(len(tmp[0])):
                    for l in range(self.n[i]):
                        for k in range(self.p[i] + 1):
                            calc[j][l] += tmp[0][j][l*(self.p[i] + 1) + k]*tmp[1][0][l*(self.p[i] + 1) + k]
                A = np.dot(calc.T, calc) 
                b = ((np.dot(calc.T, self.y[index].T)).T)/3
                set_a[i] = [0.]*self.n[i]
                x01 = np.ones((1, self.n[i]))
                b01 = []
                b01.append(b)
                b1 = np.array(b01)
                set_a[i] = Conjugate_Grad(A, b1, x01, 10, self.p_flag)
            tmp_ind = 0
            a = np.zeros((1, m), "f")
            for i in range(len(set_a)):
                for j in range(len(set_a[i][0])):
                    a[0][tmp_ind + j] = set_a[i][0][j]
                tmp_ind += len(set_a[i][0])
            print('a = ', a)
            tmp = self.find_Psi3(self.p, index)
            calc_com = np.zeros((self.N, m))
            tmp_ind = 0 
            t_i = 0
            for i in range(3): 
                for j in range(len(tmp[0])):
                    for l in range(self.n[i]):
                        for k in range(self.p[i] + 1):
                            calc_com[j][l + t_i] += tmp[0][j][l*(self.p[i] + 1) + tmp_ind + k]*tmp[1][0][l*(self.p[i] + 1) + tmp_ind + k]
                t_i += self.n[i]
                tmp_ind += (self.p[i] + 1)*self.n[i]
            lambda_ = np.zeros((len(set_lambda[0])+len(set_lambda[2])+len(set_lambda[1])),"f")
            t_i = 0
            for i in range(3):
                for j in range(len(set_lambda[i])):
                    lambda_[j + t_i] = set_lambda[i][j]
                t_i += len(set_lambda[i])
            return [a, calc_com, lambda_]
        else:
            m = 0
            for i in range(3):
                m += self.n[i]
            calc_com = np.zeros((self.N, m))
            set_a = [0.]*m
            tmp = self.find_Psi3(self.p, index)
            tmp_ind = 0 
            t_i = 0 
            for i in range(3): 
                for j in range(len(tmp[0])):
                    for l in range(self.n[i]):
                        for k in range(self.p[i] + 1):
                            calc_com[j][l + t_i] += tmp[0][j][l*(self.p[i] + 1) + tmp_ind + k]*tmp[1][0][l*(self.p[i] + 1) + tmp_ind + k]
                t_i += self.n[i]
                tmp_ind += (self.p[i] + 1)*self.n[i]
            A = np.dot(calc_com.T, calc_com) 
            b = (np.dot(calc_com.T, self.y[index].T)).T
            x01 = np.ones((1, m))
            b01 = []
            b01.append(b)
            b1 = np.array(b01)
            set_a = Conjugate_Grad(A, b1, x01, 10, self.p_flag)
        return [set_a, calc_com, tmp[1][0]]
        
    def find_F(self, index):
        tmp = self.find_Ef(index)
        set_a = tmp[0]
        F = tmp[1]
        calc_com = np.zeros((self.N, 3))
        t_i = 0
        for i in range(3): 
            for j in range(self.N):
                for l in range(self.n[i]):
                    calc_com[j][i] += F[j][l + t_i]*set_a[0][l + t_i]
            t_i += self.n[i]
        A = np.dot(calc_com.T, calc_com)
        b = (np.dot(calc_com.T, self.y[index].T)).T
        x01 = np.ones((1, 3))
        b01 = []
        b01.append(b)
        b1 = np.array(b01)
        set_c = Conjugate_Grad(A, b1, x01, 10, self.p_flag, True)
        print(f'Results {index}:')
        for i in range(len(self.y[index])):
            calc_line1 = []
            calc_line1.append(calc_com[i])
            print(self.y[index][i], ' ', np.dot(calc_line1,set_c.T), ' ')
        print('Searching c =', set_c)
        self.y_cnt[index] = np.dot(calc_com, set_c.T).T
        return [set_c, set_a, tmp[2]]

    def normalization(self):
        for i in range(3):
            for j in range(self.n[i]): 
                max_ = max(self.x[i][j])
                min_ = min(self.x[i][j])
                if (np.fabs(max_ - min_) < self.eps):
                    for k in range(self.N):
                        self.x[i][j][k] = 0
                else:
                    for k in range(self.N):
                        self.x[i][j][k] = (self.x[i][j][k]-min_)/(max_-min_)
        for i in range(self.m):
            max_ = max(self.y[i])
            min_ = min(self.y[i])
            if (np.fabs(max_-min_) < self.eps):
                for j in range(self.N):
                    self.y[i][j] = 0
            else:
                for j in range(self.N):
                    self.y[i][j] = (self.y[i][j]-min_)/(max_-min_)

    def denormalization(self):
        for i in range(3):
            for j in range(self.n[i]):
                for k in range(self.N):
                    self.x[i][j][k] = self.x[i][j][k]*(self.max_x[i][j]-self.min_x[i][j]) + self.min_x[i][j]
        for i in range(self.m):
            for j in range(self.N):
                self.y[i][j] = self.y[i][j]*(self.max_y[i]-self.min_y[i]) + self.min_y[i]
                self.y_cnt[i][j] = self.y_cnt[i][j]*(self.max_y[i]-self.min_y[i]) + self.min_y[i]

    def approximate(self, filename):
        file_out = open(filename, 'w')
        output = [] 
        file_out.write(f'Functional dependency restoration for {self.m} functions and {self.N} samples\n')
        polinomtype = ['Chebyshev', 'Legender', 'Lagger', 'Hermit']
        file_out.write(f'Type of Polinoms used: {polinomtype[self.p_type]} polinoms\n')
        file_out.write(f'Degrees:')
        for i in range(3):
            file_out.write(f' {self.p[i]}')
        file_out.write('\n\n')
        for i in range(self.m):
            output.append(self.find_F(i))
            set_c = output[i][0][0]
            set_a = output[i][1]
            lambda_ = output[i][2]
            file_out.write(f'---------------------------------Function Y{i+1}----------------------------------\n\n')
            tmp_ind = 0
            t_i = 0
            for j in range(3):
                for o in range(self.n[j]):
                    file_out.write(f'Lambda for X{j+1}{o+1} is:\n')
                    for k in range(self.p[j]+1):
                        file_out.write(' '+str(lambda_[o*(self.p[j] + 1) + tmp_ind + k]))
                    file_out.write('\n')
                t_i += self.n[j] 
                tmp_ind += self.n[j]*(self.p[j]+1)
            file_out.write('\n')
            t_i = 0
            for j in range(3):
                file_out.write('A for X'+str(j+1)+' is:\n')
                for o in range(self.n[j]):
                    file_out.write(' '+str(set_a[0][o + t_i]))
                t_i += self.n[j] 
                file_out.write('\n')
            file_out.write('\n')
            for j in range(3):
                file_out.write('C['+str(j+1)+'] is:\n')
                file_out.write(' '+str(set_c[j])+'\n')
            file_out.write('\n')
            file_out.write('F = ')
            for j in range(3):
                if(set_c[j]>0 and j>0):
                    file_out.write('+')
                file_out.write(' '+str(set_c[j])+'*Φ['+str(i+1)+'](X'+str(j+1)+') \n')
            file_out.write('\n')
            t_i = 0
            for j in range(3):
                file_out.write('Φ['+str(j)+'] = ')
                for o in range(self.n[j]):
                    if(set_a[0][o + t_i]>0 and o>0):
                        file_out.write('+')
                    file_out.write(' '+str(set_a[0][o + t_i])+'*Ψ['+str(j+1)+''+str(o+1)+'] ')
                    file_out.write('\n')
                t_i += self.n[j] 
                file_out.write('\n ')
            tmp_ind = 0
            t_i = 0
            for j in range(3):
                for o in range(self.n[j]):
                    file_out.write('Ψ['+str(j+1)+str(o+1)+'] = ')
                    for k in range(self.p[j]+1):
                        if(lambda_[o*(self.p[j] + 1) + tmp_ind + k]>0 and k>0):
                            file_out.write('+')
                        file_out.write(' '+str(lambda_[o*(self.p[j] + 1) + tmp_ind + k])+'*T'+str(k)+'(x'+str(j+1)+''+str(o+1)+') ')
                        file_out.write('\n')
                t_i += self.n[j] 
                tmp_ind += self.n[j]*(self.p[j]+1)
            file_out.write('\n')
            tmp_ind = 0
            t_i = 0
            file_out.write('Full expression from polinoms:\n')
            file_out.write('F = ')
            for j in range(3):
                for o in range(self.n[j]):
                    for k in range(self.p[j]+1):
                        if(lambda_[o*(self.p[j] + 1) + tmp_ind + k]*set_a[0][o + t_i]*set_c[j]>0 and j+o+k>0):
                            file_out.write('+')
                        file_out.write(' '+str(lambda_[o*(self.p[j] + 1) + tmp_ind + k]*set_a[0][o + t_i]*set_c[j])+'*T'+str(k)+'(x'+str(j+1)+''+str(o+1)+') ')
                file_out.write('\n')
                t_i += self.n[j] 
                tmp_ind += self.n[j]*(self.p[j]+1)
            file_out.write('\n')
            file_out.write('Function samples and their approximation:\n\n')
            file_out.write('Y \t Y_app \t Residual \n')
            for j in range(self.N):
                file_out.write(str(self.y[i][j]) + ' \t ' + str(self.y_cnt[i][j]) + ' \t ' + str(np.fabs(self.y[i][j] - self.y_cnt[i][j])) +'\n')
            file_out.write('\n')
        file_out.close()
