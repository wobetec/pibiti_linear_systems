from matrix import *
from math import sqrt, fabs
from copy import deepcopy
from functools import reduce
from random import randint

sis_dir = "./sistemas/"

class Sistema():

    def __init__(self, A, b):
        self.A = A
        self.b = b


    def open(file_name):
        filename = sis_dir + file_name

        A, b = [], []
        with open(filename, "r") as f:

            for line in f:
                linha = [int(x) for x in line.split()]
                n = len(linha)
                A.append(linha[:n-1])
                b.append(linha[n-1])

        return Sistema(A, b)


    def create(n, lim=100):
        A, b = [], []

        for i in range(n):
            eq = [randint(0, lim) for j in range(n)]
            A.append(eq)
            b.append(randint(0, lim))

        return Sistema(A, b)


    def save(self, file_name):
        filename = sis_dir + file_name
        print(filename)
        with open(filename, "w") as f:

            for i in range(len(self.b)):
                f.write(" ".join(str(a) for a in self.A[i]))
                f.write(f" {self.b[i]}")
                f.write("\n")
                print(i)


    def creation(file_name, n, lim=100):
        S = Sistema.create(n, lim=lim)
        S.show()
        try:
            S.save(file_name)
        except:
            print("Generation Erro")


    def show(self):
        for i in range(len(self.A)):
            for j in self.A[i]:
                print(f"{j}".rjust(5, " "), end="")
            print(f"{self.b[i]}".rjust(5, " "))


    def gauss_jordan(self):
        i, p = 0, 1
        
        n, A, b = len(self.b), deepcopy(self.A), self.b.copy()
        x=[]

        while i<n and p!=0:
            #PIVOT4 --> garante que o pivo n seja zero e facilita as contas

            c = A[i][i]
            l = i
            for k in range(i+1, n): 
                if fabs(c) < fabs(A[k][i]):
                    c = A[k][i]
                    l = k
            
            if l != i:
                for j in range(0, n):
                    A[i][j], A[l][j] = A[l][j], A[i][j]
                b[i], b[l] = b[l], b[i]

            #continação
            p = A[i][i]
            if p != 0: #verifica se o pivo n é zero
                for j in range(i, n):
                    A[i][j] = A[i][j]/p
                b[i] = b[i]/p
                for k in range(0,  n):
                    if k != i:
                        q = A[k][i]
                        for j in range(i, n):
                            A[k][j] = A[k][j] - q * A[i][j]
                        b[k] = b[k] - q * b[i]
                i+=1
        if p != 0:
            x = b.copy()
            return x
        else:
            print("Sistema impossível ou indeterminado")
            return False


    def gauss_jacobi(self, nmi, erro):
        d, k = 1, 0

        n, A, b = len(self.b), deepcopy(self.A), self.b.copy()

        xa = []
        xb = [b[i]/A[i][i] for i in range(n)]

        while d>0 and nmi >=k:
            d=0
            xa = xb.copy()
            
            for i in range(n):
                soma = b[i]
                for j in range(n):
                    if j != i:
                        soma -= A[i][j] * xa[j]
                xb[i] = soma/A[i][i]

            x = [xb[i] - xa[i] for i in range(n)]
            mod_x = sqrt(reduce(lambda a, b: a+b, [a*a for a in x]))
            mod_xb = sqrt(reduce(lambda a, b: a+b, [a*a for a in xb]))
            if mod_x/mod_xb > erro:
                d = 1

            k+=1

        print(k)

        if k > nmi:
            print("limite de repetições")
            return False
        else:
            return xb


    def gauss_seidel(self, nmi, erro):
        d, k = 1, 0

        n, A, b = len(self.b), deepcopy(self.A), self.b.copy()

        xa = []
        xb = [b[i]/A[i][i] for i in range(n)]

        while d>0 and nmi >= k:
            d=0
            xa = xb.copy()

            for i in range(len(self.b)):
                a = self.b[i]
                for j in range(len(self.b)):
                    if j > i:
                        a -= self.A[i][j] * xa[j]
                    elif j < i:
                        a -= self.A[i][j] * xb[j]

                xb[i] = a/A[i][i]

            x = [xb[i] - xa[i] for i in range(n)]
            mod_x = sqrt(reduce(lambda a, b: a+b, [a*a for a in x]))
            mod_xb = sqrt(reduce(lambda a, b: a+b, [a*a for a in xb]))
            if mod_x/mod_xb > erro:
                d = 1

            k+=1

        print(k)

        if k > nmi:
            print("limite de repetições")
            return False
        else:
            return xb


if __name__ == "__main__":
    pass

