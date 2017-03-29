import random
from math import ceil,log

def generateMatrix(dimension, kind):
    # choose from 0,1
    if kind ==0:
        randomMatrix = [[random.randint(0,1) for e in range(dimension)] for e in range(dimension)]
    # choose from 0,1,2
    if kind ==1:
        randomMatrix = [[random.randint(0,2) for e in range(dimension)] for e in range(dimension)]
    # choose from -1,0,1
    if kind ==2:
        randomMatrix = [[random.randint(-1,1) for e in range(dimension)] for e in range(dimension)]
    # random real in [0,1]
    if kind ==3:
        randomMatrix = [[random.random() for e in range(dimension)] for e in range(dimension)]  
    return(randomMatrix)

@profile
def strassen(AMat,BMat,cutOff=3):
    n = len(AMat)
    if n <=cutOff:
        return(matmult(AMat,BMat))
    else:
        A,B,C,D = subset(AMat)
        E,F,G,H = subset(BMat)
        # P1 = A(F-H)
        P1 = strassen(A,subtract(F,H))
        # P2 = (A+B)*H
        P2 = strassen(add(A,B),H)
        # P3 = (C+D)*E
        P3 = strassen(add(C,D),E)
        #P4 = D*(G-E)
        P4 = strassen(D,subtract(G,E))
        #P5 = (A+D)(E+H)
        P5 = strassen(add(A,D),add(E,H))
        #P6 = (B-D)(G+H)
        P6 = strassen(subtract(B,D),add(G,H))
        #P7 = (A-C)(E+F)
        P7 = strassen(subtract(A,C),add(E,F))
        # add terms to get final result
        #C1 = AE+BG = P4+P5-P2+P6
        C1 = add(subtract(add(P4,P5),P2),P6)
        #C2 = AF+BH = P1+P2
        C2 = add(P1,P2)
        # C3 = CE+DG = P3+P4
        C3 = add(P3,P4)
        # C4 = CF+DH = P5+P1-P3-P7
        C4 = subtract(subtract(add(P5,P1),P3),P7)
    #out = [C1[0]+C2[0]]+[C1[1]+C2[1]]+[C3[0]+C4[0]]+[C3[1]+C4[1]]
    out = compile_matrix(C1,C2,C3,C4)
    return(out)

def matmult(a,b):
    zip_b = zip(*b)
    # uncomment next line if python 3 : 
    zip_b = list(zip_b)
    return [[sum(ele_a*ele_b for ele_a, ele_b in zip(row_a, col_b)) 
             for col_b in zip_b] for row_a in a]

def subset(A):
    nextPowerOfTwo = lambda n: 2**int(ceil(log(n,2)))
    n = len(A)/2
    m = nextPowerOfTwo(n)
    A1 = [[0 for i in range(m)] for j in range(m)]
    A2 = [[0 for i in range(m)] for j in range(m)]
    A3 = [[0 for i in range(m)] for j in range(m)]
    A4 = [[0 for i in range(m)] for j in range(m)]
    n = int(n)
    for i in range(n):
        A1[i][0:n]=(A[i][0:n])
    for i in range(0,n):
        A2[i][0:n]=(A[i][n:n*2])
    for i in range(n):
        A3[i][0:n]=(A[i+n][0:n])
    for i in range(n):
        A4[i][0:n]=(A[i+n][n:n*2])
    return(A1,A2,A3,A4)


def add(A,B):
    t = []
    for i in range(len(A)):
        t.append([])
        for j in range(len(A)):
            t[i].append(A[i][j]+B[i][j])
    return(t)

def subtract(A,B):
    t = []
    for i in range(len(A)):
        t.append([])
        for j in range(len(A)):
            t[i].append(A[i][j]-B[i][j])
    return(t)

def compile_matrix(C11,C12,C21,C22):
   #Combine rows
   n = len(C11)
   for i in range(n):
       C11[i].extend(C12[i])
   for j in range(n):
       C21[j].extend(C22[j])
   return(C11 + C21)

if __name__ == "__main__":
    dimension = int(input("Enter the dimension for the matrices: "))
    kind = int(input("what kind of numbers would you like: "))
    cutoff = int(input("When should we switch to the normal algorithm: "))
    A = generateMatrix(dimension, kind)
    B = generateMatrix(dimension, kind)
    strassen(A,B,cutoff)

