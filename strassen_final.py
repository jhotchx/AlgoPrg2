import random
from math import ceil,log
import sys

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

#@profile
def strassen(d,AMat,BMat,cutOff=3):
    if d <=cutOff:
        return(matmult(AMat,BMat))
    else:
        chunksize = ceil(d/2)

        A = subset_matrix(AMat,1,1,chunksize)
        B = subset_matrix(AMat,1,2,chunksize)
        C = subset_matrix(AMat,2,1,chunksize)
        D = subset_matrix(AMat,2,2,chunksize)
        
        E = subset_matrix(BMat,1,1,chunksize)
        F = subset_matrix(BMat,1,2,chunksize)
        G = subset_matrix(BMat,2,1,chunksize)
        H = subset_matrix(BMat,2,2,chunksize)

        # P1 = A(F-H)
        P1 = strassen(chunksize,A,subtract(F,H))
        # P2 = (A+B)*H
        P2 = strassen(chunksize,add(A,B),H)
        # P3 = (C+D)*E
        P3 = strassen(chunksize,add(C,D),E)
        #P4 = D*(G-E)
        P4 = strassen(chunksize,D,subtract(G,E))
        #P5 = (A+D)(E+H)
        P5 = strassen(chunksize,add(A,D),add(E,H))
        #P6 = (B-D)(G+H)
        P6 = strassen(chunksize,subtract(B,D),add(G,H))
        #P7 = (A-C)(E+F)
        P7 = strassen(chunksize,subtract(A,C),add(E,F))
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
    zip_b = list(zip_b)
    return [[sum(ele_a*ele_b for ele_a, ele_b in zip(row_a, col_b)) 
             for col_b in zip_b] for row_a in a]

def subset_matrix(M,x,y,d):
    padr = False
    padc = False
    if len(M)%d!=0:
        if x==2:
            padr = True
        if y==2:
            padc = True
        
    if x==1:
        rs = 0
        re = d
    elif x==2:
        rs = d
        if padr==True:
            re = 2*d - 1
        else:
            re = 2*d
        
    if y==1:
        cs = 0
        ce = d
    elif y==2:
        cs = d
        if padc==True:
            ce = 2*d - 1
        else:
            ce = 2*d
        
    mat = [None]*d
    j = 0
    for i in range(rs,re):
        row = M[i][cs:ce]
        if padc==True:
            row.append(0)
        mat[j] = row
        j += 1
    if padr==True:
        mat[j] = [0]*d
    return(mat)


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

def wrapstras(d,A,B,cutoff):
    finalmat = strassen(d,A,B,cutoff)
    if d%2!=0:
        for i in range(d):
            finalmat[i] = finalmat[i][:-1]
        finalmat.remove([0]*(d+1))
    return(finalmat)

if __name__ == "__main__":
    dimension = int(sys.argv[1])
    cutoff = int(sys.argv[2])
    kind = int(sys.argv[3])
    A = generateMatrix(dimension, kind)
    B = generateMatrix(dimension, kind)
    print(wrapstras(dimension,A,B,cutoff))
