import random
from math import ceil,log
import sys

def process_inputfile(d,inputfile):
    A = [[0 for e in range(d)] for e in range(d)]
    B = [[0 for e in range(d)] for e in range(d)]
    i = 0
    j = 0
    M = A
    with open(inputfile) as f:
        for line in f:
            M[i][j] = float(line)
            if j==(d-1) and i==(d-1):
                M = B
                i = 0
                j = 0
            elif j==(d-1):
                j = 0
                i += 1
            else:
                j +=1
    return(A,B)   

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

def genFile(dimension,kind):
    with open('test','w') as outFile:
        A = generateMatrix(dimension,kind)
        B = generateMatrix(dimension,kind)
        aNumbers = [item for sublist in A for item in sublist]
        bNumbers = [item for sublist in B for item in sublist]
        for num in aNumbers:
            outFile.write(str(num)+'\n')
        for num in bNumbers:
            outFile.write(str(num)+'\n')

#@profile
def strassen(d,AMat,BMat,cutoff):
    if d <=cutoff:
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
        #P1 = A(F-H)
        P1 = strassen(chunksize,A,subtract(F,H),cutoff)
        # P2 = (A+B)*H
        P2 = strassen(chunksize,add(A,B),H,cutoff)
        # P3 = (C+D)*E
        P3 = strassen(chunksize,add(C,D),E,cutoff)
        #P4 = D*(G-E)
        P4 = strassen(chunksize,D,subtract(G,E),cutoff)
        #P5 = (A+D)(E+H)
        P5 = strassen(chunksize,add(A,D),add(E,H),cutoff)
        #P6 = (B-D)(G+H)
        P6 = strassen(chunksize,subtract(B,D),add(G,H),cutoff)
        #P7 = (A-C)(E+F)
        P7 = strassen(chunksize,subtract(A,C),add(E,F),cutoff)
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

#@profile
def matmult(AMat,BMat):
    B = list(zip(*BMat)) #Make row-wise matrix become column-wise
    return [[sum(Aij*Bjk for Aij, Bjk in zip(Ai, Bk)) for Bk in B] for Ai in AMat]

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

#@profile
def wrapstras(d,A,B,cutoff=161):
    finalmat = strassen(d,A,B,cutoff)
    if d%2!=0:
        try:
            finalmat.remove([0]*(2**int(ceil(log(d,2)))))
        except:
            for i in range(d):
                finalmat[i] = finalmat[i][:-1]
    return(finalmat)

def get_diag(d,C):
    diagvals = [None]*d
    for i in range(d):
        diagvals[i] = C[i][i]
    return(diagvals)


if __name__ == "__main__":
    dimension = int(sys.argv[2])
    inputfile = sys.argv[3]
    A,B = process_inputfile(dimension,'test')
    C = wrapstras(dimension,A,B)
    Diag = get_diag(dimension,C)
    for i in Diag:
        print(i)


