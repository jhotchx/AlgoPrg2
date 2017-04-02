#Implement Strassen's Algorithm
#Implement Matrix Multiplication
#Find cross-over point analytically
#Find cross-over point experimentally

#Make array feed into a queue?

import math
import random

def generateMatrix(dimension, type):
    # choose from 0,1
    if type ==0:
        randomMatrix = [[random.randint(0,1) for e in range(dimension)] for e in range(dimension)]
    # choose from 0,1,2
    if type ==1:
        randomMatrix = [[random.randint(0,2) for e in range(dimension)] for e in range(dimension)]
    # choose from -1,0,1
    if type ==2:
        randomMatrix = [[random.randint(-1,1) for e in range(dimension)] for e in range(dimension)]
    # random real in [0,1]
    if type ==3:
        randomMatrix = [[random.random() for e in range(dimension)] for e in range(dimension)]	
    return(randomMatrix)

def make_matrix(d,A):
    mat = [None]*d
    index = 0
    for i in range(d):
        row = [None]*d
        for j in range(i,d):
            row[i]=A[index]
            index += 1
        mat[i] = row
    return(mat)

def process_inputfile(d,A):
    pass

def matmult(a,b):
   zip_b = zip(*b)
   # uncomment next line if python 3 : 
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

def add_mat(d,A,B):
    C = [[0 for j in range(d)] for i in range(d)]
    for i in range(d):
        for j in range(d):
            C[i][j] = A[i][j] + B[i][j]
    return C
    
def sub_mat(d,A,B):
    C = [[0 for j in range(d)] for i in range(d)]
    for i in range(d):
        for j in range(d):
            C[i][j] = A[i][j] - B[i][j]
    return C 

def compile_matrix(C11,C12,C21,C22):
    #Combine rows
    n = len(C11)
    for i in range(n):
        C11[i].extend(C12[i])
    for j in range(n):
        C21[j].extend(C22[j])
    return(C11 + C21)

@profile
def strassen(d,A,B,cutoff):
    # d = one dimension of matrix
    # A = dxd matrix
    # B = dxd matrix
    if d<=cutoff:
        return(matmult(A,B))
    else:
        chunksize = math.ceil(d/2)
        
        A11 = subset_matrix(A,1,1,chunksize)
        A12 = subset_matrix(A,1,2,chunksize)
        A21 = subset_matrix(A,2,1,chunksize)
        A22 = subset_matrix(A,2,2,chunksize)
        
        B11 = subset_matrix(B,1,1,chunksize)
        B12 = subset_matrix(B,1,2,chunksize)
        B21 = subset_matrix(B,2,1,chunksize)
        B22 = subset_matrix(B,2,2,chunksize)
        
        #Multiplications
        M1 = strassen(chunksize,add_mat(chunksize,A11,A22),
                                add_mat(chunksize,B11,B22),cutoff)
        M2 = strassen(chunksize,add_mat(chunksize,A21,A22),
                                B11,cutoff)
        M3 = strassen(chunksize,A11,
                                sub_mat(chunksize,B12,B22),cutoff)
        M4 = strassen(chunksize,A22,
                                sub_mat(chunksize,B21,B11),cutoff)
        M5 = strassen(chunksize,add_mat(chunksize,A11,A12),
                                B22,cutoff)
        M6 = strassen(chunksize,sub_mat(chunksize,A21,A11),
                                add_mat(chunksize,B11,B12),cutoff)
        M7 = strassen(chunksize,sub_mat(chunksize,A12,A22),
                                add_mat(chunksize,B21,B22),cutoff)
        
        #Build return matrix
        C11 = add_mat(chunksize,sub_mat(chunksize,add_mat(chunksize,M1,M4),
                                                  M5),
                                M7)
        C12 = add_mat(chunksize,M3,M5)
        C21 = add_mat(chunksize,M2,M4)
        C22 = add_mat(chunksize,add_mat(chunksize,sub_mat(chunksize,M1,M2),
                                                  M3),
                                M6)
        
    return(compile_matrix(C11,C12,C21,C22))

def wrapstras(d,A,B,cutoff):
    finalmat = strassen(d,A,B,cutoff)
    if d%2!=0:
        for i in range(d):
            finalmat[i] = finalmat[i][:-1]
        finalmat.remove([0]*(d+1))
    return(finalmat)

n = 100
TEST  = generateMatrix(n,0)
TEST2 = generateMatrix(n,0)
wrapstras(n,TEST,TEST2,3)

#AMat = [[1, 0, 0, 1], [1, 0, 0, 0], [0, 0, 0, 1], [1, 0, 1, 1]]
#BMat = [[0, 1, 1, 0], [1, 1, 1, 1], [1, 0, 0, 0], [1, 1, 0, 1]]
#strassen(4,AMat,BMat,3)

    
    