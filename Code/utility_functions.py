import numpy as np


def Matrix2String01_LG(_Matrix):
    '''Convert 1-D string representation of Logic Gate Matrix into numpy array for fast matrix manipulations'''
    OutString01 = ''
    for i in range(0, _Matrix.shape[0]):
        for j in range(0, _Matrix.shape[1]):
            OutString01 = OutString01 + str(_Matrix[i][j])
    return OutString01


def LogicGatesString2Matrix(string):
    '''Convert 1-D string representation of Logic Gate Matrix into numpy array for fast matrix manipulations'''
    outmatrix = np.random.randint(3, 4, (int(len(string)/2), 2))
    stringindex = 0
    for i in range(0, outmatrix.shape[0]):
        for j in range(0, outmatrix.shape[1]):
            outmatrix[i][j] = int(string[stringindex])
            stringindex = stringindex + 1
    return outmatrix


def GetCorespondingMatrix(List012):
    '''Convert 1-D array representation of Weighted Adjacency Matrix into numpy array for fast matrix manipulations'''
    SQRTLEN = int(np.sqrt(len(List012)))
    OutMatrix = np.random.randint(8, 9, (2, SQRTLEN, SQRTLEN))
    for i in range(0, len(List012)):
        if List012[i] == 0:
            OutMatrix[0][i//SQRTLEN][i % SQRTLEN] = 0
            OutMatrix[1][i//SQRTLEN][i % SQRTLEN] = 0
        elif List012[i] == 1:
            OutMatrix[0][i//SQRTLEN][i % SQRTLEN] = 1
            OutMatrix[1][i//SQRTLEN][i % SQRTLEN] = 0
        elif List012[i] == 2:
            OutMatrix[0][i//SQRTLEN][i % SQRTLEN] = 0
            OutMatrix[1][i//SQRTLEN][i % SQRTLEN] = 1
        else:
            raise Exception('List012 must just have 012!\n{}'.format(List012))
    return OutMatrix


def String012ToMatrix(List012):
    '''Convert 1-D string representation of Weighted Adjacency Matrix into numpy array for fast matrix manipulations'''
    # String '012' to Matrix
    SQRTLEN = int(np.sqrt(len(List012)))
    OutMatrix = np.random.randint(8, 9, (2, SQRTLEN, SQRTLEN))
    for i in range(0, len(List012)):
        if List012[i] == '0':
            OutMatrix[0][i//SQRTLEN][i % SQRTLEN] = 0
            OutMatrix[1][i//SQRTLEN][i % SQRTLEN] = 0
        elif List012[i] == '1':
            OutMatrix[0][i//SQRTLEN][i % SQRTLEN] = 1
            OutMatrix[1][i//SQRTLEN][i % SQRTLEN] = 0
        elif List012[i] == '2':
            OutMatrix[0][i//SQRTLEN][i % SQRTLEN] = 0
            OutMatrix[1][i//SQRTLEN][i % SQRTLEN] = 1
        else:
            raise Exception('List012 must just have 012!\n{}'.format(List012))
    return OutMatrix


def GetCorespondingMatrix_LG(List01):
    '''Convert 1-D array representation of Logic Gate Matrix into numpy array for fast matrix manipulations'''
    SQRTLEN = int(len(List01)/2)
    OutMatrix = np.random.randint(8, 9, (SQRTLEN, 2))
    for i in range(0, len(List01)):
        OutMatrix[i//2][i % 2] = int(List01[i])
    return OutMatrix


def ConfigurationTo012(ConfigurationMatrix):
    '''Convert configuration matrix to string'''
    OutString = ''
    for j in range(0, ConfigurationMatrix.shape[1]):
        for z in range(0, ConfigurationMatrix.shape[2]):
            if ConfigurationMatrix[0][j][z] == 0 and ConfigurationMatrix[1][j][z] == 0:
                OutString = OutString + '0'
            elif ConfigurationMatrix[0][j][z] == 1 and ConfigurationMatrix[1][j][z] == 0:
                OutString = OutString + '1'
            elif ConfigurationMatrix[0][j][z] == 0 and ConfigurationMatrix[1][j][z] == 1:
                OutString = OutString + '2'
            else:
                raise Exception('Configuration (1,1)!')
    return OutString


def maskmean(array_):
    '''mean of array_ with entries called 'mask' removed'''
    out_mean = []
    for each_array in array_.T:
        temp_mean = []
        for each in each_array[0]:
            if each == 'mask':
                pass
            else:
                temp_mean.append(float(each))
        out_mean.append(np.nanmean(temp_mean))
    return out_mean
