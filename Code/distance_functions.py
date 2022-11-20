from utility_functions import *

def GetAttractorDistance(mRNA1, WT_mRNA, transcriptionalprofilemax):
    '''Compute Hamming distances for matrices with flattening (for logic gate matrices)'''
    outdistance = 0.0
    for y in range(0, len(WT_mRNA)):
        outdistance = outdistance + abs((mRNA1[y]-WT_mRNA[y])/transcriptionalprofilemax[y])
    return outdistance

def GetHammingDistance(_MatrixA, _MatrixB):
    '''Compute Hamming distances for matrices without flattening (for adjacency matrices)'''
    HammingDistance = 0
    if _MatrixA.shape == _MatrixB.shape:
        for i in range(0, _MatrixA.shape[0]):
            for j in range(0, _MatrixA.shape[1]):
                for z in range(0, _MatrixA.shape[2]):
                    if _MatrixA[i,j,z] != _MatrixB[i,j,z]:
                        HammingDistance = HammingDistance + 1
                    else:
                        pass
    else:
        raise Exception('Shapes don\'t match!')
    return HammingDistance

def GetHammingDistance_LG(_MA, _MB):
    '''Compute Hamming distances for matrices with flattening (for logic gate matrices)'''
    HammingDis = 0
    _SA = Matrix2String01_LG(_MA)
    _SB = Matrix2String01_LG(_MB)
    for i in range(0, len(_SA)):
        if _SA[i] == _SB[i]:
            pass
        else:
            HammingDis = HammingDis + 1
    return HammingDis
