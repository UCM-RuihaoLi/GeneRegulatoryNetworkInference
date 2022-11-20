import copy
from scipy.integrate import solve_ivp
import numpy as np


def ActivatorSigmoid(x, t, k):
    '''Activating Hill function: threshold t, Hill coefficient k'''
    return (x**k / (x**k + t**k))


def RepressorSigmoid(x, t, k):
    '''Repressing Hill function: threshold t, Hill coefficient k'''
    return (t**k / (x**k + t**k))


def GetAttractorDistance(mRNA1, WT_mRNA, transcriptionalprofilemax):
    '''Compute normalized distance between an two attractors'''
    outdistance = 0.0
    for y in range(0, len(WT_mRNA)):
        outdistance = outdistance + \
            abs((mRNA1[y]-WT_mRNA[y])/transcriptionalprofilemax[y])
    return outdistance


def Run_Dynamics(GRN_instance_ori, index_i_func, index_j_func, 
                 TrainningCount_func, WTTP, mRNAstate, TRMAX, return_dict):
    '''Simulate dynamics'''
    t_ticks = [int(TrainningCount_func/250)*tick for tick in range(0, 250+1)]
    t_ticks = None
    GRN_instance = copy.deepcopy(GRN_instance_ori)
    GRN_instance.SetmRNA(
        GRN_instance_ori.mRNA[mRNAstate], WTTP[str(index_i_func)][0], TRMAX)
    GRN_instance.SetProtein(
        GRN_instance_ori.Protein[mRNAstate], WTTP[str(index_i_func)][0], TRMAX)
    scipystring = GRN_instance.Delta_mRNA(WTTP[str(index_i_func)][0], TRMAX)
    scipystring = (scipystring 
                   + '\nInitialStates = list(GRN_instance.mRNA)+list(GRN_instance.Protein)' 
                   + '\nsol = solve_ivp(update_mRNA_protein, [0, TrainningCount_func], InitialStates, '
                   +                    'args=([WTTP[str(index_i_func)][0]]), t_eval=t_ticks)')
    scipystring = scipystring + '\nGRN_instance.sol = sol.y'
    exec(scipystring)
    sol = GRN_instance.sol
    GRN_instance.sol = []
    protein_final = []
    for each in [x for x in range(len(GRN_instance.mRNA), len(GRN_instance.mRNA)*2)]:
        protein_final.append(sol[each][-1])
    GRN_instance.SetProtein(np.array(protein_final),
                            WTTP[str(index_i_func)][0], TRMAX)
    return_dict[index_j_func] = sol[:len(GRN_instance.mRNA)]
    sol = []
    GRN_instance = []
    return
