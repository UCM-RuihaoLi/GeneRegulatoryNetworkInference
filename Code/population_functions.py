import numpy as np
import math


def GetNextProportions(CurrentDistance_, CurrentProportions_):
    CurrentFitness = np.array(CurrentDistance_)  # Convert distance to fitness
    CurrentFitness = -CurrentFitness + np.max(CurrentFitness)
    ExcessFitness = CurrentFitness - np.mean(CurrentFitness)
    DeltaProportions = np.multiply(ExcessFitness, CurrentProportions_)
    NextProportions = CurrentProportions_ + DeltaProportions
    NextProportions = NextProportions + (-1*np.min(NextProportions))
    if sum(NextProportions) == 0:
        NextProportions = np.array([1/len(CurrentDistance_)
                                   for i in range(0, len(CurrentDistance_))])
    else:
        NextProportions = NextProportions / sum(NextProportions)
    return NextProportions


def NaturalSelection(GRNlist, Proportionlist, SelectionPower_):
    '''Delete the last 25% GRN and replicate the most prevailing one'''
    numbers = []
    for i in range(0, len(GRNlist)):
        namesplited = GRNlist[i].split('_')
        numbers.append(int(namesplited[1]))
    thebignumber = max(numbers)
    length = math.ceil(len(GRNlist)//(1/SelectionPower_))
    for i in range(0, length):
        deleteindex = np.argmin(Proportionlist)
        del(GRNlist[deleteindex])
        del(Proportionlist[deleteindex])
    for i in range(1, length+1):
        GRNlist.append('GRN_{}'.format(i+thebignumber))
    copyindex = np.argmax(Proportionlist)
    return GRNlist, length, copyindex


def GetLongJumpList(GRNlist, Proportionlist_, SelectionPower_):
    '''Return the indece of instances that are going to have a long jump (The last SelectionPower%)'''
    length = math.ceil(len(GRNlist)//(1/SelectionPower_))
    LongJumpIndexList = []
    for i in range(0, length):
        longjumpindex = np.argmin(Proportionlist_)
        LongJumpIndexList.append(longjumpindex)
        del(Proportionlist_[longjumpindex])
    return LongJumpIndexList


def InitiateGlobalMutationDirection(totalnumberofgenes):
    '''Helper funciton for intialization of mutation parameters (adjacency matrix)'''
    GlobalMutationDirection = {}
    for i in range(0, totalnumberofgenes**2):
        GlobalMutationDirection[str(i)] = [0, 0, 0]
    return GlobalMutationDirection


def InitiateGlobalMutationDirection_LG(totalnumberofgenes):
    '''Helper funciton for intialization of mutation parameters (logic gates)'''
    GlobalMutationDirection_LG = {}
    for i in range(0, totalnumberofgenes*4):
        GlobalMutationDirection_LG[str(i)] = [0, 0]
    return GlobalMutationDirection_LG
