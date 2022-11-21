import sys
import getopt
import os
import pandas as pd
import numpy as np
import random
import copy
import multiprocessing
import time
from scipy.integrate import solve_ivp


from dynamics import *
from distance_functions import *
from GRN import GRN
from mutation_functions import *
from population_functions import *
from utility_functions import *


if __name__ == '__main__':
    np.seterr(over='ignore')

    sys_path = os.getcwd()

    sys_input_RNAseq = ""
    sys_input_ChIP = ""
    sys_promoter_strengths = ""
    sys_protein_degradation_rate = ""
    sys_mRNA_elongation_rate = 4.8
    sys_aa_elongation_rate = 8
    sys_gene_length = []
    sys_samplesize = 100
    sys_training_count = 200
    sys_PerturbationPower = 0.1
    sys_iteration_num = 800
    sys_output_name = ""

    #############################
    ### PARSE INPUT ARGUMENTS ###
    #############################
    try:
        opts, args = getopt.getopt(sys.argv[1:],
                                   "hr:c::i:t:s::n::p::m::a::l:o:d::",
                                   ["help", "input_RNAseq", "input_ChIP",
                                    "iteration_num", "promoter_strengths",
                                    "samplesize", "training_count", "PerturbationPower",
                                    "mRNA_elongation_rate", "aminoacid_elongation_rate",
                                    "gene_length", "output_name", "protein_degradation"])
    except getopt.GetoptError:
        print('Usage: EGRNM.py [-h] -r <input_RNAseq> [-c <input_ChIP>] [-i <iteration_num>]'
              ' -t <promoter_strengths> [-s <samplesize>] [-n <training_count>]'
              ' -p <PerturbationPower> [-m <mRNA_elongation_rate>] [-a <aminoacid_elongation_rate>]'
              ' -l <gene_length> -o <output_name> -[d <protein_degradation>]')
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print('Usage: EGRNM.py [-h] -r <input_RNAseq> [-c <input_ChIP>] [-i <iteration_num>]'
                  ' -t <promoter_strengths> [-s <samplesize>] [-n <training_count>]'
                  ' -p <PerturbationPower> [-m <mRNA_elongation_rate>] [-a <aminoacid_elongation_rate>]'
                  ' -l <gene_length> -o <output_name> -[d <protein_degradation>]')
            sys.exit()
        elif opt in ("-r", "--input_RNAseq"):
            sys_input_RNAseq = sys_path + '/' + arg
            sys_input_RNAseq = np.array(pd.read_csv(
                sys_input_RNAseq, header=None, delimiter='\t'))
            sys_WTTP = {}
            for i in range(0, sys_input_RNAseq.shape[0]):
                if sys_input_RNAseq[i][0] == '-1':
                    sys_WTTP[str(i)] = [[], sys_input_RNAseq[i][1:]]
                elif len(sys_input_RNAseq[i][0]) == 1:
                    sys_WTTP[str(i)] = [[int(sys_input_RNAseq[i][0])], sys_input_RNAseq[i][1:]]
                else:
                    sys_WTTP[str(i)] = [np.array(sys_input_RNAseq[i][0].split(','), dtype = int).tolist(), sys_input_RNAseq[i][1:]]
        elif opt in ("-c", "--input_ChIP"):
            sys_input_ChIP = sys_path + '/' + arg
            sys_input_ChIP = np.array(pd.read_csv(
                sys_input_ChIP, header=None, delimiter='\t'), dtype=int)
        elif opt in ("-i", "--iteration_num"):
            sys_iteration_num = arg
        elif opt in ("-m", "--mRNA_elongation_rate"):
            sys_mRNA_elongation_rate = arg
        elif opt in ("-a", "--aminoacid_elongation_rate"):
            sys_aa_elongation_rate = arg
        elif opt in ("-l", "--gene_length"):
            sys_gene_length = sys_path + '/' + arg
            sys_gene_length = np.array(pd.read_csv(
                sys_gene_length, header=None, delimiter='\t'), dtype=int)[0]
        elif opt in ("-t", "--promoter_strengths"):
            sys_promoter_strengths = sys_path + '/' + arg
            sys_promoter_strengths = np.array(pd.read_csv(
                sys_promoter_strengths, header=None, delimiter='\t'), dtype=float)
        elif opt in ("-s", "--samplesize"):
            sys_samplesize = arg
        elif opt in ("-n", "--training_count"):
            sys_training_count = arg
        elif opt in ("-p", "--PerturbationPower"):
            sys_PerturbationPower = arg
        elif opt in ("-o", "--output_name"):
            sys_output_name = arg
        elif opt in ("-d", "--protein_degradation"):
            sys_protein_degradation_rate = arg
        else:
            pass
    if sys_protein_degradation_rate == "":
        sys_protein_degradation_rate = [
            (0.00796)*60 for i in range(0, len(sys_gene_length))]
    else:
        pass

    ###########################
    ### PREPARING MAIN LOOP ###
    ###########################
    Loopcounter = 0
    SampleSize = int(sys_samplesize)
    TotalNumberOfGenes = len(sys_gene_length)
    minDistance = []
    meanDistance = []
    TrainingCount = int(sys_training_count)
    minDistance_Single = 1.01e6
    Target_Distance = 0
    #Recombination_Frequency = 1
    SelectionPower = 0.25
    Natural_Selection_Memory = 0
    TimeToMakeALongJump = 0
    InitialGlobalMutationRate = 2
    GlobalMutationRate = InitialGlobalMutationRate
    PerturbationPower = sys_PerturbationPower

    outfile_check = open(sys_path + '/' + sys_output_name +
                         '_CheckLoopCounter.txt'.format(TotalNumberOfGenes), 'a', buffering=1)
    '''Population proportions and TranscriptionProfile/mRNAList initiation'''
    InitialProportions = 1/SampleSize
    InitialTranscriptionProfile = []
    WTTP = copy.deepcopy(sys_WTTP)
    BaseMatrixCollector = []
    TranscriptionPofileMax = []
    TranscriptionPofileMin = []
    TranscriptionPofileAve = []
    Overexpression = np.zeros((len(sys_WTTP), TotalNumberOfGenes))
    Knockout = np.zeros((len(sys_WTTP), TotalNumberOfGenes))

    for keys in sys_WTTP:
        for each in sys_WTTP[keys][0]:
            if each >= 0:
                Knockout[int(keys)][each] = sys_WTTP[keys][1][each]
                sys_WTTP[keys][1][each] = np.nan
            elif each <= -2:
                Overexpression[int(keys)][-each-2] = sys_WTTP[keys][1][-each-2]
                sys_WTTP[keys][1][-each-2] = np.nan
            else:
                pass
        BaseMatrixCollector.append(sys_WTTP[keys][1])

    BaseMatrixCollector = np.array(BaseMatrixCollector)

    for row in BaseMatrixCollector.T:
        TranscriptionPofileMax.append(np.nanmax(row))
        TranscriptionPofileMin.append(np.nanmin(row))
        TranscriptionPofileAve.append(
            0.5*(TranscriptionPofileMax[-1]+TranscriptionPofileMin[-1]))

    #### Scaling promoter Strengths ###
    dt_scaler = (max(TranscriptionPofileMax/(np.min(sys_promoter_strengths,
                 axis=0)*60*sys_mRNA_elongation_rate/sys_gene_length))/300)
    print('dt_scaler: ', dt_scaler)
    sys_promoter_strengths = sys_promoter_strengths * dt_scaler
    ### Done Scaling promoter Strengths ###

    for x in range(0, SampleSize):
        InitialTranscriptionProfile.append(
            [float(random.randrange(1, 100)) for i in range(0, TotalNumberOfGenes)])
    np.array(InitialTranscriptionProfile)

    def ParametersInitiations(MutationRate_, which_i):
        '''Parameters Initiations'''
        # Name
        # mRNA
        Protein = [0.0 for i in range(0, TotalNumberOfGenes)]

        ####################################################
        RandomStringList = np.random.randint(
            2, size=(1, TotalNumberOfGenes**2))[0]
        RandomString = ''
        for INT in RandomStringList:
            RandomString = RandomString + str(INT)
        Configuration = String012ToMatrix(RandomString)
        ####################################################

        if type(sys_input_ChIP) == str:
            pass
        else:
            for ChIP_i in range(0, sys_input_ChIP.shape[0]):
                for ChIP_j in range(0, sys_input_ChIP.shape[1]):
                    if sys_input_ChIP[ChIP_i][ChIP_j] == 0:
                        Configuration[0][ChIP_i][ChIP_j] = 0
                        Configuration[1][ChIP_i][ChIP_j] = 0
                    else:
                        pass

        # Need to make sure the activators/repressors proportion is ~10%,
        # and the activators and repressors do not present simutaneously.

        MutationRate = 2  # must be integer

        TranscriptionRate = sys_promoter_strengths[which_i] * \
            60*sys_mRNA_elongation_rate/sys_gene_length.tolist()

        TranslationRate = (sys_aa_elongation_rate*60 /
                           (sys_gene_length/3)).tolist()

        DegradationRatemRNA = [0 for i in range(0, TotalNumberOfGenes)]
        for i in range(0, len(DegradationRatemRNA)):
            DegradationRatemRNA[i] = TranscriptionRate[i] / \
                TranscriptionPofileMax[i]

        DegradationRateProtein = np.array(
            sys_protein_degradation_rate)*dt_scaler

        DilutionRate = np.array(
            [0.0 for i in range(0, TotalNumberOfGenes)])*dt_scaler

        # Commented out: random initialization of thresholds
        # TranscriptionThreshold = np.random.randint(6,size=(1, TotalNumberOfGenes, TotalNumberOfGenes))[0]

        # Get thresholds from input data instead:
        TranscriptionThreshold = np.array([[5.0 for i in range(0, TotalNumberOfGenes)]
                                           for i in range(0, TotalNumberOfGenes)])
        for i in range(0, TranscriptionThreshold.shape[0]):
            for j in range(0, TranscriptionThreshold.shape[1]):
                TranscriptionThreshold[i][j] = ((TranscriptionPofileAve[i]*TranslationRate[i])
                                                / DegradationRateProtein[i])

        Sigmoid_k_init = []
        for i in range(0, TotalNumberOfGenes):
            Sigmoid_k_init.append(np.log(10**-3/(1-10**-3))
                                  / np.log(TranscriptionThreshold[i][0]
                                           / (TranscriptionPofileMax[i]*TranslationRate[i]
                                              / DegradationRateProtein[i])))

        LogicGates = np.array([np.random.randint(2, size=2).tolist()
                              for i in range(0, TotalNumberOfGenes)])

        Leakage = []
        for i in range(0, TotalNumberOfGenes):
            Leakage.append(TranscriptionPofileMin[i]*DegradationRatemRNA[i])

        f0 = []
        for i in range(0, TotalNumberOfGenes):
            f0.append((DegradationRatemRNA[i]*TranscriptionPofileAve[i] -
                      Leakage[i])/(TranscriptionRate[i]-Leakage[i]))

        return(Protein, Configuration, MutationRate,
               TranscriptionRate, TranslationRate,
               DegradationRatemRNA, DegradationRateProtein,
               DilutionRate, TranscriptionThreshold,
               LogicGates, Sigmoid_k_init, Leakage, f0)

    GRN_List = []
    Out_List = []
    for i in range(0, SampleSize):
        GRN_List.append('GRN_{}'.format(i+1))
        Out_List.append('Out_{}'.format(i+1))

    for i in range(0, SampleSize):
        (Protein, Configuration, MutationRate,
         TranscriptionRate, TranslationRate,
         DegradationRatemRNA, DegradationRateProtein,
         DilutionRate, TranscriptionThreshold,
         LogicGates, Sigmoid_ks, Leakages, f0) = ParametersInitiations(
            GlobalMutationRate, 0)
        exec('{} = GRN(\'{}\', {},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{})'.format(
            GRN_List[i], GRN_List[i],
            'InitialTranscriptionProfile[i]',
            'Protein',
            'Configuration',
            'InitialProportions',
            'MutationRate', 'TranscriptionRate', 'TranslationRate',
            'DegradationRatemRNA', 'DegradationRateProtein',
            'DilutionRate',
            'TranscriptionThreshold',
            'LogicGates', 'Sigmoid_ks', 'Leakages', 'f0',
            'sys_input_ChIP'))

    Proportions = []
    for i in range(0, SampleSize):
        Proportions.append(eval('{}.Proportion'.format(GRN_List[i])))

    DistanceToShow = []
    ProportionToShow = []
    mRNAMonitor = [[] for i in range(0, TotalNumberOfGenes)]
    ProteinMonitor = [[] for i in range(0, TotalNumberOfGenes)]
    DistanceOnLoop = [[] for m in range(0, len(WTTP))]
    OneTimeRecord = []

    
    ###########################
    ### MAIN EVOLUTION LOOP ###
    ###########################
    while Loopcounter <= int(sys_iteration_num):
        # Loop n times
        starting_time = time.time()
        CurrentProportions = []
        CurrentDistance = []
        DistanceOnLoop = [[] for m in range(0, len(WTTP))]
        f0_matrix = [[] for m in range(0, len(WTTP))]
        Attractor_collector = {}
        for j in range(0, SampleSize):
            Attractor_collector[str(j)] = []

        for SPP_i in range(0, len(WTTP)):
            '''Update f0'''
            NewmRNA = np.array(
                [0 for i in range(0, len(WTTP[str(SPP_i)][1]))], dtype=float)
            NewProtein = np.array(
                [0 for i in range(0, len(WTTP[str(SPP_i)][1]))], dtype=float)
            for j in range(0, len(NewmRNA)):
                # The New mRNA list is generated by WTTP*Perturbation
                NewmRNA[j] = WTTP[str(SPP_i)][1][j]
            for j in range(0, len(NewmRNA)):
                # The New Protein list is calculated by the steady state of mRNA
                NewProtein[j] = eval(
                    '({}.TranslationRate[j]*NewmRNA[j])/({}.DegradationRateProtein[j]+{}.DilutionRate[j])'.format(
                        GRN_List[0], GRN_List[0], GRN_List[0]))
            for j in range(0, SampleSize):
                exec('{}.SetmRNA({},{},{})'.format(
                    GRN_List[j], 'NewmRNA', WTTP[str(SPP_i)][0], 'Overexpression[SPP_i]'))
                exec('{}.SetProtein({},{},{})'.format(
                    GRN_List[j], 'NewProtein', WTTP[str(SPP_i)][0], 'Overexpression[SPP_i]'))
            for j in range(0, SampleSize):
                f0_matrix[SPP_i].append(
                    eval('{}.Updatef0_1()'.format(GRN_List[j])))

        f0_matrix = np.array(f0_matrix)
        f0_list = []
        for j in range(0, SampleSize):
            temp_f0_list = maskmean(f0_matrix)
            f0_list.append(True not in np.isnan(temp_f0_list))
            exec('{}.Updatef0_2(temp_f0_list)'.format(GRN_List[j]))

        index_to_train = np.argwhere(f0_list).flatten()

        for i in range(0, len(WTTP)):

            '''Update promoter strength and dependent parameters'''
            for j in range(0, SampleSize):
                TR = sys_promoter_strengths[i]*60 * \
                    sys_mRNA_elongation_rate/sys_gene_length.tolist()
                GRN_instance = eval(GRN_List[j])
                GRN_instance.Update_Transcription_Rate(TR)

            '''Setting up mRNA Protein numbers'''
            NewmRNA = np.array(
                [0 for i in range(0, len(WTTP[str(i)][1]))], dtype=float)
            NewProtein = np.array(
                [0 for i in range(0, len(WTTP[str(i)][1]))], dtype=float)
            for j in range(0, len(NewmRNA)):
                # The New mRNA list is generated by WTTP*Perturbation
                NewmRNA[j] = WTTP[str(i)][1][j]*(1
                                                 + PerturbationPower
                                                 * (np.random.random_sample([1])[0]
                                                    * np.random.choice([1, -1])))
            for j in range(0, len(NewmRNA)):
                # The New Protein list is calculated by the steady state of mRNA
                NewProtein[j] = eval(
                    '({}.TranslationRate[j]*NewmRNA[j])/({}.DegradationRateProtein[j]+{}.DilutionRate[j])'.format(
                        GRN_List[0], GRN_List[0], GRN_List[0]))

            for j in range(0, SampleSize):
                exec('{}.SetmRNA_AntiSelfAct({},{},{},{},{},{})'.format(GRN_List[j], 'NewmRNA', WTTP[str(
                    i)][0], 'TranscriptionPofileMax', 'Overexpression[i]', 'WTTP', 'i'))
                exec('{}.SetProtein_AntiSelfAct({},{})'.format(
                    GRN_List[j], WTTP[str(i)][0], 'Overexpression[i]'))
                # exec('{}.Reset_Memory()').format(GRN_List[j])
            expanded_index_to_train = []
            for instance_index in index_to_train:
                for mRNA_len in range(0, len(eval('{}.mRNA'.format(GRN_List[instance_index])))):
                    expanded_index_to_train.append((instance_index, mRNA_len))
            '''Run the model for a few times'''
            manager = multiprocessing.Manager()
            return_dict = manager.dict()
            jobs = []

            for each_eitt in range(0, len(expanded_index_to_train)):
                p = multiprocessing.Process(
                    target=Run_Dynamics,
                    args=(eval('{}'.format(GRN_List[expanded_index_to_train[each_eitt][0]])),
                          i, each_eitt, TrainingCount, WTTP,
                          expanded_index_to_train[each_eitt][1], Overexpression[i], return_dict))
                jobs.append(p)
                p.start()

            for proc in jobs:
                proc.join()

            CurrentDistance = []

            for j in range(0, SampleSize):
                if j not in index_to_train:
                    CurrentDistance.append(5*TotalNumberOfGenes)
                else:
                    CurrentDistance.append('TBD')

            CurrentDistance_2 = []
            for j in range(0, len(expanded_index_to_train)):
                IsPointAttractor = True
                npmRNA = np.round(
                    np.array(return_dict[j][:TotalNumberOfGenes])).T
                npmRNA_continuous = np.array(
                    return_dict[j][:TotalNumberOfGenes]).T
                exec('{}.SetMemo({})'.format(
                    GRN_List[expanded_index_to_train[j][0]], 'npmRNA_continuous'))
                for o in range(0, TotalNumberOfGenes):
                    if np.std(npmRNA[npmRNA.shape[0]-100:npmRNA.shape[0]+1, o]) < 0.05*TranscriptionPofileMax[o]:
                        continue
                    else:
                        IsPointAttractor = False
                if IsPointAttractor:
                    Memo_mRNA = np.mean(npmRNA_continuous[-50:], axis=0)
                    if expanded_index_to_train[j][1] == 0:
                        Attractor_collector[str(
                            expanded_index_to_train[j][0])].append(Memo_mRNA)

                        # Record fitness for this timepoint
                        CurrentDistance_2.append(
                            GetAttractorDistance(Memo_mRNA,
                                                 WTTP[str(i)][1],
                                                 TranscriptionPofileMax))

                    else:
                        FindCloestAttractor = []
                        for q in range(0, len(WTTP)):
                            FindCloestAttractor.append(GetAttractorDistance(Memo_mRNA, WTTP[str(q)][1],
                                                                            TranscriptionPofileMax))
                        CurrentDistance_2.append(min(FindCloestAttractor))

                else:
                    Memo_mRNA = np.array(
                        [0.1 for g in range(0, TotalNumberOfGenes)])
                    if expanded_index_to_train[j][1] == 0:
                        Attractor_collector[str(
                            expanded_index_to_train[j][0])].append(Memo_mRNA)
                        CurrentDistance_2.append(5*TotalNumberOfGenes)
                    else:
                        CurrentDistance_2.append(5*TotalNumberOfGenes)

            IDic = {}
            for j in range(0, len(expanded_index_to_train)):
                if expanded_index_to_train[j][0] not in IDic:
                    IDic[expanded_index_to_train[j][0]] = [
                        CurrentDistance_2[j]]
                else:
                    IDic[expanded_index_to_train[j][0]].append(
                        CurrentDistance_2[j])

            for keys in IDic:
                CurrentDistance[keys] = np.mean(IDic[keys])

            if 'TBD' in CurrentDistance:
                raise Exception('TBD in CurrentDistance!')
            else:
                pass

            DistanceOnLoop[i] = copy.deepcopy(CurrentDistance)

        DistanceOnLoop = np.array(DistanceOnLoop)
        DistanceOnSample = DistanceOnLoop.sum(axis=0)

        BestIndex = np.argmin(DistanceOnSample)
        exec('TheGuy=copy.deepcopy({})'.format(GRN_List[BestIndex]))
        minDistance.append(np.min(DistanceOnSample/len(WTTP)))
        meanDistance.append(np.mean(DistanceOnSample/len(WTTP)))

        for j in range(0, SampleSize):
            exec('{}.ResetMemo_mRNA()'.format(GRN_List[j]))

        # Updating Proportions
        for j in range(0, SampleSize):
            # Record current proportions
            CurrentProportions.append(
                eval('{}.Proportion'.format(GRN_List[j])))
        NextProportions = GetNextProportions(
            DistanceOnSample, CurrentProportions)
        DistanceToShow.append(np.mean(DistanceOnSample))
        for j in range(0, SampleSize):
            exec('{}.UpdateProportion({})'.format(
                GRN_List[j], NextProportions[j]))

        # Make a Long Jump
        if TimeToMakeALongJump >= 100 and minDistance_Single <= 0.01:
            if GlobalMutationRate < (TotalNumberOfGenes*3):
                GlobalMutationRate = GlobalMutationRate + 1
            else:
                pass
            for j in range(0, SampleSize):
                #exec('{}.SetMutationRate({})'.format(GRN_List[j], 'GlobalMutationRate'))
                pass
            CurrentProportionsCopy = copy.deepcopy(CurrentProportions)
            LongJumpIndex_List = GetLongJumpList(
                GRN_List, CurrentProportionsCopy, SelectionPower)
            for j in range(0, len(LongJumpIndex_List)):
                (Protein, Configuration, MutationRate,
                 TranscriptionRate, TranslationRate,
                 DegradationRatemRNA, DegradationRateProtein,
                 DilutionRate,
                 TranscriptionThreshold,
                 LogicGates, Sigmoid_ks, Leakages, f0) = ParametersInitiations(GlobalMutationRate, i)
                exec('{} = GRN(\'{}\', {},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{})'.format(
                    GRN_List[i], GRN_List[i],
                    'InitialTranscriptionProfile[i]',
                    'Protein',
                    'Configuration',
                    'InitialProportions',
                    'MutationRate', 'TranscriptionRate', 'TranslationRate',
                    'DegradationRatemRNA', 'DegradationRateProtein',
                    'DilutionRate',
                    'TranscriptionThreshold',
                    'LogicGates', 'Sigmoid_ks', 'Leakages', 'f0',
                    'sys_input_ChIP'))
            # print('*** LongJump&IncreaseMutationPower ***')
        else:
            pass

        '''Natural Selection'''
        if Natural_Selection_Memory:
            CurrentProportionsCopy = copy.deepcopy(CurrentProportions)
            GRN_List, Length, Copyindex = NaturalSelection(
                GRN_List, CurrentProportionsCopy, SelectionPower)
            for j in range(SampleSize-Length, SampleSize):
                exec('{}=copy.deepcopy({})'.format(
                    GRN_List[j], GRN_List[Copyindex]))
            # print('*** NaturalSelection ***')
        else:
            pass

        # Examine Progress and Generate Mutation . . .
        # If the Population has converged to its optima:
        if round(meanDistance[-1], 5) < ((1+2*SelectionPower)*round(minDistance[-1], 5)):
            # Good mutations occur and the population has converged
            if minDistance[-1] < minDistance_Single:
                for j in range(0, SampleSize):
                    exec('{}=copy.deepcopy({})'.format(
                        Out_List[j], GRN_List[j]))
                minDistance_Single = copy.deepcopy(minDistance[-1])
                TheBestGuy = copy.deepcopy(TheGuy) # type: ignore (TheGuy is set by 
                                                   # exec('TheGuy=copy.deepcopy({})'.format(
                                                       # GRN_List[BestIndex])))
                                                   # above
                                                   
                # print('*** Update_minDistanceSingle ***')
                GlobalMutationRate = InitialGlobalMutationRate  # Reset GlobalMutationRate
                '''Stop Criteria'''
                if minDistance_Single < Target_Distance:
                    break
                else:
                    pass
            else:  # No good mutations occur but the population has already converged
                pass
        else:  # The population has not converged
            pass

            # print('*** Recombination ***')
        # No good mutations occur, don't care if population converge
        if minDistance[-1] >= minDistance_Single:
            if minDistance[-1] == minDistance_Single:  # At least no bad mutations occur
                GlobalMutationDirection = InitiateGlobalMutationDirection(
                    TotalNumberOfGenes)
                GlobalMutationDirection_LG = InitiateGlobalMutationDirection_LG(
                    TotalNumberOfGenes)

                for j in range(0, SampleSize):
                    # Configurations
                    Temp_Matrix = eval('{}.Configuration'.format(GRN_List[j]))
                    Temp_012String = ConfigurationTo012(Temp_Matrix)
                    for aj in range(0, len(Temp_012String)):
                        GlobalMutationDirection[str(aj)][int(Temp_012String[aj])] = (
                            GlobalMutationDirection[str(aj)][int(Temp_012String[aj])] + 1)

                    # Logic Gates
                    Temp_Matrix_LG = eval('{}.LogicGates'.format(GRN_List[j]))
                    Temp_012String_LG = Matrix2String01_LG(Temp_Matrix_LG)
                    for aj in range(0, len(Temp_012String_LG)):
                        GlobalMutationDirection_LG[str(aj)][int(Temp_012String_LG[aj])] = (
                            GlobalMutationDirection_LG[str(aj)][int(Temp_012String_LG[aj])] + 1)

                for j in range(0, SampleSize):
                    exec('{}.GenerateMutation({}, {}, {})'.format(
                        GRN_List[j],
                        GlobalMutationDirection,
                        SampleSize,
                        GlobalMutationDirection_LG))  # Generate Mutation

                Natural_Selection_Memory = 0
                TimeToMakeALongJump = TimeToMakeALongJump + 1

            # Bad mutations happen and go back to the last checkpoint (Out_List)
            elif minDistance[-1] > minDistance_Single:
                GlobalMutationDirection = InitiateGlobalMutationDirection(
                    TotalNumberOfGenes)
                GlobalMutationDirection_LG = InitiateGlobalMutationDirection_LG(
                    TotalNumberOfGenes)
                for j in range(0, SampleSize):
                    # Configurations
                    Temp_Matrix = eval('{}.Configuration'.format(GRN_List[j]))
                    Temp_012String = ConfigurationTo012(Temp_Matrix)
                    for aj in range(0, len(Temp_012String)):
                        GlobalMutationDirection[str(aj)][int(Temp_012String[aj])] = (
                            GlobalMutationDirection[str(aj)][int(Temp_012String[aj])] + 1)

                    # Logic Gates
                    Temp_Matrix_LG = eval('{}.LogicGates'.format(GRN_List[j]))
                    Temp_012String_LG = Matrix2String01_LG(Temp_Matrix_LG)
                    for aj in range(0, len(Temp_012String_LG)):
                        GlobalMutationDirection_LG[str(aj)][int(Temp_012String_LG[aj])] = (
                            GlobalMutationDirection_LG[str(aj)][int(Temp_012String_LG[aj])] + 1)

                for j in range(0, SampleSize):
                    exec('{}=copy.deepcopy({})'.format(
                        GRN_List[j], Out_List[j]))
                    exec('{}.GenerateMutation({}, {}, {})'.format(
                        GRN_List[j],
                        GlobalMutationDirection,
                        SampleSize,
                        GlobalMutationDirection_LG))  # Generate Mutation

                Natural_Selection_Memory = 0
                TimeToMakeALongJump = TimeToMakeALongJump + 1
                # print('*** ReSetPopulation ***')
            else:
                pass
            # print('*** Mutation ***')
        else:  # Good mutations occur, do nothing and let the population converge
            Natural_Selection_Memory = 1
            TimeToMakeALongJump = 0
        outfile_check.write('{}_{}\n'.format(
            Loopcounter, time.time()-starting_time))

        OutputFile3 = open(
            sys_path + '/' + sys_output_name + '_ResultFlow.txt', 'a')
        CurrentProportions = []
        for j in range(0, SampleSize):
            CurrentProportions.append(
                eval('{}.Proportion'.format(GRN_List[j])))
        FinalIndex = np.argmax(CurrentProportions)
        TheBestGuy_Configuration = eval(
            '{}.Configuration'.format(GRN_List[FinalIndex]))
        TheBestGuy_LogicGates = eval(
            '{}.LogicGates'.format(GRN_List[FinalIndex]))
        OutputFile3.write('{}\t{}\t{}\n'.format(
            ConfigurationTo012(TheBestGuy_Configuration),
            Matrix2String01_LG(TheBestGuy_LogicGates),
            minDistance_Single))

        OutputFile3.close()

        Loopcounter = Loopcounter+1
    
    ##################################
    ### END OF MAIN EVOLUTION LOOP ###
    ##################################
    
    outfile_check.close()
    
    #######################################
    ### SAVE ADDITIONAL OUTPUT TO FILES ###
    #######################################

    OutputFile = open(
        sys_path + '/' + 'BetaTest_G{}_Result.txt'.format(TotalNumberOfGenes), 'a')
    TheBestGuy_Configuration = TheBestGuy.Configuration
    TheBestGuy_LogicGates = TheBestGuy.LogicGates
    OutputFile.write('{}\t{}\t{}\n'.format(
        ConfigurationTo012(TheBestGuy_Configuration),
        Matrix2String01_LG(TheBestGuy_LogicGates),
        minDistance_Single))
    OutputFile.close()

    OutputFile = open(
        sys_path + '/' + 'BetaTest_G{}_Result_f0.txt'.format(TotalNumberOfGenes), 'a')
    TheBestGuy_Configuration = TheBestGuy.Configuration
    TheBestGuy_LogicGates = TheBestGuy.LogicGates
    TheBestGuy_f0 = TheBestGuy.f0
    TheBestGuy_sigmoid = TheBestGuy.Sigmoid_k
    TheBestGuy_TranscriptionRate = TheBestGuy.TranscriptionRate
    TheBestGuy_TranslationRate = TheBestGuy.TranslationRate
    TheBestGuy_DRmRNA = TheBestGuy.DegradationRatemRNA
    TheBestGuy_DRPro = TheBestGuy.DegradationRateProtein
    TheBestGuy_Leakage = TheBestGuy.Leakage
    TheBestGuy_Threshold = TheBestGuy.TranscriptionThreshold

    OutputFile.write('Adjacency Matrix:\t{}\nLogic Gate:\t{}\nf0:\t{}\nHill Coefficient:\t{}\nTranscription Rate:\t{}\nTranslation Rate:\t{}\nmRNA Degradation Rate:\t{}\nProtein Degradation Rate:\t{}\nLeakage Rate:\t{}\nTF Effective Threshold:\t{}\n\n'.format(
        ConfigurationTo012(TheBestGuy_Configuration),
        Matrix2String01_LG(TheBestGuy_LogicGates),
        TheBestGuy_f0,
        TheBestGuy_sigmoid,
        TheBestGuy_TranscriptionRate,
        TheBestGuy_TranslationRate,
        TheBestGuy_DRmRNA,
        TheBestGuy_DRPro,
        TheBestGuy_Leakage,
        TheBestGuy_Threshold))

    OutputFile.close()

    outfile = open(sys_path + '/' + sys_output_name + '_CheckPoint.txt', 'a')
    outfile.write('SS_0')
    outfile.write('\n')

    for j in range(0, SampleSize):
        outfile.write(ConfigurationTo012(
            eval('{}.Configuration'.format(Out_List[j]))))
        outfile.write('\t')
        outfile.write(str(eval('{}.Proportion'.format(Out_List[j]))))
        outfile.write('\t')
        outfile.write(Matrix2String01_LG(
            eval('{}.LogicGates'.format(Out_List[j]))))
        outfile.write('\t')
        outfile.write(str(eval('{}.f0'.format(Out_List[j]))))
        outfile.write('\n')

    outfile.write('SS_1')
    outfile.write('\n')
    outfile.close()

    Calculate_AttractorDistance = []
    for Global_i in range(0, len(WTTP)):

        '''Setting up mRNA Protein numbers'''
        NewmRNA = np.array(
            [0 for i in range(0, len(WTTP[str(Global_i)][1]))], dtype=float)
        NewProtein = np.array(
            [0 for i in range(0, len(WTTP[str(Global_i)][1]))], dtype=float)
        for j in range(0, len(NewmRNA)):
            # The New mRNA list is generated by WTTP*Perturbation
            NewmRNA[j] = WTTP[str(Global_i)][1][j]

        for j in range(0, len(NewmRNA)):
            NewProtein[j] = eval(
                '({}.TranslationRate[j]*NewmRNA[j])/({}.DegradationRateProtein[j]+{}.DilutionRate[j])'.format(
                    GRN_List[0], GRN_List[0], GRN_List[0]))

        #  Set overexpression or knockout
        for j in range(0, SampleSize):
            #exec('{}.SetmRNA({},{})'.format(GRN_List[j], 'NewmRNA', knocklist[Global_i]))
            exec('{}.SetmRNA({},{},{})'.format('TheBestGuy', 'NewmRNA',
                 WTTP[str(Global_i)][0], 'Overexpression[Global_i]'))
            #exec('{}.SetProtein({},{})'.format(GRN_List[j], 'NewProtein', knocklist[Global_i]))
            exec('{}.SetProtein({},{},{})'.format('TheBestGuy', 'NewProtein',
                 WTTP[str(Global_i)][0], 'Overexpression[Global_i]'))

        # Run the model for a few times
        mRNACheckList = []
        ProteinCheckList = []

        scipystring = eval('{}.Delta_mRNA(WTTP[str(Global_i)][0], Overexpression[i])'.format(
            'TheBestGuy'))  # Calculate delta_mRNA
        exec(scipystring)

        sol = solve_ivp(update_mRNA_protein, # type: ignore (update_mRNA_protein is by exec(scipystring))
                        [0, TrainingCount],
                        list(NewmRNA)+list(NewProtein),
                        args=([WTTP[str(Global_i)][0]]),
                        t_eval=[int(TrainingCount/250)*tick for tick in range(0, 250+1)]) 

        npmRNA_continuous = sol.y[:TotalNumberOfGenes].T
        npmRNA = np.round(sol.y[:TotalNumberOfGenes]).T

        Plot_Curve = [[] for i in range(0, TotalNumberOfGenes)]
        Average_level = []

        exec('{}.SetMemo({})'.format('TheBestGuy', 'npmRNA_continuous'))

        for i in range(0, TotalNumberOfGenes):
            for j in range(0, npmRNA.shape[0]):
                Plot_Curve[i].append(npmRNA_continuous[j][i])

        for i in range(0, TotalNumberOfGenes):
            Average_level.append(np.mean(Plot_Curve[i][-50:]))

        Calculate_AttractorDistance.append(
            np.mean(TheBestGuy.Memo_mRNA[-50:], axis=0))

    CurrentDistance = []
    for G1 in range(0, len(WTTP)):
        CurrentDistance.append(GetAttractorDistance(
            Calculate_AttractorDistance[G1], WTTP[str(G1)][1], TranscriptionPofileMax))
    print(np.sum(CurrentDistance))

    Attractor_Distance_matrix = np.zeros((len(WTTP), len(WTTP)))
    for novel_i in range(0, len(Calculate_AttractorDistance)):
        for novel_j in range(0, len(Calculate_AttractorDistance)):
            Attractor_Distance_matrix[novel_i][novel_j] = GetAttractorDistance(
                Calculate_AttractorDistance[novel_i],
                Calculate_AttractorDistance[novel_j],
                TranscriptionPofileMax)
    Unique_attractor_N = (
        1+np.count_nonzero(np.mean(Attractor_Distance_matrix, axis=0) > 1))

    print(Unique_attractor_N)

    CurrentDistance = np.array(CurrentDistance)
    print(np.mean(CurrentDistance)/Unique_attractor_N)
