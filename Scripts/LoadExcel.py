import xlrd
import numpy as np

def exceltomatrix(path, ModelPeriod):

    data  = xlrd.open_workbook(path)
    table = data.sheets()[0]
    nrows = table.nrows
    ncols = table.ncols

    datamatrix = np.zeros((nrows-3, ncols))

    for i in range(ncols):
        cols = table.col_values(i)
        del cols[0:3]
        if i > 0 and ModelPeriod == 'ModelPeriod_2426_2480_1BR':
            cols[46]=0
            cols[47]=0
        datamatrix[:,i] = cols
    return datamatrix


def LoadModelFromAslam(ModelPeriod):

    if ModelPeriod == 'ModelPeriod_2426_2480_1BR':
        pathX = '/home/bo791269/Software/AntiprotonAnalysis/ReferenceFiles/AntiprotonToProtonRatio/Antiproton-Proton-Ratio_AslamModel_2426_2480.xlsx'
    elif ModelPeriod == 'ModelPeriod_2480_2506_6BR':
        pathX = '/home/bo791269/Software/AntiprotonAnalysis/ReferenceFiles/AntiprotonToProtonRatio/Antiproton-Proton-Ratio_AslamModel_2480_2506.xlsx'
    elif ModelPeriod == 'ModelPeriod_2481_2520_1BR':
        pathX = '/home/bo791269/Software/AntiprotonAnalysis/ReferenceFiles/AntiprotonToProtonRatio/Antiproton-Proton-Ratio_AslamModel_2481_2520.xlsx'

    X = exceltomatrix(pathX, ModelPeriod)

    X=np.array(X)
    print(X.shape)
    return X
 



