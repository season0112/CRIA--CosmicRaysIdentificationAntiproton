

def CovNumber(parametername, parameterindex):

    covnumbers = 2
    
    if parametername == 'ParameterError_TOFAntiproton':
        if parameterindex == 1:
            covnumbers = 0 
        if parameterindex == 2:
            covnumbers = 0

    elif parametername == 'ParameterError_TOFElectron':
        if parameterindex == 0:
            covnumbers = 5
        if parameterindex == 1:
            covnumbers = 4
        if parameterindex == 2:
            covnumbers = 4
        if parameterindex == 3:
            covnumbers = 4

    elif parametername == 'ParameterError_TOFPion':
        if parameterindex == 0:
            covnumbers = 3
        if parameterindex == 1:
            covnumbers = 0
        if parameterindex == 2:
            covnumbers = 0
        if parameterindex == 3:
            covnumbers = 0

    
    elif parametername == 'ParameterError_TRDAntiproton':
        if parameterindex == 0:
            covnumbers = 6
        if parameterindex == 1:
            covnumbers = 4
        if parameterindex == 2:
            covnumbers = 4
     
    elif parametername == 'ParameterError_TRDElectron':
        if parameterindex == 0:
            covnumbers = 6
        if parameterindex == 1:
            covnumbers = 6
        if parameterindex == 2:
            covnumbers = 6
        if parameterindex == 3:
            covnumbers = 6

    elif parametername == 'ParameterError_TRDPion':
        if parameterindex == 0:
            covnumbers = 16
        if parameterindex == 1:
            covnumbers = 16
        if parameterindex == 2:
            covnumbers = 16
        if parameterindex == 3:
            covnumbers = 16
        if parameterindex == 4:
            covnumbers = 16
        if parameterindex == 5:
            covnumbers = 16
        if parameterindex == 6:
            covnumbers = 16

    return covnumbers
