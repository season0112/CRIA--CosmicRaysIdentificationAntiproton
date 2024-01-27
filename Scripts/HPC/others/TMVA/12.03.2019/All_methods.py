import ROOT
from ROOT import TMVA, TFile, TTree, TCut, TString
import matplotlib.pyplot as plt
import numpy as np


'''
# Select Theano as backend for Keras
from os import environ
#environ['KERAS_BACKEND'] = 'theano'
environ['KERAS_BACKEND'] = 'tensorflow'
# Set architecture of system (AVX instruction set is not supported on SWAN)
environ['THEANO_FLAGS'] = 'gcc.cxxflags=-march=corei7'
from keras.models import Sequential
from keras.layers.core import Dense, Dropout
from keras.optimizers import Adam

###############################################################################################
learningrate = 0.001
epochsnumber = 1000
batchsize = 200
###############################################################################################
'''

B1042_1or9_p2 = TFile.Open("data/B1042_1or9_p2_train.root")
B1042_1and9_n2 = TFile.Open("data/B1042_1and9_n2_80.root")

# Get signal and background trees from file
background1 = B1042_1and9_n2.Get('tree')
signal1 = B1042_1or9_p2.Get('tree')


# Add variables to dataloader
dataloader = ROOT.TMVA.DataLoader('dataset_pymva')
numVariables = len(signal1.GetListOfBranches())
for branch in signal1.GetListOfBranches():
    dataloader.AddVariable(branch.GetName())
    
# Add trees to dataloader
dataloader.AddSignalTree(signal1, 1.0)
dataloader.AddBackgroundTree(background1, 1.0)
trainTestSplit = 0.7
dataloader.PrepareTrainingAndTestTree(ROOT.TCut(''),
        'TrainTestSplit_Signal={}:'.format(trainTestSplit)+\
        'TrainTestSplit_Background={}:'.format(trainTestSplit)+\
        'SplitMode=Random')

# Setup TMVA
ROOT.TMVA.Tools.Instance()
#ROOT.TMVA.PyMethodBase.PyInitialize()

outputFile = ROOT.TFile.Open('TMVAOutputPyMVA.root', 'RECREATE')
factory = ROOT.TMVA.Factory('TMVAClassification', outputFile,
        '!V:!Silent:Color:DrawProgressBar:Transformations=I,G:'+\
        'AnalysisType=Classification')
'''
# Define model
model = Sequential()
model.add(Dense(32,   activation='relu',
        input_dim=numVariables))
model.add(Dropout(0.5))
model.add(Dense(200,  activation='relu'))
model.add(Dropout(0.5))
model.add(Dense(500,  activation='relu'))
model.add(Dropout(0.5))
model.add(Dense(200,  activation='relu'))
model.add(Dropout(0.5))
model.add(Dense(100,  activation='relu'))
 model.add(Dropout(0.5))
model.add(Dense(2,    activation='softmax'))

# Set loss and optimizer
opt3 = Adam(lr=learningrate, beta_1=0.9, beta_2=0.999, epsilon=1e-08)
model.compile(loss='categorical_crossentropy', optimizer=opt3,
        metrics=['accuracy',])

# Store model to file
model.save('model.h5')

# Print summary of model
model.summary()
'''

# Keras interface with previously defined model
'''
factory.BookMethod(dataloader, ROOT.TMVA.Types.kPyKeras, 'PyKeras',
        'H:!V:VarTransform=None:FilenameModel=model.h5:'+\
        'NumEpochs=epochsnumber:BatchSize=batchsize:'+\
        'TriesEarlyStopping=5')
'''

# Gradient tree boosting from scikit-learn package
#factory.BookMethod(dataloader, ROOT.TMVA.Types.kPyGTB, 'GTB',
#        'H:!V:VarTransform=None:'+\
#        'NEstimators=100:LearningRate=0.1:MaxDepth=3')

# other
#factory.BookMethod(dataloader, TMVA.Types.kBDT, "BDT","!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" )
#factory.BookMethod(dataloader, TMVA.Types.kBDT, "BDT","!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad:GradBaggingFraction=0.5:nCuts=20:NNodesMax=5" )
#factory.BookMethod(dataloader, TMVA.Types.kBDT, "BDTG","!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" )
#factory.BookMethod(dataloader, TMVA.Types.kBDT, "BDTB","!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20" )
#factory.BookMethod(dataloader, TMVA.Types.kBDT, "BDTD","!H:!V:NTrees=400:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate" )
## Used
factory.BookMethod(dataloader, TMVA.Types.kBDT, "BDTG","!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:GradBaggingFraction=0.5:nCuts=50:MaxDepth=3:IgnoreNegWeightsInTraining")
#factory.BookMethod(dataloader, TMVA.Types.kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm" )
factory.BookMethod(dataloader, TMVA.Types.kBDT, "BDT","!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" )
#factory.BookMethod(dataloader, TMVA.Types.kMLP, "MLPBNN", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator" )
#factory.BookMethod( dataloader, TMVA.Types.kLikelihood, "Likelihood","H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50" )
#factory.BookMethod( dataloader, TMVA.Types.kLikelihood, "LikelihoodKDE","!H:!V:!TransformOutput:PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:KDEFineFactor=0.3:KDEborder=None:NAvEvtPerBin=50" )
#factory.BookMethod( dataloader, TMVA.Types.kFisher, "Fisher", "H:!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" )

## Fabian
#factory.BookMethod(dataloader,TMVA.Types.kBDT, "BDT", "!H:!V:NTrees=300:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:IgnoreNegWeightsInTraining")
#factory.BookMethod(dataloader, TMVA.Types.kBDT, "BDTG", "!H:!V:NTrees=200:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:GradBaggingFraction=0.5:nCuts=20:MaxDepth=2:IgnoreNegWeightsInTraining")



factory.TrainAllMethods()

factory.TestAllMethods()

factory.EvaluateAllMethods()

# Enable Javascript for ROOT so that we can draw the canvas
#%jsroot on

# Print ROC
canvas = factory.GetROCCurve(dataloader)
canvas.Draw()
canvas.Print("All_methods.png")






'''

# fill histograms for signal and background from the test sample tree
ROOT.TestTree.Draw("BDT>>hSig(22,-1.1,1.1)","classID == 0","goff")  # signal
ROOT.TestTree.Draw("BDT>>hBg(22,-1.1,1.1)","classID == 1", "goff")  # background
 
ROOT.hSig.SetLineColor(ROOT.kRed); ROOT.hSig.SetLineWidth(2)  # signal histogram
ROOT.hBg.SetLineColor(ROOT.kBlue); ROOT.hBg.SetLineWidth(2)   # background histogram
 
# use a THStack to show both histograms
hs = ROOT.THStack("hs","")
hs.Add(ROOT.hSig)
hs.Add(ROOT.hBg)
 
# show the histograms
gcSaver.append(ROOT.TCanvas())
hs.Draw()

'''










'''
###########################################################################################
################## Plot Loss and Accuracy Prediction ###################
plt.figure(figsize=(18,18))
plt.subplot(221)
loss = model.fit.history['loss']
val_loss = model.fit.history['val_loss']
epochs = range(1, len(loss) + 1)
plt.plot(epochs, loss, 'bo', label='Training loss')
plt.plot(epochs, val_loss, 'b', label='Validation loss')
plt.title('Training and validation loss',fontsize=26)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlabel('Epochs',fontsize=22)
plt.ylabel('Loss',fontsize=22)
plt.legend(loc='best',fontsize=20)
plt.subplot(222)
acc = model.fit.history['acc']
val_acc = model.fit.history['val_acc']
plt.plot(epochs, acc, 'bo', label='Training acc')
plt.plot(epochs, val_acc, 'b', label='Validation acc')
plt.title('Training and validation accuracy',fontsize=26)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlabel('Epochs',fontsize=22)
plt.ylabel('Accuracy',fontsize=22)
plt.legend(loc='best',fontsize=20)
plt.subplot(223)
plt.hist(y_pred[0:testset1.shape[0],1],bins=binnumber,range=(0,1),log=True )
plt.title("Proton",fontsize=26)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlabel('Estimator $_{CC}$',fontsize=22)
plt.ylabel('Count',fontsize=22)
plt.subplot(224)
plt.hist(y_pred[testset1.shape[0]:y_pred.shape[0],1],bins=binnumber,range=(0,1),log=True )
plt.title("CCProton",fontsize=26)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlabel('Estimator $_{CC}$',fontsize=22)
plt.ylabel('Count',fontsize=22)
plt.savefig('ML.png')
'''










