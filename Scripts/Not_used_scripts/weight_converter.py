from __future__ import division
import numpy as np
import math
import json
import collections
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from keras.layers import Activation, Dropout, Flatten, Dense
from keras.layers import Conv1D, GlobalAveragePooling1D, MaxPooling1D,ZeroPadding1D,Convolution1D
from keras.models import Sequential, Model, model_from_json
from keras.optimizers import SGD, RMSprop, Adam
from keras.preprocessing.image import ImageDataGenerator
from keras.utils import np_utils
from keras import initializers
from keras.applications.resnet50 import ResNet50
from PIL import ImageFont
from PIL import Image
from PIL import ImageDraw
plt.switch_backend('agg')

features = 16

import CNN_models
model = CNN_models.VGG16(features)
model.load_weights('VGG16.h5')

model.save('VGG16_archi.h5')

