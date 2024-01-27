from keras.layers import Activation, Dropout, Flatten, Dense
from keras.layers import Dense, Conv1D, GlobalAveragePooling1D, MaxPooling1D,ZeroPadding1D,Convolution1D,GlobalMaxPooling1D
from keras.models import Sequential, Model, model_from_json
from keras import initializers

def LeNet5(features):
    model = Sequential()
    model.add(Conv1D(filters=6,kernel_size=5,padding='valid',activation='tanh',input_shape=(features,1)))
    model.add(MaxPooling1D(pool_size=2))
    model.add(Conv1D(filters=16,kernel_size=5,padding='valid',activation='tanh'))
    model.add(MaxPooling1D(pool_size=2))
    model.add(Flatten())
    model.add(Dense(120,activation='tanh'))
    model.add(Dense(84,activation='tanh'))
    model.add(Dense(2,activation='softmax'))
    return model

def MLP(features):
    model = Sequential()
    model.add(Dense(1000, activation='relu',input_dim=features))
    model.add(Dropout(0.5))
#   model.add(Dense(2096, activation='relu'))
#   model.add(Dropout(0.5))
#   model.add(Dense(500, activation='relu'))
#   model.add(Dropout(0.5))
    model.add(Dense(100, activation='relu'))
    model.add(Dropout(0.5))
#   model.add(Dense(1, activation='sigmoid'))
    model.add(Dense(2,activation='softmax'))
    return model

'''
    model.add(Dense(300, activation='relu', kernel_initializer=initializers.random_normal(stddev=0.01), input_dim=features))
    model.add(Dropout(0.5))
    model.add(Dense(1500, activation='relu',kernel_initializer=initializers.random_normal(stddev=0.01)))
    model.add(Dropout(0.5))
    model.add(Dense(200, activation='relu',kernel_initializer=initializers.random_normal(stddev=0.01)))
    model.add(Dropout(0.5))
    model.add(Dense(200, activation='relu',kernel_initializer=initializers.random_normal(stddev=0.01)))
    model.add(Dropout(0.5))
    model.add(Dense(2, activation='softmax'))
'''

def test(features):
    model = Sequential()
    model.add(Conv1D(15, 3, activation='relu', input_shape=(features, 1)))
    model.add(Conv1D(15, 3, activation='relu',))
    model.add(MaxPooling1D(pool_size=3, padding='valid'))
    model.add(ZeroPadding1D((1,1)))
    model.add(Conv1D(9, 3, activation='relu',))
    #model.add(Conv1D(9, 3, activation='relu',))
    model.add(GlobalAveragePooling1D())
    model.add(Dropout(0.5))
    model.add(Dense(2, activation='softmax'))
    return model


def VGG16(features):
    model = Sequential()
    model.add(Convolution1D(64, 3, activation='relu',padding='same',input_shape=(features,1)))
    model.add(Convolution1D(64, 3, activation='relu',padding='same'))
    model.add(MaxPooling1D(2, strides=1,padding='same'))
    model.add(Convolution1D(128, 3, activation='relu',padding='same'))
    model.add(Convolution1D(128, 3, activation='relu',padding='same'))
    model.add(MaxPooling1D(2, strides=1,padding='same'))
    model.add(Convolution1D(256, 3, activation='relu',padding='same'))
    model.add(Convolution1D(256, 3, activation='relu',padding='same'))
    model.add(Convolution1D(256, 3, activation='relu',padding='same'))
    model.add(MaxPooling1D(2, strides=1,padding='same'))
    model.add(Convolution1D(512, 3, activation='relu',padding='same'))
    model.add(Convolution1D(512, 3, activation='relu',padding='same'))
    model.add(Convolution1D(512, 3, activation='relu',padding='same'))
    model.add(MaxPooling1D(2, strides=1,padding='same'))
    model.add(Convolution1D(512, 3, activation='relu',padding='same'))
    model.add(Convolution1D(512, 3, activation='relu',padding='same'))
    model.add(Convolution1D(512, 3, activation='relu',padding='same'))
    model.add(MaxPooling1D(2, strides=1,padding='same'))
    #model.add(GlobalAveragePooling1D())
    model.add(Flatten())
    model.add(Dense(4096, activation='relu'))
    model.add(Dropout(0.5))
    model.add(Dense(4096, activation='relu'))
    model.add(Dropout(0.5))
    model.add(Dense(1000, activation='relu'))
    model.add(Dropout(0.5))
    model.add(Dense(2, activation='softmax'))
    return model

