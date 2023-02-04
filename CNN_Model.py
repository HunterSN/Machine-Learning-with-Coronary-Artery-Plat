import numpy as np
import pandas as pd
from keras.optimizers import SGD
from keras.models import Sequential
from keras.utils import np_utils as utils
from keras.layers import Dropout, Dense, Flatten, Activation
from keras.layers.convolutional import Conv2D, MaxPooling2D
from sklearn.model_selection import train_test_split
from keras.layers import Convolution2D, Dense

zgxfVcouple = np.loadtxt('CHD_CAQ_425.csv',skiprows = 1)
ngxyVcouple = np.loadtxt('CHD_NCAQ_425.csv',skiprows = 1)


yVal = np.loadtxt('CHD_CAQO_850.csv')

groupedZgxfVcouple = zgxfVcouple.reshape(425,216)
groupedNgxyVcouple = ngxyVcouple.reshape(425,216)


zgxf = []
for m in range(0,54,1):
    zgxf.append(groupedZgxfVcouple[0][m*4:m*4+4])


for n in range(1,425):
    catch = []
    for i in range(0,54,1):
        catch.append(groupedZgxfVcouple[0][i*4:i*4+4])
    zgxf = np.vstack((zgxf, catch))        
zgxf = zgxf.reshape(425,54,4)




ngxy = []
for m in range(0,54,1):
    ngxy.append(groupedNgxyVcouple[0][m*4:m*4+4])


for n in range(1,425):
    catch = []
    for i in range(0,54,1):
        catch.append(groupedNgxyVcouple[0][i*4:i*4+4])
    ngxy = np.vstack((ngxy, catch))        
ngxy = ngxy.reshape(425,54,4)




xTrain = np.vstack((zgxf, ngxy))
xTrain = xTrain.reshape(850,54,4,1)


target = yVal

X_train, X_test, y_train, y_test = train_test_split(xTrain, yVal, test_size=0.30, random_state=0)




X_train = X_train.astype('float32')
X_test = X_test.astype('float32')




Y_train = utils.to_categorical(y_train, 10)
Y_test = utils.to_categorical(y_test, 10)


model = Sequential()


model.add(Conv2D(32, 3, 3, input_shape=(54, 4, 1), padding='same', activation='relu'))
# Dropout
model.add(Dropout(0.2))
# 添加另一个卷积层 padding ='valid'表示输出尺寸可以采用任何形式
model.add(Conv2D(32, 3, 3,activation ='relu',padding ='same'))
# 添加一个最大池化层
model.add(MaxPooling2D(pool_size =(2,2),padding ='same'))
# 展平
model.add(Flatten())
# Dense层 隐藏单元数为521
model.add(Dense(54, activation='relu'))
# Dropout
model.add(Dropout(0.3))
#output 
model.add(Dense(10, activation='softmax'))
# 编译模型 激活器选择SGD
model.compile(loss='categorical_crossentropy',optimizer=SGD(momentum=0.5, decay=0.0004), metrics=['accuracy'])


model.fit(X_train, Y_train, validation_data=(X_test, Y_test), epochs=20,batch_size=1,validation_split=0.1)











