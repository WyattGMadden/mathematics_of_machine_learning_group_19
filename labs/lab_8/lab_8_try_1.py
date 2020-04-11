#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 15:46:53 2020

@author: dancrowley
"""
imort seaborn as sns

import numpy as np
#training data
XTrain = 2*np.pi*np.random.random(1000)  # x, sampled uniformly from [0,2pi]
YTrain = np.sin(XTrain)+0.1*np.random.normal(0, 1, 1000)    #; % corresponding y w/ noise


# validation data
XV = 2*np.pi*np.random.random(100)
YV = np.sin(XV)+0.1*np.random.normal(0, 1, 100) 

from keras.models import Sequential
from keras.layers import Dense, Activation


model = Sequential()
model.add(Dense(512, activation='relu', input_shape=(1,)))
model.add(Dense(512, activation='relu'))
model.add(Dense(512, activation='relu'))
model.add(Dense(1, activation='selu'))

#model.add(Activation('relu'))

model.summary()

# For a mean squared error regression problem
model.compile(optimizer='rmsprop',
              loss='mse')

model.fit(XTrain, YTrain,
          epochs=20,
          batch_size=128)

score = model.evaluate(XV, YV, batch_size=128)
print(score)

predict = model.predict(XV)


sns.scatterplot(np.reshape(predict, 100), np.reshape(YV, 100))