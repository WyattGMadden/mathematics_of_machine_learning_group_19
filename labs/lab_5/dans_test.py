#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 12:40:34 2020

@author: dancrowley
"""

import scipy.io as scipy
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import numpy.linalg as lg
from make_cloud import *
#from boosteval import *
from weakeval import *
from weaklearn import *
from boostlearn import *

temp = make_cloud()
dat = pd.DataFrame(np.transpose(temp[0]), columns = ("X0", "X1"))
dat['t'] = temp[1]

sns.scatterplot(x="X0", y="X1", hue="t",data=dat, legend = False)

weak_params = weaklearn(temp[0], temp[1], v = None)
dat['t_weak_pred'] = weakeval(temp[0], weak_params)
dat['correct_class'] = dat['t'] == dat['t_weak_pred']

sns.scatterplot(x="X0", y="X1", hue="t", style = "correct_class", data=dat, legend = False)

temp_boostlearn = boostlearn(X = temp[0], t = temp[1], M = 5)

dat['t_boost_5'] = boosteval(X = temp[0], 
                             params = temp_boostlearn[0], 
                             alpha = temp_boostlearn[1])
dat['correct_class_boost_5'] = dat['t'] == dat['t_boost_5']
sns.scatterplot(x="X0", y="X1", hue="t", style = "correct_class_boost_5", data=dat, legend = False)




