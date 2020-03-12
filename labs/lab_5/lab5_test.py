#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 18:27:34 2020

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

X = temp[0]
t = temp[1]
M = 10

bl_results = boostlearn(X = X, t = t, M = M)

