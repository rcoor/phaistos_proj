import pickle
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from sklearn import linear_model
import scipy

import plotly.plotly as py
from plotly.graph_objs import Scatter

df = pickle.load(open("mumuDataframe2.p", "rb"))
df = pd.concat(df)

# Let's initiate the model
model = linear_model.LinearRegression()

# ddG = []
# for i in df.keys():
#     print df[i]['ddG']
#     print df[i]['logP']
#
#     plt.plot(df[i]['ddG'],df[i]['logP'], 'ro')
#     plt.show()

maximum = 20
minimum = -10

# Filter the data.
df = df[(df.logP > minimum) & (df.logP < maximum)]

# Replace Nan with 0.
df.logP.fillna(0)
df.ddG.fillna(0)

y = df.ddG.convert_objects(convert_numeric=True).dropna()
X = df.logP.convert_objects(convert_numeric=True).dropna()

# Fit the model.
model.fit(X[:,np.newaxis],y)

# And plot.
fig = plt.figure()
ax = fig.add_subplot(311)
ax.hist(X)
ax.set_title('dlogP')
ax = fig.add_subplot(312)
ax.hist(y)
ax.set_title('ddG')
ax = fig.add_subplot(313)

ax.scatter(X,y)
ax.plot(X, model.predict(X[:,np.newaxis]))
ax.set_xlabel('dlogP'), ax.set_ylabel('ddG')

plt.show()

print(scipy.stats.pearsonr(X,y))
