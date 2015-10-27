import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pylab as plt
import matplotlib.mlab as mlab
import pandas as pd
import os
from sklearn import linear_model
import scipy

import plotly.plotly as py
from plotly.graph_objs import Scatter


class analytic(object):
    def __init__(self):

        #Plot settings.
        # Setting plot options
        font = {'family' : 'normal',
            'weight' : 'light',
            'size'   : 12}

        matplotlib.rc('font', **font)

        df = pickle.load(open("mumuDataframe2.p", "rb"))
        self.df = pd.concat(df)

    def filter_data(self, minimum, maximum):
        self.df = self.df[(self.df.logP > minimum) & (self.df.logP < maximum)]

    def filter_data_ddG(self, minimum, maximum):
        self.df = self.df[(self.df.ddG > minimum) & (self.df.ddG < maximum)]


    def clean_data(self):
        # Replace Nan with 0.
        # self.df  = self.df.dd.dropna(0)
        # self.df  = self.df.logP.dropna(0)
        # self.df.ddG = filter(lambda x: x != "Non", self.df.ddG)
        # self.df.logP = filter(lambda x: x != "Non", self.df.logP)
        # self.df = self.df[self.df.logP.str.contains("Non") == False]
        # self.df[self.df['ddG'].apply(lambda x: str(x).isdigit())]
        self.df = self.df[self.df['logP'].apply(lambda x: x != "Non")]


        X = self.df.logP.fillna(0)
        y = self.df.ddG.fillna(0)


        return X,y



analytic_obj = analytic()
analytic_obj.filter_data(-10,20)
# analytic().filter_data_ddG(-15,15)
X, y = analytic_obj.clean_data()
print len(X), len(y)

model = linear_model.LinearRegression()
model.fit(X[:,np.newaxis] ,y)

#
# print df.Mutation
# # The mean square error
# print("Residual sum of squares: %.2f"% np.mean((model.predict(X[:,np.newaxis]) - y) ** 2))
# print X[:,np.newaxis]
# print y[:, np.newaxis].T
#
# # Explained variance score: 1 is perfect prediction. A good fit would show the same variance.
# print('Variance score: %.2f' % model.score(X[:,np.newaxis], y.T))
#
# And plot.
fig = plt.figure()
ax = fig.add_subplot(311)
ax.hist(X, bins = 40)
ax.set_title('dlogP')
ax = fig.add_subplot(312)
ax.hist(y, bins = 40)
ax.set_title('ddG')
ax = fig.add_subplot(313)

ax.scatter(X,y)
ax.plot(X, model.predict(X[:,np.newaxis]))

# Pearson correlation.
r, p = scipy.stats.pearsonr(X,y)

# Add pearson to title.
ax.set_title('Linear Regression, r = %s, p = %s'%(r,p))
ax.set_xlabel('dlogP'), ax.set_ylabel('ddG')

plt.show()
#
# # ddG = []
# # for i in df.keys():
# #     print df[i]['ddG']
# #     print df[i]['logP']
# #
# #     plt.plot(df[i]['ddG'],df[i]['logP'], 'ro')
# #     plt.show()
