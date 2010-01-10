#!/usr/bin/env,python
# a stacked bar plot with errorbars
from pylab import *

fullData=(85.94, 92.59, 88.83, 90.39, 91, 93.73)
fullData_std=(4, 3.04, 3.68, 3.33, 3.11, 2.52)

cfsData=(86.78, 89.14, 84.77, 90.26, 89.25, 85.98)
cfsData_std=(3.62, 3.6, 4.04, 3.28, 3.36, 4.03)

pcaData=(82.48, 92.11, 88.67, 89.13, 93, 89.94)
pcaData_std=(4.32, 3.2, 3.88, 3.81, 2.82, 3.53)

ind = arange(len(fullData))    # the x locations for the groups
width = 0.25       # the width of the bars: can also be len(x) sequence

p1 = bar(ind, fullData,   width, color='r', yerr=fullData_std)
p2 = bar(ind+width, cfsData,   width, color='b', yerr=cfsData_std)
p3 = bar(ind+2*width, pcaData,   width, color='g', yerr=pcaData_std)

ylabel('Percent correct classification')
title('Results by model and representation')
xticks(ind+width, ('DT', 'LR',\
        'NB', 'BN',
        'SVM', 'SVM-rbf') )
yticks(arange(50,100,10))
legend( (p1[0], p2[0], p3[0]), ('FULL', 'CFS', 'PCA') )
xlim(-width,len(ind))
ylim(50,100)

show()
