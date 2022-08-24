#########################
### Imported packages ###
#########################
import numpy as np
import pandas as pd

from vertex import *
# Class: Vertex
# Functions: subVertex, addVertex, lambdaVertex, negVertex, dot, norm, xangle, fetchY, fetchX, distVertex

from segment import * 
# Class: Segment
##-methods: checkInput, length, dim
# Functions: dotDist

from polygons import * 
# Class: Polygon
##-method: sort, angles, plot, scatter
# Functions: minkowskiSum, lambdaPolygon, dHausdorff, hausdorff

###############################
###   Defined Structures:   ###
###############################  

class Options:
    def __init__(self, MC_iterations, seed, rng, conf_level):
        self.MC_iterations = MC_iterations
        self.seed = seed
        self.rng = rng
        self.conf_level = conf_level
        
class testResults:
    def __init__(self, testStat, criticalVal, ConfidenceInterval):
        self.testStat = testStat 
        self.criticalVal = criticalVal
        self.ConfidenceInterval = ConfidenceInterval
        
class Results:
    def __init__(self, bound, Htest, dHtest):
        self.bound = bound
        self.Htest = Htest
        self.dHtest = dHtest  

#####################
###   Constants   ###
#####################

default_options = Options(2000, 15217, np.random.MT19937(15217), 0.95)

#####################
###  Functions:   ###
#####################

def plus(x):
    return max(0.0, x)

def minus(x):
    return max(0.0, -x)

def HdistInterval(v1:list, v2:list):
    v = np.array(v1) - np.array(v2)
    return max(abs(v))

def dHdistInterval(v1:list, v2:list):
    v = np.array(v1) - np.array(v2)
    return max(plus(v[0]), minus(v[1]))

def EY(yl:list, yu:list, H0:list, options:Options=default_options):
    LB = np.mean(yl)
    UB = np.mean(yu)
    bound = [LB, UB]
    
    # test Statistic
    n = len(yl)
    sqrt_n = np.sqrt(n)
    testStat_H = sqrt_n*HdistInterval(bound, H0)
    testStat_dH = sqrt_n * dHdistInterval(bound, H0)
    
    # critical value based on Hausdorff distance
    Pi = np.cov(yl, yu) #covariance matrix for yl yu
    
    B = options.MC_iterations #number of MC iterations to compute the critical value
    alpha = options.conf_level #confidence level for the critical value1
    
    ## Following Algorithm on page 780 in BM2008:
    rr = np.random.multivariate_normal([0,0], Pi, B) #drawing B pairs from a bivariate-normal distribution.
    
    ## test based on Hausdorff distance:
    r_H = np.amax(abs(rr),axis=1) #row max
    r_H.sort()
    c_H = r_H[np.floor(alpha*len(rr)).astype(int)]
    CI_H = [LB - c_H/sqrt_n, UB+c_H/sqrt_n]
    Htest = testResults(testStat_H,c_H,CI_H) 
    
    ## test based on directed Hausdorff distance:
    r_dH = np.amax(np.array([list(map(plus,rr[:,0])), list(map(minus,rr[:,1]))]), axis = 0)
    r_dH.sort()
    c_dH = r_dH[np.floor(alpha*len(rr)).astype(int)]
    CI_dH = [LB - c_dH/sqrt_n, UB + c_dH/sqrt_n]
    dHtest = testResults(testStat_dH,c_dH,CI_dH)
    
    results = Result(bound, Htest, dHtest)
    return results
