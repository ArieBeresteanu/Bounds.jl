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

## Plan: add DataFrame capabilities
def EY(yl:list,yu:list,H0:list,options:Options=default_options,method="Asymptotic"):
    #THis is the shell function that calls either the asymtotic distribution version or the bootstrap version of EY
    if method =="Asymptotic":
        return EYasy(yl,yu,H0,options)
    else:
        return EYboot(yl,yu,H0,options)

def EYboot(yl:list, yu:list, H0:list, options:Options=default_options):
    #This function uses a bootstrap test. This option is not in BM(2008) for EY but it is proved for BLP in section 4
    LB = np.mean(yl)
    UB = np.mean(yu)
    bound = [LB, UB]
    
    # test Statistic
    n = len(yl)
    sqrt_n = np.sqrt(n)
    testStat_H = sqrt_n*HdistInterval(bound, H0)
    testStat_dH = sqrt_n*dHdistInterval(bound, H0)
    
    B = options.MC_iterations #number of MC iterations to compute the critical value
    alpha = options.conf_level #confidence level for the critical value1
    
    r_H = []
    r_dH = []
    rng = np.random.Generator(options.rng)
    for i in range(B):
        indx = rng.integers(low=0, high=n, size=n)
        yl_b = yl[indx]
        yu_b = yu[indx]
        bound_b = [np.mean(yl_b), np.mean(yu_b)]
        r_H.append(sqrt_n*HdistInterval(bound_b, bound))
        r_dH.append(sqrt_n*dHdistInterval(bound_b, bound))
    
    r_H.sort()
    c_H = r_H[np.floor(alpha*B).astype(int)]
    CI_H = [LB - c_H/sqrt_n, UB+c_H/sqrt_n]
    Htest = testResults(testStat_H,c_H,CI_H) 
    
    r_dH.sort()
    c_dH = r_dH[np.floor(alpha*B).astype(int)]
    CI_dH = [LB - c_dH/sqrt_n, UB + c_dH/sqrt_n]
    dHtest = testResults(testStat_dH,c_dH,CI_dH)
       
    results = Results(bound, Htest, dHtest)
    return results            

def EYasy(yl:list, yu:list, H0:list, options:Options=default_options):
    #This function uses the test based on the asymptotic distributin as developed in BM(2008) pp. 778-779
    LB = np.mean(yl)
    UB = np.mean(yu)
    bound = [LB, UB]
    
    # test Statistic
    n = len(yl)
    sqrt_n = np.sqrt(n)
    testStat_H = sqrt_n*HdistInterval(bound, H0)
    testStat_dH = sqrt_n*dHdistInterval(bound, H0)
    
    #Simulating the asy. distribution using a MC method to establish a critical value (quantile):
    
    # critical value based on Hausdorff distance
    Pi = np.cov(yl, yu) #covariance matrix for yl yu
    
    B = options.MC_iterations #number of MC iterations to compute the critical value
    alpha = options.conf_level #confidence level for the critical value1
    
    ## Following Algorithm on page 780 in BM2008:
    rr = np.random.multivariate_normal([0,0], Pi, B) #drawing B pairs from a bivariate-normal distribution.
    
    ## test based on Hausdorff distance:
    r_H = np.amax(abs(rr),axis=1) #row max
    r_H.sort()
    c_H = r_H[np.floor(alpha*B).astype(int)]
    CI_H = [LB - c_H/sqrt_n, UB+c_H/sqrt_n]
    Htest = testResults(testStat_H,c_H,CI_H) 
    
    ## test based on directed Hausdorff distance:
    r_dH = np.amax(np.array([list(map(plus,rr[:,0])), list(map(minus,rr[:,1]))]), axis = 0)
    r_dH.sort()
    c_dH = r_dH[np.floor(alpha*B).astype(int)]
    CI_dH = [LB - c_dH/sqrt_n, UB + c_dH/sqrt_n]
    dHtest = testResults(testStat_dH,c_dH,CI_dH)
    
    results = Results(bound, Htest, dHtest)
    return results

def oneDproj(yl:list, yu:list, x:list):
    M1 = np.multiply(x, yl)
    M2 = np.multiply(x, yu)
    s = np.dot(x,x)
    bound = [sum(np.minimum(M1, M2))/s, sum(np.maximum(M1, M2))/s]
    return bound

def CI1d(yl:list, yu:list, H0:list, x:list, options:Options=default_options):
    ## computes the 1D projection of the identification set on a specific dinesion of the explanatory variable
    
    #step 1: demean x 
    x = np.array(x)-np.mean(x)
    
    #step 2: Compute the formula on page 787 in BM2008
    bound = oneDproj(yl,yu,x)
    
    n = len(yl)
    sqrt_n = sqrt(n)
    testStat_H = sqrt_n*HdistInterval(bound,H0)
    testStat_dH = sqrt_n*dHdistInterval(bound,H0)

    B = options.MC_iterations #number of MC iterations to compute the critical value
    alpha = options.conf_level  #confidence level for the critical value1

    r_H = []
    r_dH = []
    
    for i in range(B):
        indx = rng.integers(low=0, high=n, size=n)
        yl_b = yl[indx]
        yu_b = yu[indx]
        x_b = x[indx]
        bound_b = oneDproj(yl_b,yu_b,x_b)
        r_H.append(sqrt_n*HdistInterval(bound_b, bound))
        r_dH.append(sqrt_n*dHdistInterval(bound_b, bound))
    
    r_H.sort()
    c_H = r_H[np.floor(alpha*B).astype(int)]
    CI_H = [LB - c_H/sqrt_n, UB+c_H/sqrt_n]
    Htest = testResults(testStat_H,c_H,CI_H) 
    
    r_dH.sort()
    c_dH = r_dH[np.floor(alpha*B).astype(int)]
    CI_dH = [LB - c_dH/sqrt_n, UB + c_dH/sqrt_n]
    dHtest = testResults(testStat_dH,c_dH,CI_dH)
       
    results = Results(bound, Htest, dHtest)
    return results  
