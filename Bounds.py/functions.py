#########################
### Imported packages ###
#########################
import numpy as np

#####################
## kernel fuctions ##
#####################

import smoothingKernels

##################################
## exported functions and types ##
##################################
#export SimpleBound,missingObs, Assumptions, Results, default_options

######################
## type definitions ##
######################
class Assumptions:
    def __init__(self, tol:float=0.000001, Yl:float=0.0, Yu:float=1.0, Bootstrap_iterations:int=100, kernel = smoothingKernels.epanechnikov):
        self.tol = tol
        self.Yl = Yl
        self.Yu = Yu
        self.Bootstrap_iterations=Bootstrap_iterations
        self.kernel = kernel
    
### Worst Case Scenario Bounds 

class Results:
    def __init__(self, prob0, prob1, yhat0, yhat1, bound0L, bound0U, bound1L, bound1U, treatL, treatU, model):
        self.prob0 = prob0
        self.prob1 = prob1
        self.yhat0 = yhat0
        self.yhat1 = yhat1
        self.bound0L = bound0l
        self.bound0U = bound0U
        self.bound1L = bound1L
        self.bound1U = bound1U
        self.treatL = treatL
        self.treatU = treatU
        self.model = model

# default options
default_options = Assumptions()

#############################################
###      M A I N   F U N C T I O N S      ###
#############################################

#############################
# 1. Misceleneous functions #
#############################

### Silverman Rule of Thumb for bandwidth selection ##
def silverman(x:list):
    # Silverman rule of thumb
    return 1.06*(len(x)**(-0.2))*np.std(x)

#########################
# 2. Kernel Regression ##
#########################

### Kernel density estimation ###
def kdens(X:list,x0:float=0.0,h:float=1.0, assumptions:Assumptions = default_options):
    f = sum([assumptions.kernel((x-x0)/h) for x in X])
    return f/(len(X)*h)

### Kernel regression estimation ###
def kreg(Y:list, X:list, x0:float=0.0, h:float=1.0, assumptions:Assumptions = default_options):
    f=0.0
    m=0.0
    nx=len(X)
    ny=len(Y)
    if nx != ny:
        raise ValueError("length of Y and X do not match")
    else:
        df = [assumptions.kernel((x-x0)/h) for x in X]
        f = sum(df)
        m = sum([y*d for y,d in zip(Y, df)])
    if f>assumptions.tol:
        return m/f
    else:
        return 0.0

def kreg_w(Y:list,X:list,x0:float,h:float,w:list, assumptions:Assumptions = default_options):
    # This is a version of kreg where we can introduce weights. 
    # It is done mostly so we can conditon on binary weights, i.e. exclude some observations

    f = 0.0
    m = 0.0
    nx=len(X)
    ny=len(Y)
    nw=len(w)
    
    if nx != ny or nx !=nw:
        raise ValueError("length of Y,X, and w do not match")
    else:
        df = [weight*assumptions.kernel((x-x0)/h) for x,weight in zip(X,w)]
        f = sum(df)
        m = sum([weight*y*d for weight,y,d in zip(w,Y,df)])
    if f > assumptions.tol:
        return m/f
    else:
        return 0.0

############################
##  Missing Observations  ##
############################

def missingObs(y:list,
                z:list,
                x:list,
                x0:float,
                cont:bool,
                h:float=1.0, 
                K0:float=0.0, 
                K1:float=1.0, 
                assumptions:Assumptions = default_options):
    # if z=0,  y is missing
    ny = len(y)
    nx = len(x)
    nz = len(z)
    if ny !=nx or ny != nz:
        raise ValueError("vector length do not match")
    res = Results(0,0,'NA',0,'NA', 'NA',0,0,'NA', 'NA',"Missing Y observations")
    if cont:
        res.prob1 = kreg(z,x,x0,h,assumptions)
        res.yhat1 = kreg_w(y,x,x0,h,z,assumptions)
    else:
        nx0 = x.count(x0)
        temp1 = [Z for (Z,X) in zip(z,x) if X==x0]
        nz1 = sum(temp1)
        
        res.prob1 = nz1/nx0 
        res.yhat1 = sum([Z for (Z,X,Y) in zip(z,x,y) if X==x0 and Y != 0.0])/nz1
    res.prob0 = 1- res.prob1
    res.bound1L = K0 * res.prob0 + res.yhat1 * res.prob1
    res.bound1U = K1 * res.prob0 + res.yhat1 * res.prob1
    return res

########################
##  Treatment Effect  ##
########################

def treatmentEffect(y:list,
                    z:list,
                    x:list,
                    x0:float,
                    cont:bool,
                    h:float=1.0,
                    assumptions:Assumptions = default_options):
    ny = len(y)
    nx = len(x)
    nz = len(z)
    if ny !=nx or ny != nz:
        raise ValueError("vector length do not match")
    res = Results(0,0,0,0,'NA','NA','NA','NA','NA','NA',"Treatment Effect, no assumptions")
    if cont:
        res.prob1 = kreg(z,x,x0,h,assumptions)
        res.yhat1 = kreg_w(y,x,x0,h,z,assumptions)
        res.yhat0 = kreg_w(y,x,x0,h,[1-Z for Z in z],assumptions) 
    else:
        nx0 = x.count(x0)
        temp1 = [Z for (Z,X) in zip(z,x) if X==x0]
        nz1 = sum(temp1)

        res.prob1 = nz1/nx0
        res.yhat1 = sum([Z for (Z,X,Y) in zip(z,x,y) if X==x0 and Y != 0.0])/nz1
        res.yhat0 = sum([1-Z for (Z,X,Y) in zip(z,x,y) if X==x0 and Y != 0.0])/nz1
    res.prob0 = 1-res.prob1
    return res

class SimpleBound:
    def __init__(self, LB, UB, method):
        self.LB = LB
        self.UB = UB
        self.method = method
