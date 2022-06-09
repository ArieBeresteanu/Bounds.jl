import numpy as np

### Silverman Rule of Thumb for bandwidth selection ##
def silverman(x:list):
    # Silverman rule of thumb
    return 1.06*(len(x)**(-0.2))*np.std(x)

### Epanechnikov Kernel function ###
#     Version 1: x, x0, and h are scalars
def epa(x:float,x0:float=0.0,h:float=1.0):
    temp = min(abs((x-x0)/h),1.0)
    return 0.75*(1.0-temp**2)

# #     Version 2: x is a vector, x0 and h are scalars
# def epa(x:list,x0:float=0.0,h:float=1.0):
#     temp = [min((x1-x0)/h, 1.0) for x1 in x]
#     return [0.75*(1.0-a**2) for a in temp]

### Kernel density estimation ###
def kdens(X:list,x0:float=0.0,h:float=1.0):
    f = sum([epa(x, x0, h) for x in X])
    return f/(len(X)*h)

### Kernel regression estimation ###
def kreg(Y:list, X:list, x0:float=0.0, h:float=1.0):
    f=0.0
    m=0.0
    nx=len(X)
    ny=len(Y)
    if nx != ny:
        raise ValueError("length of Y and X do not match")
    else:
        df = [epa(x, x0, h) for x in X]
        f = sum(df)
        m = sum([y*d for y,d in zip(Y, df)])
    if f>0:
        return m/f
    else:
        return 0.0

def kreg_w(Y:list,X:list,x0:float,h:float,w:list):
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
        df = [weight*epa(x, x0, h) for x,weight in zip(X,w)]
        f = sum(df)
        m = sum([weight*y*d for weight,y,d in zip(w,Y,df)])
    if f > 0.000001:
        return m/f
    else:
        return 0.0

# mutable struct assumptions
#     tol   :: Float64
#     Yₗ     :: Float64
#     Yᵤ    :: Float64
#     model :: String
# end

### Worst Case Scenario Bounds 
# mutable struct Results{T} 
#     prob0 :: T
#     prob1 :: T
#     yhat0 :: T
#     yhat1 :: T 
#     bound0L :: T 
#     bound0U :: T 
#     bound1L :: T 
#     bount1U :: T 
#     treatL :: T 
#     treatU :: T 
#     model  :: String
# end

def missingObs(y:list,z:list,x:list,x0:float,cont:bool,h:float=1.0, K0:float=0.0, K1:float=1.0):
    # if z=0,  y is missing
    ny = len(y)
    nx = len(x)
    nz = len(z)
    if ny !=nx or ny != nz:
        raise ValueError("vector length do not match")
    res= {}
    if cont:
        res['prob1'] = kreg(z,x,x0,h)
        res['prob0']= 1-res['prob1']
        res['yhat1'] = kreg_w(y,x,x0,h,z)
        res['yhat0'] = 'NA'
    else:
        nx0 = x.count(x0)
        temp1 = [Z for (Z,X) in zip(z,x) if X==x0]
        nz1 = sum(temp1)

        res['prob1'] = nz1/nx0
        res['prob0'] = 1-res['prob1']   
        res['yhat1'] = sum([Z for (Z,X,Y) in zip(z,x,y) if X==x0 and Y != 0.0])/nz1
        res['yhat0'] = 'NA' 
    res['bound0L'] = 'NA'
    res['bound0U'] = 'NA'
    res['bound1L'] = K0 * res['prob0'] + res['yhat1'] * res['prob1']
    res['bound1U'] = K1 * res['prob0'] + res['yhat1'] * res['prob1']
    res['treatL'] = 'NA'
    res['treatU'] = 'NA'
    res['model']  = "Missing Y observations"
    return res

def treatmentEffect(y:list,z:list,x:list,x0:float,cont:bool,h:float=1.0):
    ny = len(y)
    nx = len(x)
    nz = len(z)
    if ny !=nx or ny != nz:
        raise ValueError("vector length do not match")
    res={}
    if cont:
        res['prob1'] = kreg(z,x,x0,h)
        res['prob0']= 1-res['prob1']
        res['yhat1'] = kreg_w(y,x,x0,h,z)
        res['yhat0'] = 'NA'
        res.yhat0 = kreg(y,x,x0,h,[1-Z for Z in z]) 
    else:
        nx0 = x.count(x0)
        temp1 = [Z for (Z,X) in zip(z,x) if X==x0]
        nz1 = sum(temp1)

        res['prob1'] = nz1/nx0
        res['prob0'] = 1-res['prob1']   
        res['yhat1'] = sum([Z for (Z,X,Y) in zip(z,x,y) if X==x0 and Y != 0.0])/nz1
        res['yhat0'] = sum([1-Z for (Z,X,Y) in zip(z,x,y) if X==x0 and Y != 0.0])/nz1
    res['model']  = "Treatment Effect, no assumptions"
    return res

