import numpy as np
from numpy import mean, abs 
from numpy.linalg import norm

def corr_comps(x, xest):
    idx = np.where(abs(x)>0)
    x_aux = x[idx]
    xest_aux = xest[idx]
    cor = (
        abs(sum((x_aux-mean(x_aux)) * (xest_aux-mean(xest_aux)))) 
        /(norm(x_aux-mean(x_aux)) * norm(xest_aux-mean(xest_aux))+ 1e-15)
            )
    return cor

def mse(x, xest):
    assert len(x) == len(xest), 'Should be of equal length.'
    idx = np.where(abs(x)>0)
    x_aux = x[idx]
    xest_aux = xest[idx]
    error = np.mean((x_aux-xest_aux)**2)
    return error

def compute_qrf(x, x_hat, tmin=None,tmax=None):
    """
    Quality reconstruction factor
    """
    if tmin is None:
        tmin = 0
    
    if tmax is None:
        tmax = len(x)

    x = x[tmin:tmax]
    x_hat = x_hat[tmin:tmax]
    qrf = 10*np.log10(np.sum(x**2)/np.sum((x_hat-x)**2))
    return qrf


# def order_components(Xest, X, metric = corr_comps):
#     order = []
#     values = np.array([[metric(x,xest) for x in X] for xest in Xest])
#     # order = np.argmax(values, axis=0)
#     for i in range(values.shape[1]):
#         col = values[:,i]
#         sort_col = np.sort(col)[-1::-1]
#         for j in range(len(sort_col)):
#             row = values[np.where(col == sort_col[j])]
#             if sort_col[j] == np.max(row):
#                 order.append(np.where(col == sort_col[j])[0][0])
#                 break
#     return order

def order_components(Xest, X, minormax = 'max', metric = corr_comps):
    order = [[] for aaa in range(len(X))]
    values = np.array([[metric(x,xest) for x in X] for xest in Xest], dtype=object)
    if minormax=='max':
        fun = np.argmax
        factor = -1
    if minormax == 'min':
        fun = np.argmin
        factor = 1

    while np.any([k == [] for k in order]):
        ind = np.unravel_index(fun(values, axis=None), values.shape)
        if (ind[0] not in order) and (order[ind[1]] == []):
            order[ind[1]] = ind[0]
        values[ind] = factor*np.inf
    return order    


def compare_qrf_block(signal, method_output, tmin=None, tmax=None):
    X = signal.comps
    output = []
    for Xest in method_output:
        order = order_components(Xest, X)
        Xaux = Xest[order]
        qrfs = []
        for x,xaux in zip(X,Xaux):
            indx = np.where(np.abs(x)>0)
            qrfs.append(compute_qrf(x[indx], xaux[indx],tmin=tmin,tmax=tmax))
        output.append(qrfs)    
    output = np.array(output, dtype=object)
    dict_output = {'Comp.{}'.format(i):output[:,i] for i in range(output.shape[1])}
    return dict_output

    
def compare_instf_block(signal, method_output, tmin=None, tmax=None):
    X = signal.instf
    output = []
    for Xest in method_output:
        order = order_components(Xest, X, minormax = 'min', metric = mse)
        Xaux = Xest[order]
        qrfs = []
        for x,xaux in zip(X,Xaux):
            indx = np.where(np.abs(x)>0)
            # qrfs.append(compute_qrf(x[indx], xaux[indx],tmin=tmin,tmax=tmax))
            qrfs.append(mse(x[indx], xaux[indx]))
        output.append(qrfs)    
    output = np.array(output, dtype=object)
    dict_output = {'Comp.{}'.format(i):output[:,i] for i in range(output.shape[1])}
    return dict_output    
