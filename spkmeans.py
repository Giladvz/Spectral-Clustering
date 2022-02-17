from enum import Enum
import sys
import numpy as np
import pandas as pd
import spkmeansmodule

# Calculates probability for kmeans++ algo
def prob (dataframe,centroids,j):
    # Destination function
    dist = lambda x,y: sum(f(x,y) for f in (lambda x,y: ((x.to_list()[i]-y[i])**2) for i in range(0,len(y))))
    if (j==0):
        # Creates dis on first run
        dataframe['dis'] = dataframe.apply(dist,args=([centroids[0]]),axis=1)
        
    else:
        # Saves the minimum between all runs
        func = lambda row,centroids, j: min(dist(row,centroids[j]),row['dis'])
        dataframe['dis'] = dataframe.apply(func,args=(centroids,j),axis=1)
    total = dataframe['dis'].sum()
    dataframe['prob'] = (dataframe['dis']/total)

# Choses the wanted k centroids
def choseCent (inputpan,k):
    centroids = [[0 for x in range(d)] for y in range(k)] 
    cent = np.random.choice(list(inputpan.index.values))
    chosen = []
    chosen.append(cent)
    centroids[0] = inputpan.loc[cent,:].to_list()
    # Chooses all the other centroids to pass to c
    for j in range(1,k):
        prob(inputpan,centroids,j-1)
        cent = np.random.choice(list(inputpan.index.values),
            p=inputpan['prob'].to_list())
        centroids[j] = inputpan.loc[cent,:].to_list()[:-2]
        chosen.append(cent)
    return centroids,chosen

if __name__ == '__main__':
    k = int(sys.argv[1])
    np.random.seed(0)
    goal = sys.argv[2]
    inputnum = np.genfromtxt(sys.argv[3], delimiter=",")
    if (inputnum.ndim == 1):
        inputnum = inputnum.reshape(len(inputnum),1)
    inputpan = pd.DataFrame(data=inputnum[0:,0:])
    d = len(inputnum[0])

    ret = spkmeansmodule.fit(inputnum.tolist(),goal,k,d)
    ret = pd.DataFrame(ret)
    print(ret)
    if (sys.argv[2] == "spk"):
        k = len(ret.iloc[0])
        centroids, chosen = choseCent(ret,k)
        print(*chosen,sep=",")
        if (k!=1):
            ret = ret.drop('dis',1)
            ret = ret.drop('prob',1)
        ret = ret.to_numpy()
        spkmeansmodule.Kfit(ret.tolist(),centroids)
