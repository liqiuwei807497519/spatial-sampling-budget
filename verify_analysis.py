"""Independent small-frame enumeration and seeded Monte Carlo verification."""
import itertools,math
from decimal import Decimal,localcontext
import numpy as np,pandas as pd
from pathlib import Path
from multi_target_and_economics import fixed_request_probability,evaluate,k_threshold,resources,repeat_count
def enumerate_fixed(L,Y,k):
    distributions=[]
    for l,y in zip(L,Y):
        mol=[j for j,n in enumerate(y) for _ in range(n)]+[-1]*(l-sum(y))
        distributions.append([set(mol[i] for i in draw) for draw in itertools.combinations(range(l),k)])
    combinations=list(itertools.product(*distributions))
    return sum(set(range(len(Y[0])))<=set.union(*choice) for choice in combinations)/len(combinations)
def main():
    L=np.array([4,5,6,3]);Y=np.array([[1,1],[1,3],[2,1],[1,0]])
    for k in range(1,5):
        raw=[]
        for choose in itertools.combinations(range(4),2):
            i=list(choose);l=L[i];y=Y[i]
            if min(l-y.sum(1))<k:raw.append(0);continue
            exact=enumerate_fixed(l,y,k)
            assert abs(fixed_request_probability(l,y,k)-exact)<1e-12;raw.append(exact)
        assert abs(evaluate(L,Y,2,k)[2]-np.mean(raw))<1e-12
    for tau in [.05,.12,.25,.4,.8]:
        expected=next((k for k in range(1,7) if evaluate(L,Y,2,k)[2]>=tau),None)
        assert k_threshold(L,Y,2,tau)[0]==expected
    assert fixed_request_probability([0,2],[[0,0],[1,1]],2)==0
    assert fixed_request_probability([0,2],[[0,0],[1,1]],2,False)==1
    assert repeat_count(0)[1]=='no finite count' and repeat_count(1)==(1,'finite')
    rng=np.random.default_rng(42);L=rng.integers(80,300,size=12)
    Y=np.column_stack([rng.binomial(L,p) for p in [.03,.02,.015,.01]]);Y=np.minimum(Y,(L//5)[:,None])
    exact=fixed_request_probability(L,Y,50);acc=np.zeros((100000,4),bool)
    for l,y in zip(L,Y):acc|=(rng.multivariate_hypergeometric([*y,l-sum(y)],50,size=100000)[:,:4]>0)
    simulated=acc.all(1).mean();assert abs(simulated-exact)<5*np.sqrt(exact*(1-exact)/100000)
    print(f'M=4: exact {exact:.9f}; Monte Carlo {simulated:.9f} (100000 requests).')
if __name__=='__main__':main()
