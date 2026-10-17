"""Minimal continuous-time colonization and floral-competence model.

State order for each region: P (pigmented), Wp (white, competent),
Wm (white, degraded). Rates are per arbitrary time unit; all numerical
examples are synthetic. Occupancy/lineage branching is outside this toy model.
"""
import numpy as np

STATES = ("P", "Wp", "Wm")

def generator(regions, dispersal, loss, reexpression, degradation, repair=0.0):
    """Row generator Q; p(t)=p(0)@expm(Q*t).

    dispersal: dict of (origin,destination)->rate.
    loss, reexpression, degradation, repair: scalar or region-keyed dict.
    Each event is a state transition; no population-level colonization is inferred.
    """
    n = len(regions)
    Q = np.zeros((3*n, 3*n), dtype=float)
    def rate(v, region):
        return float(v[region] if isinstance(v, dict) else v)
    def add(a,b,v):
        if v < 0 or not np.isfinite(v):
            raise ValueError("Rates must be finite and nonnegative")
        if a != b:
            Q[a,b] += v
            Q[a,a] -= v
    for i,r in enumerate(regions):
        add(3*i, 3*i+1, rate(loss,r))
        add(3*i+1,3*i,rate(reexpression,r))
        add(3*i+1,3*i+2,rate(degradation,r))
        add(3*i+2,3*i+1,rate(repair,r))
    for (a,b),v in dispersal.items():
        if a == b:
            raise ValueError("Dispersal needs distinct regions")
        i,j=regions.index(a),regions.index(b)
        for k in range(3):
            add(3*i+k,3*j+k,float(v))
    assert np.allclose(Q.sum(axis=1),0)
    return Q

def distribution(Q, p0, time):
    from scipy.linalg import expm
    p0=np.asarray(p0,dtype=float)
    if time < 0 or not np.isclose(p0.sum(),1) or np.any(p0 < 0):
        raise ValueError("Time >= 0 and p0 must be a probability vector")
    return p0 @ expm(Q*time)

def observable_colour(p, regions):
    return {r:{"pigmented":float(p[3*i]),
               "white":float(p[3*i+1]+p[3*i+2]),
               "competent_white":float(p[3*i+1])}
            for i,r in enumerate(regions)}

def conditional_reexpression(t,H,r,delta):
    """Successful reexpression by horizon H given still white at t;
    initial W+, no repair, immigration or new white input."""
    if min(t,H,r,delta)<0: raise ValueError("Nonnegative arguments required")
    import math
    z=math.exp(-(r+delta)*t)
    den=delta+r*z
    if den==0: return 0.0
    return (r*z/den)*(-math.expm1(-(r+delta)*H))

if __name__=="__main__":
    regions=["source","island"]
    Q=generator(regions,{("source","island"):0.3,("island","source"):0.05},
                loss=0.2,reexpression=0.1,
                degradation={"source":0.01,"island":0.25})
    p0=np.array([1.,0,0,0,0,0])
    for t in [0,2,8,20]:
        print(t,observable_colour(distribution(Q,p0,t),regions))
