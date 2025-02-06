from simeloss import GetDeltaE
import numpy as np
import time


A, Z, Q, energy, Target = 238, 92, 86, 250, "Au"
# A, Z, Q, energy, Target = 238, 92, 86, 345, "Au"
thickness = 0.1

start = time.time()
dEtotal, dEcol, dEcc, charges, histories = GetDeltaE(A, Z, Q, energy, Target, thickness)
print("Elapsed time:", time.time()-start)

print(np.mean(dEtotal))
print(np.mean(dEcol))
print(np.mean(dEcc))
print(np.mean(charges[Q]))

if (A, Z, Q, energy, Target, thickness) == (238, 92, 86, 345, "Au", 0.1):
    assert np.mean(dEtotal) ==135.6391976008193
    assert np.mean(dEcol) ==135.6392022497982
    assert np.mean(dEcc) ==135.6370945870273
    assert np.mean(charges[Q]) ==0.04914285714285714
    print("Passed")