import numpy as np
def scatter_probability(E_i, E_j, alpha) :
    p = (1.0/E_j/(1.0-alpha)) * 1.0*((E_i >= alpha*E_j))
    return p 
def compute_spectrum(E, Sigma_t, Sigma_s, alpha) :
    N_E = len(E)
    phi = np.zeros(N_E)
    phi[N_E-1] = 1.0
    for i in range(N_E-2, -1, -1) :
        Q_i = 0.0
        for j in range(N_E-1, i, -1) :
            dE = E[j] - E[j-1]
            E_bar = np.sqrt(E[j]*E[j-1])
            Q_i += phi[j] * Sigma_s * scatter_probability(E[i], E[j], alpha) * dE
        phi[i] = Q_i / Sigma_t
        print i, phi[i]
    return phi
E = np.logspace(-1, 2, 1000)
phi = compute_spectrum(E, 1.0, 1.0, 0.0) 