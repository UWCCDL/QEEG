# Alpha

import numpy as np
import scipy.stats as stat

samples = 100

def generate_matrix(num, rho, samples=samples):
    """Produces a list of variables with a mean correlation of rho"""
    r = np.sqrt(rho)
    x1 = stat.uniform.rvs(size = samples)
    x2 = stat.uniform.rvs(size = samples)
    m = np.zeros((samples, num))
    for j in range(num):
        m[:,j] = r * x1 + np.sqrt(1 - r**2) * x2
        x2 = stat.uniform.rvs(size = samples)
    return m

def mean_r(mat):
    """Mean correlation across columns of a matrix"""
    n = mat.shape[1]
    val = 0.0
    m = np.zeros((n,n))
    for i in range(n):
        for j in range(i):
            val += stat.pearsonr(mat[:,i], mat[:,j])[0]
            m[i,j] = stat.pearsonr(mat[:,i], mat[:,j])[0]
    return val / (n * (n - 1) / 2)

def test_sims(samples=1000, pval=0.05):
    """Tests that the math is correct. Should produce values ~ pval"""
    x1 = stat.uniform.rvs(size = samples)
    counter = 0
    for i in range(samples):
        x2 = stat.uniform.rvs(size = samples)
        r, p = stat.pearsonr(x1, x2)
        if p <= pval:
            counter += 1

    return (counter * 1.0) / samples

def test_matrix(num=14, rho=0.5, pval=0.05, samples=100, n=10000):
    res = 0.0
    c2 = 0.0
    c5 = 0.0
    c10 = 0.0
    for i in range(n):
        y = stat.uniform.rvs(size = samples)
        m = generate_matrix(num = num, rho = rho, samples = samples)
        count = 0
        for j in range(num):
            r, p = stat.pearsonr(y, m[:,j])
            if p <= pval:
                count += 1
        res += count
        if count >= 2:
            c2 += 1
        if count >= 5:
            c5 += 1
        if count >= 10:
            c10 += 1
    return np.array([res, c2, c5, c10])/n


if __name__ == "__main__":
    num_channels = int(sys.argv[1])
    rho = float(sys.argv[2])
    run_sims(num_channels, rho)
