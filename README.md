# Spectral-Clustering
Normalized Spectral Clustering We present the Normalized Spectral Clustering algorithm. Given a set of n points X = x1, x2, . . . , xN ∈ R^d

The Normalized Spectral Clustering Algorithm:

1: Form the weighted adjacency matrix W from X

2: Compute the normalized graph Laplacian Lnorm

3: Determine k and obtain the first k eigenvectors u1, . . . , uk of Lnorm

4: Let U ∈ R^(n×k) be the matrix containing the vectors u1, . . . , uk as columns

5: Form the matrix T ∈ R^(n×k) from U by renormalizing each of U’s rows to have unit length, that is
set t(ij) = u(ij)/(Σ(u^2(ij))^0.5

6: Treating each row of T as a point in R^k, cluster them into k clusters via the K-means algorithm

7: Assign the original point xi to cluster j if and only if row i of the matrix T was assigned to cluster j

The project returns step 6

side note:u^2(ij) = (u_(ij))^2
