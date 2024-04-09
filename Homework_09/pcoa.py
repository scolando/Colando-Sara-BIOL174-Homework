import pandas as pd
import numpy as np
from scipy.linalg import eigh
import scipy.spatial

def calc_distance(df, distance="euclidean"):
    """Function to calculate a distance matrix dataframe from a dataframe

    Args:
        df (dataframe): DataFrame containing the variables in columns and observations in rows.
        distance (str, optional): Metric to be used to calcualte pair-wise distances. Uses scipy.spatial.distance.pdist to calculate distance. The distance function can be ‘braycurtis’, ‘canberra’, ‘chebyshev’, ‘cityblock’, ‘correlation’, ‘cosine’, ‘dice’, ‘euclidean’, ‘hamming’, ‘jaccard’, ‘jensenshannon’, ‘kulsinski’, ‘kulczynski1’, ‘mahalanobis’, ‘matching’, ‘minkowski’, ‘rogerstanimoto’, ‘russellrao’, ‘seuclidean’, ‘sokalmichener’, ‘sokalsneath’, ‘sqeuclidean’, ‘yule’. Defaults to "euclidean".
    """


    cols = list(df.columns)
    data = df.T.copy()


    dist = scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(data,metric = distance))

    dist = pd.DataFrame(dist,index=cols, columns = cols)
    return dist


def pcoa(df):
    res = {}
    names = list(df.columns)
    d_matrix = df.values
    d_matrix = d_matrix * d_matrix / -2

    row_means = d_matrix.mean(axis=1, keepdims=True)
    col_means = d_matrix.mean(axis=0, keepdims=True)
    matrix_mean = d_matrix.mean()
    result = d_matrix - row_means - col_means + matrix_mean
    result

    eigvals, eigvecs = eigh(result)
    idxs_descending = eigvals.argsort()[::-1]
    eigvals = eigvals[idxs_descending]
    eigvecs = eigvecs[:, idxs_descending]

    num_positive = (eigvals >= 0).sum()
    eigvecs[:, num_positive:] = np.zeros(eigvecs[:, num_positive:].shape)
    eigvals[num_positive:] = np.zeros(eigvals[num_positive:].shape)

    sum_eigenvalues = np.sum(eigvals)
    proportion_explained = eigvals / sum_eigenvalues

    number_of_dimensions = len(eigvals)
    eigvecs = eigvecs[:, :number_of_dimensions]
    eigvals = eigvals[:number_of_dimensions]
    explain = proportion_explained[:number_of_dimensions]

    coords = eigvecs * np.sqrt(eigvals)

    axis_labels = ["PC%d" % i for i in range(1, number_of_dimensions + 1)]
    res['eigvals'] = pd.Series(eigvals, index=axis_labels),
    res['loadings'] = pd.DataFrame(coords, index = names,columns=axis_labels)
    res['eigvecs'] = pd.DataFrame(eigvecs, index = names,columns=axis_labels)
    res['explain'] = pd.Series(explain,index=axis_labels)
    res['labels'] = names
    res['data'] = df


    return(res)