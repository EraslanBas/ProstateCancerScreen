import numpy as np
from sklearn.linear_model import LinearRegression

def multivariate_r2(Y, X):
    """
    Computes the proportion of total variance in Y explained by covariates X
    (multivariate R^2 via linear projection).
    
    Parameters:
    - Y: (n_samples, n_features) array-like response matrix
    - X: (n_samples,) or (n_samples, n_covariates) covariate(s)
    
    Returns:
    - explained_variance: float, proportion of variance in Y explained by X
    """
    Y = np.asarray(Y)
    X = np.asarray(X)

    # Ensure X is 2D
    if X.ndim == 1:
        X = X.reshape(-1, 1)

    # Fit linear regression model of Y ~ X
    model = LinearRegression()
    model.fit(X, Y)
    Y_hat = model.predict(X)

    # Center Y
    Y_mean = Y.mean(axis=0)
    Y_centered = Y - Y_mean
    
    # Residuals and total variance
    residuals = Y - Y_hat
    ss_res = np.sum(residuals ** 2) 
    ss_total = np.sum(Y_centered ** 2) 

    return 1 - (ss_res / ss_total)