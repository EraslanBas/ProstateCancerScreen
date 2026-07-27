import libraries

which = lambda lst:list(np.where(lst)[0])

# def ot_distance(exp1, exp2, N=None):
#     if N is not None:
#         s1 = np.random.choice(exp1.shape[0], N, replace=True)  
#         exp1 = exp1[s1,:]
#         s2 = np.random.choice(exp2.shape[0], N, replace=True)
#         exp2 = exp2[s2,:]
#     M = ot.dist(exp1, exp2)
#     #M /= M.max()
    
#     a = np.ones(len(exp1)) / len(exp1)
#     b = np.ones(len(exp2)) / len(exp2)
    
#     ot_distance = ot.emd(a, b, M)
#     return np.sum(ot_distance*M)

def violinPlotObsValue(columnName, groupName, adata, testedPairs):
    
    # Extract the gene expression and group
    df = pd.DataFrame({
        'expression': adataSub.obs[columnName],  # or .X.flatten() if dense
        'group': adataSub.obs[groupName].values
    })

    # Plot
    plt.figure(figsize=(8,6))
    ax = sns.violinplot(x="group", y="expression", data=df, inner="point")
    ax.set_ylabel("NE score")


    annotator = Annotator(ax, testedPairs, data=df, x="group", y="expression")
    annotator.configure(test='Mann-Whitney', text_format='star', loc='outside')
    annotator.apply_and_annotate()

    plt.show()

def violinPlotGeneWithSignificance(geneName, groupName, adata, testedPairs):
    
    # Extract the gene expression and group
    df = pd.DataFrame({
        'expression': adata[:, geneName].X.toarray().flatten(),  # or .X.flatten() if dense
        'group': adata.obs[groupName].values
    })

    # Plot
    plt.figure(figsize=(8,6))
    ax = sns.violinplot(x="group", y="expression", data=df, inner="point")
    ax.set_ylabel("ASCL1 expression")

    # Add significance test


    annotator = Annotator(ax, testedPairs, data=df, x="group", y="expression")
    annotator.configure(test='Mann-Whitney', text_format='star', loc='outside')
    annotator.apply_and_annotate()

    plt.show()

def plot_all_observations(y_true, y_pred, y_max, titleStr=""):
    
    y_pred_censored = np.where(y_pred > y_max, y_max, y_pred)
    
    plt.figure(figsize=(15, 15))
    plt.scatter(y_pred_censored, y_test, 
                c="blue", edgecolor="black", 
                label=f'R^2={round(r2_score(y_test, y_pred_censored),2)}')
    
    plt.title(titleStr, size=22)
    plt.xlabel("Predicted values", size=16)
    plt.ylabel("Observed true values", size=16)
    plt.ticklabel_format(style='plain')
    plt.legend()
    
    plt.show()
    
    
    
## Plot the coefficients in a beta matrix
def draw_coefficients(beta_df):
    figure, axis = plt.subplots(5,2, figsize=(15,15))


    for key, ax in zip(beta_df.columns, axis.ravel()):
        ax.set_title(key)
        sns.distplot(beta_df[key], 
                     ax=ax, 
                     bins=100, 
                     color="blue", 
                     kde=True, 
                     axlabel=False, 
                     hist_kws=dict(edgecolor="black"))

    plt.subplots_adjust(hspace=0.5)
    plt.show()
  
#
# x_scaler = StandardScaler()
# X_scaled = my_x_scaler.fit_transform(X)
# reconstruct the beta matrix for unscaled data if the beta_df
# comes from a model where X_scaled is used instead of X

def create_beta_df(beta_df, x_scaler, feature_names):        
    
    i = 0
    # for each column except the intercept and the additional sigma term
    for col in beta_df:
        if (col != 'beta_intercept'):
            if ('beta_' in col):

                # subtract the appropriate value from the intercept (review intercept final expression)
                beta_df['beta_intercept'] -= (beta_df[col] * x_scaler.mean_[i])/x_scaler.scale_[i]     
                
                # scale the coefficient (review each coefficient final expression)
                beta_df[col] /= x_scaler.scale_[i]
                i += 1

    return beta_df

def predict_linear_combination(beta_df, X):
    
    # Don't grab the last column, that is our estimate of the error standard deviation, "sigma"
    coefficients = beta_df.iloc[:, :-1].mean()

    # Find our linear combination again
    linear_combination = X.dot(coefficients[1:]) + coefficients.iloc[0]
    
    return linear_combination