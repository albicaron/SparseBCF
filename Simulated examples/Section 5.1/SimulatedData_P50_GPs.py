# Importing packages
import numpy as np
import pandas as pd
import time
import gc
from copulae import GaussianCopula

from sklearn.neighbors import KNeighborsRegressor
from scipy import stats as sts

# Load local package
from causal_models import CMGP


# # # FUNCTION DEFINITION
# PEHE
def PEHE(T_true, T_est):
    return np.sqrt(np.mean((T_true.reshape((-1, 1)) - T_est.reshape((-1, 1))) ** 2))


# Monte Carlo Standard Error
def MC_se(x, B):
    return sts.t.ppf(0.975, B - 1) * np.std(np.array(x)) / np.sqrt(B)


# Data Generating Process
def get_data(N, P, rng):
    # Set seed
    np.random.seed(rng+100)

    # Generate X with a little bit of correlation b/w the continuous variables (i.e. 0.3 Pearson coeffincient c.a.)
    X = np.zeros((N, P))

    # Generate Gaussian copula correlated uniforms
    g_cop = GaussianCopula(dim=P)

    mycova = np.ones((P, P))

    for i in range(P):
        for j in range(P):

            mycova[i, j] = 0.3**(np.abs(i- j)) + np.where(i == j, 0, 0.1)

    g_cop[:] = mycova

    rv = g_cop.random(N)

    # Generate correlated covariates
    X[:, 0:5] = np.asarray([sts.norm.ppf(rv[:, i]) for i in range(5)]).T
    X[:, 5:P] = np.asarray([sts.binom.ppf(rv[:, i], 1, 0.3) for i in range(5, P)]).T


    # Generate Z
    und_lin = -0.5 + 0.2*X[:, 0] + 0.1*X[:, 1] + 0.4*X[:, 20] + sts.uniform.rvs(size=N)/10
    pscore = sts.norm.cdf(und_lin)

    Z = sts.binom.rvs(1, pscore)

    # Generate Y
    mu = 3 + 1.5*np.sin(np.pi*X[:, 0]) + 0.5*(X[:, 1] - 0.5)**2 + 1.5*(2 - np.abs(X[:, 2])) + 1.5*X[:, 3]*(X[:, 20] + 1)
    ITE = 0.1 + 1*np.abs(X[:, 0] - 1)*(X[:, 20] + 2)

    sigma = np.std(mu)/2
    Y = mu + ITE*Z + sts.norm.rvs(0, sigma, N)

    return Y, X, Z, ITE


# Train-Test Split
def tr_te_split(X, split):
    train = np.array(X[split])
    test = np.array(X[~split])

    return train, test


# T-KNN
def T_KNN(x_train, z_train, y_train, x_test, k=7):
    x_train_1 = x_train[z_train == 1]
    x_train_0 = x_train[z_train == 0]

    y_train_1 = y_train[z_train == 1]
    y_train_0 = y_train[z_train == 0]

    myKNN_1 = KNeighborsRegressor(n_neighbors=k)
    myKNN_1.fit(x_train_1, y_train_1)

    myKNN_0 = KNeighborsRegressor(n_neighbors=k)
    myKNN_0.fit(x_train_0, y_train_0)

    Y1_train = myKNN_1.predict(x_train)
    Y0_train = myKNN_0.predict(x_train)

    Y1_test = myKNN_1.predict(x_test)
    Y0_test = myKNN_0.predict(x_test)

    ITE_est_train = Y1_train - Y0_train
    ITE_est_test = Y1_test - Y0_test

    return ITE_est_train, ITE_est_test


# Options
N, P, B = [1000, 50, 1000]

# Results storage
Train_PEHE_KNN = list(); Test_PEHE_KNN = list()
Train_PEHE_CMGP = list(); Test_PEHE_CMGP = list()
Train_PEHE_NSGP = list(); Test_PEHE_NSGP = list()


# Simulation Study
start = time.time()

for i in range(B):

    print("\n*** Iteration", i+1, "\n")

    # Generate Data
    myY, myX, myZ, ITE = get_data(N, P, i)

    # Train-Test Split (70-30%)
    split = np.random.choice(np.array([True, False]), N, replace=True, p=np.array([0.7, 0.3]))

    x_train, x_test = tr_te_split(myX, split)
    y_train, y_test = tr_te_split(myY, split)
    z_train, z_test = tr_te_split(myZ, split)
    ITE_train, ITE_test = tr_te_split(ITE, split)


    # 1) KNN (T-KNN)
    ITE_est_train, ITE_est_test = T_KNN(x_train, z_train, y_train, x_test, k=7)

    Train_PEHE_KNN.append(PEHE(ITE_train, ITE_est_train))
    Test_PEHE_KNN.append(PEHE(ITE_test, ITE_est_test))


    # 2) CMGP
    myCMGP = CMGP(dim=P, mode="CMGP", mod="Multitask", mkern="LCM")
    myCMGP.fit(X=x_train, Y=y_train, W=z_train)

    ITE_est_train = myCMGP.predict(x_train)[0]
    ITE_est_test = myCMGP.predict(x_test)[0]

    Train_PEHE_CMGP.append(PEHE(ITE_train, np.array(ITE_est_train)))
    Test_PEHE_CMGP.append(PEHE(ITE_test, np.array(ITE_est_test)))

    # 3) NSGP
    myNSGP = CMGP(dim=P, mode="NSGP", mod='Multitask', kern='Matern')
    myNSGP.fit(X=x_train, Y=y_train, W=z_train)

    ITE_est_train = myNSGP.predict(x_train)[0]
    ITE_est_test = myNSGP.predict(x_test)[0]

    Train_PEHE_NSGP.append(PEHE(ITE_train, np.array(ITE_est_train)))
    Test_PEHE_NSGP.append(PEHE(ITE_test, np.array(ITE_est_test)))

    # Garbage Collection
    gc.collect()



elapsed = time.time() - start
print("\n\nElapsed time (in h) is", round(elapsed/3600, 2))


# PEHE results
PEHE_final = {'KNN': [np.mean(np.array(Train_PEHE_KNN)), MC_se(Train_PEHE_KNN, B),
                      np.mean(np.array(Test_PEHE_KNN)), MC_se(Test_PEHE_KNN, B)],
              'CMGP': [np.mean(np.array(Train_PEHE_CMGP)), MC_se(Train_PEHE_CMGP, B),
                       np.mean(np.array(Test_PEHE_CMGP)), MC_se(Test_PEHE_CMGP, B)],
              'NSGP': [np.mean(np.array(Train_PEHE_NSGP)), MC_se(Train_PEHE_NSGP, B),
                       np.mean(np.array(Test_PEHE_NSGP)), MC_se(Test_PEHE_NSGP, B)]}

PEHE_final = pd.DataFrame(PEHE_final, index=["Train", "SE_Train", "Test", "SE_Test"])


PEHE_Single = {'KNN_Train': np.array(Train_PEHE_KNN),
               'KNN_Test': np.array(Test_PEHE_KNN),
               'CMGP_Train': np.array(Train_PEHE_CMGP),
               'CMGP_Test': np.array(Test_PEHE_CMGP),
               'NSGP_Train': np.array(Train_PEHE_NSGP),
               'NSGP_Test': np.array(Test_PEHE_NSGP)}
PEHE_Single = pd.DataFrame(PEHE_Single)


# Saving
PEHE_Single.to_csv("YourDirectory/PEHE_All_P%i_B%i_GPs.csv" % (P, B),
                  index=False, header=True)


PEHE_final.to_csv("YourDirectory/PEHE_Final_P%i_B%i_GPs.csv" % (P, B),
                  index=True, header=True)


print("\n\n++++++++  FINISHED  +++++++++")
