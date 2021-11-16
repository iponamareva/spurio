# This is a set of functions used by gpm_train, gpm_test or both.

import numpy as np
import glob
import os
from sklearn.metrics import (roc_auc_score, accuracy_score, mean_squared_error,
                             roc_curve, auc, precision_recall_curve,
                             average_precision_score)
from scipy import interp

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


# def load_data(name_prefix, filedir='../output/np/', n_samples=None, userpath=0):
# changed on 25 february
def load_data(name_prefix, filedir='../output_feb/np/', n_samples=None, userpath=0):
    '''
    This function is used to collect the .npy output from the homology
    search and return it as one big matrix [n_proteins x n_features].
    '''

    print('filedir')
    print(filedir)
    if userpath==1:
        feature_files = glob.glob('{}*.npy'.format(filedir))[:n_samples]
    else:
        feature_files = glob.glob('{}*.npy'.format(os.path.join(filedir, name_prefix)))[:n_samples]

    print(len(feature_files))
    count = 0
    accessions = []

    feature_mat = np.zeros((len(feature_files), 4))

    for n, feature_file in enumerate(feature_files):
        protein_summary = np.load(feature_file)
        acc = feature_file.split('.')[-2].split('_')[-1]
        accessions.append(acc)
        feature_mat[n] = protein_summary
        count += 1
        if count % 200 == 0:
            print(count)
    notnan_idx = ~np.isnan(feature_mat[:, -1])

    return accessions, feature_mat[notnan_idx]

def prepare_train_test_data(n_train='all'):
    '''
    For training/testing, we use Swissprot/Antifam data.
    This function is used to load, merge and preprocess
    them.
    '''

    spnames, spdata = load_data('sp', n_samples=n_train)
    afnames, afdata = load_data('af', n_samples=n_train)

    if n_train == 'all':
        n_train = np.min((len(spdata), len(afdata)))
    spnames = spnames[:n_train]
    afnames = afnames[:n_train]
    spdata = spdata[:n_train]
    afdata = afdata[:n_train]

    X = np.vstack((spdata, afdata))
    y = np.concatenate((np.zeros(len(spdata)), np.ones(len(afdata))))
    names = np.concatenate((spnames, afnames))

    return X, y, names


def plot_pr_cv(ytest, ypred, aps=[], recalls=[], precisions=[]):
    #TODO Docstring
    plt.figure(2)
    ax = plt.subplot(111)
    mean_precision = np.linspace(0, 1, 100)
    precision, recall, _ = precision_recall_curve(ytest, ypred)
    recalls.append(interp(mean_precision, recall, precision))
    average_precision = average_precision_score(ytest, ypred)
    aps.append(average_precision)
    precisions.append(interp(mean_precision, precision, recall))
    ax.step(recall[:-1], precision[:-1], alpha=0.2, where='post')
    return aps, recalls, precisions


def plot_roc_cv(ytest, ypred, tprs=[], fprs=[], aucs=[]):
    #TODO Docstring
    plt.figure(1)
    ax1 = plt.subplot(111)
    mean_fpr = np.linspace(0, 1, 100)
    fpr, tpr, thresholds = roc_curve(ytest, ypred)
    tprs.append(interp(mean_fpr, fpr, tpr))
    tprs[-1][0] = 0.0
    roc_auc = auc(fpr, tpr)
    aucs.append(roc_auc)
    ax1.plot(fpr, tpr, lw=1, alpha=0.3)

    print(fprs, tprs, aucs)
    return tprs, fprs, aucs


def finish_ap_plot(tprs, fprs, aucs):
    #TODO Docstring
    plt.figure(2)
    ax = plt.subplot(1, 1, 1)
    mean_tpr = np.linspace(0, 1, 100)
    mean_fpr = np.mean(fprs, axis=0)
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    ax.step(mean_tpr[:-1], mean_fpr[:-1], color='b',
            label=r'Mean PR (AP = %0.3f $\pm$ %0.3f)' % (mean_auc, std_auc),
            lw=2, alpha=0.8)
    std_fpr = np.std(fprs, axis=0)
    fprs_upper = np.minimum(mean_fpr + std_fpr, 1)
    fprs_lower = np.maximum(mean_fpr - std_fpr, 0)
    plt.fill_between(mean_tpr[:-1], fprs_lower[:-1], fprs_upper[:-1],
                     color='grey', alpha=0.2, label=r'$\pm$ 1 std. dev.')

    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('Recovery')
    plt.ylabel('Precision')
    plt.legend(loc="lower left")
    plt.savefig('RPcurve.png')


def finish_roc_plot(tprs, fprs, aucs):
    #TODO Docstring
    plt.figure(1)
    ax1 = plt.subplot(1, 1, 1)
    ax1.plot([0, 1], [0, 1], linestyle='--', lw=2, color='grey', alpha=0.8)
    mean_fpr = np.linspace(0, 1, 100)
    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    ax1.plot(mean_fpr, mean_tpr, color='b',
             label=r'Mean ROC (AUC = %0.3f $\pm$ %0.3f)' % (mean_auc, std_auc),
             lw=2, alpha=0.8)

    std_tpr = np.std(tprs, axis=0)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    plt.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=0.2,
                     label=r'$\pm$ 1 std. dev.')

    # Hide the right and top spines
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)

    # Only show ticks on the left and bottom spines
    ax1.yaxis.set_ticks_position('left')
    ax1.xaxis.set_ticks_position('bottom')

    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.legend(loc="lower right")
    plt.savefig('ROCcurve.png')
    plt.close()


def preprocess_data(X, scaler):
    #TODO Docstring
    X = add_pseudocount(X)
    X = keep_desired_features(X)
    X = logit(X)
    X = scaler.transform(X)

    return X



# HMMMMMMM
# def preprocess_data(X):
#     #TODO Docstring
#     X = add_pseudocount(X)
#     X = keep_desired_features(X)
#     X = logit(X)
#     scaler = RobustScaler()
#     scaler.fit(X)

#     X = scaler.transrom(X)

#     return X, scaler


def add_pseudocount(X, pc=0.5):
    #TODO Docstring
    if np.shape(X) == (4,):
        X = X.reshape(1, 4)
        
    if X[0, 0] == 0:
        X = np.array(X, dtype=float)
    
    X[: ,0] += (pc / X[:, 3])

    if np.shape(X)[1] == 5:
        X[: ,4] += (0.5)
    return X


def keep_desired_features(X):
    #TODO Docstring
    return X[:, [0, 1, 2]]


# W's function
def normalise_(X, min_all, max_all):
    '''Normalise all features to 0-1'''
    if np.shape(X) == (3,):
        X = X.reshape((1, 3))
    for i in range(3):
        X[:, i] = (X[:, i] - min_all[i]) / (max_all[i] - min_all[i])
    return X


def logit(X):
    '''Apply logarithm to all features'''
    if np.shape(X) == (3,):
        X = X.reshape((1, 3))

    X = np.log10(X)
    return X


def eval_result(p_predicted, y_true):
    '''
    Given y_true and y_predicted, this function
    returns the area under the curve, accuracy
    and mean squared error
    '''

    mse = mean_squared_error(y_true, p_predicted)
    auc = roc_auc_score(y_true, p_predicted)
    acc = accuracy_score(y_true, p_predicted.round())

    return auc, acc, mse


def error_analysis(p_predicted, y_true):
    idx_FP = np.intersect1d(np.nonzero(p_predicted > 0.8), np.nonzero(y_true == 0))
    idx_FN = np.intersect1d(np.nonzero(p_predicted < 0.2), np.nonzero(y_true == 1))
    # FP_data = p_predicted[idx_FP]

    return idx_FP, idx_FN

def err_plot(X1, X2, figname):
    print('x1 shape ', X1.shape)
    print('x2 shape ', X2.shape)
    ax = plt.subplot(111)
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    for axis in ['bottom', 'left']:
            ax.spines[axis].set_linewidth(2)
    
    plt.xlabel('log(Stop codons per MSA residue), normalised')
    plt.ylabel("log(length in amino acidcs), normalised")
    
    ax.scatter(X2[:, 0], X2[:, 1], s=7, label='FN mispredictions', alpha=0.5)
    ax.scatter(X1[:, 0], X1[:, 1], s=7, label='FP mispredictions', alpha=0.5)
    
    plt.legend()
    # plt.show()

    plt.savefig(figname)
    plt.close()



def eval_binned_result(p_predicted, y_true, lengths, bins=12):

    ''' 
    Binned
    returns 3 x bins matrix 
    '''
    
    result = np.zeros((bins + 1, 3))
    conds = [40, 50, 60, 70, 80, 90, 100, 110, 120, 140, 160, 180, 200]

    for i in range(bins):
        down_cond = conds[i]
        up_cond = conds[i + 1]

        print('conditions: ' + str(down_cond) + ' ' + str(up_cond))

        mask = (lengths >= down_cond) & (lengths < up_cond)

        #print(mask)
        
        # if all elements are masked
        if np.count_nonzero(mask) == len(p_predicted):
            print('encountered empty bin while running with restrictions {}-{}'.format(down_cond, up_cond))
            result[i] = [None, None, None]
            continue
        
        masked_pred = p_predicted[mask]
        masked_true = y_true[mask]

        print('with this condition found sequences: ' + (str(len(masked_true))))
        spur_n = np.count_nonzero(masked_true)
        print('spurious out of them:' + str(spur_n))

        if ((spur_n/len(masked_true) > 0.9) | (spur_n/len(masked_true) < 0.1)):
            print('highly unbalanced class')
            continue
        
        if i == 1:
            for j in range(len(masked_pred)):
                print(str(masked_pred[i]) + ' ' + str(masked_true[i]))
                if j == 20:
                    break



        # masked_true = ma.masked_where((lengths >= down_cond and lengths < up_cond), y_true)

        mse = mean_squared_error(masked_true, masked_pred)
        auc = roc_auc_score(masked_true, masked_pred)
        acc = accuracy_score(masked_true, masked_pred.round())

        result[i] = [auc, acc, mse]

    mask = (lengths >= 200)
    if np.count_nonzero(mask) == len(p_predicted):
        print('encountered empty bin while running with last restriction')
        result[len(result) - 1] = [None, None, None]
    else:
        masked_pred = p_predicted[mask]
        masked_true = y_true[mask]
        mse = mean_squared_error(masked_true, masked_pred)
        auc = roc_auc_score(masked_true, masked_pred)
        acc = accuracy_score(masked_true, masked_pred.round())
        print(auc, acc, mse)
        
        for i in range(len(masked_pred)):
            print(str(masked_pred[i]) + ' ' + str(masked_true[i]))
            if i == 20:
                break


        result[len(result) - 1] = [auc, acc, mse]
    # result[len(result) - 1] = [1.0, 1.0, 1.0]



    return result


def binned_mean(data, bins=12, axis=0):

    ''' data is an array  splits x bins x 3 '''

    mask = (data == None)

    masked_data = np.ma.array(data, mask=mask)
    result = np.mean(masked_data, axis=0)

    return result


def plot_binned_auc(data, name='default_name', bins=12):
    name = name + '.png'

    x = np.arange(bins + 1)

    plt.figure(facecolor='white')
    plt.plot(x, data)
    matplotlib.rcParams.update({'font.size': 8})
    # plt.xticks(x, ('40-50', '50-60', '60-70', '70-80', '80-100', '100-120', '120-140', '140-160', '160-180', '180-200', '>200'))
    plt.xticks(x, ('40-50', '50-60', '60-70', '70-80', '80-90', '100-110', '110-120', '120-140', '140-160', '160-180', '180-200', '>200'))
    plt.savefig(name, dpi=200)









