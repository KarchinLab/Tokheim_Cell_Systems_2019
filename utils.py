import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_curve, auc, average_precision_score,precision_recall_curve
import glob
import os

def process_cgc(path):
    """Get the list of CGC genes with small somatic variants."""
    # read in data
    df = pd.read_table(path)

    # keep small somatic variants
    s = df['Mutation Types']
    #is_small = s.str.contains('Mis|F|N|S').fillna(False)
    is_small = s.str.contains('Mis').fillna(False)
    is_somatic = ~df['Tumour Types(Somatic)'].isnull()
    df = df[is_small & is_somatic].copy()

    # get gene names
    cgc_genes = df['Gene Symbol'].tolist()

    return cgc_genes


############################
# Function to read CHASM2 result
############################
def read_result(cancer_type, only_significant=False):
    """Read CHASM2 results"""
    # read mutations
    mut_df = pd.read_table('CHASM2/data/aggregated_results/{0}.maf'.format(cancer_type))
    mut_df['UID'] = range(len(mut_df))

    # read CHASM2 result
    useful_cols = ['UID', 'ID', 'CHASM2', 'CHASM2_genome', 'CHASM2_pval', 'CHASM2_genome_pval', 'CHASM2_qval', 'CHASM2_genome_qval']
    result_df = pd.read_table('CHASM2/data/aggregated_results/{0}.txt'.format(cancer_type), usecols=useful_cols)

    # merge mutation information into CHASM2 result
    mut_cols = ['UID', 'Hugo_Symbol', 'Transcript_ID', 'Protein_position', 'HGVSp_Short', 'CODE']
    result_df = pd.merge(result_df, mut_df[mut_cols], on='UID', how='left')

    if only_significant:
        result_df = result_df.rename(columns={'CHASM2_genome_qval': 'gwCHASM2 qvalue'})
        # merge full CHASM2 result into MAF dataframe
        result_df = pd.merge(mut_df, result_df, on=['Hugo_Symbol', 'Transcript_ID', 'HGVSp_Short'], how='left')
        result_df = result_df[result_df['gwCHASM2 qvalue']<=0.01]
    else:
        result_df = result_df.rename(columns={'CHASM2_genome_qval': cancer_type})

    return result_df

def read_all_results(base_dir='CHASM2/data/aggregated_results'):
    """Reads all of the results"""
    cancer_types = [os.path.basename(f)[:-4] for f in glob.glob('{0}/*.txt'.format(base_dir)) if 'PANCAN' not in f]
    merged_df = read_result('PANCAN')
    for c in cancer_types:
        tmp = read_result(c)
        tmp[c+'_flag'] = (tmp[c]<=.01).astype(int)
        merged_df = pd.merge(merged_df, tmp[['Hugo_Symbol', 'Transcript_ID', 'HGVSp_Short', c, c+'_flag']],
                             on=['Hugo_Symbol', 'Transcript_ID', 'HGVSp_Short'], how='left')
    return merged_df

############################
# Function to read rarity result
############################
def read_all_rarity_results(base_dir='CHASM2/data/rarity_analysis/'):
    """Reads all of the results"""
    cancer_types = [os.path.basename(f)[:-4] for f in glob.glob('{0}/*.txt'.format(base_dir))]
    result_list = []
    for c in cancer_types:
        tmp_path = os.path.join(base_dir, c+'.txt')
        tmp = pd.read_table(tmp_path)
        counts = tmp.groupby('category')['number of mutations'].sum()
        frac = counts / counts.sum()
        frac.name = c
        result_list.append(frac)
    concat_df = pd.concat(result_list, axis=1)
    concat_df = concat_df.fillna(0)
    return concat_df


############################
# Functions for ATM analysis
############################
def read_atm_result(chasm2_path, maf_path):
    """Read the CHASM2 result for ATM and merge information with mutations"""
    # read chasm2
    useful_cols = ['gene', 'UID', 'ID', 'driver score', 'CHASM2', 'CHASM2_genome', 'CHASM2_pval',
                   'CHASM2_genome_pval', 'CHASM2_qval', 'CHASM2_genome_qval']
    df = pd.read_table(chasm2_path, usecols=useful_cols)

    # read mutations
    mut_df = pd.read_table(maf_path)
    mut_df['UID'] = range(len(mut_df))

    # calculate mutation recurrence
    counts = mut_df.groupby(['Hugo_Symbol', 'HGVSp_Short'])['Tumor_Sample_Barcode'].nunique().reset_index(name='recurrence')
    mut_df = pd.merge(mut_df, counts, on=['Hugo_Symbol', 'HGVSp_Short'], how='left')

    # merge in the mutation data
    df = pd.merge(df, mut_df[['UID', 'Hugo_Symbol', 'Transcript_ID', 'HGVSp_Short', 'Protein_position', 'recurrence']], on='UID', how='left')
    chasm_cols = ['CHASM2', 'CHASM2_genome', 'CHASM2_pval', 'CHASM2_genome_pval', 'CHASM2_qval', 'CHASM2_genome_qval']
    mut_df = pd.merge(mut_df, df, on=['Hugo_Symbol', 'Transcript_ID', 'HGVSp_Short', 'recurrence'], how='left')

    # get all of the significant mutations
    df['Protein_change'] = df['Hugo_Symbol'] + '_' + df['Transcript_ID'] + '_' + df['HGVSp_Short']
    mut_df['Protein_change'] = mut_df['Hugo_Symbol'] + '_' + mut_df['Transcript_ID'] + '_' + mut_df['HGVSp_Short']
    is_signif = df['CHASM2_genome_qval']<=0.01
    is_atm = df['Hugo_Symbol'] == 'ATM'
    signif_atm_df = mut_df[mut_df.Protein_change.isin(df[is_signif & is_atm]['Protein_change'])].drop_duplicates(['Hugo_Symbol', 'HGVSp_Short'])
    is_atm2 = mut_df['Hugo_Symbol'] == 'ATM'
    is_missense = mut_df['Variant_Classification'] == 'Missense_Mutation'
    all_atm_df = mut_df[is_atm2 & is_missense]

    # fix x/y suffixes
    rename_dict = {'Protein_position_x': 'Protein_position'}
    all_atm_df = all_atm_df.rename(columns=rename_dict)
    signif_atm_df = signif_atm_df.rename(columns=rename_dict)

    return signif_atm_df, all_atm_df


############################
# Functions for ROC plot
############################
def roc_plot(data, methods, other_methods):
    """Create a receiver operating characteristic curve of methods."""
    y = data['y']
    # plot the low performing methods
    for method in other_methods:
        if method == 'CanDrA.plus': method = 'CanDrA plus'
        pred = data[method].astype(float).dropna()
        fpr, tpr, thresholds = roc_curve(y[pred.index], pred)
        myauc = auc(fpr, tpr)
        plt.plot(fpr, tpr,
                 color='lightgray')
        plt.xlabel('False positive rate')
        plt.ylabel('True positive rate')

    # plot the top methods
    for method in methods:
        if method == 'CanDrA.plus': method = 'CanDrA plus'
        pred = data[method].astype(float).dropna()
        fpr, tpr, thresholds = roc_curve(y[pred.index], pred)
        myauc = auc(fpr, tpr)
        plt.plot(fpr, tpr,
                 label='{0} (area = {1:0.3f})'.format(method, myauc))
        plt.legend(loc='best')
        plt.xlabel('False positive rate')
        plt.ylabel('True positive rate')

    # format axis
    plt.gca().xaxis.set_ticks_position('bottom')
    plt.gca().yaxis.set_ticks_position('left')

def top5(data):
    """Figure out the top 5 best beforming methods."""
    replace_dict = {'CHASM2_genome': 'CHASM2', 'CanDrA.plus': 'CanDrA',
                    'Polyphen2_hdiv': 'Polyphen2', 'Polyphen2_hvar': 'Polyphen2',
                    }
    data['method'] = data.method.replace(replace_dict)
    tmp_top_methods = data.groupby('method')['auc'].max().sort_values(ascending=False).head(5).index.tolist()
    top_methods = data[data.method.isin(tmp_top_methods)].sort_values('auc', ascending=False).index.tolist()
    return top_methods


############################
# Functions for a p-value QQ plot
############################

def fix_formatting(ax, title, remove_xlab=True, remove_ylab=False):
    """simple function to set title and remove xlabel"""
    ax.set_title(title)
    if remove_xlab:
        ax.set_xlabel('')
    if remove_ylab:
        ax.set_ylabel('')

def set_axes_label(fig, xlab, ylab,
                   ylab_yoffset=.55, ylab_xoffset=0.04,
                   xlab_yoffset=.04, xlab_xoffset=0.5):
    txt1 = fig.text(xlab_xoffset, xlab_yoffset, xlab, ha='center', size=22)
    txt2 = fig.text(ylab_xoffset, ylab_yoffset, ylab, ha='center', size=22, rotation=90)
    return txt1, txt2

def qqplot(data,
           ax=None,
           log=False, title=None,
           use_xlabel=True, use_ylabel=True,
           **kwargs):
    """qq-plot with uniform distribution"""
    tmp = data.copy()
    tmp.sort_values(inplace=True)
    dist_quant = np.arange(1, len(tmp)+1)/float(len(tmp)+1)
    if log:
        log_quant = -np.log10(dist_quant)
        if ax is None:
            plt.plot(log_quant, -np.log10(tmp),'o', markersize=3, **kwargs)
            plt.plot([0, log_quant[0]], [0, log_quant[0]], ls="-", color='red')
        else:
            ax.plot(log_quant, -np.log10(tmp),'o', markersize=3, **kwargs)
            ax.plot([0, log_quant[0]], [0, log_quant[0]], ls="-", color='red')
        # set axis labels
        if use_xlabel:
            if ax is None: plt.xlabel('Theoretical (-log10(p))')
            else: ax.set_xlabel('Theoretical (-log10(p))')
        if use_ylabel:
            if ax is None: plt.ylabel('Observed (-log10(p))')
            else: ax.set_ylabel('Observed (-log10(p))')
    else:
        if ax is None:
            plt.plot(dist_quant, tmp,'o', markersize=3, **kwargs)
            plt.plot([0, 1], [0, 1], ls="-", color='red')
        else:
            ax.plot(dist_quant, tmp,'o', markersize=3, **kwargs)
            ax.plot([0, 1], [0, 1], ls="-", color='red')
            ax.set_ylabel('p-value')
        if use_xlabel:
            if ax is None: plt.xlabel('Theoretical p-value')
            else: ax.set_xlabel('Theoretical p-value')
        if use_ylabel:
            if ax is None: plt.ylabel('Observed p-value')
            else: ax.set_ylabel('Observed p-value')
    if title:
        ax.set_title(title)
    sns.despine()


def mean_log_fold_change(data):
    """Mean log fold change function

    Parameters
    ----------
    data : pd.Series
        a series of p-values

    Returns
    -------
    mlfc : float
        mean log fold change.
    """
    tmp = data.copy()
    tmp.sort_values(ascending=True, inplace=True)
    tmp[tmp==0] = tmp[tmp>0].min()  # avoid infinity in log by avoiding zero pvals
    dist_quant = np.arange(1, len(tmp)+1)/float(len(tmp))
    mlfc = np.mean(np.abs(np.log2(tmp/dist_quant)))
    return mlfc
