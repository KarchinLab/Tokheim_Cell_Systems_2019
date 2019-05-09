import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_curve, auc, average_precision_score,precision_recall_curve
from matplotlib_venn import venn2, venn2_circles, venn3, venn3_circles
import glob
import os
import bisect

# list of methods
methods = ['VEST', 'CADD', 'Polyphen2_hdiv', 'Polyphen2_hvar', 'SIFT', 'MutationAssessor', 'REVEL', 'MCAP',
           'ParsSNP', 'CHASM', 'raw CHASMplus', 'gwCHASMplus', 'CanDrA', 'CanDrA plus', 'TransFIC', 'FATHMM']

def process_cgc(path, return_dataframe=False):
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
    if not return_dataframe:
        cgc_genes = df['Gene Symbol'].tolist()
    else:
        cgc_genes = df

    return cgc_genes


############################
# Function to read CHASM2 result
############################
def read_result(cancer_type, 
                only_significant=False,
                change_col_name=True,
                base_dir='CHASMplus/data/aggregated_results'):
    """Read CHASM2 results"""
    # read mutations
    mut_df = pd.read_table('{0}/{1}.maf'.format(base_dir, cancer_type))
    mut_df['UID'] = range(len(mut_df))
    
    # some minor renaming of columns
    rename_dict = {'Protein_Change': 'HGVSp_Short'}
    mut_df = mut_df.rename(columns=rename_dict)

    # read CHASM2 result
    useful_cols = ['UID', 'ID', 'CHASM2', 'CHASM2_genome', 'CHASM2_pval', 'CHASM2_genome_pval', 'CHASM2_qval', 'CHASM2_genome_qval']
    result_df = pd.read_table('{0}/{1}.txt'.format(base_dir, cancer_type), usecols=useful_cols)

    # merge mutation information into CHASM2 result
    mut_cols = list(set(['UID', 'Hugo_Symbol', 'Transcript_ID', 'Protein_position', 'HGVSp_Short', 'CODE']) & set(mut_df.columns))
    result_df = pd.merge(result_df, mut_df[mut_cols], on='UID', how='left')

    if only_significant:
        if change_col_name:
            result_df = result_df.rename(columns={'CHASM2_genome_qval': 'gwCHASM2 qvalue'})
        # merge full CHASM2 result into MAF dataframe
        result_df = pd.merge(mut_df, result_df, on=['Hugo_Symbol', 'Transcript_ID', 'HGVSp_Short'], how='left')
        result_df = result_df[result_df['gwCHASM2 qvalue']<=0.01]
    else:
        if change_col_name:
            result_df = result_df.rename(columns={'CHASM2_genome_qval': cancer_type})

    return result_df

def read_all_results(base_dir='CHASMplus/data/aggregated_results'):
    """Reads all of the results"""
    # read in pan-cancer 
    merged_df = read_result('PANCAN')
    merged_df['PANCAN_flag'] = (merged_df['PANCAN']<=0.01).astype(int)
    
    # add cancer types
    cancer_types = [os.path.basename(f)[:-4] for f in glob.glob('{0}/*.txt'.format(base_dir)) if 'PANCAN' not in f]
    for c in cancer_types:
        tmp = read_result(c)
        tmp[c+'_flag'] = (tmp[c]<=.01).astype(int)
        merged_df = pd.merge(merged_df, tmp[['Hugo_Symbol', 'Transcript_ID', 'HGVSp_Short', c, c+'_flag']],
                             on=['Hugo_Symbol', 'Transcript_ID', 'HGVSp_Short'], how='left')
    
    # create a flag for significant in any cancer type specific analysis
    cancer_type_flags = [
        'ACC_flag','BLCA_flag', 'BRCA_flag', 'CESC_flag', 'CHOL_flag', 
        'COAD_flag', 'DLBC_flag', 'ESCA_flag', 'GBM_flag', 'HNSC_flag',
        'KICH_flag', 'KIRC_flag', 'KIRP_flag', 'LAML_flag', 'LGG_flag', 
        'LIHC_flag', 'LUAD_flag', 'LUSC_flag', 'MESO_flag', 'OV_flag',
        'PAAD_flag', 'PCPG_flag', 'PRAD_flag', 'READ_flag', 'SARC_flag', 
        'STAD_flag', 'TGCT_flag', 'THCA_flag', 'THYM_flag', 'UCEC_flag', 
        'UCS_flag', 'UVM_flag' 
    ]
    is_signif_cancer_type = (merged_df[cancer_type_flags].sum(axis=1)>=1).astype(int)
    merged_df['Any_cancer_type_flag'] = is_signif_cancer_type
    
    return merged_df

############################
# Function to read rarity result
############################
def read_all_rarity_results(base_dir='CHASMplus/data/rarity_analysis/'):
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

def read_all_rarity_count(base_dir='CHASMplus/data/rarity_analysis/'):
    """Reads all of the results"""
    cancer_types = [os.path.basename(f)[:-4] for f in glob.glob('{0}/*.txt'.format(base_dir))]
    result_list = []
    for c in cancer_types:
        tmp_path = os.path.join(base_dir, c+'.txt')
        tmp = pd.read_table(tmp_path)
        counts = tmp.groupby('category')['number of mutations'].sum()
        counts.name = c
        result_list.append(counts)
    concat_df = pd.concat(result_list, axis=1)
    concat_df = concat_df.fillna(0)
    return concat_df

############################
# Read significant mutations from MSK-IMPACT
############################
def read_msk_impact(chasm2_path, maf_path):
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
    signif_df = mut_df[mut_df.Protein_change.isin(df[is_signif]['Protein_change'])].drop_duplicates(['Hugo_Symbol', 'HGVSp_Short'])
    #signif_df = mut_df.drop_duplicates(['Hugo_Symbol', 'HGVSp_Short'])

    # fix x/y suffixes
    rename_dict = {'Protein_position_x': 'Protein_position'}
    signif_df = signif_df.rename(columns=rename_dict)
    
    return signif_df

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
    zorder = 10
    for method in methods:
        if method == 'CanDrA.plus': method = 'CanDrA plus'
        pred = data[method].astype(float).dropna()
        fpr, tpr, thresholds = roc_curve(y[pred.index], pred)
        myauc = auc(fpr, tpr)
        plt.plot(fpr, tpr,
                 label='{0} (area = {1:0.3f})'.format(method, myauc),
                 zorder=zorder)
        plt.legend(loc='best')
        plt.xlabel('False positive rate')
        plt.ylabel('True positive rate')
        zorder -= 1

    # format axis
    plt.gca().xaxis.set_ticks_position('bottom')
    plt.gca().yaxis.set_ticks_position('left')

def top5(data):
    """Figure out the top 5 best beforming methods."""
    replace_dict = {'gwCHASMplus': 'raw CHASMplus', 'CanDrA.plus': 'CanDrA',
                    'Polyphen2_hdiv': 'Polyphen2', 'Polyphen2_hvar': 'Polyphen2',
                     #'CHASM2': 'raw CHASMplus',
                    #'gwCHASM2': 'gwCHASMplus'
                    }
    replace_dict2 = {#'CHASM2_genome': 'raw CHASMplus', 
                     'CanDrA.plus': 'CanDrA plus',
                    '1-SIFT': 'SIFT', '1-CHASM': 'CHASM', 
                    #'CHASM2': 'raw CHASMplus',
                    #'gwCHASM2': 'gwCHASMplus'
                    }
    data['method'] = data.method.replace(replace_dict)
    tmp_top_methods = data.groupby('method')['auc'].max().sort_values(ascending=False).head(5).index.to_series().replace(replace_dict2).tolist()
    top_methods = data[data.method.isin(tmp_top_methods)].sort_values('auc', ascending=False).index.to_series().replace(replace_dict2).tolist()
    return top_methods

def fetch_methods(path):
    comp_df = pd.read_table(path)
    top_methods = top5(comp_df)
    other_methods = list(set(methods) - set(top_methods))
    return top_methods, other_methods
    

def pr_curve(data, methods, other_methods):
    """Create a precision recall curve of methods."""
    y = data['y']
    # plot the low performing methods
    for method in other_methods:
        if method == 'CanDrA.plus': method = 'CanDrA plus'
        pred = data[method].astype(float).dropna()
        prec, recall, thresholds = precision_recall_curve(y[pred.index], pred)
        myauc = average_precision_score(y[pred.index], pred)
        #myauc = auc(recall[:-1], prec[:-1])
        plt.plot(recall[:-1], prec[:-1],
                 color='lightgray')
        plt.xlabel('Recall')
        plt.ylabel('Precision')

    # plot the top methods
    zorder = 10
    for method in methods:
        if method == 'CanDrA.plus': method = 'CanDrA plus'
        pred = data[method].astype(float).dropna()
        prec, recall, thresholds = precision_recall_curve(y[pred.index], pred)
        myauc = average_precision_score(y[pred.index], pred)
        #myauc = auc(recall[:-1], prec[:-1])
        plt.plot(recall[:-1], prec[:-1], 
                 label='{0} (area = {1:0.3f})'.format(method, myauc),
                 zorder=zorder)
        plt.legend(loc='best')
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        zorder -= 1

    # format axis
    plt.gca().xaxis.set_ticks_position('bottom')
    plt.gca().yaxis.set_ticks_position('left')

    
def box_plot_with_significance(x, y, **kwargs):
    """Creates a Facet grid of boxplots with annotations on whether differes 
    are significant.
    """
    # take in data
    data = kwargs['data']
    signif = kwargs['signif']
    facet_var = data.loc[x.index, 'var'].iloc[0]
    
    # establish significance
    miss_pvalue = signif.loc[signif.variable==facet_var, 'wt vs missense p-value'].iloc[0]
    lof_pvalue = signif.loc[signif.variable==facet_var, 'wt vs lof p-value'].iloc[0]
    miss_text = 'ns'
    if miss_pvalue <= 0.05: miss_text = '*'
    if miss_pvalue <= 0.01: miss_text += '*'
    if miss_pvalue <= 0.001: miss_text += '*'
    lof_text = 'ns'
    if lof_pvalue <= 0.05: lof_text = '*'
    if lof_pvalue <= 0.01: lof_text += '*'
    if lof_pvalue <= 0.001: lof_text += '*'
    
    # set up max y value
    max_val = y.max() + (y.max() - y.min()) * .05
    
    # do missense vs control
    h = (y.max() - y.min()) * .05
    x1, x2 = 0, 1
    plt.plot([x1, x1, x2, x2], [max_val, max_val+h, max_val+h, max_val], lw=1.5, color='black')
    if '*' in miss_text:
        plt.text((x1+x2)*.5, max_val, miss_text, ha='center', va='bottom', color='black', fontsize=16)
    else:
        plt.text((x1+x2)*.5, max_val+1.2*h, miss_text, ha='center', va='bottom', color='black')
    
    # do lof vs control
    max_val = max_val + 3*h
    x1, x2 = 0, 2
    plt.plot([x1, x1, x2, x2], [max_val, max_val+h, max_val+h, max_val], lw=1.5, color='black')
    if '*' in lof_text:
        plt.text((x1+x2)*.5, max_val, lof_text, ha='center', va='bottom', color='black', fontsize=16)
    else:
        plt.text((x1+x2)*.5, max_val+1.2*h, lof_text, ha='center', va='bottom', color='black')

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

######################
# Venn diagram plot
######################

def venn_diagram(set1, set2, name1, name2, title=''):
    len_intersect = len(set(set1) & set(set2))
    len_set1 = len(set1)
    len_set2 = len(set2)
    overlap = (len_set1 - len_intersect, len_set2 - len_intersect, len_intersect)

    # plot venn diagram
    with sns.plotting_context('paper', font_scale=1.4):
        venn2(subsets=overlap, set_labels=(name1, name2))
        venn2_circles(subsets=overlap, linestyle='solid', linewidth=.75)
        plt.title(title, size=16)
        
    return overlap

def venn_diagram3(set1, set2, set3, name1, name2, name3, title='', ax=None):
    # make sure to convert to set object
    set1 = set(set1)
    set2 = set(set2)
    set3 = set(set3)
    
    # do all possible intersections
    intersect_12 = set(set1) & set(set2)
    intersect_13 = set(set1) & set(set3)
    intersect_23 = set(set2) & set(set3)
    full_intersect = set(set1) & set(set2) & set(set3)
    
    # figure out number only specific to one set
    len_set1_specific = len(set1 - set2 - set3)
    len_set2_specific = len(set2 - set1 - set3)
    len_set3_specific = len(set3 - set1 - set2)
    
    # Figure out length of full intersect
    len_full_intersect = len(full_intersect)
    
    # figure out length of two-set specific overlaps
    len_set12 = len(intersect_12 - full_intersect)
    len_set13 = len(intersect_13 - full_intersect)
    len_set23 = len(intersect_23 - full_intersect)

    # create overlap object
    #overlap = (len_set1 - len_intersect, len_set2 - len_intersect, len_intersect)
    overlap = {'100': len_set1_specific, '010': len_set2_specific, '001': len_set3_specific,
               '110': len_set12, '101': len_set13, '011': len_set23,
               '111': len_full_intersect}

    # plot venn diagram
    with sns.plotting_context('notebook', font_scale=1.0):
        if ax:
            venn3(subsets=overlap, set_labels=(name1, name2, name3), ax=ax)
            venn3_circles(subsets=overlap, linestyle='solid', linewidth=.75, ax=ax)
        else:
            venn3(subsets=overlap, set_labels=(name1, name2, name3))
            venn3_circles(subsets=overlap, linestyle='solid', linewidth=.75)            
        plt.title(title, size=16)
        
    return overlap

########################
# Read in OncoKB
########################

def read_oncokb(path='CHASMplus/data/misc/oncokb_4_3_2017.txt'):
    """Read in the OncoKB mutations"""
    oncokb = pd.read_table(path)
    oncokb['HGVSp_Short'] = 'p.' + oncokb['Alteration']
    oncokb = oncokb.rename(columns={'Gene': 'Hugo_Symbol'})
    oncokb['OncoKB'] = oncokb['Oncogenicity'].isin(['Oncogenic', 'Likely Oncogenic']).astype(int)
    return oncokb


########################
# stat funcs
########################
def cummin(x):
    """A python implementation of the cummin function in R"""
    for i in range(1, len(x)):
        if x[i-1] < x[i]:
            x[i] = x[i-1]
    return x


def bh_fdr(pval):
    """A python implementation of the Benjamani-Hochberg FDR method.
    This code should always give precisely the same answer as using
    p.adjust(pval, method="BH") in R.
    Parameters
    ----------
    pval : list or array
        list/array of p-values
    Returns
    -------
    pval_adj : np.array
        adjusted p-values according the benjamani-hochberg method
    """
    pval_array = np.array(pval)
    sorted_order = np.argsort(pval_array)
    original_order = np.argsort(sorted_order)
    pval_array = pval_array[sorted_order]
    
    # calculate the needed alpha
    n = float(len(pval))
    pval_adj = np.zeros(int(n))
    i = np.arange(1, int(n)+1, dtype=float)[::-1]  # largest to smallest
    pval_adj = np.minimum(1, cummin(n/i * pval_array[::-1]))[::-1]
    return pval_adj[original_order]

def compute_p_value(scores, null_p_values):
    """Get the p-value for each score by examining the list null distribution
    where scores are obtained by a certain probability.
    NOTE: uses score2pval function
    Parameters
    ----------
    scores : pd.Series
        series of observed scores
    null_p_values: pd.Series
        Empirical null distribution, index are scores and values are p values
    Returns
    -------
    pvals : pd.Series
        Series of p values for scores
    """
    num_scores = len(scores)
    pvals = pd.Series(np.zeros(num_scores))
    null_p_val_scores = list(reversed(null_p_values.index.tolist()))
    #null_p_values = null_p_values.ix[null_p_val_scores].copy()
    null_p_values.sort_values(inplace=True, ascending=False)
    pvals = scores.apply(lambda x: score2pval(x, null_p_val_scores, null_p_values))
    return pvals


def score2pval(score, null_scores, null_pvals):
    """Looks up the P value from the empirical null distribution based on the provided
    score.
    NOTE: null_scores and null_pvals should be sorted in ascending order.
    Parameters
    ----------
    score : float
        score to look up P value for
    null_scores : list
        list of scores that have a non-NA value
    null_pvals : pd.Series
        a series object with the P value for the scores found in null_scores
    Returns
    -------
    pval : float
        P value for requested score
    """
    # find position in simulated null distribution
    pos = bisect.bisect_right(null_scores, score)

    # if the score is beyond any simulated values, then report
    # a p-value of zero
    if pos == null_pvals.size and score > null_scores[-1]:
        return 0
    # condition needed to prevent an error
    # simply get last value, if it equals the last value
    elif pos == null_pvals.size:
        return null_pvals.iloc[pos-1]
    # normal case, just report the corresponding p-val from simulations
    else:
        return null_pvals.iloc[pos]