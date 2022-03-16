import argparse
import sys
import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import glob
from sklearn import metrics
from scipy.stats import pearsonr, spearmanr
from scipy.optimize import curve_fit
from collections import Counter
import pickle
import pdb

parser = argparse.ArgumentParser(description = '''Visualize and analyze the DockQ scores.''')
#Bench4
parser.add_argument('--bench4_dockq_aa', nargs=1, type= str, default=sys.stdin, help = 'Path to dockq scores for benchmark 4 AF in csv.')
parser.add_argument('--bench4_dockq_RF', nargs=1, type= str, default=sys.stdin, help = 'Path to dockq scores for benchmark 4 from RF in csv.')
parser.add_argument('--plDDT_bench4', nargs=1, type= str, default=sys.stdin, help = 'Path to plDDT metrics in csv.')
#parser.add_argument('--pconsdock_bench4', nargs=1, type= str, default=sys.stdin, help = 'Path to pconsdock metrics in csv.')
#parser.add_argument('--pconsdock_marks', nargs=1, type= str, default=sys.stdin, help = 'Path to pconsdock metrics in csv.')
parser.add_argument('--bench4_kingdom', nargs=1, type= str, default=sys.stdin, help = 'Path to kingdoms for bench4 in csv.')
parser.add_argument('--dssp_bench4', nargs=1, type= str, default=sys.stdin, help = 'Path to dssp annotations for bench4 in csv.')
parser.add_argument('--afdefault_neff_bench4', nargs=1, type= str, default=sys.stdin, help = 'Path to default AF alignments Neff in csv.')
parser.add_argument('--tophits_neff_bench4', nargs=1, type= str, default=sys.stdin, help = 'Path to top hits alignments Neff in csv.')
#Marks positivef
parser.add_argument('--marks_dockq_RF', nargs=1, type= str, default=sys.stdin, help = 'Path to dockq scores for marks set RF in csv.')
parser.add_argument('--marks_dockq_AF_bb', nargs=1, type= str, default=sys.stdin, help = 'Path to dockq scores for marks set AF back bone atoms in csv.')
parser.add_argument('--marks_dockq_AF_aa', nargs=1, type= str, default=sys.stdin, help = 'Path to dockq scores for marks set AF all atoms in csv.')
parser.add_argument('--marks_dockq_GRAMM', nargs=1, type= str, default=sys.stdin, help = 'Path to dockq scores for marks set GRAMM in csv.')
parser.add_argument('--marks_dockq_TMfull', nargs=1, type= str, default=sys.stdin, help = 'Path to dockq scores for marks set TMdock in csv.')
parser.add_argument('--marks_dockq_TMint', nargs=1, type= str, default=sys.stdin, help = 'Path to dockq scores for marks set interface TMdock in csv.')
parser.add_argument('--marks_dockq_mdockpp', nargs=1, type= str, default=sys.stdin, help = 'Path to dockq scores for marks set MdockPP in csv.')

parser.add_argument('--plDDT_marks_af', nargs=1, type= str, default=sys.stdin, help = 'Path to plDDT metrics in csv.')
parser.add_argument('--plDDT_marks_fused', nargs=1, type= str, default=sys.stdin, help = 'Path to plDDT metrics in csv.')
parser.add_argument('--dssp_marks', nargs=1, type= str, default=sys.stdin, help = 'Path to dssp metrics in csv.')
parser.add_argument('--ifstats_marks', nargs=1, type= str, default=sys.stdin, help = 'Path to if metrics in csv.')
parser.add_argument('--aln_scores_marks', nargs=1, type= str, default=sys.stdin, help = 'Path to aln scores in csv.')
parser.add_argument('--oxstats_marks', nargs=1, type= str, default=sys.stdin, help = 'Path to statistics over organisms in csv.')
parser.add_argument('--afdefault_neff_marks', nargs=1, type= str, default=sys.stdin, help = 'Path to default AF alignments Neff in csv.')
parser.add_argument('--tophits_neff_marks', nargs=1, type= str, default=sys.stdin, help = 'Path to top hits alignments Neff in csv.')
parser.add_argument('--af_chain_overlap_marks', nargs=1, type= str, default=sys.stdin, help = 'Path to chain overlap for AF a3m in csv.')
#Marks negative
parser.add_argument('--plDDT_marks_negative_af', nargs=1, type= str, default=sys.stdin, help = 'Path to plDDT metrics in csv.')
#Negatome
parser.add_argument('--plDDT_negatome_af', nargs=1, type= str, default=sys.stdin, help = 'Path to plDDT metrics in csv.')
#New set
parser.add_argument('--newset_dockq_AF', nargs=1, type= str, default=sys.stdin, help = 'Path to dockq scores for new set AF in csv.')
parser.add_argument('--plDDT_newset', nargs=1, type= str, default=sys.stdin, help = 'Path to plDDT metrics in csv for newset.')
#Output directory
parser.add_argument('--outdir', nargs=1, type= str, default=sys.stdin, help = 'Path to output directory. Include /in end')


################FUNCTIONS#################
def dockq_box(bench4_dockq, outdir):
    '''Plot a boxplot of the dockq score for the different modes
    '''

    #Plot
    fig,ax = plt.subplots(figsize=(24/2.54,12/2.54))
    modes = bench4_dockq.columns[1:]
    all_modes = []
    all_scores = []
    all_msas = []
    all_model_options = []
    accuracies = {}
    for mode in modes:
        #Frac correct and avg score
        fraq_correct = np.argwhere(bench4_dockq[mode].values>=0.23).shape[0]/len(bench4_dockq)
        accuracies[mode]=fraq_correct
        av = np.average(bench4_dockq[mode].values)
        print(mode, np.round(fraq_correct,3),np.round(av,3))
        #Save scores
        all_scores.extend([*bench4_dockq[mode].values])
        mode = '_'.join(mode.split('_')[4:])
        mode = mode.split('_')
        msa = mode[0]
        model = '_'.join(mode[1:-1])
        option = mode[-1]
        #save
        all_modes.extend([msa+'\n'+model+'\n'+option]*len(bench4_dockq))
        all_msas.extend([msa]*len(bench4_dockq))
        all_model_options.extend([model+' '+option]*len(bench4_dockq))






def correlate_scores(bench4_dockq, outdir):
    '''Correlate the scores for all different modeling strategies
    '''

    modes = ['DockQ_dockqstats_bench4_af2_hhblits_model_1_ens8',
     'DockQ_dockqstats_bench4_af2_hhblits_model_1_rec10',
     'DockQ_dockqstats_bench4_af2_af2stdmsa_model_1_ens8',
     'DockQ_dockqstats_bench4_af2_af2stdmsa_model_1_rec10',
     'DockQ_dockqstats_bench4_af2_af2andhhblitsmsa_model_1_ens8',
     'DockQ_dockqstats_bench4_af2_af2andhhblitsmsa_model_1_rec10',
     'DockQ_dockqstats_bench4_af2_hhblits_model_1_ptm_ens8',
      'DockQ_dockqstats_bench4_af2_hhblits_model_1_ptm_rec10',
      'DockQ_dockqstats_bench4_af2_af2stdmsa_model_1_ptm_ens8',
      'DockQ_dockqstats_bench4_af2_af2stdmsa_model_1_ptm_rec10',
      'DockQ_dockqstats_bench4_af2_af2andhhblitsmsa_model_1_ptm_ens8',
      'DockQ_dockqstats_bench4_af2_af2andhhblitsmsa_model_1_ptm_rec10']


    corr_matrix = np.zeros((len(modes),len(modes)))
    for i in range(len(modes)):
        scores_i = bench4_dockq[modes[i]].values
        for j in range(i+1,len(modes)):
            scores_j = bench4_dockq[modes[j]].values
            #Correlate
            R,p = pearsonr(scores_i,scores_j)
            corr_matrix[i,j]=np.round(R,2)
            corr_matrix[j,i]=np.round(R,2)

    print(modes)
    print(corr_matrix)
    #Create df
    corr_df = pd.DataFrame()
    modes = ['_'.join(x.split('_')[4:]) for x in modes]
    corr_df['Comparison'] = modes
    for i in range(len(modes)):
        corr_df[modes[i]]=corr_matrix[i,:]

    corr_df.to_csv(outdir+'model_correlations.csv')


def fetch_missing_dockq(marks_dockq_AF_bb,marks_dockq_AF_aa):
    '''Fetch missing DockQ scores
    '''

    ids = ['_'.join(x.split('-')) for x in marks_dockq_AF_aa.complex_id.values]
    #Get mising scores
    missing = marks_dockq_AF_bb[~marks_dockq_AF_bb.complex_id.isin(ids)]
    ids = [x[:6]+'-'+x[7:] for x in missing.complex_id.values]
    missing['complex_id']=ids
    marks_dockq_AF_aa = pd.concat([marks_dockq_AF_aa,missing[marks_dockq_AF_aa.columns]])
    return marks_dockq_AF_aa

def pdockq(if_plddt_contacts, dockq_scores, outdir):

    #pdockq
    fig,ax = plt.subplots(figsize=(12/2.54,12/2.54))
    #Create RA
    x_ra = []
    y_ra = []
    y_std = []
    y_av_err = []
    step = 20
    for t in np.arange(0,max(if_plddt_contacts)-step,step):
        inds = np.argwhere((if_plddt_contacts>=t)&(if_plddt_contacts<t+step))[:,0]
        x_ra.append(t+step/2)
        y_ra.append(np.average(dockq_scores[inds]))
        y_std.append(np.std(dockq_scores[inds]))
        y_av_err.append(np.average(np.absolute(dockq_scores[inds]-y_ra[-1])))

    #Do a simple sigmoid fit
    def sigmoid(x, L ,x0, k, b):
        y = L / (1 + np.exp(-k*(x-x0)))+b
        return (y)

    xdata = if_plddt_contacts[np.argsort(if_plddt_contacts)]
    ydata = dockq_scores[np.argsort(if_plddt_contacts)]
    p0 = [max(ydata), np.median(xdata),1,min(ydata)] # this is an mandatory initial guess

    popt, pcov = curve_fit(sigmoid, xdata, ydata,p0, method='dogbox')
    y = sigmoid(xdata, *popt)
    plt.plot(xdata,y,color='r',label='Sigmoidal fit')
    #Calc error
    print('Sigmoid params:',*popt)
    plt.scatter(if_plddt_contacts,dockq_scores,s=1)
    #plt.plot(x_ra,y_ra,label='Running average', color='tab:blue')
    #plt.fill_between(x_ra,np.array(y_ra)-np.array(y_av_err),np.array(y_ra)+np.array(y_av_err),color='tab:blue',alpha=0.25, label='Average error')
    plt.title('pDockQ')
    plt.xlabel('IF plDDT⋅log(IF contacts)')
    plt.ylabel('DockQ')
    plt.legend()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(outdir+'pDockQ.svg',format='svg',dpi=300)
    plt.close()
    print('Average error for sigmoidal fit:',np.average(np.absolute(y-ydata)))
    print('L=',np.round(popt[0],3),'x0=',np.round(popt[1],3) ,'k=',np.round(popt[2],3), 'b=',np.round(popt[3],3))

    return popt


def ROC_pred_marks(marks_dockq_AF, plDDT_marks, outdir):
    '''Compare the separation in the marks dataset for AF using metrics from the
    predicted structures
    '''
    #Merge dfs
    plDDT_marks['complex_id']=plDDT_marks.id1+'-'+plDDT_marks.id2
    merged = pd.merge(marks_dockq_AF,plDDT_marks,on=['complex_id'],how='inner')


    #Get min of chains
    single_chain_plddt = np.min(merged[['ch1_plddt_av_1', 'ch2_plddt_av_1']].values,axis=1)
    merged['min_chain_plddt_av_1'] = single_chain_plddt
    #Analyze ROC as a function of
    plDDT_metrics = ['if_plddt_av',  'min_chain_plddt_av',
                     'plddt_av', 'num_atoms_in_interface', 'num_res_in_interface']

    plDDT_nice_names = {'if_plddt_av':'IF_plDDT',  'min_chain_plddt_av':'Min plDDT per chain',
                    'plddt_av':'Average plDDT', 'num_atoms_in_interface':'IF_contacts',
                        'num_res_in_interface':'IF_residues'}
    run='1'
    dockq_scores = merged['DockQ_dockqstats_marks_af2_af2andhhblitsmsa_model_1_rec10_run'+run].values
    correct = np.zeros(len(dockq_scores))
    correct[np.argwhere(dockq_scores>=0.23)]=1
    #Plot
    fig,ax = plt.subplots(figsize=(12/2.54,12/2.54))
    colors = {0:'darkblue',1:'magenta',2:'orange',3:'darkgreen',4:'tab:blue',5:'tab:yellow',6:'tab:black'}

    for i in range(len(plDDT_metrics)):
        plDDT_metric_vals = merged[plDDT_metrics[i]+'_'+run].values
        #Create ROC
        fpr, tpr, threshold = metrics.roc_curve(correct, plDDT_metric_vals, pos_label=1)
        roc_auc = metrics.auc(fpr, tpr)
        label = plDDT_metrics[i]
        plt.plot(fpr, tpr, label = plDDT_nice_names[label]+': AUC = %0.2f' % roc_auc,color=colors[i])

    #Add log(if contacts)*if_plddt_av
    if_plddt_contacts = merged['if_plddt_av_1'].values*np.log10(merged['num_atoms_in_interface_1'].values+1)
    #Create ROC
    fpr, tpr, threshold = metrics.roc_curve(correct, if_plddt_contacts, pos_label=1)
    roc_auc = metrics.auc(fpr, tpr)
    plt.plot(fpr, tpr, label = 'IF_plDDT⋅log(IF_contacts)'+': AUC = %0.2f' % roc_auc,color='tab:cyan')

    #Get pDockQ
    def sigmoid(x, L ,x0, k, b):
        y = L / (1 + np.exp(-k*(x-x0)))+b
        return (y)
    sigmoid_params = pdockq(if_plddt_contacts, dockq_scores, outdir)
    #Create ROC
    fpr, tpr, threshold = metrics.roc_curve(correct, sigmoid(if_plddt_contacts,*sigmoid_params), pos_label=1)
    roc_auc = metrics.auc(fpr, tpr)
    plt.plot(fpr, tpr, label = 'pDockQ'+': AUC = %0.2f' % roc_auc,color='k',linestyle='--')
    plt.plot([0,1],[0,1],linewidth=1,linestyle='--',color='grey')
    plt.legend(fontsize=9)
    plt.title('ROC as a function of different metrics')
    plt.xlabel('FPR')
    plt.ylabel('TPR')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(outdir+'ROC_marks.svg',format='svg',dpi=300)
    plt.close()

    #pDockQ vs DockQ
    fig,ax = plt.subplots(figsize=(12/2.54,12/2.54))
    plt.scatter(sigmoid(if_plddt_contacts,*sigmoid_params),dockq_scores,s=1)
    plt.title('pDockQ vs DockQ')
    plt.xlabel('pDockQ')
    plt.ylabel('DockQ')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(outdir+'pdockq_vs_dockq.svg',format='svg',dpi=300)
    plt.close()


    #plot if plddt vs log contacts and color by dockq
    fig,ax = plt.subplots(figsize=(12/2.54,12/2.54))
    plt.scatter(merged['num_atoms_in_interface_1'].values+1, merged['if_plddt_av_1'].values,c=dockq_scores,s=2)
    cbar = plt.colorbar()
    cbar.set_label('DockQ')
    plt.xscale('log')
    plt.ylim([40,100])
    plt.title('Interface contacts, plDDT and DockQ')
    plt.xlabel('Interface contacts')
    plt.ylabel('Average interface plDDT')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(outdir+'if_conctacts_vs_plddt.svg',format='svg',dpi=300)
    plt.close()




    return sigmoid_params



def score_marks_5runs_paired_af(marks_dockq_AF, plDDT_marks, sigmoid_params, outdir):
    '''Analyze the variation in DockQ using 5 identical runs of the same settings
    '''

    plDDT_marks['complex_id'] = plDDT_marks.id1+'-'+plDDT_marks.id2
    merged = pd.merge(marks_dockq_AF,plDDT_marks,on='complex_id',how='inner')

    #Get separator
    separator1 = merged[['if_plddt_av_1', 'if_plddt_av_2','if_plddt_av_3','if_plddt_av_4','if_plddt_av_5']].values
    separator2 = merged[['num_atoms_in_interface_1', 'num_atoms_in_interface_2','num_atoms_in_interface_3','num_atoms_in_interface_4','num_atoms_in_interface_5']].values
    separator = separator1*np.log10(separator2+1)
    def sigmoid(x, L ,x0, k, b):
        y = L / (1 + np.exp(-k*(x-x0)))+b
        return (y)
    separator = sigmoid(separator, *sigmoid_params)

    scores = merged[['DockQ_dockqstats_marks_af2_af2andhhblitsmsa_model_1_rec10_run1','DockQ_dockqstats_marks_af2_af2andhhblitsmsa_model_1_rec10_run2',
    'DockQ_dockqstats_marks_af2_af2andhhblitsmsa_model_1_rec10_run3','DockQ_dockqstats_marks_af2_af2andhhblitsmsa_model_1_rec10_run4',
    'DockQ_dockqstats_marks_af2_af2andhhblitsmsa_model_1_rec10_run5']].values
    #Get max and min scores
    max_scores = np.max(scores,axis=1)
    min_scores = np.min(scores,axis=1)
    #Get success rates per initializations
    srs = []
    for i in range(scores.shape[1]):
        srs.append(np.argwhere(scores[:,i]>=0.23).shape[0]/len(scores))
    print('Test set scoring:')
    print('Success rate 5 runs top scores',np.argwhere(max_scores>=0.23).shape[0]/len(max_scores))
    print('Av diff',np.average(max_scores-min_scores))
    print('Std diff',np.std(max_scores-min_scores))
    print('Avg and std success rate', np.average(srs),np.std(srs))
    #Separate the models using the number of contacts in the interface
    max_inds = np.argmax(separator,axis=1)
    first_ranked_scores = []
    first_ranked_separators = []
    #Get max separator scores
    for i in range(len(max_inds)):
        first_ranked_scores.append(scores[i,max_inds[i]])
        first_ranked_separators.append(separator[i,max_inds[i]])
    #Convert to array
    first_ranked_scores = np.array(first_ranked_scores)
    first_ranked_separators = np.array(first_ranked_separators)
    #Get success rate
    print('Ranking test set success rate using av plDDT*log(if_contacts) in interface',np.argwhere(first_ranked_scores>=0.23).shape[0]/len(first_ranked_scores))
    #Get AUC using that success rate
    correct = np.zeros(len(first_ranked_scores))
    correct[np.argwhere(first_ranked_scores>=0.23)]=1
    fpr, tpr, threshold = metrics.roc_curve(correct, first_ranked_separators, pos_label=1)
    roc_auc = metrics.auc(fpr, tpr)
    print('AUC using the same ranking', roc_auc)

    #Plot
    fig,ax = plt.subplots(figsize=(12/2.54,12/2.54))
    plt.scatter(first_ranked_scores, max_scores,s=3,color='tab:blue',label='Max')
    plt.scatter(first_ranked_scores, min_scores,s=3,color='mediumseagreen',label='Min')
    plt.title('Model ranking on the test set')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.plot([0,1],[0,1],color='k',linewidth=1,linestyle='--')
    plt.xlabel('DockQ first ranked model')
    plt.ylabel('DockQ')
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir+'DockQ_marks_5runs.svg',format='svg',dpi=300)
    plt.close()

    #Assign top ranked scores and origin
    marks_dockq_AF['top_ranked_model_DockQ_af2']=first_ranked_scores
    marks_dockq_AF['top_ranked_pDockQ']=first_ranked_separators
    marks_dockq_AF['top_ranked_model_run_af2']=max_inds+1

    #Create and save df
    roc_df = pd.DataFrame()
    roc_df['FPR']=fpr
    roc_df['TPR']=tpr
    roc_df['FPR*SR']=fpr*(np.argwhere(first_ranked_scores>=0.23).shape[0]/len(first_ranked_scores))
    roc_df['TPR*SR']=tpr*np.argwhere(first_ranked_scores>=0.23).shape[0]/len(first_ranked_scores)
    roc_df['PPV']=tpr/(tpr+fpr)
    roc_df['pDockQ']=threshold
    roc_df.to_csv(outdir+'roc_df_af2_marks.csv')
    #Select a reduced version of the df
    roc_df.loc[np.arange(0,len(roc_df),10)].to_csv(outdir+'roc_df_af2_marks_reduced.csv')

    return marks_dockq_AF

def score_marks_5runs_paired_fused(marks_dockq_AF, plDDT_marks, sigmoid_params, outdir):
    '''Analyze the variation in DockQ using 5 identical runs of the same settings
    '''

    plDDT_marks['complex_id'] = plDDT_marks.id1+'-'+plDDT_marks.id2
    merged = pd.merge(marks_dockq_AF,plDDT_marks,on='complex_id',how='inner')

    #Get separator
    separator1 = merged[['if_plddt_av_1', 'if_plddt_av_2','if_plddt_av_3','if_plddt_av_4','if_plddt_av_5']].values
    separator2 = merged[['num_atoms_in_interface_1', 'num_atoms_in_interface_2','num_atoms_in_interface_3','num_atoms_in_interface_4','num_atoms_in_interface_5']].values
    separator = separator1*np.log10(separator2+1)

    def sigmoid(x, L ,x0, k, b):
        y = L / (1 + np.exp(-k*(x-x0)))+b
        return (y)
    separator = sigmoid(separator, *sigmoid_params)


    scores = merged[['DockQ_dockqstats_marks_af2_pairedandfused_model_1_rec10_run1','DockQ_dockqstats_marks_af2_pairedandfused_model_1_rec10_run2',
      'DockQ_dockqstats_marks_af2_pairedandfused_model_1_rec10_run3','DockQ_dockqstats_marks_af2_pairedandfused_model_1_rec10_run4',
      'DockQ_dockqstats_marks_af2_pairedandfused_model_1_rec10_run5']].values
    #Get max and min scores
    max_scores = np.max(scores,axis=1)
    min_scores = np.min(scores,axis=1)
    #Get success rates per initializations
    srs = []
    for i in range(scores.shape[1]):
        srs.append(np.argwhere(scores[:,i]>=0.23).shape[0]/len(scores))
    print('FUSED test set scoring:')
    print('Success rate 5 runs top scores',np.argwhere(max_scores>=0.23).shape[0]/len(max_scores))
    print('Av diff',np.average(max_scores-min_scores))
    print('Std diff',np.std(max_scores-min_scores))
    print('Avg and std success rate', np.average(srs),np.std(srs))

    #Separate the models using the number of contacts in the interface
    max_inds = np.argmax(separator,axis=1)
    first_ranked_scores = []
    first_ranked_separators = []
    #Get max separator scores
    for i in range(len(max_inds)):
        first_ranked_scores.append(scores[i,max_inds[i]])
        first_ranked_separators.append(separator[i,max_inds[i]])
    #Convert to array
    first_ranked_scores = np.array(first_ranked_scores)
    first_ranked_separators = np.array(first_ranked_separators)
    #Get success rate
    print('Ranking test set success rate using if_plddt_av and num contacts in interface',np.argwhere(first_ranked_scores>=0.23).shape[0]/len(first_ranked_scores))
    #Get AUC using that success rate
    correct = np.zeros(len(first_ranked_scores))
    correct[np.argwhere(first_ranked_scores>=0.23)]=1
    fpr, tpr, threshold = metrics.roc_curve(correct, first_ranked_separators, pos_label=1)
    roc_auc = metrics.auc(fpr, tpr)
    print('FUSED AUC using the same ranking', roc_auc)

    #Assign top ranked scores and origin
    marks_dockq_AF['top_ranked_model_DockQ_fused']=first_ranked_scores
    marks_dockq_AF['top_ranked_model_run_fused']=max_inds+1

    #Create and save df
    roc_df = pd.DataFrame()
    roc_df['FPR']=fpr
    roc_df['TPR']=tpr
    roc_df['FPR*SR']=fpr*(np.argwhere(first_ranked_scores>=0.23).shape[0]/len(first_ranked_scores))
    roc_df['TPR*SR']=tpr*np.argwhere(first_ranked_scores>=0.23).shape[0]/len(first_ranked_scores)
    roc_df['PPV']=tpr/(tpr+fpr)
    roc_df['pDockQ']=threshold
    roc_df.to_csv(outdir+'roc_df_fused_marks.csv')
    #Select a reduced version of the df
    roc_df.loc[np.arange(0,len(roc_df),10)].to_csv(outdir+'roc_df_fused_marks_reduced.csv')

    return marks_dockq_AF

def marks_box(marks_dockq_AF, marks_dockq_GRAMM, marks_dockq_mdockpp, marks_dockq_TMfull, marks_dockq_TMint, marks_dockq_RF,outdir):
    '''Box df of Marks set
    '''
    marks_dockq_TMint = marks_dockq_TMint.dropna()
    marks_dockq_TMfull = marks_dockq_TMfull.dropna()

    #Get data
    rf_scores = marks_dockq_RF.DockQ_dockqstats_marks_RF.values
    gramm_scores = marks_dockq_GRAMM[1].values
    mdockpp_scores = marks_dockq_mdockpp.DockQ.values
    TMfull_scores = marks_dockq_TMfull.dockq.values
    TMint_scores = marks_dockq_TMint.dockq.values
    paired_scores = marks_dockq_AF.DockQ_dockqstats_marks_af2_hhblitsn2_model_1_rec10.values
    af2_std_scores = marks_dockq_AF.DockQ_dockqstats_marks_af2_af2stdmsa_model_1_rec10.values
    run1_both_scores= marks_dockq_AF.DockQ_dockqstats_marks_af2_af2andhhblitsmsa_model_1_rec10_run1.values
    run1_fused_scores = marks_dockq_AF.DockQ_dockqstats_marks_af2_pairedandfused_model_1_rec10_run1.values
    top_paired_af_scores = marks_dockq_AF.top_ranked_model_DockQ_af2.values
    top_paired_fused_scores = marks_dockq_AF.top_ranked_model_DockQ_fused.values


    data1 = [rf_scores, gramm_scores, mdockpp_scores, TMint_scores, af2_std_scores,  paired_scores, top_paired_af_scores, top_paired_fused_scores]
    data2 = [run1_both_scores, run1_fused_scores, top_paired_af_scores,top_paired_fused_scores]
    all_data = [data1,data2]
    xlabels1 = ['RF','GRAMM', 'MDockPP', 'TMdock\nInterfaces', 'AF2', 'Paired', 'AF2+Paired\ntop ranked','Block+Paired\ntop ranked']
    xlabels2 = ['AF2+Paired', 'Block+Paired', 'AF2+Paired\ntop ranked', 'Block+Paired\ntop ranked']
    all_xlabels = [xlabels1, xlabels2]

    #Color
    colors = sns.husl_palette(len(xlabels1)+2)
    all_colors = [colors[:len(xlabels1)],colors[-len(xlabels2):]]

    for i in range(len(all_data)):
        #Boxplot
        fig,ax = plt.subplots(figsize=(24/2.54,12/2.54))
        data = all_data[i] #Get data and xlabel variation
        xlabels = all_xlabels[i]
        colors = all_colors[i]
        #Success rates
        srs = []
        for j in range(len(data)):
            sr = np.argwhere(data[j]>=0.23).shape[0]/len(data[j])
            median = np.median(data[j])
            print(xlabels[j],'sr:',np.round(sr,3),len(data[j]),median)
            #xlabels[j]+='\nSR: '+str(np.round(100*sr,1))+'%'
            #xlabels[j]+='\nM: '+str(np.round(median,3))

        # Creating plot
        #ax.violinplot(data)
        bp = ax.boxplot(data, patch_artist = True, notch=True, showfliers=False)

        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.75)

        # changing color and linewidth of
        # medians
        for median in bp['medians']:
            median.set(color ='k',linewidth = 3)

        # #Add swarm
        # for i in range(len(data)):
        #     # Add some random "jitter" to the x-axis
        #     x = np.random.normal(i, 0.04, size=len(data[i]))
        #     plt.plot(x+1, data[i], 'r.', alpha=0.2)

        # changing color and linewidth of
        # whiskers
        for whisker in bp['whiskers']:
            whisker.set(color ='grey',
                        linewidth = 1)

        # changing color and linewidth of
        # caps
        for cap in bp['caps']:
            cap.set(color ='grey',
                    linewidth = 1)

        plt.title('DockQ scores for the test set',fontsize=20)
        plt.xticks(np.arange(1,len(xlabels)+1),xlabels,fontsize=12)
        plt.ylabel('DockQ')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.tight_layout()
        plt.savefig(outdir+'DockQ_box_test'+str(i)+'.svg',format='svg',dpi=300)
        plt.close()




def AF_vs_RF_marks(marks_dockq_RF,marks_dockq_AF, outdir):
    '''Compare the scores for RF vs AF
    '''

    merged = pd.merge(marks_dockq_RF,marks_dockq_AF,on='complex_id',how='inner')
    print('Number of complexes in merged Marks RF and AF', len(merged))
    #Plot
    fig,ax = plt.subplots(figsize=(12/2.54,12/2.54))
    plt.scatter(merged['DockQ_dockqstats_marks_RF'],merged['DockQ_dockqstats_marks_af2_af2andhhblitsmsa_model_1_rec10_run1'],s=1)

    plt.plot([0,1],[0,1],linewidth=1,linestyle='--',color='grey')
    #Plot correct cutoff
    plt.plot([0.23,0.23],[0,0.23],linewidth=1,linestyle='--',color='k')
    plt.plot([0,0.23],[0.23,0.23],linewidth=1,linestyle='--',color='k',label='Success cutoff')
    plt.title('RF vs AF2 performance on the test set')
    plt.xlabel('RF DockQ')
    plt.ylabel('AF DockQ')
    plt.legend()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(outdir+'RF_vs_AF_test.svg',format='svg',dpi=300)

    #Get num correct
    num_correct_RF = np.argwhere(merged['DockQ_dockqstats_marks_RF'].values>=0.23).shape[0]
    num_correct_AF = np.argwhere(merged['DockQ_dockqstats_marks_af2_af2andhhblitsmsa_model_1_rec10_run1'].values>=0.23).shape[0]
    num_total = len(merged)
    print('Success rate RF:',num_correct_RF,'out of',num_total,'|',np.round(100*num_correct_RF/num_total,2),'%')
    print('Success rate AF:',num_correct_AF,'out of',num_total,'|',np.round(100*num_correct_AF/num_total,2),'%')

    #Get where RF outperforms AF
    scores = merged[['DockQ_dockqstats_marks_RF','DockQ_dockqstats_marks_af2_af2andhhblitsmsa_model_1_rec10_run1']].values
    rf_pos = scores[np.argwhere(scores[:,0]>=0.23)[:,0],:]
    max_scores = np.argmax(rf_pos,axis=1)
    print('RF outperform AF', np.argwhere(max_scores==0).shape[0], 'out of',len(rf_pos),'times|',np.argwhere(max_scores==0).shape[0]/len(rf_pos))


def AF_vs_GRAMM_marks(marks_dockq_GRAMM, marks_dockq_AF, outdir):
    '''Compare the scores for GRAMM vs AF
    '''

    marks_dockq_GRAMM = marks_dockq_GRAMM.rename(columns={1: 'DockQ GRAMM'})
    marks_dockq_GRAMM['complex_id'] = ['_'.join(x.split('-')) for x in marks_dockq_GRAMM[0]]
    merged = pd.merge(marks_dockq_GRAMM,marks_dockq_AF,on='complex_id',how='inner')
    print('Number of complexes in merged Marks GRAMM and AF', len(merged))

    #Plot
    fig,ax = plt.subplots(figsize=(12/2.54,12/2.54))
    plt.scatter(merged['DockQ GRAMM'],merged['DockQ_dockqstats_marks_af2_af2andhhblitsmsa_model_1_rec10_run1'],s=1)

    plt.plot([0,1],[0,1],linewidth=1,linestyle='--',color='grey')
    #Plot correct cutoff
    plt.plot([0.23,0.23],[0,0.23],linewidth=1,linestyle='--',color='k')
    plt.plot([0,0.23],[0.23,0.23],linewidth=1,linestyle='--',color='k',label='Success cutoff')
    plt.title('GRAMM vs AF2 performance on the test set')
    plt.xlabel('GRAMM DockQ')
    plt.ylabel('AF DockQ')
    plt.legend()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(outdir+'GRAMM_vs_AF_test.svg',format='svg',dpi=300)

    #Get num correct
    num_correct_GRAMM = np.argwhere(merged['DockQ GRAMM'].values>=0.23).shape[0]
    num_correct_AF = np.argwhere(merged['DockQ_dockqstats_marks_af2_af2andhhblitsmsa_model_1_rec10_run1'].values>=0.23).shape[0]
    num_total = len(merged)
    print('Success rate GRAMM:',num_correct_GRAMM,'out of',num_total,'|',np.round(100*num_correct_GRAMM/num_total,2),'%')
    print('Success rate AF:',num_correct_AF,'out of',num_total,'|',np.round(100*num_correct_AF/num_total,2),'%')

    #Get where GRAMM outperforms AF
    scores = merged[['DockQ GRAMM','DockQ_dockqstats_marks_af2_af2andhhblitsmsa_model_1_rec10_run1']].values
    GRAMM_pos = scores[np.argwhere(scores[:,0]>=0.23)[:,0],:]
    max_scores = np.argmax(GRAMM_pos,axis=1)
    print('GRAMM outperform AF', np.argwhere(max_scores==0).shape[0], 'out of',len(GRAMM_pos),'times|',np.argwhere(max_scores==0).shape[0]/len(GRAMM_pos))

def AF_vs_TMint_marks(marks_dockq_TMint, marks_dockq_AF, outdir):
    '''Compare the scores for GRAMM vs AF
    '''

    marks_dockq_TMint = marks_dockq_TMint.rename(columns={'dockq': 'DockQ TMint'})
    merged = pd.merge(marks_dockq_TMint,marks_dockq_AF,on='complex_id',how='inner')
    print('Number of complexes in merged Marks TMint and AF', len(merged))

    #Plot
    fig,ax = plt.subplots(figsize=(12/2.54,12/2.54))
    plt.scatter(merged['DockQ TMint'],merged['DockQ_dockqstats_marks_af2_af2andhhblitsmsa_model_1_rec10_run1'],s=1)

    plt.plot([0,1],[0,1],linewidth=1,linestyle='--',color='grey')
    #Plot correct cutoff
    plt.plot([0.23,0.23],[0,0.23],linewidth=1,linestyle='--',color='k')
    plt.plot([0,0.23],[0.23,0.23],linewidth=1,linestyle='--',color='k',label='Success cutoff')
    plt.title('TMint vs AF2 performance on the test set')
    plt.xlabel('TMint DockQ')
    plt.ylabel('AF DockQ')
    plt.legend()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(outdir+'TMint_vs_AF_test.svg',format='svg',dpi=300)

    #Get num correct
    num_correct_TMint = np.argwhere(merged['DockQ TMint'].values>=0.23).shape[0]
    num_correct_AF = np.argwhere(merged['DockQ_dockqstats_marks_af2_af2andhhblitsmsa_model_1_rec10_run1'].values>=0.23).shape[0]
    num_total = len(merged)
    print('Success rate TMint:',num_correct_TMint,'out of',num_total,'|',np.round(100*num_correct_TMint/num_total,2),'%')
    print('Success rate AF:',num_correct_AF,'out of',num_total,'|',np.round(100*num_correct_AF/num_total,2),'%')

    #Get where GRAMM outperforms AF
    scores = merged[['DockQ TMint','DockQ_dockqstats_marks_af2_af2andhhblitsmsa_model_1_rec10_run1']].values
    TMint_pos = scores[np.argwhere(scores[:,0]>=0.23)[:,0],:]
    max_scores = np.argmax(TMint_pos,axis=1)
    print('TMint outperform AF', np.argwhere(max_scores==0).shape[0], 'out of',len(TMint_pos),'times|',np.argwhere(max_scores==0).shape[0]/len(TMint_pos))



def real_features_marks(marks_dockq_AF, dssp_marks, ifstats_marks, aln_scores_marks, AFneffs_marks, topneffs_marks, outdir):
    '''Compare the separation in the marks dataset for AF using metrics from the
    real structures
    '''
    #Change DSSP df
    dssp_marks['Helix']=dssp_marks.G+dssp_marks.H+dssp_marks.I
    dssp_marks['Sheet']=dssp_marks.E+dssp_marks.B
    dssp_marks['Loop']=dssp_marks[' '].values
    ss = dssp_marks[['Helix','Sheet','Loop']].values #0,1,2
    dssp_marks['ss_class']=np.argmax(dssp_marks[['Helix','Sheet','Loop']].values,axis=1)
    dssp_marks = dssp_marks[['id1','id2','ss_class']]
    #Merge dfs
    dssp_marks['complex_id']=dssp_marks.id1+'-'+dssp_marks.id2
    ifstats_marks['complex_id']=ifstats_marks.id1+'-'+ifstats_marks.id2
    aln_scores_marks['complex_id']=aln_scores_marks.id1+'-'+aln_scores_marks.id2
    aln_scores_marks = aln_scores_marks[['complex_id','aln_score']]

    merged_dssp = pd.merge(marks_dockq_AF,dssp_marks,on=['complex_id'],how='inner')

    merged_if = pd.merge(marks_dockq_AF,ifstats_marks,on=['complex_id'],how='inner')
    merged_if = pd.merge(merged_if,aln_scores_marks,on=['complex_id'],how='inner')

    #AFneffs_marks['complex_id']=[code.replace('-', '_') for code in AFneffs_marks['complex_id']]
    #topneffs_marks['complex_id']=[code.replace('-', '_') for code in topneffs_marks['complex_id']]
    merged_if = pd.merge(merged_if,AFneffs_marks,on=['complex_id'],how='inner')
    merged_if = pd.merge(merged_if,topneffs_marks,on=['complex_id'],how='inner')
    '''
    G = 3-turn helix (310 helix). Min length 3 residues.
    H = 4-turn helix (α helix). Minimum length 4 residues.
    I = 5-turn helix (π helix). Minimum length 5 residues.
    T = hydrogen bonded turn (3, 4 or 5 turn)
    E = extended strand in parallel and/or anti-parallel β-sheet conformation. Min length 2 residues.
    B = residue in isolated β-bridge (single pair β-sheet hydrogen bond formation)
    S = bend (the only non-hydrogen-bond based assignment).
    C = coil (residues which are not in any of the above conformations).
    '''

    print('Num complexes in DSSP feature analysis',len(merged_dssp))
    #Plot success rate per ss class
    ss_classes = {0:'Helix',1:'Sheet',2:'Loop'}
    fig,ax = plt.subplots(figsize=(12/2.54,12/2.54))
    for i in range(3):
        sel = merged_dssp[merged_dssp.ss_class==i]
        success=np.argwhere(sel.top_ranked_model_DockQ_af2.values>=0.23).shape[0]/len(sel)
        print(ss_classes[i],'success rate',np.round(success,3),'over',len(sel),'structures')
        #
        sns.distplot(sel.top_ranked_model_DockQ_af2,label=ss_classes[i]+' : '+str(np.round(100*success,1))+' % successful',hist=False)
    plt.title('DockQ and SS for the test set')
    plt.xlabel('DockQ')
    plt.ylabel('Density')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir+'DockQ_per_SS_marks.svg',format='svg',dpi=300)
    plt.close()

    #Plot feature vs DockQ
    #Get min chain len
    merged_if['smallest chain length'] = np.min(merged_if[['l1','l2']].values,axis=1)
    #Get max chain len
    merged_if['biggest chain length'] = np.max(merged_if[['l1','l2']].values,axis=1)
    vars = ['num_if_contacts_total','smallest chain length', 'biggest chain length', 'aln_score', 'AFdefault_Neff', 'tophit_Neff']
    nicer_names = {'num_if_contacts_total':'number of interface contacts','smallest chain length':'smallest chain length', 'biggest chain length':'biggest chain length',
    'aln_score':'alignment score', 'AFdefault_Neff':'AF Neff', 'tophit_Neff':'Paired Neff'}
    print('Num complexes in real feature analysis',len(merged_if))
    #Plot each third and the distribution vs vars
    for var in vars:
        fig,ax = plt.subplots(figsize=(12/2.54,12/2.54))
        fig,ax = plt.subplots(figsize=(12/2.54,12/2.54))
        print (np.quantile(merged_if[var],0.5,axis=0))
        l=[np.min(merged_if[var])]
        l+=[np.quantile(merged_if[var],0.33,axis=0)]
        l+=[np.quantile(merged_if[var],0.67,axis=0)]
        l+=[np.max(merged_if[var])]
        print (l)
        j=0
        for i in l[0:3]:
            j+=1
            #print ("test: ",i,j,l[j])
            sel = merged_if.loc[ (merged_if[var] > i) &  (merged_if[var] < l[j])  ]
            success=np.argwhere(sel.top_ranked_model_DockQ_af2.values>=0.23).shape[0]/len(sel)
            print(j,str(i)+" - "+ str(l[j])+":",'success rate',np.round(success,3),'over',len(sel),'structures')
            #
            sns.kdeplot(sel.top_ranked_model_DockQ_af2,label=str(round(i,0))+"-"+str(round(l[j],0))+' : '+str(np.round(100*success,1))+' % successful')
        plt.title('DockQ and ' + nicer_names[var] + '\nfor the test set')
        plt.xlabel('DockQ')
        plt.ylabel('Density')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.legend()
        plt.tight_layout()
        plt.savefig(outdir+'DockQ_per_'+var+'_marks.svg',format='svg',dpi=300)
        plt.close()


def marks_dockq_per_org(marks_dockq_AF, oxstats_marks, ifstats_marks, aln_scores_marks, AFneffs_marks, topneffs_marks, outdir):
    '''Analyze the dockq per organism
    '''
    #Merge
    oxstats_marks['complex_id'] = oxstats_marks.id1+'-'+oxstats_marks.id2
    ifstats_marks['complex_id']=ifstats_marks.id1+'-'+ifstats_marks.id2
    #AFneffs_marks['complex_id']=[code.replace('-', '_') for code in AFneffs_marks['complex_id']]
    #topneffs_marks['complex_id']=[code.replace('-', '_') for code in topneffs_marks['complex_id']]
    aln_scores_marks['complex_id']=aln_scores_marks.id1+'-'+aln_scores_marks.id2
    aln_scores_marks = aln_scores_marks[['complex_id','aln_score']]
    merged = pd.merge(marks_dockq_AF,oxstats_marks,on='complex_id',how='left')
    merged = pd.merge(merged,ifstats_marks,on=['complex_id'],how='inner')
    merged = pd.merge(merged,aln_scores_marks,on=['complex_id'],how='inner')
    merged = pd.merge(merged,AFneffs_marks,on=['complex_id'],how='inner')
    merged = pd.merge(merged,topneffs_marks,on=['complex_id'],how='inner')

    #Get min chain len
    merged['smallest chain length'] = np.min(merged[['l1','l2']].values,axis=1)
    #Get max chain len
    merged['biggest chain length'] = np.max(merged[['l1','l2']].values,axis=1)

    organisms = ['Homo sapiens','Saccharomyces cerevisiae', 'Escherichia coli']
    vars = ['num_if_contacts_total','smallest chain length', 'biggest chain length', 'aln_score','AFdefault_Neff', 'tophit_Neff']
    #Save
    orgs = []
    dockq_scores = []
    fig,ax = plt.subplots(figsize=(12/2.54,12/2.54))

    for org in organisms:
        sel = merged[merged.Org1==org]
        sel = sel[sel.Org2==org]
        print('Number of complexes for',org,len(sel))
        #Successs rate
        sel_scores = sel.top_ranked_model_DockQ_af2.values
        sr = np.argwhere(sel_scores>=0.23).shape[0]/len(sel_scores)
        print('Success rate',sr)
        #correlation
        for var in vars:
            R,p = spearmanr(sel[var].values,sel['top_ranked_model_DockQ_af2'].values)
            print(var, np.round(R,2))
        if org =='Saccharomyces cerevisiae':
            org = 'S.cerevisiae'
        if org =='Escherichia coli':
            org = 'E.coli'
        sns.distplot(sel_scores,label=org+' : '+str(np.round(sr*100,1))+' % successful',hist=False)


    plt.title('DockQ per organism for the test set')
    plt.xlabel('DockQ')
    plt.ylabel('Density')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir+'DockQ_per_org_marks.svg',format='svg',dpi=300)
    plt.close()

def marks_dockq_per_kingdom(marks_dockq_AF, oxstats_marks, AFneffs_marks, topneffs_marks, outdir):
    '''Analyze the dockq per organism
    '''
    #Merge
    oxstats_marks['complex_id'] = oxstats_marks.id1+'-'+oxstats_marks.id2
    #AFneffs_marks['complex_id']=['_'.join(x.split('-')) for x in AFneffs_marks.complex_id]
    #topneffs_marks['complex_id']=['_'.join(x.split('-')) for x in topneffs_marks.complex_id]
    merged = pd.merge(marks_dockq_AF,oxstats_marks,on='complex_id',how='left')
    merged = pd.merge(merged,AFneffs_marks,on=['complex_id'],how='inner')
    merged = pd.merge(merged,topneffs_marks,on=['complex_id'],how='inner')
    kingdoms = ['E', 'B', 'A', 'V']
    nice_labels = {'top_ranked_model_DockQ_af2':'DockQ', 'AFdefault_Neff':'AF Neff', 'tophit_Neff':'Paired Neff'}
    for var in ['top_ranked_model_DockQ_af2', 'AFdefault_Neff', 'tophit_Neff']:
        fig,ax = plt.subplots(figsize=(12/2.54,12/2.54))
        for kd in kingdoms:
            sel = merged[merged.kingdom1==kd]
            sel = sel[sel.kingdom2==kd]

            #Successs rate
            sel_scores = sel[var].values
            if var=='top_ranked_model_DockQ_af2':
                sr = np.argwhere(sel_scores>=0.23).shape[0]/len(sel_scores)
                print('Success rate for',kd,sr,len(sel_scores))
                sns.distplot(sel_scores,label=kd+' : '+str(np.round(sr*100,1))+' % successful',hist=False)
            else:
                sns.distplot(sel_scores,label=kd,hist=False)


        plt.title(nice_labels[var]+' per kingdom for the test set')
        plt.xlabel(nice_labels[var])
        plt.ylabel('Density')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.legend()
        plt.tight_layout()
        plt.savefig(outdir+var+'_per_kd_marks.svg',format='svg',dpi=300)
        plt.close()


def marks_dockq_vs_aln_overlap(marks_dockq_AF, af_chain_overlap_marks, outdir):
    '''Analyze the dockq vs chain overlap
    '''
    #Merge
    cid = ['_'.join(x.split('-')) for x in af_chain_overlap_marks.complex_id.values]
    af_chain_overlap_marks['complex_id']=cid
    merged = pd.merge(marks_dockq_AF,af_chain_overlap_marks,on='complex_id',how='inner')
    #Plot tertiles
    fig,ax = plt.subplots(figsize=(12/2.54,12/2.54))

    l=[np.min(merged.Overlap)]
    l+=[np.quantile(merged.Overlap,0.33,axis=0)]
    l+=[np.quantile(merged.Overlap,0.67,axis=0)]
    l+=[np.max(merged.Overlap)]

    j=0
    for i in l[0:3]:
        j+=1
        sel = merged.loc[ (merged['Overlap'] > i) &  (merged['Overlap'] < l[j])  ]
        success=np.argwhere(sel.DockQ_dockqstats_marks_af2_af2stdmsa_model_1_rec10.values>=0.23).shape[0]/len(sel)
        print(j,str(i)+" - "+ str(l[j])+":",'success rate',np.round(success,3),'over',len(sel),'structures')
        #
        sns.kdeplot(sel.DockQ_dockqstats_marks_af2_af2stdmsa_model_1_rec10,label=str(round(i,2))+"-"+str(round(l[j],2))+' : '+str(np.round(100*success,1))+' % successful')


    plt.title('DockQ vs chain overlap in AF2 msas')
    plt.xlabel('DockQ')
    plt.ylabel('Density')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir+'dockq_vs_overlap.svg',format='svg',dpi=300)
    plt.close()
    #Plot overlap distribution
    fig,ax = plt.subplots(figsize=(12/2.54,12/2.54))
    sns.distplot(merged.Overlap)
    plt.title('Chain overlap distribution in AF2 msas')
    plt.xlabel('Overlap fraction')
    plt.ylabel('Density')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(outdir+'overlap_distr.svg',format='svg',dpi=300)
    plt.close()


def score_newset_5runs(newset_dockq_AF, plDDT_newset, sigmoid_params, outdir):
    '''Compare the separation in the newset dataset for AF using metrics from the
    predicted structures
    '''
    #Merge dfs
    plDDT_newset['complex_id'] = plDDT_newset.id1+'-'+plDDT_newset.id2
    merged = pd.merge(newset_dockq_AF,plDDT_newset,on=['complex_id'],how='inner')
    #Get num res in interface
    separator1 = merged[['num_atoms_in_interface_1', 'num_atoms_in_interface_2','num_atoms_in_interface_3','num_atoms_in_interface_4','num_atoms_in_interface_5']].values
    separator2 = merged[['num_atoms_in_interface_1', 'num_atoms_in_interface_2','num_atoms_in_interface_3','num_atoms_in_interface_4','num_atoms_in_interface_5']].values
    separator = separator1*np.log10(separator2+1)
    def sigmoid(x, L ,x0, k, b):
        y = L / (1 + np.exp(-k*(x-x0)))+b
        return (y)
    separator = sigmoid(separator, *sigmoid_params)


    scores = merged[newset_dockq_AF.columns[1:]].values
    #Get max and min scores
    max_scores = np.max(scores,axis=1)
    min_scores = np.min(scores,axis=1)

    #Get success rates per initializations
    srs = []
    for i in range(scores.shape[1]):
        srs.append(np.argwhere(scores[:,i]>=0.23).shape[0]/len(scores))
    print('New set scoring:')
    print('Success rate 5 runs top scores',np.argwhere(max_scores>=0.23).shape[0]/len(max_scores))
    print('Av diff',np.average(max_scores-min_scores))
    print('Std diff',np.std(max_scores-min_scores))
    print('Avg and std success rate', np.average(srs),np.std(srs))
    #Separate the models using the number of contacts in the interface
    max_inds = np.argmax(separator,axis=1)
    first_ranked_scores = []
    first_ranked_separators = []
    #Get max separator scores
    for i in range(len(max_inds)):
        first_ranked_scores.append(scores[i,max_inds[i]])
        first_ranked_separators.append(separator[i,max_inds[i]])
    #Convert to array
    first_ranked_scores = np.array(first_ranked_scores)
    first_ranked_separators = np.array(first_ranked_separators)
    #Get success rate
    print('Ranking test set success rate using num contacts in interface',np.argwhere(first_ranked_scores>=0.23).shape[0]/len(first_ranked_scores))

    #Plot
    fig,ax = plt.subplots(figsize=(12/2.54,12/2.54))
    plt.scatter(first_ranked_scores, max_scores,s=10,color='darkblue',label='Max')
    plt.scatter(first_ranked_scores, min_scores,s=10,color='mediumseagreen',label='Min')
    plt.title('Model ranking from 5 runs on the new dimer set\n(both MSAs, model 1, 10 recycles)')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.plot([0,1],[0,1],color='k',linewidth=1,linestyle='--')
    plt.xlabel('DockQ first ranked model')
    plt.ylabel('DockQ')
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir+'DockQ_newset_5runs.svg',format='svg',dpi=300)
    plt.close()




def dev_vs_test(marks_dockq_AF, oxstats_marks, ifstats_marks, aln_scores_marks, AFneffs_marks,
                topneffs_marks, bench4_kingdom,  dssp_bench4, AFneffs_bench4, topneffs_bench4, outdir):
    '''Analyze the distributions of different features for the dev vs the test set
    Neff
    Kingdom
    SS in interface
    Number of interface contacts
    Chain length (biggest and smallest)
    '''

    #Merge bench4
    bench4_merged = pd.merge(bench4_kingdom,dssp_bench4,on=['id1','id2'],how='inner')
    bench4_merged['complex_id'] = bench4_merged.PDB+'_u1-'+bench4_merged.PDB+'_u2'
    bench4_merged = pd.merge(bench4_merged,AFneffs_bench4,on='complex_id',how='inner')
    bench4_merged = pd.merge(bench4_merged,topneffs_bench4,on='complex_id',how='inner')
    bench4_merged['min_chain_len'] = np.min(bench4_merged[['Sequence Length 1','Sequence Length 2']].values,axis=1)
    bench4_merged['max_chain_len'] = np.max(bench4_merged[['Sequence Length 1','Sequence Length 2']].values,axis=1)
    bench4_merged['sum_chain_len'] = np.sum(bench4_merged[['Sequence Length 1','Sequence Length 2']].values,axis=1)
    bench4_kingdom['Kingdom'] = bench4_kingdom['Kingdom'].replace({' Bacteria ':'B', ' Eukaryota ':'E','Virus':'V'})
    bench4_merged['if_fraction'] = np.divide(bench4_merged['num_if_contacts_total'],bench4_merged['sum_chain_len'])


    #Merge Marks
    marks_merged = pd.merge(oxstats_marks, ifstats_marks, on=['id1','id2'],how='inner')
    marks_merged['complex_id'] = marks_merged.id1+'_'+marks_merged.id2
    marks_merged = pd.merge(marks_merged,AFneffs_marks,on='complex_id',how='inner')
    marks_merged = pd.merge(marks_merged,topneffs_marks,on='complex_id',how='inner')
    marks_merged['min_chain_len'] = np.min(marks_merged[['l1','l2']].values,axis=1)
    marks_merged['max_chain_len'] = np.max(marks_merged[['l1','l2']].values,axis=1)
    marks_merged['sum_chain_len'] = np.sum(marks_merged[['l1','l2']].values,axis=1)
    marks_merged['if_fraction'] = np.divide(marks_merged['num_if_contacts_total'],marks_merged['sum_chain_len'])

    #Get kingdom fractions
    kingdoms = ['E', 'B', 'A', 'V']
    print('KD','Bench4','Marks')
    for kd in kingdoms:
        sel_bench4 = bench4_kingdom[bench4_kingdom.Kingdom==kd]
        sel_marks = marks_merged[(marks_merged.kingdom1==kd)&(marks_merged.kingdom2==kd)]
        print(kd,len(sel_bench4)/len(bench4_kingdom),len(sel_marks)/len(marks_merged))

    #Plot vars
    vars = ['num_if_contacts_total', 'min_chain_len','max_chain_len', 'sum_chain_len', 'AFdefault_Neff' ,'tophit_Neff','if_fraction']
    for var in vars:
        fig,ax = plt.subplots(figsize=(12/2.54,12/2.54))
        sns.distplot(bench4_merged[var],label='Dev. set',hist=True,kde=True,norm_hist=True)
        sns.distplot(marks_merged[var],label='Test. set',hist=True,kde=True,norm_hist=True)
        plt.legend()
        plt.title('Dev. vs Test '+var)
        plt.xlabel(var)
        plt.ylabel('Density')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.tight_layout()
        plt.savefig(outdir+'dev_vs_test_'+var+'.svg',format='svg',dpi=300)
        plt.close()

def neg_vs_pos(plDDT_marks_af, plDDT_marks_negative_af, plDDT_negatome_af, sigmoid_params):
    '''Compare the interfaces of positive and neg marks set + negatome
    '''

    #Filter out the homodimers from the negatome
    keep_inds = []

    for i in range(len(plDDT_negatome_af)):
        row = plDDT_negatome_af.loc[i]
        if row.id1!=row.id2:
            keep_inds.append(i)
    print('Num homodimers:',len(plDDT_negatome_af)-len(keep_inds))
    plDDT_negatome_af = plDDT_negatome_af.loc[keep_inds]

    #Get AUC using the different metrics
    #Get min of chains
    #Pos
    single_chain_plddt = np.min(plDDT_marks_af[['ch1_plddt_av_1', 'ch2_plddt_av_1']].values,axis=1)
    plDDT_marks_af['min_chain_plddt_av_1'] = single_chain_plddt
    #Neg marks
    single_chain_plddt = np.min(plDDT_marks_negative_af[['ch1_plddt_av', 'ch2_plddt_av']].values,axis=1)
    plDDT_marks_negative_af['min_chain_plddt_av'] = single_chain_plddt
    #Negatome
    single_chain_plddt = np.min(plDDT_negatome_af[['ch1_plddt_av', 'ch2_plddt_av']].values,axis=1)
    plDDT_negatome_af['min_chain_plddt_av'] = single_chain_plddt


    #Analyze ROC as a function of
    feature_nice_names = {'if_plddt_av':'IF_plDDT',  'min_chain_plddt_av':'Min plDDT per chain',
                    'plddt_av':'Average plDDT', 'num_atoms_in_interface':'IF_contacts',
                        'num_res_in_interface':'IF_residues'}
    colors = {'if_plddt_av':'darkblue','min_chain_plddt_av':'magenta','plddt_av':'orange',
    'num_atoms_in_interface':'darkgreen','num_res_in_interface':'tab:blue', 'IF_cp':'cyan', 'pDockQ':'k'}


    fig,ax = plt.subplots(figsize=(12/2.54,12/2.54))
    for key in feature_nice_names:
        pos_features = plDDT_marks_af[key+'_1'].values
        neg_features = np.concatenate([plDDT_marks_negative_af[key].values, plDDT_negatome_af[key].values])
        #ROC
        correct = np.zeros(len(pos_features)+len(neg_features))
        correct[:len(pos_features)]=1
        all_features = np.concatenate([pos_features, neg_features])
        fpr, tpr, threshold = metrics.roc_curve(correct, all_features, pos_label=1)
        roc_auc = metrics.auc(fpr, tpr)
        #Plot ROC
        plt.plot(fpr, tpr, label = feature_nice_names[key]+': AUC = %0.2f' % roc_auc, color=colors[key])
        #TPRs
        print(key,'TPR at FPR 1%=',np.round(100*tpr[np.argwhere(np.round(fpr,2)<=0.01)[-1][0]]))
        print(key,'TPR at FPR 5%=',np.round(100*tpr[np.argwhere(np.round(fpr,2)<=0.05)[-1][0]]))

    #Add log(if contacts)*if_plddt_av
    pos_features_if_cp = plDDT_marks_af['if_plddt_av_1'].values*np.log10(plDDT_marks_af['num_atoms_in_interface_1'].values+1)
    neg_features_marks_if_cp = plDDT_marks_negative_af['if_plddt_av'].values*np.log10(plDDT_marks_negative_af['num_atoms_in_interface'].values+1)
    neg_features_negatome_if_cp = plDDT_negatome_af['if_plddt_av'].values*np.log10(plDDT_negatome_af['num_atoms_in_interface'].values+1)
    neg_features_if_cp = np.concatenate([neg_features_marks_if_cp, neg_features_negatome_if_cp])
    correct = np.zeros(len(pos_features_if_cp)+len(neg_features_if_cp))
    correct[:len(pos_features_if_cp)]=1
    all_features = np.concatenate([pos_features_if_cp, neg_features_if_cp])
    #Create ROC
    fpr, tpr, threshold = metrics.roc_curve(correct, all_features, pos_label=1)
    roc_auc = metrics.auc(fpr, tpr)
    #plt.plot(fpr, tpr, label = 'IF_plDDT⋅log(IF_contacts)'+': AUC = %0.2f' % roc_auc,color='tab:cyan')

    #Do the same with pDockQ
    def sigmoid(x, L ,x0, k, b):
        y = L / (1 + np.exp(-k*(x-x0)))+b
        return (y)
    pos_features_pdockq = sigmoid(pos_features_if_cp, *sigmoid_params)
    neg_features_pdockq = sigmoid(neg_features_if_cp, *sigmoid_params)
    all_features = np.concatenate([pos_features_pdockq, neg_features_pdockq])
    #Create ROC
    fpr, tpr, threshold = metrics.roc_curve(correct, all_features, pos_label=1)
    roc_auc = metrics.auc(fpr, tpr)
    plt.plot(fpr, tpr, label = 'pDockQ'+': AUC = %0.2f' % roc_auc,color='k',linestyle='--')
    #TPRs
    print('pDockQ TPR at FPR 1%=',np.round(100*tpr[np.argwhere(np.round(fpr,2)<=0.01)[-1][0]]))
    print('pDockQ TPR at FPR 5%=',np.round(100*tpr[np.argwhere(np.round(fpr,2)<=0.05)[-1][0]]))
    #Plot formatting
    plt.plot([0,1],[0,1],linewidth=1,linestyle='--',color='grey')
    plt.legend(fontsize=9)
    plt.title('Identifying interacting proteins\nROC as a function of different metrics')
    plt.xlabel('FPR')
    plt.ylabel('TPR')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(outdir+'ROC_pos_neg.svg',format='svg',dpi=300)

    #Marks comparison
    print('Only marks negative')
    neg_features = sigmoid(neg_features_marks_if_cp, *sigmoid_params)
    #ROC
    correct = np.zeros(len(pos_features_pdockq)+len(neg_features))
    correct[:len(pos_features)]=1
    all_features = np.concatenate([pos_features_pdockq, neg_features])
    fpr, tpr, threshold = metrics.roc_curve(correct, all_features, pos_label=1)
    roc_auc = metrics.auc(fpr, tpr)

    t=0.01
    print('Average interface pDockQ TPR at FPR '+str(t*100)+'%=',np.round(100*tpr[np.argwhere(np.round(fpr,2)<=t)[-1][0]],2),'\nAUC:',roc_auc, 'FPR:',np.round(100*fpr[np.argwhere(np.round(fpr,2)<=t)[-1][0]],2))
    t=0.05
    print('Average interface pDockQ TPR at FPR '+str(t*100)+'%=',np.round(100*tpr[np.argwhere(np.round(fpr,2)<=t)[-1][0]],2),'\nAUC:',roc_auc, 'FPR:',np.round(100*fpr[np.argwhere(np.round(fpr,2)<=t)[-1][0]],2))

    #Plot distribution of separators
    feature_nice_names = {'if_plddt_av':'IF_plDDT', 'num_atoms_in_interface':'IF_contacts',
                         'pDockQ':'pDockQ'}
    xlims = {'if_plddt_av':[-20,120],  'num_atoms_in_interface':[-100,500],
            'IF_cp':[0,250], 'pDockQ':[0,1]}
    bins = {'if_plddt_av':20,  'num_atoms_in_interface':50,
            'IF_cp':20, 'pDockQ':20}
    matplotlib.rcParams.update({'font.size': 9})
    for key in feature_nice_names:
        fig,ax = plt.subplots(figsize=(6/2.54,6/2.54))
        if key not in ['IF_cp','pDockQ']:
            pos_features = plDDT_marks_af[key+'_1'].values
            neg_features = np.concatenate([plDDT_marks_negative_af[key].values, plDDT_negatome_af[key].values])
        if key=='IF_cp':
            pos_features = pos_features_if_cp
            neg_features = neg_features_if_cp
        if key=='pDockQ':
            pos_features = pos_features_pdockq
            neg_features = neg_features_pdockq

        plt.hist(pos_features,label='pos',color=colors[key],alpha=0.75,bins=bins[key],density=True)
        plt.hist(neg_features,label='neg',color='gray',alpha=0.75,bins=bins[key],density=True)
        plt.legend(fontsize=9)
        #plt.title('Distribution of '+feature_nice_names[key],fontsize=9)
        plt.xlim(xlims[key])
        plt.xlabel(feature_nice_names[key])
        plt.ylabel('Density')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.tight_layout()
        plt.savefig(outdir+key+'_pos_vs_neg_distr.svg',format='svg',dpi=300)
        plt.close()
        # #Create df of fpr vs tpr
        # roc_df = pd.DataFrame()
        # roc_df['FPR']=fpr
        # roc_df['TPR']=tpr
        # roc_df['Number of interface contacts'] = threshold

def ppv_vs_dockq_marks(marks_dockq_AF,ifstats_marks,outdir):
    '''Analysis of the relationship between if stats and DockQ
    '''

    ifstats_marks['complex_id']=ifstats_marks.id1+'-'+ifstats_marks.id2

    merged_if = pd.merge(marks_dockq_AF,ifstats_marks,on=['complex_id'],how='inner')
    #Calc PPV
    merged_if['PPV'] = merged_if['num_accurate_if']/merged_if['num_if_contacts_total']
    #Plot
    fig,ax = plt.subplots(figsize=(12/2.54,12/2.54))
    plt.scatter(merged_if.PPV,merged_if.DockQ_dockqstats_marks_af2_af2andhhblitsmsa_model_1_rec10_run1,
    s=1,label='AF2+Paired',c='tab:blue')
    plt.scatter(merged_if.PPV,merged_if.DockQ_dockqstats_marks_af2_hhblitsn2_model_1_rec10,
    s=1,label='Paired',c='tab:orange')

    #RA
    ra_x = []
    ra_y_paired = []
    ra_y_both = []
    sr_paired = []
    sr_both = []
    step = 0.05
    for i in np.arange(0,0.5,step):
        sel = merged_if[(merged_if.PPV>=i)&(merged_if.PPV<i+0.1)]
        if len(sel)<1:
            continue
        ra_x.append(i+step/2)
        #ra_y_paired.append(np.average(sel.DockQ_dockqstats_marks_af2_hhblitsn2_model_1_rec10))
        #ra_y_both.append(np.average(sel.DockQ_dockqstats_marks_af2_af2andhhblitsmsa_model_1_rec10_run1))
        #SR
        sr_paired.append(np.argwhere(sel.DockQ_dockqstats_marks_af2_hhblitsn2_model_1_rec10.values>=0.23).shape[0]/len(sel))
        sr_both.append(np.argwhere(sel.DockQ_dockqstats_marks_af2_af2andhhblitsmsa_model_1_rec10_run1.values>=0.23).shape[0]/len(sel))
    #RA
    # plt.plot(ra_x,ra_y_paired,label='RA Paired',c='tab:orange')
    # plt.plot(ra_x,ra_y_both,label='RA AF2+Paired',c='tab:blue')
    #SR
    plt.plot(ra_x,sr_paired,label='SR Paired',c='tab:orange')
    plt.plot(ra_x,sr_both,label='SR AF2+Paired',c='tab:blue')
    plt.legend()
    plt.title('Interface PPV vs DockQ')
    plt.xlabel('Interface PPV')
    plt.ylabel('DockQ')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(outdir+'ppv_vs_dockq.svg',format='svg',dpi=300)
    plt.close()



#################MAIN####################

#Parse args
args = parser.parse_args()
#Data
#bench4
bench4_dockq_aa = pd.read_csv(args.bench4_dockq_aa[0])
bench4_dockq_RF = pd.read_csv(args.bench4_dockq_RF[0])
plDDT_bench4 = pd.read_csv(args.plDDT_bench4[0])
#pconsdock_bench4 = pd.read_csv(args.pconsdock_bench4[0])
#pconsdock_marks = pd.read_csv(args.pconsdock_marks[0])
bench4_kingdom = pd.read_csv(args.bench4_kingdom[0])
dssp_bench4 = pd.read_csv(args.dssp_bench4[0])
AFneffs_bench4 = pd.read_csv(args.afdefault_neff_bench4[0])
topneffs_bench4 = pd.read_csv(args.tophits_neff_bench4[0])
#Marks positive
marks_dockq_RF = pd.read_csv(args.marks_dockq_RF[0])
marks_dockq_AF_bb = pd.read_csv(args.marks_dockq_AF_bb[0])
marks_dockq_AF_aa = pd.read_csv(args.marks_dockq_AF_aa[0])
marks_dockq_GRAMM = pd.read_csv(args.marks_dockq_GRAMM[0],header=None)
marks_dockq_TMfull = pd.read_csv(args.marks_dockq_TMfull[0])
marks_dockq_TMint = pd.read_csv(args.marks_dockq_TMint[0])
marks_dockq_mdockpp = pd.read_csv(args.marks_dockq_mdockpp[0])

plDDT_marks_af = pd.read_csv(args.plDDT_marks_af[0])
plDDT_marks_fused = pd.read_csv(args.plDDT_marks_fused[0])
dssp_marks = pd.read_csv(args.dssp_marks[0])
ifstats_marks = pd.read_csv(args.ifstats_marks[0])
aln_scores_marks = pd.read_csv(args.aln_scores_marks[0])
oxstats_marks = pd.read_csv(args.oxstats_marks[0])
AFneffs_marks = pd.read_csv(args.afdefault_neff_marks[0])
topneffs_marks = pd.read_csv(args.tophits_neff_marks[0])
af_chain_overlap_marks = pd.read_csv(args.af_chain_overlap_marks[0])
#Marks negative
plDDT_marks_negative_af = pd.read_csv(args.plDDT_marks_negative_af[0])
#Negatome
plDDT_negatome_af = pd.read_csv(args.plDDT_negatome_af[0])
#New set
newset_dockq_AF = pd.read_csv(args.newset_dockq_AF[0])
plDDT_newset = pd.read_csv(args.plDDT_newset[0])
#Outdir
outdir = args.outdir[0]

#Visualize score distribution
#Dev set
#dockq_box(bench4_dockq_aa, outdir)
#correlate_scores(bench4_dockq_aa, outdir)
# AF_vs_RF_bench4(bench4_dockq_aa, bench4_dockq_RF, outdir)
#
# # #Test set
#Get marks_dockq_AF_aa for the missing in marks_dockq_AF_aa
marks_dockq_AF = fetch_missing_dockq(marks_dockq_AF_bb,marks_dockq_AF_aa)
sigmoid_params = ROC_pred_marks(marks_dockq_AF, plDDT_marks_af, outdir)

marks_dockq_AF = score_marks_5runs_paired_af(marks_dockq_AF, plDDT_marks_af, sigmoid_params, outdir)
marks_dockq_AF = score_marks_5runs_paired_fused(marks_dockq_AF, plDDT_marks_fused, sigmoid_params, outdir)
marks_box(marks_dockq_AF, marks_dockq_GRAMM, marks_dockq_mdockpp, marks_dockq_TMfull, marks_dockq_TMint, marks_dockq_RF,outdir)

#AF_vs_RF_marks(marks_dockq_RF,marks_dockq_AF, outdir)
#AF_vs_GRAMM_marks(marks_dockq_GRAMM,marks_dockq_AF, outdir)
#AF_vs_TMint_marks(marks_dockq_TMint, marks_dockq_AF, outdir)

real_features_marks(marks_dockq_AF, dssp_marks, ifstats_marks, aln_scores_marks, AFneffs_marks, topneffs_marks, outdir)
marks_dockq_per_org(marks_dockq_AF, oxstats_marks, ifstats_marks, aln_scores_marks, AFneffs_marks, topneffs_marks, outdir)
marks_dockq_per_kingdom(marks_dockq_AF, oxstats_marks, AFneffs_marks, topneffs_marks, outdir)

#marks_dockq_vs_aln_overlap(marks_dockq_AF, af_chain_overlap_marks, outdir)

# New set
#score_newset_5runs(newset_dockq_AF, plDDT_newset, sigmoid_params, outdir)

# #Dev set vs test set
# dev_vs_test(marks_dockq_AF, oxstats_marks, ifstats_marks, aln_scores_marks, AFneffs_marks,
# topneffs_marks, bench4_kingdom, dssp_bench4, AFneffs_bench4, topneffs_bench4, outdir)

#Interaction/not: Marks negative and positive
neg_vs_pos(plDDT_marks_fused, plDDT_marks_negative_af, plDDT_negatome_af, sigmoid_params)

#PPV vs DockQ
ppv_vs_dockq_marks(marks_dockq_AF,ifstats_marks,outdir)
