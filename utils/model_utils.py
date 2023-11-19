#!/usr/bin/python3.9
from sklearn.preprocessing import LabelEncoder
from sklearn.multioutput import MultiOutputClassifier
from sklearn.metrics import auc, roc_curve
from sklearn.preprocessing import label_binarize
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
from sklearn.multiclass import OneVsRestClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
from sklearn.multiclass import OneVsRestClassifier
from scipy import interp
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns

def rf_3class_model(data1, data2, data3, names, train_size, test_size, color_curve=['blue', 'green', 'red'], shared_roc=False):
    # Define lists and arrays
    list_features = []
    n_iter = 100
    arr_pred = np.zeros(100)
    imp = []
    mean_fpr = np.linspace(0, 1, 100)
    roc_auc_overall = []
    tprs_list = []
    feat_imp_list = []
    for i in range(n_iter):
        print(f'Iteration number: {str(i+1)}')
        # Generate a training set made by $train_size features each dataset
        train_set1 = (data1[['peptide_name', 'is_pos_cntrl', 'is_neg_cntrl', 'is_rand_cntrl', 'is_auto', 'is_infect',
               'is_EBV', 'is_toxin', 'is_EM', 'is_MPA', 'is_patho', 'is_IgA',
               'is_probio', 'is_flagellum', 'diamond_mmseqs_intersec_toxin', 'is_topgraph_new_&_old', 'signalp6_slow']]).fillna(0).sample(train_size)*1
        train_set1['dataset'] = names[0]
        train_set2 = (data2[['peptide_name', 'is_pos_cntrl', 'is_neg_cntrl', 'is_rand_cntrl', 'is_auto', 'is_infect',
               'is_EBV', 'is_toxin', 'is_EM', 'is_MPA', 'is_patho', 'is_IgA',
               'is_probio', 'is_flagellum', 'diamond_mmseqs_intersec_toxin', 'is_topgraph_new_&_old', 'signalp6_slow']]).fillna(0).sample(train_size)*1
        train_set2['dataset'] = names[1]

        #agilent_total_clean = df_agilent[df_agilent['peptide'].isin(df_agilent_final['peptide'].tolist() + df_not_HLA_agilent['peptide'].tolist()) == False]

        train_set3 = (data3[['peptide_name', 'is_pos_cntrl', 'is_neg_cntrl', 'is_rand_cntrl', 'is_auto', 'is_infect',
               'is_EBV', 'is_toxin', 'is_EM', 'is_MPA', 'is_patho', 'is_IgA',
               'is_probio', 'is_flagellum', 'diamond_mmseqs_intersec_toxin', 'is_topgraph_new_&_old', 'signalp6_slow']]).fillna(0).sample(train_size)*1
        train_set3['dataset'] = names[2]


        test_set1 = data1[data1['peptide_name'].isin(train_set1['peptide_name']) == False][['peptide_name', 'is_pos_cntrl', 'is_neg_cntrl', 'is_rand_cntrl', 'is_auto', 'is_infect',
               'is_EBV', 'is_toxin', 'is_EM', 'is_MPA', 'is_patho', 'is_IgA',
               'is_probio', 'is_flagellum', 'diamond_mmseqs_intersec_toxin', 'is_topgraph_new_&_old', 'signalp6_slow']].fillna(0).sample(test_size)*1
        test_set1['dataset'] = names[0]


        test_set2 = data2[data2['peptide_name'].isin(train_set2['peptide_name']) == False][['peptide_name', 'is_pos_cntrl', 'is_neg_cntrl', 'is_rand_cntrl', 'is_auto', 'is_infect',
               'is_EBV', 'is_toxin', 'is_EM', 'is_MPA', 'is_patho', 'is_IgA',
               'is_probio', 'is_flagellum', 'diamond_mmseqs_intersec_toxin', 'is_topgraph_new_&_old', 'signalp6_slow']].fillna(0).sample(test_size)*1
        test_set2['dataset'] = names[1]


        test_set3 = data3[data3['peptide_name'].isin(train_set3['peptide_name']) == False][['peptide_name', 'is_pos_cntrl', 'is_neg_cntrl', 'is_rand_cntrl', 'is_auto', 'is_infect',
               'is_EBV', 'is_toxin', 'is_EM', 'is_MPA', 'is_patho', 'is_IgA', 'is_PNP', 'is_nonPNP_strains',
               'is_probio', 'is_flagellum', 'diamond_mmseqs_intersec_toxin', 'is_topgraph_new_&_old', 'signalp6_slow']].fillna(0).sample(test_size)*1
        test_set3['dataset'] = names[2]

        # Concat test and training sets
        training_set = pd.concat([train_set1, train_set2, train_set3])
        test_set = pd.concat([test_set1, test_set2, test_set3])
        training_set

        training_set_x = training_set[['is_pos_cntrl', 'is_neg_cntrl', 'is_rand_cntrl', 'is_auto', 'is_infect',
               'is_EBV', 'is_toxin', 'is_EM', 'is_MPA', 'is_patho', 'is_IgA',
               'is_probio', 'is_flagellum', 'diamond_mmseqs_intersec_toxin', 'is_topgraph_new_&_old', 'signalp6_slow']]
        training_set_labels = training_set['dataset']
        training_set_labels = training_set['dataset']
        label_encoder = LabelEncoder()
        label_encoder.fit(training_set['dataset'])
        training_set_labels = label_encoder.transform(training_set_labels)
        classes = label_encoder.classes_
        n_classes = len(classes)

        test_set_x = test_set [['is_pos_cntrl', 'is_neg_cntrl', 'is_rand_cntrl', 'is_auto', 'is_infect',
               'is_EBV', 'is_toxin', 'is_EM', 'is_MPA', 'is_patho', 'is_IgA',
               'is_probio', 'is_flagellum', 'diamond_mmseqs_intersec_toxin', 'is_topgraph_new_&_old', 'signalp6_slow']]
        test_set_labels = test_set['dataset']
        label_encoder2 = LabelEncoder()
        label_encoder2.fit(test_set['dataset'])
        test_set_labels = label_encoder2.transform(test_set_labels)

        # Build Random forest
        clf = OneVsRestClassifier(RandomForestClassifier(n_estimators=1000, criterion='gini'))
        clf.fit(training_set_x, training_set_labels)
        y_pred=clf.predict(test_set_x)
        y_score = clf.predict_proba(test_set_x)

        # Create confusion matrix
        #metrics.confusion_matrix(test_set_labels, y_pred)
        #print(metrics.accuracy_score(test_set_labels, y_pred)) # Display accuracy score
        #print(metrics.f1_score(test_set_labels, y_pred, average='micro')) # Display F1 score

        #print(f"Accuracy: {metrics.accuracy_score(test_set_labels, y_pred)}\nF1 score: {metrics.f1_score(test_set_labels, y_pred, average='micro')}")

        test_labels_binarized = label_binarize(test_set_labels, classes=np.unique(test_set_labels))

        fpr = {}
        tpr = {}
        tprs = {}
        roc_auc = {}
        dict_importances = {}

        for c in range(n_classes):
            # Make dictionary of importances for each class
            dict_importances[c] = clf.estimators_[c].feature_importances_
            # Make 3 dictionaries for each class representing fpr tpr and _
            fpr[c], tpr[c], _ = roc_curve(test_labels_binarized[:,c], y_score[:,c])
            tprs[c] = np.interp(mean_fpr, fpr[c], tpr[c])
            roc_auc[c] = auc(fpr[c], tpr[c])
        # Append dictionaries of every iter to list
        roc_auc_overall.append(roc_auc)
        tprs_list.append(tprs)
        feat_imp_list.append(dict_importances)
    #print(dict_importances)
    features = ['is_pos_cntrl', 'is_neg_cntrl', 'is_rand_cntrl', 'is_auto', 'is_infect',
               'is_EBV', 'is_toxin', 'is_EM', 'is_MPA', 'is_patho', 'is_IgA',
               'is_probio', 'is_flagellum', 'diamond_mmseqs_intersec_toxin', 'is_topgraph_new_&_old', 'signalp6_slow']
    # Need to get one list per class containing only its class tprs arrays
    list_class_1 = []
    list_class_2 = []
    list_class_3 = []
    roc_auc1 = []
    roc_auc2 = []
    roc_auc3 = []
    feat_imp1 = []
    feat_imp2 = []
    feat_imp3 = []
    for i in range(n_iter):
        list_class_1.append(tprs_list[i][0])
        list_class_2.append(tprs_list[i][1])
        list_class_3.append(tprs_list[i][2])
        roc_auc1.append(roc_auc_overall[i][0])
        roc_auc2.append(roc_auc_overall[i][1])
        roc_auc3.append(roc_auc_overall[i][2])
        feat_imp1.append(feat_imp_list[i][0])
        feat_imp2.append(feat_imp_list[i][1])
        feat_imp3.append(feat_imp_list[i][2])

    feat_df = pd.DataFrame(feat_imp1).T
    feat_df.index = features

    roc_auc_avg = {}   
    roc_auc_avg[0] = np.mean(roc_auc1)
    roc_auc_avg[1] = np.mean(roc_auc2)
    roc_auc_avg[2] = np.mean(roc_auc3)

    # Plot ROC curve on the average TP FN values and add std dev
    colors = color_curve
    title = 'ROC_3class_model_100boot_agilent_' + names[0].replace(' ', '') + '_' + names[1].replace(' ', '') + names[2].replace(' ', '') + '.png'
    if shared_roc:
        plt.figure(figsize=(8,6))
        for i, list_class in zip(range(n_classes), [list_class_1, list_class_2, list_class_3]):
            array_class = np.array(list_class)
            class_avg = np.average(array_class, axis=0)
            class_std = np.std(array_class, axis=0)

            tprs_upper = class_avg + class_std
            tprs_lower = class_avg - class_std
            plt.plot(mean_fpr, class_avg, colors[i],
                    label='%s vs Rest (AUC average 100 boot=%0.2f)'%(classes[i], roc_auc_avg[i]))
            plt.fill_between(mean_fpr, tprs_lower, tprs_upper, color=colors[i], alpha=0.1);

        plt.plot([0, 1], [0, 1],'--')
        plt.xlim([-0.01, 1.01])
        plt.ylim([-0.01, 1.01])
        plt.ylabel('True Positive Rate')
        plt.xlabel('False Positive Rate')
        plt.legend(loc='lower right')


        plt.savefig(title,
                    facecolor='white', 
                    transparent=False,
                   bbox_inches = 'tight')

        plt.close()

    # Also plot the 3 plots in 3 panels
    else:
        fig,axes = plt.subplots(nrows=3, figsize=(6,12))
        plt.rcParams.update({'font.size': 10})
        for i, list_class, ax in zip(range(n_classes), [list_class_1, list_class_2, list_class_3], axes):
            array_class = np.array(list_class)
            class_avg = np.average(array_class, axis=0)
            class_std = np.std(array_class, axis=0)

            tprs_upper = class_avg + class_std
            tprs_lower = class_avg - class_std
            ax.plot(mean_fpr, class_avg, colors[i],
                    label='%s vs Rest (AUC average 100 boot=%0.2f)'%(classes[i], roc_auc_avg[i]))
            ax.fill_between(mean_fpr, tprs_lower, tprs_upper, color=colors[i], alpha=0.1);
            ax.plot([0, 1], [0, 1],'--')
            ax.set_xlim([-0.01, 1.01])
            ax.set_ylim([-0.01, 1.01])
            ax.set_ylabel('True Positive Rate', size=10)
            ax.set_xlabel('False Positive Rate', size=10)
            #ax.set_xticklabels(mean_fpr, size=10)
            ax.legend(loc='lower right')
        plt.savefig(title,
                    dpi=300,
                    facecolor='white', 
                    transparent=False,
                   bbox_inches = 'tight')
    return feat_imp1, feat_imp2, feat_imp3


def rf_2class_model(data1, data2, names, train_size, test_size, features, color_curve='orange'):
    '''
    Runs 2 class model random forest with 100 bootstraps with random sampling of train and test sets.
    Returns feature importance matrix for further plot analysis.
    '''
    from sklearn.preprocessing import LabelEncoder
    from sklearn.multioutput import MultiOutputClassifier
    from sklearn.metrics import auc, roc_curve
    from sklearn.preprocessing import label_binarize
    import numpy as np
    # Define lists and arrays
    list_features = []
    n_iter = 100
    arr_pred = np.zeros(100)
    imp = []
    mean_fpr = np.linspace(0, 1, 100)
    roc_auc_overall = []
    tprs_list = []
    feat_imp_list = []
    id_features = ['peptide_name'] + features
    for i in range(n_iter):
        train_set1 = (data1[['peptide_name', 'is_pos_cntrl', 'is_neg_cntrl', 'is_rand_cntrl', 'is_auto', 'is_infect',
               'is_EBV', 'is_toxin', 'is_EM', 'is_MPA', 'is_patho', 'is_IgA',
               'is_probio', 'is_flagellum', 'diamond_mmseqs_intersec_toxin', 'is_topgraph_new_&_old', 'signalp6_slow']]).fillna(0).sample(train_size)*1
        train_set1['dataset'] = names[0]
        train_set2 = (data2[['peptide_name', 'is_pos_cntrl', 'is_neg_cntrl', 'is_rand_cntrl', 'is_auto', 'is_infect',
               'is_EBV', 'is_toxin', 'is_EM', 'is_MPA', 'is_patho', 'is_IgA',
               'is_probio', 'is_flagellum', 'diamond_mmseqs_intersec_toxin', 'is_topgraph_new_&_old', 'signalp6_slow']]).fillna(0).sample(train_size)*1
        train_set2['dataset'] = names[1]

        test_set1 = data1[data1['peptide_name'].isin(train_set1['peptide_name']) == False][id_features].fillna(0).sample(test_size)*1
        test_set1['dataset'] = names[0]


        test_set2 = data2[data2['peptide_name'].isin(train_set2['peptide_name']) == False][id_features].fillna(0).sample(test_size)*1
        test_set2['dataset'] = names[1]


        comparisons = names[0] + ' vs.\n' + names[1]

        # Using 2 labels 
        training_set = pd.concat([train_set1, train_set2])
        test_set = pd.concat([test_set1, test_set2])

        training_set_x = training_set[features]
        training_set_labels = pd.factorize(training_set['dataset'])[0]


        test_set_x = test_set [features]
        test_set_labels = pd.factorize(test_set['dataset'])[0]

        # Run random forest classifier
        clf = RandomForestClassifier(n_estimators=100)
        clf.fit(training_set_x, training_set_labels)
        test_set_labels
        y_pred = clf.predict(test_set_x)
        y_score = clf.predict_proba(test_set_x)[:,1]

        # Make dictionary of importances for each class
        dict_importances = clf.feature_importances_
        # Make 3 dictionaries for each class representing fpr tpr and _
        fpr, tpr, _ = roc_curve(test_set_labels,y_score)
        #print(mean_fpr, fpr, tpr)
        tprs = np.interp(mean_fpr, fpr, tpr)
        tprs[0] = 0.
        #print('Resulting:', tprs)
        roc_auc = auc(fpr, tpr)
        # Append dictionaries of every iter to list
        #print(roc_auc)
        roc_auc_overall.append(roc_auc)
        tprs_list.append(tprs)
        feat_imp_list.append(dict_importances)

    # Need to get one list per class containing only its class tprs arrays
    list_class_1 = []
    roc_auc1 = []
    feat_imp = []
    for i in range(n_iter):
        list_class_1.append(tprs_list[i])
        roc_auc1.append(roc_auc_overall[i])
        feat_imp.append(feat_imp_list[i])

    roc_auc_avg = np.mean(roc_auc1)

    array_class = np.array(list_class_1)
    class_avg = np.average(array_class, axis=0)
    class_std = np.std(array_class, axis=0)

    tprs_upper = class_avg + class_std
    tprs_lower = class_avg - class_std
    fig, ax = plt.subplots(1, figsize=(3,3))
    # Set axis formatting 
    # Increase the size of the x-axis ticks
    ax.tick_params(axis='both', which='major', labelsize=11)
    ax.plot(mean_fpr, class_avg, color=color_curve,
            label='%s\n(AUC average=%0.2f)'%(comparisons, roc_auc_avg))
    ax.fill_between(mean_fpr, tprs_lower, tprs_upper, color=color_curve, alpha=0.1);
    ax.plot([0, 1], [0, 1],'--')
    ax.set_xlim([-0.01, 1.01])
    ax.set_ylim([-0.01, 1.01])
    ax.set_ylabel('True Positive Rate', size=11)
    ax.set_xlabel('False Positive Rate', size=11)
    #ax.set_xticklabels(mean_fpr, size=10)
    ax.legend(loc='lower right', fontsize=13)
    return feat_imp



def featureImp_3class_boxplot(feat_imp1, feat_imp2, feat_imp3, classes, plot_title, grouped=False):
    feature_df1 = pd.DataFrame(feat_imp1)
    feature_df2 = pd.DataFrame(feat_imp2)
    feature_df3 = pd.DataFrame(feat_imp3)
    #print(feature_df1)
    title = classes    
    # Rename variables
    renamed = {'is_EM' : 'Microbiota genes',  'is_infect' : 'Infectious diseases (IEDB)', 
               'is_topgraph_new_&_old' : 'Membrane proteins',
               'is_toxin' : 'Virulence Factor DB', 'is_MPA' : 'Microbiota strains',
               'is_patho' : 'Gut Pathogens', 'signalp6_slow' : 'Secreted proteins',
               'is_auto' : 'Autoantigens', 'is_probio' : 'Probiotic strains',
               'is_IgA' : 'IgA coated strains', 'is_pos_cntrl' : 'Positive control',
               'is_EBV' : 'EBV', 'diamond_mmseqs_intersec_toxin' : 'Predicted toxins', 
               'is_flagellum' : 'Flagellum proteins (annotations)', 
               'is_neg_cntrl' : 'Negative control', 'is_rand_cntrl' : 'Random control'}
    
    # Define palette
    palette = {'Positive control': 'rosybrown', 
               'Negative control': 'blue', 
               'Random control': 'lightblue', 
               'Autoantigens': 'lightgreen', 
               'Infectious diseases (IEDB)': 'yellow', 
               'EBV': 'brown', 
               'Virulence Factor DB': 'pink', 
               'Microbiota genes': 'purple', 
               'Microbiota strains': 'teal', 
               'Gut Pathogens': 'darkgrey', 
               'IgA coated strains': 'red', 
               'Probiotic strains': 'khaki', 
               'Flagellum proteins (annotations)': 'chocolate', 
               'Predicted toxins': 'lawngreen', 
               'Membrane proteins': 'sienna', 
               'Secreted proteins': 'tab:orange'}
    
    renamed_vars = [renamed[x] for x in features]
    #fig, axes = plt.subplots(3, figsize=(18,12))
    for i,name in zip([feature_df1, feature_df2, feature_df3], title):
        fig, ax = plt.subplots(figsize=(16,4))
        i.columns = renamed_vars
        order=i.mean(axis=0).sort_values(ascending=False).index
        melt_features = pd.melt(i.loc[:,order])
        sns.boxplot(data=melt_features, x='variable', y='value', ax=ax, palette=palette);
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=15, ha="right")
        ax.set_xlabel('')
        ax.set_title('%s vs Rest'%(name))
        ax.set_ylabel('Feature importance\n(Mean decrease in impurity)')
        ax.set_xlabel('')
        #ax.grid()
        plot_title = plot_title + '_' + '%s vs Rest'%(name)
        plt.savefig(plot_title,
                     facecolor='white', 
                     transparent=False,
                     bbox_inches = 'tight')

# Plot feature importances using mean decrease in impurity
def grouped_featureImp_box(feat_imp1, feat_imp2, feat_imp3, df_agilent_final, df_tested_agilent_clean, agilent_total_clean2, features, comparisons):
    # Define renamed variables
    renamed = {'is_EM' : 'Microbiota genes',  'is_infect' : 'Infectious diseases (Immune Epitope Database)', 
               'is_topgraph_new_&_old' : 'Membrane proteins',
               'is_toxin' : 'Virulence Factor Database', 'is_MPA' : 'Microbiota strains',
               'is_patho' : 'Gut pathogens', 'signalp6_slow' : 'Secreted proteins',
               'is_auto' : 'Autoantigens', 'is_probio' : 'Probiotic strains',
               'is_IgA' : 'IgA coated strains', 'is_pos_cntrl' : 'Positive control',
               'is_EBV' : 'EBV', 'diamond_mmseqs_intersec_toxin' : 'Predicted toxins', 
               'is_flagellum' : 'Flagellum proteins', 
               'is_neg_cntrl' : 'Negative control', 'is_rand_cntrl' : 'Random control'}
    renamed_vars = [renamed[x] for x in features]
    
    # Concatenate dataframes
    feature_df = pd.concat([pd.DataFrame(feat_imp1), pd.DataFrame(feat_imp2), pd.DataFrame(feat_imp3)])

    # Prepare indexes:
    import itertools
    # feature_df.columns = renamed_vars
    feature_df.columns = features
    # This commands for sorting based on means:
    # order = feature_df.loc[:,renamed_vars].mean(axis=0).sort_values(ascending=False).index
    feature_df['comparison'] = list(itertools.repeat(comparisons[0], 100)) + \
     list(itertools.repeat(comparisons[1], 100)) +  list(itertools.repeat(comparisons[2], 100))

    
    ordered_means = feature_df.groupby('comparison').mean().loc[:,features].mean(axis=0).sort_values(ascending=False).index

    ordered_means = list(ordered_means)

    ordered_means = ordered_means + ['comparison']
    # Order values based on mean considering all 3 comparisons
    feature_df = feature_df.loc[:,ordered_means]

    # This command to sort based on HLA vs GWAS category
    # order = feature_df.groupby('comparison').mean().transpose().sort_values(by='HLA associated vs. GWAS tested', ascending=False).index
    # order = list(order) + ['comparison']
    # print(order)
    # feature_df = feature_df.loc[:,order]
    # feature_df = feature_df.groupby('comparison').sum().sort_values(by='HLA associated vs. GWAS tested', axis=1, ascending=False)
    # feature_df['comparison'] = list(itertools.repeat(comparisons[0], 100)) + \
    # list(itertools.repeat(comparisons[1], 100)) +  list(itertools.repeat(comparisons[2], 100))

    feature_df_melt = pd.melt(feature_df, id_vars='comparison')
    #print(feature_df_melt)
    enrichment_column = []
    for comparison, type_prot in feature_df_melt[['comparison', 'variable']].values:
        if comparison == 'HLA associated vs. GWAS tested':
            if len(df_agilent_final[df_agilent_final[type_prot] == 1])/len(df_agilent_final) > len(df_tested_agilent_clean[df_tested_agilent_clean[type_prot] == 1])/len(df_tested_agilent_clean):
                enrichment_column.append('enriched in HLA associated')
            else:
                enrichment_column.append('enriched in GWAS tested')
        if comparison == 'HLA associated vs. Entire library':
            if len(df_agilent_final[df_agilent_final[type_prot] == 1])/len(df_agilent_final) > len(agilent_total_clean2[agilent_total_clean2[type_prot] == 1])/len(agilent_total_clean2):
                enrichment_column.append('enriched in HLA associated')
            else:
                enrichment_column.append('enriched in Entire library')
        if comparison == 'GWAS tested vs. Entire library':
            if len(df_tested_agilent_clean[df_tested_agilent_clean[type_prot] == 1])/len(df_tested_agilent_clean) > len(agilent_total_clean2[agilent_total_clean2[type_prot] == 1])/len(agilent_total_clean2):
                enrichment_column.append('enriched in GWAS tested')
            else:
                enrichment_column.append('enriched in Entire library')
    print(feature_df_melt)
    feature_df_melt['enrichment'] = enrichment_column
    feature_df_melt['variable'].replace(renamed, inplace=True)
    # Plotting: some parameters
    fig, axes = plt.subplots(1,3, figsize=(14,8), sharex=True)
    for comp,ax in zip(comparisons, axes): 

        ax.tick_params(axis='both', which='major', labelsize=13)
        palette= {'enriched in HLA associated' : 'tab:red', 'enriched in GWAS tested' : 'lightskyblue',
                  'enriched in Entire library' : 'lightgrey'}
        # Set the colors of the boxplots:
        PROPS = {
    'boxprops':{'edgecolor':'black'},
    'medianprops':{'color':'black'},
    'whiskerprops':{'color':'black'},
    'capprops':{'color':'black'}
}
        flierprops = dict(markeredgecolor='black')

        # Plot all features together based on ordered means
        sns.boxplot(data=feature_df_melt[feature_df_melt['comparison'] == comp], x='value', y='variable', ax=ax, hue='enrichment', palette=palette,
                    dodge=False, flierprops=flierprops, **PROPS);
        #plt.setp(ax.xaxis.get_majorticklabels(), rotation=35, ha="right") #, fontsize=15)
        ax.set_ylabel(ax.get_ylabel(), fontsize=13)
        ax.get_legend().remove()
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_xlabel('Feature Importance\n(Mean decrease in impurity)', fontsize=13)
        # adding transparency to colors
        # for patch in ax.artists:
        #     r, g, b, a = patch.get_facecolor()
        #     patch.set_facecolor((r, g, b, .4))
        # # Add trasparency to legend
        # leg = plt.legend()
        # for lh in leg.legendHandles: 
        #     lh.set_alpha(.4)
    axes[1].set_yticks([])
    axes[2].set_yticks([])
    return feature_df
    #return palette
    #plt.grid()


    
def single_featureImp_boxplot(feat_imp, comparison):
    #### Feature importance
    feature_df = pd.DataFrame(feat_imp)
    title = comparison
    # Rename variables
    renamed = {'is_EM' : 'Microbiota genes',  'is_infect' : 'Infectious pathogens', 
               'is_topgraph_new_&_old' : 'Membrane protein',
               'is_toxin' : 'Toxin', 'is_MPA' : 'Microbiota strains',
               'is_patho' : 'Pathogens', 'signalp6_slow' : 'Secreted proteins',
               'is_auto' : 'Human Autoantigens', 'is_probio' : 'Probiotic strains',
               'is_IgA' : 'IgA coated strains', 'is_pos_cntrl' : 'Positive control',
               'is_EBV' : 'EBV', 'diamond_mmseqs_intersec_toxin' : 'Predicted toxin (diamond)', 
               'is_flagellum' : 'Predicted flagellum protein (mapped)', 
               'is_neg_cntrl' : 'Negative control', 'is_rand_cntrl' : 'Random control'}
    
    # Define palette
    palette = {'Positive control': 'rosybrown', 
               'Negative control': 'blue', 
               'Random control': 'lightblue', 
               'Human Autoantigens': 'lightgreen', 
               'Infectious pathogens': 'yellow', 
               'EBV': 'brown', 
               'Toxin': 'pink', 
               'Microbiota genes': 'purple', 
               'Microbiota strains': 'teal', 
               'Pathogens': 'darkgrey', 
               'IgA coated strains': 'red', 
               'Probiotic strains': 'khaki', 
               'Predicted flagellum protein (mapped)': 'chocolate', 
               'Predicted toxin (diamond)': 'lawngreen', 
               'Membrane protein': 'sienna', 
               'Secreted proteins': 'tab:orange'}
    
    renamed_vars = [renamed[x] for x in features]
    # First melt dataframe and order according to mean values for each feature
    fig, ax = plt.subplots(1, figsize=(8,5))
    feature_df.columns = renamed_vars
    order = feature_df.mean(axis=0).sort_values(ascending=False).index
    melt_features = pd.melt(feature_df.loc[:,order])
    
    # Plot boxplots
    sns.boxplot(data = melt_features, x='variable', y='value', ax=ax, palette=palette);
    plt.setp(ax.xaxis.get_majorticklabels(), fontsize=14, rotation=35, ha="right")
    ax.set_xlabel('')
    ax.set_title(title)
    #ax.set_yticklabels()
    plt.grid()  
