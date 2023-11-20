#!/bin/python3.9

import matplotlib.pyplot as plt
from statsmodels.stats import multitest
from scipy.stats import fisher_exact
import pandas as pd
import numpy as np
from collections import defaultdict
import seaborn as sns
from statsmodels.stats.multitest import multipletests
from scipy.stats import fisher_exact
from scipy.stats import chi2_contingency

def multitest_BH(dictio, list_of_categories, list_of_datasets, names):
    from statsmodels.stats import multitest
    from scipy.stats import fisher_exact
    if len(list_of_datasets) > 2:
        # Reshape the pvalues to use the multitest correction
        pvalues_tot = np.concatenate(tuple(dictio.values()), axis=0)
        pvalues_tot = pvalues_tot.reshape(len(list_of_categories), len(names))
        # Compute correction
        corrected = []
        for n in range(len(list_of_datasets)):
            #print(pvalues_tot_func[:,n])
            is_sign,pcorr = multitest.fdrcorrection(pvalues_tot[:,n])
            #print(pcorr)
            corrected.append(pcorr)

        correctedpvalues = np.array(tuple(corrected)).transpose()
        # generate a dictionary
        corr_dict = {}
        for cat,corr in zip(list_of_categories, correctedpvalues):
            corr_dict[cat] = corr
        print('Correction DONE')
        return corr_dict
    else:
        pvalues = tuple([x[0] for x in dictio.values()])
        print(pvalues)
        is_sign,pcorr = multitest.fdrcorrection(pvalues)        
        correctedpvalues = pcorr
        # generate a dictionary
        corr_dict = {}
        for cat,corr in zip(list_of_categories, correctedpvalues):
            corr_dict[cat] = [corr]
        return corr_dict
        
def compute_pvalues(dict_counts, classes, renamed): # first have to be the total dataset, second the HLA third the non HLA
    #print(relevant)
    if len(classes) == 3:
        comb = [(0,1), (0,2), (1,2)]
    elif len(classes) == 2:
        comb = [(0,1)]
    
    pvalues_dict = {}
    tot_values = [len(x) for x in classes]
    #print(tot_values)
    for i in range(len(renamed)):
        #print(dict_counts)
        values = [v[i] for k, v in dict_counts.items()]
        #print(values)
        pvalues = []
        for n in comb:
            cont1 = [(values[n[0]], tot_values[n[0]] - values[n[0]]), (values[n[1]], tot_values[n[1]] - values[n[1]])]
            #print(cont1)
            odds, pvalue = fisher_exact(cont1)
            pvalues.append(pvalue)
        pvalues_dict[renamed[i]] = pvalues
    #print(pvalues_dict)
    return pvalues_dict

def annot_stat(star, x1, x2, y, h, col='k', ax=None):
            ax = plt.gca() if ax is None else ax 
            ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
            ax.text((x1+x2)*.5, y+h, star, ha='center', va='bottom', color=col, fontsize=16)

def create_relative_abund_plots(list_of_subsets, 
                                relevant_categories, 
                                renamed, 
                                multitest_correction = False, 
                                classes = False,
                               color = None,
                               log = False,
                               count_on_bars = True,
                               custom_size=None): # can be agilent sets (already processed) or twist
    # Define classes
    if classes:
        pass
    else:
        if len(list_of_subsets) == 3:
            classes = ['input', 'HLA', 'non HLA']
        elif len(list_of_subsets) == 2:
            classes = ['input', 'associations']
        
    # Initalize count and percentage dictionary
    dict_count = {} # Dict in which every key is the library name (or subset) and value is list of counts for every category
    dict_perc = {} # Dict in which every key is the library name (or subset) and value is list of percentages
    
    # Iterate over all the classes and then (nested) iterate on each feature
    for names_lib,lib in zip(classes,list_of_subsets):
        list_perc = []
        list_counts = []
        # if functional categories will be 4
        for cat in relevant_categories:
            try:
                perc = (sum(lib[lib[cat].isna()==False][cat].values) / len(lib)) * 100
            except:
                perc = 0
            count = sum(lib[lib[cat].isna()==False][cat].values)
            list_perc.append(perc)
            list_counts.append(count)
        dict_count[names_lib] = (list_counts)
        dict_perc[names_lib] = (list_perc)
    # Compute pvalues
    pvalues = compute_pvalues(dict_count, list_of_subsets, renamed)
    if multitest_correction:
        pvalues = multitest_BH(pvalues, renamed, list_of_subsets, classes)
        #print(corr_pvalues)
    print(pvalues)
    # Create a dataframe for plotting percentages
    perc_df = pd.DataFrame(dict_perc)
    perc_df['categories'] = relevant_categories
    melted_perc = pd.melt(perc_df, id_vars='categories')

    # Plot
    if not custom_size:
        if len(relevant_categories) == 10:
            fig, axes= plt.subplots(ncols= 5, nrows=2, figsize=(28,12))
        elif len(relevant_categories) == 12:
            fig, axes= plt.subplots(ncols= 4, nrows=3, figsize=(24,21))# figsize=(32,12))   
        elif len(relevant_categories) < 10:
            fig, axes = plt.subplots(ncols= 4, nrows=3, figsize=(24,21))
    else:
        fig, axes = plt.subplots(ncols= 2, nrows=2, figsize=custom_size)
    axes = axes.flatten()
    plots = []
    for i,ax,title in zip(relevant_categories, axes, renamed):
        plot = sns.barplot(data=melted_perc[melted_perc['categories'] == i], 
                    x='variable', 
                    y= 'value',
                    edgecolor='black',
                    ax=ax,
                    color=color)
        plots.append(plot)
        ax.set_title(title, fontweight='bold')
        # Take ticks and change orientation
        # labels = ax.get_xticklabels
        ax.set_ylabel('Relative abundance (%)', fontsize=12)
        ax.set_xlabel('')
        if log==True:
            ax.set_yscale('log')

    # Create annotations and plot number of peptides on each bar
    dict_count2 = defaultdict(list)
    for i in dict_count.values():
        for n1,n2 in zip(i, relevant_categories):
            dict_count2[n2].append(n1) # associating each pvalue
            
    list_heights = []
    for ax, values in zip(plots, dict_count2.values()):
        count = 0
        maximum = 0
        #print(ax)
        #print(f"ax patches are: {[x for x in ax.patches]}")
        for rect in ax.patches:
            #print(rect)
            if np.isnan(rect.get_height()) == False:
                #print(count)
                #print(rect.get_height(), values[count], sep='\t')
                if count_on_bars:
                    ax.annotate(f"{values[count]:,}", (rect.get_x() + rect.get_width() / 2., rect.get_height()),
                     ha='center', va='center', fontsize=20, color='black', xytext=(0, 10),
                     textcoords='offset points')
                count += 1
                if rect.get_height() > maximum:
                    maximum = rect.get_height()
                    maximum = maximum 
        list_heights.append(maximum)
        if not count_on_bars:
            newlabs = []
            for label, count in zip(classes, values):
                newlabel = '\n'.join(label.split(' ')) + '\n' + str(count)
                newlabs.append(newlabel)
            ax.set_xticklabels(newlabs, fontsize=14)
     
    
    for ax, maximum in zip(plots,list_heights):
        pvalue = pvalues[ax.get_title()]
        c = 0
        if maximum <= 2:
            levels = (int(maximum)+0.5, int(maximum)+1, int(maximum)+1.5)
        else:
            levels = (int(maximum)+2, int(maximum)+5, int(maximum)+8)
        for i, comparison, level in zip(pvalue,[(0,1),(0,2),(1,2)], levels):
            #print(level)
            if maximum <= 2:
                wide = 0.1
            elif maximum <= 1:
                wide = 0.1
            elif maximum > 2 and maximum < 50:
                wide = 1
            elif maximum >= 50:
                wide = 2
            if i == 'NaN':
                pass
            else:
                if i < 0.05 and i >= 0.01:
                    star = '*'
                if i < 0.01 and i >= 0.001:
                    star = '**'
                if i < 0.001:
                    star = '***'
                if i >= 0.05:
                    star = 'ns'
                annot_stat(star, comparison[0], comparison[1], level, wide, ax=ax)
        

    
    return pvalues, dict_count, perc_df

def clean_df_total(df_total, df_1, df_2):
    df_total_clean = df_total[~df_total['peptide_name'].isin(df_1['peptide_name'].tolist()+df_2['peptide_name'].tolist())]
    
    return df_total_clean


def plot_microorganism_counts(df, top='all', figsize=None, barv=False, ascending=True):
    if barv:
        if figsize:
            df[df.Organism_name.str.contains(' ')]\
            .groupby('Organism_name')['Organism_name']\
            .count()\
            .sort_values(ascending=ascending)\
            .plot\
            .bar(figsize=figsize)
        else:
            df.groupby('Organism_name')['Organism_name']\
            .count()\
            .sort_values(ascending=ascending)\
            .plot\
            .barv()
    else:
        if figsize:
            if top == 'all':
                df[df.Organism_name.str.contains(' ')]\
                .groupby('Organism_name')['Organism_name']\
                .count()\
                .sort_values(ascending=ascending)\
                .plot\
                .barh(figsize=figsize)
            else:
                df[df.Organism_name.str.contains(' ')]\
                .groupby('Organism_name')['Organism_name']\
                .count()\
                .sort_values(ascending=ascending)[:top]\
                .plot\
                .barh(figsize=figsize)
        else:
            df.groupby('Organism_name')['Organism_name']\
            .count()\
            .sort_values(ascending=ascending)\
            .plot\
            .barh()
    return


def plot_peptide_dist(df, relevant_categories, newnames,color=None):
    df_selected = df.loc[:,relevant_categories]
    fig, ax = plt.subplots(figsize=(10,10))
    values_df = [df_selected[x].sum() for x in relevant_categories]

    values_df = pd.DataFrame({'Category': newnames, 'number of peptides' : values_df})
    values_df = values_df.sort_values(by='number of peptides')
    width=.75
    if color == None:
        bars = ax.barh(values_df.Category, values_df['number of peptides'],
                width, color='lightblue')
    else:
        bars = ax.barh(values_df.Category, values_df['number of peptides'],
                width, color=color)
    for bars in ax.containers:
        ax.bar_label(bars)
    ax.set_xlim(1,800000)
    ax.set_xscale('log')
    ax.set_xlabel('Number of peptides (log)')
    return


def add_organism_column(df_agilent):
    import ast
    from collections import namedtuple
    #Creare funzione per aggiungere microrganismi tutti su una colonna
    index_IEDB_NA_uniref = df_agilent[(df_agilent['IEDB_organism_name'].isna()) &
                                      (df_agilent['uniref_func'].isna() == False) &
                                     (df_agilent['uniref_func'].str.contains('Tax'))].index
    list_names_uniref = df_agilent[(df_agilent['IEDB_organism_name'].isna()) &
                                        (df_agilent['uniref_func'].isna() == False) &
                           (df_agilent['uniref_func'].str.contains('Tax'))]['uniref_func'].apply(lambda x: x.split('Tax=')[1].split(' TaxID')[0]).tolist()
    #print(len(index_NA_IEDB_uniref), len(list_names))
    index_IEDB = df_agilent[(df_agilent['IEDB_organism_name'].isna() == False)].index
    dict_organisms_IEDB = () # tuple contains 2 values, the first is the first organism entry, the second the whole list
    list_names_IEDB = df_agilent[df_agilent['IEDB_organism_name'].isna() == False]['IEDB_organism_name'].apply(
            lambda x: set(ast.literal_eval(x.replace('] & [', ',').replace(';',',')))).tolist()
    list_whole = []
    for x in list_names_IEDB:
        dict_organisms_IEDB = (list(x)[0], [i for i in list(x)])
        list_whole.append(dict_organisms_IEDB)
    print(len(list_whole), len(index_IEDB))

    df_agilent['Organism_name'] = np.nan
    df_agilent['Organism_name'].iloc[index_IEDB] = [x[0] for x in list_whole]
    df_agilent['Organism_name'].iloc[index_IEDB_NA_uniref] = list_names_uniref
    df_agilent['Organism_name'] = df_agilent['Organism_name'].fillna('Not Known')
    return df_agilent, list_whole


def create_relative_abund_plots_single(list_of_subsets, 
                                relevant_categories, 
                                renamed, 
                                multitest_correction = False, 
                                classes = False,
                               color = None,
                               log = False,
                               count_on_bars = True,
                               custom_size=None): # can be agilent sets (already processed) or twist
    # Define classes
    if classes:
        pass
    else:
        if len(list_of_subsets) == 3:
            classes = ['input', 'HLA', 'non HLA']
        elif len(list_of_subsets) == 2:
            classes = ['input', 'associations']
        
    # Initalize count and percentage dictionary
    dict_count = {} # Dict in which every key is the library name (or subset) and value is list of counts for every category
    dict_perc = {} # Dict in which every key is the library name (or subset) and value is list of percentages
    
    # Iterate over all the classes and then (nested) iterate on each feature
    for names_lib,lib in zip(classes,list_of_subsets):
        list_perc = []
        list_counts = []
        # if functional categories will be 4
        for cat in relevant_categories:
            try:
                perc = (sum(lib[lib[cat].isna()==False][cat].values) / len(lib)) * 100
            except:
                perc = 0
            count = sum(lib[lib[cat].isna()==False][cat].values)
            list_perc.append(perc)
            list_counts.append(count)
        dict_count[names_lib] = (list_counts)
        dict_perc[names_lib] = (list_perc)
    # Compute pvalues
    pvalues = compute_pvalues(dict_count, list_of_subsets, renamed)
    if multitest_correction:
        pvalues = multitest_BH(pvalues, renamed, list_of_subsets, classes)
        #print(corr_pvalues)
    print(pvalues)
    # Create a dataframe for plotting percentages
    perc_df = pd.DataFrame(dict_perc)
    perc_df['categories'] = relevant_categories
    melted_perc = pd.melt(perc_df, id_vars='categories')


    dict_count2 = defaultdict(list)
    for i in dict_count.values():
        for n1,n2 in zip(i, relevant_categories):
            dict_count2[n2].append(n1)
    print(dict_count2)
    # Plot
    
    #plots = []
    for i,title in zip(relevant_categories, renamed):
        fig, ax = plt.subplots(figsize=(4,4))
        plot = sns.barplot(data=melted_perc[melted_perc['categories'] == i], 
                    x='variable', 
                    y= 'value',
                    edgecolor='black',
                    palette=color)
        #plots.append(plot)
        #ax.set_title(title, fontweight='bold')
        # Take ticks and change orientation
        # labels = ax.get_xticklabels
        ax.set_ylabel('Relative abundance (%)', fontsize=12)
        ax.set_xlabel('')
        ax.set_xticklabels(classes, fontsize=12)
        if log==True:
            ax.set_yscale('log')
        #print(dict_count)
        # Create annotations and plot number of peptides on each bar
        # associating each pvalue
        #print(dict_count2)     
        list_heights = []
        #for ax, values in zip(plots, dict_count2.values()):
        count = 0
        maximum = 0
        #print(ax)
        #print(f"ax patches are: {[x for x in ax.patches]}")
        for rect in ax.patches:
            #print(rect)
            if np.isnan(rect.get_height()) == False:
                #print(count)
                #print(rect.get_height(), values[count], sep='\t')
                if count_on_bars:
                    ax.annotate(f"{dict_count2[i]:,}", (rect.get_x() + rect.get_width() / 2., rect.get_height()),
                        ha='center', va='center', fontsize=20, color='black', xytext=(0, 10),
                        textcoords='offset points')
                count += 1
                if rect.get_height() > maximum:
                    maximum = rect.get_height()
                    maximum = maximum 
        list_heights.append(maximum)
        if not count_on_bars:
            newlabs = []
            
            for label, count in zip(classes, dict_count2[i]):
                #print(count)
                newlabel = '\n'.join(label.split(' ')) + '\n' + f"n={count:,}" # str(count)
                newlabs.append(newlabel)
            
            ax.set_xticklabels(newlabs, fontsize=14)
        pvalue = pvalues[title]
        c = 0
        if maximum <= 2:
            levels = (int(maximum)+0.5, int(maximum)+1, int(maximum)+1.5)
        else:
            if maximum > 2 and maximum <= 20:
                levels = (int(maximum)+2, int(maximum)+5, int(maximum)+8)
            else:
                levels = (int(maximum)+5, int(maximum)+12, int(maximum)+20)
        for i, comparison, level in zip(pvalue,[(0,1),(0,2),(1,2)], levels):
            #print(level)
            if maximum <= 2:
                wide = 0.1
            elif maximum <= 1:
                wide = 0.1
            elif maximum > 2:
                wide = 0.9
            if i == 'NaN':
                pass
            else:
                if i < 0.05 and i >= 0.01:
                    star = '*'
                if i < 0.01 and i >= 0.001:
                    star = '**'
                if i < 0.001:
                    star = '***'
                if i >= 0.05:
                    star = 'ns'
                annot_stat(star, comparison[0], comparison[1], level, wide, ax=ax)
        #print(ax.get_ylim())
        if maximum <=2:
            ax.set_ylim(0,3)
        else:
            ax.set_ylim(0, ax.get_ylim()[1] + 2)
        savefile = title + '.png'
        #plt.grid(color='gainsboro', linestyle='-', linewidth=1)
        plt.savefig(savefile,
            facecolor='white',
            dpi=300,
            transparent=False,
           bbox_inches = 'tight')



def parse_and_plot_species(whole_df, groups, sources=True, thr=2):
    from matplotlib.cm import get_cmap
    # List species and subsets
    list_species = []
    list_subsets = []
    # Use a colormap and map to it all the bacteria
    # Choose a colormap
    #fig, axes = plt.subplots(5, figsize=(4,20))
    for i,title in zip(groups, ['ImmunoEpitope DB (infectious diseases)', 'EBV', 'Virulence Factor DB', 'Microbiota genes', 'Microbiota strains', 'Gut pathogens']):
        #print(i)
        # First we filter by removing root bacterium, then we count the species and we change the names of 
        # species with prevalence <10% to 'others'
        if sources:
            whole_df_filt = whole_df[(whole_df['reformatted_name'] != 'root bacterium') & (whole_df['class'] == i) &\
            (whole_df['reformatted_name'] != 'Bacteria bacterium') ]
        else:
            # In case we use predictions instead of sources
            whole_df_filt = whole_df[(whole_df['reformatted_name'] != 'root bacterium') & (whole_df['prediction'] == i) &\
            (whole_df['reformatted_name'] != 'Bacteria bacterium') ]

        #human_IEBD = whole_df_filt[whole_df_filt['reformatted_name'] == 'Homo']['peptide_name'].tolist()
        count_species = whole_df_filt.groupby(['reformatted_name']).count().iloc[:,0]

        # Check the total number of the counts
        total_count = count_species.sum()
        
        # # Calculate the series of percentages
        percentage_series = (count_species / total_count) * 100
        #percentage_series.plot(kind='bar')

        # # Change to others every species with prevalence < 2%
        # Can change the threshold from the function
        idx_to_change = list(percentage_series.index[percentage_series < thr])
        print(f'The values with counts < {thr} % will be converted to Others, and they are {len(list(idx_to_change))} taxa')
        
        # Create dictionary for changing values
        dict_for_change = {}
        for x,y in zip(idx_to_change, np.repeat('others',len(idx_to_change))):
            dict_for_change[x] = y
        whole_df_filt['reformatted_name'].replace(dict_for_change, inplace=True)
        list_subsets.append(whole_df_filt)
        list_species.append(set(whole_df_filt['reformatted_name']))
    whole_df_filt_concat = pd.concat(list_subsets)
    # Obtain colormap
    list_species_unique = set().union(*list_species)
    cmap = get_cmap('gist_ncar')
    print(cmap)
    # try Paired
    #cmap = get_cmap('Paired')
    # try define your own colors
    # List of distinct RGB colors
    
    colors = [cmap(i / len(list_species_unique)) for i in range(len(list_species_unique))]
    # Create a dictionary to map categories to colors
    category_color_map = dict(zip(sorted(list_species_unique), colors))
    list_subsets = plot_species(whole_df_filt_concat, list_subsets, category_color_map, groups, titlelist=['ImmunoEpitope DB (infectious diseases)',\
     'EBV',\
         'Virulence Factor DB',\
            'Microbiota genes',\
                 'Microbiota strains',\
                     'Gut pathogens'])
    return list_subsets, category_color_map

def plot_species(subsets_concat, list_subsets, colormap, groups, titlelist=None):
    from matplotlib.colors import ListedColormap
    list_wide= []
    # UNCOMMENT FROM HERE
    # Generating a custom color palette using matplotlib's 'tab20' colormap
    # You can choose any other colormap from matplotlib or create your own custom list of colors.
    list_plots = []
    #fig,axes= plt.subplots(5,figsize=(4,22))
    for whole_df_filt,title,i in zip(list_subsets, titlelist, groups):
        cimap = ListedColormap([colormap[c] for c in set(whole_df_filt['reformatted_name'])])
        
        #num_colors = 20
        #base_colors = plt.cm.tab20(np.linspace(0, 1, num_colors))
        # Creating a custom colormap using ListedColormap
        # Now recalculate again the counts with the new names
        subset_wide = whole_df_filt.groupby(['group', 'reformatted_name']).count().iloc[:,0].unstack().fillna(0).astype(int)
        idx = ['Entire lib', 'GWAS tested', 'HLA']
        # Plot using percentages
        subset_wide = subset_wide.reindex(idx)
        rowsums = subset_wide.sum(axis=1)
        #subset_wide.index[0] = 'whole'
        subset_wide['Updated_index'] = subset_wide.index + '\nn=' + [f'{x:,}' for x in rowsums]
        subset_wide = subset_wide.set_index('Updated_index')
        
        row_sum = subset_wide.sum(axis=1)
        subset_wide_perc = (subset_wide.div(row_sum, axis=0) * 100).round(3)
        cimap = ListedColormap([colormap[c] for c in subset_wide_perc.columns])
        #print({c:colormap[c] for c in set(whole_df_filt['reformatted_name'])})
        # New index order
        # Name each file:
        filename = 'Species_differences_stckd_bar' + i + '.png'
        fig, ax = plt.subplots(1, figsize=(4,3))
        subset_wide_perc.plot(kind='bar', stacked=True,\
        ax=ax, edgecolor='black', cmap=cimap, label=title, alpha=0.5)
        #ax.legend(None)
        ax.set_xlabel('')
        ax.set_xticklabels(ax.get_xticklabels(), rotation=0, fontsize=12)
        ax.set_yticklabels(range(0,100,20), fontsize=12)
        #ax.set_title(title)
        ax.set_ylim(0,100)
        ax.get_legend().remove()
        #if i == 'is_patho':
        ## PLOT
        plt.savefig(filename,
        facecolor='white',
        dpi=300, 
        bbox_inches = 'tight')
        list_wide.append(subset_wide)
        list_plots.append(ax)
    
    return list_wide



def barplot_annotate_brackets(num1, num2, data, center, height, yerr=None, dh=.05, barh=.05, fs=None, maxasterix=None):
    """ 
    Annotate barplot with p-values.

    :param num1: number of left bar to put bracket over
    :param num2: number of right bar to put bracket over
    :param data: string to write or number for generating asterixes
    :param center: centers of all bars (like plt.bar() input)
    :param height: heights of all bars (like plt.bar() input)
    :param yerr: yerrs of all bars (like plt.bar() input)
    :param dh: height offset over bar / bar + yerr in axes coordinates (0 to 1)
    :param barh: bar height in axes coordinates (0 to 1)
    :param fs: font size
    :param maxasterix: maximum number of asterixes to write (for very small p-values)
    """

    if type(data) is str:
        text = data
    else:
        # * is p < 0.05
        # ** is p < 0.005
        # *** is p < 0.0005
        # etc.
        text = ''
        p = .05

        while data < p:
            text += '*'
            p /= 10.

            if maxasterix and len(text) == maxasterix:
                break

        if len(text) == 0:
            text = 'n. s.'

    lx, ly = center[num1], height[num1]
    rx, ry = center[num2], height[num2]

    if yerr:
        ly += yerr[num1]
        ry += yerr[num2]

    ax_y0, ax_y1 = plt.gca().get_ylim()
    dh *= (ax_y1 - ax_y0)
    barh *= (ax_y1 - ax_y0)

    y = max(ly, ry) + dh

    barx = [lx, lx, rx, rx]
    bary = [y, y+barh, y+barh, y]
    mid = ((lx+rx)/2, y+barh)

    plt.plot(barx, bary, c='black')

    kwargs = dict(ha='center', va='bottom')
    if fs is not None:
        kwargs['fontsize'] = fs

    plt.text(*mid, text, **kwargs)


def chi_squared_test_groups_alleles(df_merged_assoc_gene, groups, by='allele', save=False, statistic=True):
    '''This functions calculates statistics (Chi2) for the table that is given as input:
    df_merged_assoc_gene is the associations table which needs 2 mandatory columns (allele and peptide).
    by argument is set to calculate the statistics for the alleles by default, but can be used for categories
    setting by==\'peptide\'. The table can be saved using 'Save' '''
    if by=='allele':
        group_by = ['peptide', 'allele']
        variable = 'allele'
    if by=='peptide':
        group_by = ['allele', 'peptide']
        variable = 'peptide'
    group_list = []
    for i in groups:
        if by=='allele':
            group_list.append(df_merged_assoc_gene[df_merged_assoc_gene[i] == 1].groupby(group_by).count().iloc[:,0].unstack().fillna(0).sum(axis=0))
        elif by=='peptide':
            group_list.append(df_merged_assoc_gene[df_merged_assoc_gene[i] == 1].groupby(group_by).count().iloc[:,0].unstack().fillna(0).sum(axis=1))
    grouped_concat = pd.concat(group_list, axis=1).fillna(0)
    grouped_concat.columns=groups
    #return grouped_concat
    row_sums = grouped_concat.sum(axis=1)
    grouped_concat = grouped_concat.loc[:, (grouped_concat.sum(axis=0) != 0)]
    #print(grouped_concat)
    #print(grouped_concat)
    #print(row_sums)
    grouped_concat_perc = (grouped_concat.div(row_sums, axis=0) * 100).round(3)
    grouped_concat_perc = grouped_concat_perc.loc[:, (grouped_concat_perc.sum(axis=0) != 0)]
    grouped_concat_perc.index = [x+'\n(n='+str(int(y))+')' for x,y in zip(grouped_concat_perc.index,
                                                    row_sums)]
    variable = 'allele'
    if by=='peptide':
        variable = 'category'
        #grouped_concat_perc = grouped_concat_perc.transpose()
        grouped_concat = grouped_concat.transpose()
        grouped_concat['total'] = grouped_concat.sum(axis=0)
        row_sums = grouped_concat.sum(axis=1)
        
        grouped_concat.index.name = 'category'
        grouped_concat = grouped_concat.sort_values(by=['category', 'total'], ascending=[True, False])
        grouped_concat = grouped_concat.drop('total', axis=1)
        grouped_concat_perc = (grouped_concat.div(row_sums, axis=0) * 100).round(3)
        grouped_concat_perc.index = [x+'\n(n='+str(int(y))+')' for x,y in zip(grouped_concat_perc.index,
                                                    row_sums)]

    #print(grouped_concat_perc)
    grouped_concat2 = grouped_concat.reset_index()
    #print(grouped_concat2)
    if statistic:
        long_contingency_count = pd.melt(grouped_concat2, id_vars=[variable])
        pvalues = {}
        alpha = 0.05
        # Compare each Category against the rest in the whole dataset of associations
        for e,i in enumerate(set(grouped_concat2[variable])):
            val1 = tuple(long_contingency_count[long_contingency_count[variable] == i]['value'])
            if by == 'allele':
                val_rest = tuple(long_contingency_count[long_contingency_count[variable] != i].groupby('variable').sum(numeric_only=True)['value'])
            if by == 'peptide':
                val_rest = tuple(long_contingency_count[long_contingency_count[variable] != i].groupby('allele').sum(numeric_only=True)['value'])
                #print(val_rest)
            cntg_table = np.array([val1, val_rest])
            #print(cntg_table)
            # Do chi2 test on each category against the expected values
            try:
                chi2, p_value, dof, expected = chi2_contingency(cntg_table)
                pvalues[i] = p_value
            except:
                print(f'Internally computed table present some unexpected values (e.g. frequency of expected too low),\n \
                    skipping Chi2 test for {i}')
        #print(pvalues)
        #print(list(pvalues.values()))

        # Bonferroni correction
        bonferroni_corr_pvalues = multipletests(list(pvalues.values()), method='bonferroni')[1]


        # Benjamini-Hochberg procedure (FDR correction)
        bh_corr_pvalues = multipletests(list(pvalues.values()), method='fdr_bh')[1]

        #print("Bonferroni corrected p-values:", bonferroni_corr_pvalues)
        #print("Benjamini-Hochberg corrected p-values:", bh_corr_pvalues)

        Chi_squared_df = pd.DataFrame({'Category' : [c for c in pvalues.keys()],
                                        'pvalues' : [c for c in pvalues.values()],
                                    'Bonferroni_corr' : bonferroni_corr_pvalues,
                                    'BH corrected' : bh_corr_pvalues}).sort_values(by='Bonferroni_corr')
        #print(Chi_squared_df)
        Chi_squared_df#.to_csv('Analysis_flagellum_peptides_CFS/Chi_test_single_alleles.csv',
        # sep='\t')
        #print(Chi_squared_df)
        return grouped_concat, grouped_concat_perc, Chi_squared_df
    else:
        return grouped_concat, grouped_concat_perc

# Define function for parsing dataframe
def plot_upset(df, filename, save=False):
    import upsetplot
    dict_peptides = {}
    cols = df.columns
    for i in cols:
        #print(i)
        #print(df_twist_final[i].tolist())
        list_peptides = df[df[i] == 1]['peptide'].tolist()
        #print(list_peptides)
        dict_peptides[i] = set(list_peptides)

    plt.figure(figsize=(10,10));
    upset_data_sub = upsetplot.from_contents({k: v for k, v in dict_peptides.items() if len(v) != 0})
    upsetplot.plot(upset_data_sub)
    if save:
        plt.savefig(filename)
    return