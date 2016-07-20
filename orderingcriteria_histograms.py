# -*- coding: utf-8 -*-
"""
Created on Mon Jun 20 16:36:09 2016

@author: Isavannah Reyes
"""

import pandas
import seaborn as sb
import matplotlib.pyplot as plt

 
#### SET data_native TYPES ####
#### convert all columns to numeric except for last description column
def set_numeric(dataframe):
    for term in return_list_of_numeric_columns(dataframe):
        dataframe[str(term)]=pandas.to_numeric(dataframe[str(term)])
        
###convert description columns to strings       
def set_description_string(dataframe):
    dataframe['description']=dataframe['description'].astype("str")

### define function that returns column names of numeric columns only
def return_list_of_numeric_columns(datascores):
    numeric_column_list=[column for column in datascores.columns if datascores[column].dtype== 'float64' or datascores[column].dtype== 'int64' ]
    #print(numeric_column_list)   
    return numeric_column_list
#return_list_of_numeric_columns(data_flex_bb_score)    
    
    
#makes histograms from fourscourfiles for each numeric column and outputs to a pdf    
def make_hists_from_four_scorefiles(dataframe1, dataframe2, dataframe3, dataframe4, outputname):#dataframe4,
    ### this is only true for Rosetta output score file, otherwise pay attention to included terms in list slice aka remove [1:-1]
    plt.figure(figsize = [9,150])
    for i, term in enumerate(return_list_of_numeric_columns(dataframe1)):
        ### subplot position 1. # rows, 2. # columns, 3. # of plots
        plt.subplot((len(dataframe1.columns)+1)/2, 2 , i+1)
        sb.distplot(dataframe1[term].dropna(), kde=True, norm_hist=True, bins=10, label="Flex")
        sb.distplot(dataframe2[term].dropna(), kde=True, norm_hist=True, bins=10, label="Fixed")
        #sb.distplot(dataframe3[term].dropna(), kde=True, norm_hist=True, bins=10, label="0.1")
        #sb.distplot(dataframe4[term].dropna(), kde=True, norm_hist=True, bins=10, label="0.001")
        #plt.xlabel(str(term))
        plt.ylabel("Frequency [%]")
        plt.title( str(term)+" Hist")
        plt.legend()
    plt.subplots_adjust(wspace=1,hspace=1)
    plt.tight_layout()
    plt.rc('font',size=3, )
    #plt.show()
    plt.savefig(str(outputname))



            
def main():        
 
    ########DATA_PATHS#########
    flexbb_monomer_score='/Users/philipke/Dropbox/Promotion/Projekte/Baker lab/PPID/Ctla4/bb_unconst/combined_enhanced_unconst.sc'
    fixedbb_monomer_score='/Users/philipke/Dropbox/Promotion/Projekte/Baker lab/PPID/Ctla4/bb_const/combined_enhanced_const.sc'
    # coordev01bb_monomer_score='/Users/ireyes2/Documents/Monomer/score_enhanced_01.sc'
    # coordev0001bb_monomer_score='/Users/ireyes2/Documents/Monomer/score_enhanced_0001.sc'
#    flexbb_complex_score='/Users/philipke/Dropbox/Promotion/Projekte/Baker lab/PPID/trka/bb_const/score_enhanced.sc'
#    fixedbb_complex_score='/Users/philipke/Dropbox/Promotion/Projekte/Baker lab/PPID/trka/bb_const/score_enhanced.sc'
    # coordev01bb_complex_score='/Users/ireyes2/Documents/Complex/coordev0.1.sc'
    # coordev0001bb_complex_score='/Users/ireyes2/Documents/Complex/coordev0.001.sc'
    
    data_flexbb_monomer = pandas.read_csv(flexbb_monomer_score, low_memory=False, header=0, skiprows=1, sep = '\s+', skipinitialspace=True, warn_bad_lines=True)
    data_fixedbb_monomer = pandas.read_csv(fixedbb_monomer_score, low_memory=False, header=0, skiprows=1, sep = '\s+', skipinitialspace=True, warn_bad_lines=True)
    # data_coordev01bb_monomer = pandas.read_csv(coordev01bb_monomer_score, low_memory=False, header=0, skiprows=1, sep = '\s+', skipinitialspace=True, warn_bad_lines=True)
    # data_coordev0001bb_monomer = pandas.read_csv(coordev0001bb_monomer_score, low_memory=False, header=0, skiprows=1, sep = '\s+', skipinitialspace=True, warn_bad_lines=True)
    
    # data_coordev01bb_complex = pandas.read_csv(coordev01bb_complex_score, low_memory=False, header=0, skiprows=1, sep = '\s+', skipinitialspace=True, warn_bad_lines=True)
    # data_coordev0001bb_complex = pandas.read_csv(coordev0001bb_complex_score, low_memory=False, header=0, skiprows=1, sep = '\s+', skipinitialspace=True, warn_bad_lines=True)
    
    #constraint_data_list = ['data_flex_bb_score','data_fixed_bb_score','data_coordev01_flex_bb_score','data_coordev0001_flex_bb_score'] 
    
    make_hists_from_four_scorefiles(data_flexbb_monomer, data_fixedbb_monomer, data_fixedbb_monomer, data_fixedbb_monomer, 'MonomerScore_Ctla4.pdf')  
    # make_hists_from_four_scorefiles(data_fixedbb_complex, data_flexbb_complex, data_coordev01bb_complex, data_coordev01bb_complex, 'ComplexScore.pdf')  

    #print (len(data_flexbb_complex['p_aa_pp'].describe()))
    #print (len(data_flexbb_monomer['p_aa_pp'].describe()))

    df_order_area1_noconstraint= data_flexbb_monomer[(data_flexbb_monomer["ddg_fa_atr_norepack_per_1000sasa"] <= -25) & (data_flexbb_monomer["interface_buried_sasa"] >= 1100) & (data_flexbb_monomer["worstfrag"] <= 1.5) & (data_flexbb_monomer["interface_unsat_hbond2"] <= 3) & (data_flexbb_monomer["interface_sc"] >= 0.65)]

    df_order_area1_defaultconstraint= data_fixedbb_monomer[(data_fixedbb_monomer["ddg_fa_atr_norepack_per_1000sasa"] <= -25) & (data_fixedbb_monomer["interface_buried_sasa"] >= 1100) & (data_fixedbb_monomer["worstfrag"] <= 1.5) & (data_fixedbb_monomer["interface_unsat_hbond2"] <= 3) & (data_fixedbb_monomer["interface_sc"] >= 0.65)]

    print ("Percentage of designs fulfilling these ORDER criteria in area1_noconstraint: "+ str(float(len(df_order_area1_noconstraint["description"].values))/len(data_flexbb_monomer["description"].values)) + ' NUMBER OF DESIGNS:'+str(len(df_order_area1_noconstraint["description"].values)))
    print ("Percentage of designs fulfilling these ORDER criteria in area1_defaultconstraint: "+ str(float(len(df_order_area1_defaultconstraint["description"].values))/len(data_fixedbb_monomer["description"].values)) + ' NUMBER OF DESIGNS:'+str(len(df_order_area1_defaultconstraint["description"].values)))
    
    ########MERGING###########
    
    
    # #### alter description to merge on:
    # data_flexbb_complex['merge_description']=[str(x)[:-5] for x in data_flexbb_complex['description'].values]
    # data_fixedbb_complex['merge_description']=[str(x)[:-5] for x in data_fixedbb_complex['description'].values]
    # data_coordev01bb_complex['merge_description']=[str(x)[:-5] for x in data_coordev01bb_complex['description'].values]
    # data_coordev0001bb_complex['merge_description']=[str(x)[:-5] for x in data_coordev0001bb_complex['description'].values]
    
    # #print (df_IL13_all_methods_interface["merge_description"].values[5000])
     
    # data_flexbb_monomer['merge_description']=[str(x)[:-24] for x in data_flexbb_monomer['description'].values]
    # data_fixedbb_monomer['merge_description']=[str(x)[:-24] for x in data_fixedbb_monomer['description'].values]
    # data_coordev01bb_monomer['merge_description']=[str(x)[:-24] for x in data_coordev01bb_monomer['description'].values]
    # data_coordev0001bb_monomer['merge_description']=[str(x)[:-24] for x in data_coordev0001bb_monomer['description'].values]
    
    
    
    
    # #df_hopu1_giant.to_csv(path_or_buf="giant_test.csv", sep=' ', na_rep='--', float_format=None, columns=None, header=True, index=True, index_label=None, mode='w', encoding=None, compression=None, quoting=None, quotechar='"', line_terminator='\n', chunksize=None, tupleize_cols=False, date_format=None, doublequote=True, escapechar=None, decimal='.')
    # #merge complex and monomer    
    # df_flex_giant=pandas.merge(data_flexbb_complex, data_fixedbb_monomer,suffixes=['_complex','_monomer'], how='outer', on=['merge_description'])
    # df_fix_giant=pandas.merge(data_fixedbb_complex, data_fixedbb_monomer,suffixes=['_complex','_monomer'], how='outer', on=['merge_description'])
    # df_01_giant=pandas.merge(data_coordev01bb_complex, data_coordev01bb_monomer,suffixes=['_complex','_monomer'], how='outer', on=['merge_description'])
    # df_0001_giant=pandas.merge(data_coordev0001bb_complex, data_coordev0001bb_monomer,suffixes=['_complex','_monomer'], how='outer', on=['merge_description'])
    
    # #df_flex_giant.to_csv("/Users/ireyes2/Desktop/test.csv",na_rep='--')
    
    # #return_list_of_numeric_columns(df_flex_giant)
    # #print (len(df_flex_giant['p_aa_pp_x'].describe()))
    # #print (len(df_flex_giant['p_aa_pp_y'].describe()))
    
    # #add a new column to data fra,es
    # df_flex_giant['backbone_constraint']='flex'
    # df_fix_giant['backbone_constraint']='fixed'
    # df_01_giant['backbone_constraint']='0.1'
    # df_0001_giant['backbone_constraint']='0.001'
    
    
    
    # #Dataframe w/ all data
    # data_mega = pandas.concat([df_flex_giant, df_fix_giant, df_01_giant, df_0001_giant])

    
    # #data_mega.to_csv("/Users/ireyes2/Desktop/test.csv",na_rep='--')
    
    # #########ORDER CRITERIA#############   
    
    # df_order_mega= data_mega[(data_mega["ddg_fa_atr_norepack_per_1000sasa"] <= -25) & (data_mega["interface_buried_sasa"] >= 1100) & (data_mega["worstfrag"] <= 1.5)
    # & (data_mega["interface_unsat_hbond2"] <= 3) & (data_mega["interface_sc"] >= 0.65)]
    
    # df_order_fix = df_order_mega[df_order_mega['backbone_constraint']=='fix']
    # df_order_flex = df_order_mega[df_order_mega['backbone_constraint']=='flex']
    # df_order_01 = df_order_mega[df_order_mega['backbone_constraint']=='0.1']
    # df_order_0001 = df_order_mega[df_order_mega['backbone_constraint']=='0.001']
    # print ("TOTAL DESIGNS: "+str(len(data_mega['merge_description'].values)))
    # print ("Number of designs fulfilling these ORDER criteria: "+ str(len(df_order_mega["merge_description"].values)))
    # print ("Percentage of designs fulfilling these ORDER criteria with FLEXIBLE BB: "+ str(float(len(df_order_flex["merge_description"].values))/len(df_order_mega["merge_description"].values)))
    # print ("Percentage of designs fulfilling these ORDER criteria with FIXED BB: "+ str(float(len(df_order_fix["merge_description"].values))/len(df_order_mega["merge_description"].values)))
    # print ("Percentage of designs fulfilling these ORDER criteria with 0.1 BB: "+ str(float(len(df_order_01["merge_description"].values))/len(df_order_mega["merge_description"].values)))
    # print ("Percentage of designs fulfilling these ORDER criteria with 0.001 BB: "+ str(float(len(df_order_0001["merge_description"].values))/len(df_order_mega["merge_description"].values)))
    # y="* \n".join(df_order_mega["merge_description"].values)
    # with open("test.list", 'w') as output:
    #     output.write(y)

main()
