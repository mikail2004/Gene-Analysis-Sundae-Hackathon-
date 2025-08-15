import pandas as pd
import streamlit as st
import gseapy as gp
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

st.title("Gene Analysis Web App")
st.write("Compare treatment vs control wells from LINCS Level 5 data.")

# File uploaders
l5_file = st.file_uploader("Upload LINCS Level 5 subset CSV", type="csv")
gene_info_file = st.file_uploader("Upload LINCS gene info CSV", type="csv")

if l5_file and gene_info_file:
    L5_subset = pd.read_csv(l5_file)
    LG_gene = pd.read_csv(gene_info_file, sep="\t")
    
    # Let user select wells
    all_wells = [c[-3:] for c in L5_subset.columns if ":" in c]
    treatment_wells = st.multiselect("Select treatment wells", sorted(set(all_wells)))
    control_wells = ['A09', 'A17', 'B03', 'B14', 'C03', 'C06', 
                       'C12', 'C15', 'C18', 'C21', 'D14', 'D15', 
                       'D21', 'D22', 'E07', 'E14', 'E17', 'E24',
                       'F01', 'F03', 'F06', 'F09', 'G04', 'G09',
                       'G10', 'G13', 'G18', 'H04', 'H10', 'H13',
                       'H17', 'H18', 'I02', 'I09', 'I10', 'I21',
                       'I22', 'J14', 'J19', 'J23', 'K05', 'K06',
                       'K13', 'K16', 'K19', 'L02', 'L05', 'L06',
                       'L15', 'L18', 'M01', 'M04', 'M17', 'M22',
                       'M23', 'N07', 'N08', 'N17', 'N22', 'O07',
                       'O18', 'P05', 'P10']
    
    if treatment_wells and control_wells:
        treatment_cols = [col for col in L5_subset.columns if col[-3:] in treatment_wells]
        control_cols = [col for col in L5_subset.columns if col[-3:] in control_wells]
        
        # Data extraction
        L5_target = L5_subset[['rid'] + treatment_cols + control_cols]
        result = pd.merge(L5_target, LG_gene, left_on='rid', right_on='pr_gene_id', how='inner')
        result = result[result['pr_is_lm'] == 1]
        
        treated_mean = result[treatment_cols].mean(axis=1)
        control_mean = result[control_cols].mean(axis=1)
        result['diff_expr'] = treated_mean - control_mean
        
        result_sorted = result.sort_values(by='diff_expr', ascending=False)
        
        st.subheader("Top Differentially Expressed Genes")
        st.dataframe(result_sorted[['pr_gene_symbol', 'diff_expr']].head(20))
        
        # --- Volcano Plot ---
        st.subheader("Volcano Plot")
        plt.figure(figsize=(10,6))
        sns.scatterplot(
            x=result_sorted['diff_expr'],
            y=np.abs(result_sorted['diff_expr']),
            hue=(result_sorted['diff_expr']>0),
            palette={True:'red', False:'blue'},
            legend=False
        )
        plt.xlabel("Differential Expression (Treated - Control)")
        plt.ylabel("|Differential Expression|")
        plt.title("Volcano Plot of Gene Expression")
        st.pyplot(plt.gcf())
        plt.clf()

        # GSEA
        with st.spinner("Running pathway enrichment..."):
            gsea_df = gp.enrichr(
                gene_list=result_sorted['pr_gene_symbol'].dropna().tolist(),
                gene_sets=['KEGG_2019_Human', 'GO_Biological_Process_2021'],
                organism='Human',
                cutoff=0.05
            ).results
        
        st.subheader("Top Pathway Enrichment Results")
        st.dataframe(gsea_df.head(10))

        # --- Pathway Bar Plot ---
        st.subheader("Top Pathways Bar Plot")
        top_pathways = gsea_df.head(10)
        plt.figure(figsize=(10,6))
        sns.barplot(
            x='Combined Score',
            y='Term',
            data=top_pathways,
            palette="viridis"
        )
        plt.xlabel("Combined Score")
        plt.ylabel("Pathway")
        plt.title("Top Enriched Pathways")
        st.pyplot(plt.gcf())
        plt.clf()
