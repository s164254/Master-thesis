import ecm_analysis_base as ecm

protein_groups = (
    ('Collagen', False),
    ('Fibronectin', False),
    ('Laminin', False),
    ('Elastin', False),
    ('Proteoglycan', False),
)

title_dict = {
    'Collagen': 'Abundance of selected collagen types for the samples',
    'Fibronectin': 'Abundance of fibronectin for the samples',
    'Laminin': 'Abundance of all laminin proteins for the samples',
}

yaxis_overload = { 
    'ecm_foldchange.m1_p1.Fibronectin': (1e8,1e10),
    'ecm_foldchange.m1_m2.Fibronectin': (1e8,1e10),
    'ecm_foldchange.fa1a_fa2a.Fibronectin': (5e6,5e7),
}

def analysis_func(exp, samples, protein_group, use_ratio_filter, filter_func):
    remove_non_existing = True # =False does not work at the moment, see comments in fold_analysis function
    exp.fold_analysis(samples,
                      protein_group,
                      use_ratio_filter,
                      filter_func,
                      title = title_dict.get(protein_group, ''),
                      xlabel='Gene names',
                      ylabel='Abundance',
                      remove_non_existing=remove_non_existing, yaxis_overload = yaxis_overload)


ecm.run_analysis(protein_groups, analysis_func)
