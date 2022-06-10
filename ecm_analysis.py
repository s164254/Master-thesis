import ecm_analysis_base as ecm

protein_groups = (
    ('Collagen', False),
    ('Fibronectin', False),
    ('Laminin', False),
    ('Elastin', False),
    ('Proteoglycan', False),
)


def analysis_func(exp, samples, protein_group, use_ratio_filter, filter_func):
    remove_non_existing = True # =False does not work at the moment, see comments in fold_analysis function
    exp.fold_analysis(samples,
                      protein_group,
                      use_ratio_filter,
                      filter_func,
                      title = '',
                      xlabel='Gene names',
                      ylabel='Abundance',
                      remove_non_existing=remove_non_existing)


ecm.run_analysis(protein_groups, analysis_func)

#Abundance of selected collagen types for the samples
#Abundance of all laminin proteins for the samples
#Abundance of fibronectin for the samples