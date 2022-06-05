import ecm_analysis_base as ecm

protein_groups = (
    ('Collagen', False),
    ('Laminin', False),
    ('Fibronectin', False),
    ('Elastin', False),
    ('Proteoglycan', False),
)

def analysis_func(exp, samples, protein_group, use_filter, filter_func):
    exp.fold_analysis(samples, protein_group, use_filter, filter_func, normalize=True)

ecm.run_analysis(protein_groups, analysis_func)