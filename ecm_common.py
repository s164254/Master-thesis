import ecm_analysis_base as ecm

protein_groups = (
    ('Collagen', False),
    ('Laminin', False),
    ('Fibronectin', False),
)

def analysis_func(exp, samples, protein_group, use_filter, filter_func):
    exp.ecm_common(samples, protein_group, use_filter, filter_func)

ecm.run_analysis(protein_groups, analysis_func)
