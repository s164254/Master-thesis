import math


def parse_csv_cell(cell_value, not_valid_value):
    if cell_value.find(';') >= 0:
        return tuple([
            parse_csv_cell(v, not_valid_value)
            for v in cell_value.split(';')
        ])

    temp = cell_value.strip()
    if not temp:
        return not_valid_value

    if cell_value[-1] == '%':
        return float(cell_value[:-1])

    return float(cell_value)


def unique_peptides_valid_value(value):
    if not value:
        return False

    return value and (value != '1') and (value[0].isdigit() or value[0] == '.')


def unique_peptides_valid_row(row_data):
    return any((1 for v in row_data if unique_peptides_valid_value(v)))


def label_free_quant_transform(value, unique_peptide_value):
    return unique_peptide_value[0].isdigit(
    ) and unique_peptide_value != '1' and math.log2(float(value)) or '0'


#def csv_row_criteria =
