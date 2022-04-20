from os import path
import re
import csv

ATTR_RE = '_([a-zA-Z0-9]+)\.raw\.PG\.(.+)'

#c = '20220401_KMJ_LB_PB_140min_15cm_M2.raw.PG.PSMs'
#m = re.search(ATTR_RE,c)
def parse(cell_value):
    if cell_value.find(';')>=0:
        return tuple([parse(v) for v in cell_value.split(';')])

    temp = cell_value.strip()
    if not temp:
        return -1
    
    if cell_value[-1] == '%':
        return float(cell_value[:-1])
    
    return float(cell_value)

def write_csv(csv_file_name,header,data):
    with open(csv_file_name, 'w') as csvfile:
        wrtr = csv.writer(csvfile, delimiter=',')
        wrtr.writerow(header)
        for row in data:
            wrtr.writerow(row)

def common_proteins_criteria(v1,v2):
    return v1>0 and v2>0


def is_fold(v1,v2,n):
    ma = max(v1,v2)
    if ma == 0:
        return False

    mi = min(v1,v2)
    if mi == 0:
        return True

    return ma/mi >= 2


def is_two_fold(v1,v2):
    return is_fold(v1,v2,2)

class MsAnalysis:
    def __init__(self, csv_file_name):
        with open(path.join(script_dir, csv_file_name)) as csvfile:
            rdr = csv.reader(csvfile, delimiter=',')
            rows = [row for row in rdr]

        self.header = rows[0]
        data = rows[1:]
        item_attrs = [(i,m.groups()[0],m.groups()[1]) for i,m in [(i,re.search(ATTR_RE,c)) for i,c in enumerate(self.header)] if m and len(m.groups())==2]
        self.attrs = tuple(set([attr_name for i,item_name,attr_name in item_attrs]))
        self.item_names = tuple(set([item_name for i,item_name,attr_name in item_attrs]))

        self.UniprotIds = tuple([row[0] for row in data])
        self.Genes = tuple([row[1] for row in data])
        self.ProteinDescriptions = tuple([row[1] for row in data])

        item_data = {}
        for i,item_name,attr_name in item_attrs:
            if not item_name in item_data:
                item_data[item_name] = {}
            item_data[item_name][attr_name] = tuple([parse(row[i]) for row in data])
        self.item_data = item_data

    def helper(self,item_names,attr_name,criteria):
        item1,item2 = item_names
        v1 = self.item_data[item1][attr_name]
        v2 = self.item_data[item2][attr_name]
        temp = [(i,v1[i],v2[i]) for i in range(len(self.Genes)) if criteria(v1[i],v2[i])]
        return [(self.UniprotIds[i],self.Genes[i],self.ProteinDescriptions[i],v1,v2) for i,v1,v2 in temp]

    def common_proteins(self,item_names,attr_name):
        return self.helper(item_names,attr_name,common_proteins_criteria)

    def fold_change(self,item_names,attr_name):
        return self.helper(item_names,attr_name,is_two_fold)

script_dir = path.dirname(path.realpath(__file__))
fname = 'Batch2_data.output.csv'

m = MsAnalysis(path.join(script_dir, fname))
item_names = ['M2','P2']
common_proteins = m.common_proteins(item_names,'Label-Free Quant')
fold_change = m.fold_change(item_names,'Label-Free Quant')
write_csv(path.join(script_dir, 'csv.csv'),m.header[:3]+item_names,common_proteins)
