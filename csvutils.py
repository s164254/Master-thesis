from os import path
import csv


def write_csv(csv_file_name, header, data):
    with open(csv_file_name, 'w') as csvfile:
        wrtr = csv.writer(csvfile, delimiter=',')
        wrtr.writerow(header)
        for row in data:
            wrtr.writerow(row)


def load_csv(fname, delim, with_header):
    with open(fname, 'r') as csvfile:
        rdr = csv.reader(csvfile, delimiter=delim)
        rows = [row for row in rdr]
        if len(rows) > 0:
            if with_header:
                return (rows[0], rows[1:])
            else:
                return (None, rows)
        return None
