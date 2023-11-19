#!/usr/bin/python3.9

def replace_vals(row):
        if type(row['uniref_org']) != float:
            if ' ' not in row['uniref_org'] and type(row['IEDB_org']) != float:
                return row['IEDB_org']
            elif ' ' not in row['uniref_org'] and type(row['VFDB_org']) != float:
                return row['VFDB_org']
            elif ' ' not in row['uniref_org'] and type(row['bac_source_org']) != float:
                return row['bac_source_org']
            else:
                return row['uniref_org']
        else:
            return row['uniref_org']