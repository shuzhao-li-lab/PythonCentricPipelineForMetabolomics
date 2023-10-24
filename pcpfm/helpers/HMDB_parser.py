# parse HMDB ver 4, single XML files for all metabolites
# Retrieved from https://hmdb.ca/downloads
# Using Python 3 and lxml
# lxml is not in current Python 3 package. 
# It has similar API as ElementTree but faster. To install:
# pip install lxml

'''
# This took ~5 seconds on my i5 laptop, infile is 1.2 GB.
# Larger files will have to be partitioned, depending on memeory limit.

In [1]: from parse_hmdb4 import *                                                              

In [2]: myList = parse_a_file(infile, wanted)                                                  
Found 25411 entries.

In [3]: myList[0]                                                                              
Out[3]: 
{'accession': 'HMDB0000001',
 'name': '1-Methylhistidine',
 'chemical_formula': 'C7H11N3O2',
 'monisotopic_moleculate_weight': '',
 'iupac_name': '(2S)-2-amino-3-(1-methyl-1H-imidazol-4-yl)propanoic acid',
 'traditional_iupac': '1 methylhistidine',
 'cas_registry_number': '332-80-9',
 'smiles': 'CN1C=NC(C[C@H](N)C(O)=O)=C1',
 'inchi': 'InChI=1S/C7H11N3O2/c1-10-3-5(9-4-10)2-6(8)7(11)12/h3-4,6H,2,8H2,1H3,(H,11,12)/t6-/m0/s1',
 'inchikey': 'BRMWTNUJHUMWMS-LURJTMIESA-N',
 'pathways': '',
 'normal_concentrations': '',
 'abnormal_concentrations': '',
 'diseases': '',
 'drugbank_id': 'DB04151',
 'drugbank_metabolite_id': '',
 'phenol_explorer_compound_id': '',
 'phenol_explorer_metabolite_id': '',
 'foodb_id': 'FDB093588',
 'knapsack_id': '',
 'chemspider_id': '83153',
 'kegg_id': 'C01152',
 'biocyc_id': '',
 'bigg_id': '',
 'wikipidia': '',
 'nugowiki': '',
 'metagene': '',
 'metlin_id': '3741',
 'pubchem_compound_id': '92105',
 'het_id': '',
 'chebi_id': '50599',
 'protein_associations': ''}

In [40]: with open('HMDB0000031.xml', 'w') as O: 
    ...:     O.write(str(ET.tostring(x), 'utf-8')) 
    ...:                                                                                       

# Alternatively, cut the big file into 5 parts, then parse:

In [2]: myResults = []                                                                         

In [3]: wanted = ['accession', 'name', 'chemical_formula', 'monisotopic_molecular_weight', 'keg
   ...: g_id', 'pubchem_compound_id', ]                                                        

In [4]: import os                                                                              

In [5]: files = os.listdir('parts/')                                                           

In [6]: files                                                                                  
Out[6]: 
['README.txt',
 'HMDB4_part5.xml',
 'HMDB4_part3.xml',
 'HMDB4_part2.xml',
 'HMDB4_part1.xml',
 'HMDB4_part4.xml']

In [7]: for f in files[1:]: 
   ...:     myResults += parse_a_file('parts/' + f, wanted) 
   ...:                                                                                        
Found 30770 entries.
Found 21801 entries.
Found 17478 entries.
Found 22506 entries.
Found 21667 entries.

In [8]: write_tsv(myResults, wanted, "HMDB4_compounds.tsv")                

'''

from lxml import etree as ET 
import json
from mass2chem.formula import calculate_formula_mass

infile = 'serum_metabolites.xml'

# fields to retrieve
wanted = ['accession', 'name', 'chemical_formula', 'monisotopic_molecular_weight', 'iupac_name', 
    'traditional_iupac', 'cas_registry_number', 'smiles', 'inchi', 'inchikey']

def extract(obj_metabolite, x):
    try:
        y = obj_metabolite.find(x).text.strip()
        if y:
            return y
        else: return ''
    except AttributeError: 
        return ''

def extract_dict(obj_metabolite, wanted, prefix):
    result = {}
    for x in wanted:
        result[x] = extract(obj_metabolite, prefix+x)
        if x == "chemical_formula":
                result["neutral_formula"] = result[x]
                result["neutral_formula_mass"] = calculate_formula_mass(result[x])

                #result["neutral_formula"] = None
                #result["neutral_mass"] = None
        elif x == "accession":
            result["primary_id"] = result[x]
    result["primary_db"] = "HMDBv5"
    return result

def parse_a_file(f, wanted, prefix='{http://www.hmdb.ca}'):
    '''
    Input
    =====
    f: HMDB 4 XML file
    wanted: fields to retrieve
    prefix: prefix that is used in the XML file

    Return
    ======
    return as a list of dictionaries. 
    '''
    tree = ET.parse(f)
    root = tree.getroot()
    print("Found %d entries." %len(root))

    results = []
    for child in root:
        result = extract_dict(child, wanted, prefix=prefix)
        if result:
            if "neutral_formula" in result and "neutral_formula_mass" in result and "primary_id":
                if result["neutral_formula_mass"] and result["neutral_formula"] and result["primary_id"]:
                    results.append(extract_dict(child, wanted, prefix=prefix))
    return results

def write_json(results, wanted, outfile):
    json.dump(results, open(outfile, 'w'), indent=4)

def write_tsv(results, wanted, outfile):
    s = '\t'.join(wanted) + '\n'
    for R in results:
        s += '\t'.join([R[x] for x in wanted]) + '\n'
    with open(outfile, 'w') as O:
        O.write(s)

if __name__ == '__main__':
    # wanted = ['accession', 'name', 'chemical_formula', 'monisotopic_molecular_weight']
    infile = '/Users/mitchjo/Downloads/hmdb_metabolites.xml'
    # infile= 'urine_metabolites.xml'
    # infile= 'feces_metabolites.xml'

    write_tsv(
        parse_a_file(infile, wanted),
        wanted, 
        infile.replace(".xml", ".tsv")
    )
    write_json(
        parse_a_file(infile, wanted),
        wanted,
        infile.replace(".xml", ".json")
    )
