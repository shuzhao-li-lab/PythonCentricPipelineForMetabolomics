import logging
import os 
import sys
import csv
import itertools
import json

log = logging.getLogger(__name__)
def validate_path(path_to_file, fatal=False):
    log.info("Determining if " + path_to_file + " is a valid path")
    is_valid = os.path.isfile(path_to_file)
    if fatal and not is_valid:
        log.exception(path_to_file + " is not a valid file path, does it exist?")
        exit()
    return is_valid

def retrieve_abs_path(path_to_file):
    log.info("Determing the absolute path for " + path_to_file)
    abs_path = os.path.abspath(path_to_file)
    return abs_path

def create_directory(path_to_directory, overwrite=False):
    try:
        log.info("Creating directory " + path_to_directory)
        os.makedirs(path_to_directory)
    except:
        log.exception("Creating directory " + path_to_directory + " failed!")

def adductify_standards(standards_csv, adducts_csv):
    isotope_mass_table = {
        "12C": 12.0,
        "13C": 13.00335483507,
        "1H": 1.00782503223,
        "2H": 2.01410177812,
        "16O": 15.99491461957,
        "14N": 14.00307400443,
        "15N": 15.00010889888,
        "32S": 31.9720711744,
        "31P": 30.97376199842,
        "23Na": 22.9897692820,
        "e": 0.00054858
    }
    standards = []
    with open(standards_csv, encoding='utf-8-sig') as standards_fh:
        for standard in csv.DictReader(standards_fh):
            standards.append(standard)
    adducts = []
    with open(adducts_csv, encoding='utf-8-sig') as adducts_fh:
        for adduct in csv.DictReader(adducts_fh):
            adducts.append(adduct)
    adducted_standards = []
    for standard, adduct in itertools.product(standards, adducts):
        new_name = standard["Name"] + '[' + adduct["Name"] + '],z=' + str(adduct["Charge"])
        new_formula = {}
        for isotope, count in json.loads(standard['Isotope Dictionary']).items():
            count = int(count)
            if isotope not in new_formula:
                new_formula[isotope] = count
            else:
                new_formula[isotope] += count
        for isotope, count in json.loads(adduct['Isotope Dictionary']).items():
            count = int(count)
            if isotope not in new_formula:
                new_formula[isotope] = count
            else:
                new_formula[isotope] += count
        new_formula['e'] = -1 * int(adduct['Charge'])
        uncharged_mass = sum([isotope_mass_table[isotope] * isotope_count for isotope, isotope_count in new_formula.items()])
        mz = uncharged_mass / abs(int(adduct['Charge']))
        adducted_standards.append({
            "Name": new_name,
            "Isotope Formula Dict": new_formula,
            "Search Mass": mz,
        })
    return adducted_standards