import logging
import os 
import sys

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