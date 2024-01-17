import sys
from lib2to3.main import main as l2to3_main


def main():
    sys.exit(l2to3_main("lib2to3.fixes"))
