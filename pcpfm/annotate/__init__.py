"""
``pcpfm annotate``: annotation workflows.

Each workflow lives in its own prefixed modules (``lc*.py`` for LC-MS, ``gc*.py`` for GC-MS).
Only a workflow's ``run(args, record)`` touches the project record and leaves inside the master module of that annotation workflow.
The master modules are lc.py, gc.py and lcmsms.py which form the equivalent of process.py
The other modules take explicit data, references and parameters and have the helper functions
"""
