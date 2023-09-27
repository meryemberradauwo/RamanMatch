"""
This is a setup.py script generated by py2applet

Usage:
    python setup.py py2app
"""
from setuptools import setup

APP = ['RamanSept21.py']
DATA_FILES = ['RRUFFRaman_database.db', 'RRUFFRaman_databaseSEARCH.db']
OPTIONS = {
    'argv_emulation': True,
    'packages': ["openpyxl", "scipy", "matplotlib", "Pillow", "pandas", "IPython"],
}

setup(
    app=APP,
    data_files=DATA_FILES,
    options={'py2app': OPTIONS},
    setup_requires=['py2app'],  # Add this line
)