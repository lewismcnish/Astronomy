import os
import sys

for file in os.listdir('/Users/lewismcnish/Documents/vscode/Astronomy/Binary_black_holes/inv_psd_data'):
    if file.endswith('strain.txt'):
        os.replace(file, '/Users/lewismcnish/Documents/vscode/Astronomy/Binary_black_holes/strain_data/'+file)