'''
This file contains useful functions for the a345 practical projects.

Feel free to add your own functions, and ensure you add a brief description of the function below, and an example
of usage in a345_utilities_examples.ipynb:

mem_usage: function to display the memory used by variables of type ndarray visible in the local scope. It is easier to 
read than %whos, and also works in normal python scripts (%whos only works in ipython) 
    mem_usage(symbols=locals(), [minsize = 1000000], print_summary=True)
        locals(): pointer to the variables in local scope. can also use globals()
        minsize: only list variables about this size (default = 1000000 bytes)
        eg:  mem_usage(locals())                         # print all local variables above 2MB
             a = mem_usage(locals(), minsize=5000000)              # print all local variables above 5MB, total bytes returned as a
             mem_usage(globals(), minsize=2000000)       # print all global variables above 2MB
             tot = mem_usage(locals(), print_summary=False)       # get total number of bytes in ndarrays    
          



'''

import numpy as np


#-------------------------------------------------------------------------------------------------------------
def mem_usage(symbols, minsize = 1000000, print_summary=True):
    tot_bytes = 0
    for k in symbols.keys():
        if isinstance(symbols[k], np.ndarray):
            v = symbols[k]
            tot_bytes += v.nbytes
            if v.nbytes > minsize:
                sh=symbols[k].shape
                shapestr = str(sh[0])
                if len(sh) > 1:
                    for s in sh[1:]:
                        shapestr += ' x ' + str(s) 
                else:
                    shapestr += ' (1-D)'

                if print_summary:
                    # note: results printed in MB (10^6bytes, not MiB 2^20 bytes)
                    print('  {:30} {:>20s}   {:10s}  {:8.1f} MB'.format(k,
                                                                         shapestr,
                                                                         str(v.dtype),
                                                                         v.nbytes/(1000000)
                                                                         ))
    if print_summary:
        print('-'*90)
        print('{:>71}{:0.1f} MB'.format('TOTAL (inc unlisted small vars): ',tot_bytes/(1000000)))
        print('')
    return tot_bytes 


#-------------------------------------------------------------------------------------------------------------
def print_header(header, exclude_list=['COMMENT', 'HISTORY','SWOWNER','OBSERVER','']):
    # Prints the important parts of the header in two columns
   
    # first find how many keys are not in the exclude_list (default ['COMMENT', 'HISTORY']): we wont print these 
    # for some reason we seem to end up with blank keys sometimes - exculde those too
    keys = [x for x in header.keys() if not x in exclude_list]
  
    if len(keys) > 0:
        if exclude_list:
            print('Header info (excluding keys: {} '.format(exclude_list[0]),end='')
            for excl in exclude_list[1:]:
                print(', {}'.format(excl),end='')
            print('):')
                
        # now split list into two and print in two columns
        for k1,k2 in zip(keys[1:len(keys)//2], keys[len(keys)//2:]):
            print('    {:12s} : {:<50}  {:12s} : {:<50}'.format(k1, header[k1], k2, header[k2]))
        # if there was an odd number of keys print the last one
        if (len(keys) % 2) > 0:
            print('    {:12s} : {:<50}' .format(keys[-1], header[keys[-1]]))
    else:
        print('ERROR: no keys found in header')

        
        
#-------------------------------------------------------------------------------------------------------------
