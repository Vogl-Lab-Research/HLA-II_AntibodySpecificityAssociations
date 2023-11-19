def replacer_code(code_str, elem1, replacer): # code str has to be given with the 3 commas ''' ''' to include all hidden symbols (newlines, tabs etc)
    str_format = str(code_str)
    str_out = str_format.replace(elem1, replacer)
    print(str_out)