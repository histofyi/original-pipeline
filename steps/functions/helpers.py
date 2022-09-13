from typing import List, Dict, Tuple

def slugify(string:str) -> str:
    slug_char = '_'
    to_replace = [' ','-','.',',','[',']','{','}','(',')','/','\\','*',':']
    for replace_me in to_replace:
        if replace_me in string:
            string = string.replace(replace_me, slug_char)
    return string.lower()


def deduplicate_list(pdb_code_list:List) -> List:
    deduplicated_list = [pdb_code for pdb_code in set(pdb_code_list)]
    sorted_list = sorted(deduplicated_list)
    return sorted_list