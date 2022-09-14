from typing import List, Dict, Tuple

import datetime

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


def parse_date_to_isoformat(datestring:str) -> str:
    """
        This function takes an 8 digit datestring from the PDBe return and delivers an ISO formated date string
    """
    try:
        date = datetime.date(int(datestring[0:4]), int(datestring[4:6]), int(datestring[6:8]))
        return date.isoformat()
    except:
        return datestring

def parse_date_to_year(datestring:str) -> str:
    """
        This function takes an 8 digit datestring from the PDBe return and delivers a four digit year as a string
    """
    return datestring[0:4]