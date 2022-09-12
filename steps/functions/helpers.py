

def slugify(string):
    slug_char = '_'
    to_replace = [' ','-','.',',','[',']','{','}','(',')','/','\\','*',':']
    for replace_me in to_replace:
        if replace_me in string:
            string = string.replace(replace_me, slug_char)
    return string.lower()