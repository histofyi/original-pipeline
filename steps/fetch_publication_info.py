from typing import Dict, List, Tuple


from common.providers import PDBeProvider, PMCeProvider


from functions.pdb import load_pdb_lists
from functions.cli import load_config, print_spacer
from functions.files import load_constants, load_facet, write_facet

import doi

from rich.progress import Progress
from rich.console import Console
from rich import print_json

console = Console()
config = load_config()


mhc_class = 'class_i'

pdb_codes = load_pdb_lists('class_i', config['WAREHOUSE_PATH'], console)

#pdb_codes = ['3tbw']


pdbe_downloaded = []
pdbe_failed = []
pmce_failed = []
already_stored = []
missing_journal_info = []
no_doi = []
no_pubmed = []
no_bibjson = []

def convert_pdbe_authors(author: Dict) -> Dict:
    converted_author = {
        "lastname": author['last_name'],
        "name": author['full_name'],
        "initials": author["initials"]
    }
    return converted_author


def pdbe_to_bibjson(pdb_code, pdbe_publication_info:Dict) -> Tuple[Dict, bool, str]:
    has_doi = False
    article_doi = None
    bibjson = {
        'title':pdbe_publication_info['title'],
        'author':[convert_pdbe_authors(author) for author in  pdbe_publication_info['author_list']],
        'type':'article',
        'url':'',
        'identifier':[]
    }
    if pdbe_publication_info['journal_info']['pdb_abbreviation'] == 'To be published':
        missing_journal_info.append(pdb_code)
        print ('To be published')
    else:    
        bibjson['year'] = pdbe_publication_info['journal_info']['year'],
        bibjson['journal'] = {
            'name':'', 
            'iso_abbreviation':pdbe_publication_info['journal_info']['ISO_abbreviation']
        }

        bibjson['volume'] = pdbe_publication_info['journal_info']['volume'],
        bibjson['issue'] = pdbe_publication_info['journal_info']['issue'],
        bibjson['pages'] = pdbe_publication_info['journal_info']['pages'],
        
        if pdbe_publication_info['doi']:
            article_doi = pdbe_publication_info['doi']
            print (article_doi)
            has_doi = True
            bibjson['identifier'].append({'type':'doi', 'id':pdbe_publication_info['doi']})
        else:
            no_doi.append(pdb_code)
        if pdbe_publication_info['pubmed_id']:
            bibjson['identifier'].append({'type':'pubmed', 'id':pdbe_publication_info['pubmed_id']})
        else:
            no_pubmed.append(pdb_code)
    return bibjson, has_doi, article_doi



def fetch_pmce_data(doi:str, info:Dict) -> Dict:
    pmce_publication_info, success, errors = PMCeProvider().fetch(doi)
    if pmce_publication_info:
        yes_no = lambda x : True if x=='y' else False
        exists = lambda x,y : x[y] if y in x else None
        
        fields = [('open_access','isOpenAccess'),('in_pmc','inPMC'),('in_pmce','inEPMC'),('abstract','abstractText')]
        for field in fields:
            info[field[0]] = exists(pmce_publication_info, field[1])

    return info



with Progress() as progress:
        
    task = progress.add_task("[white]Processing...", total=len(pdb_codes))

    for pdb_code in pdb_codes:

        print (pdb_code)

        publication_info = load_facet(pdb_code, 'publication_details')

        if not publication_info:
            print ('NO PUBLICATION_INFO ALREADY')
            publication_info = {}

            pdbe_publication_info, success, errors = PDBeProvider(pdb_code).fetch_publications()

            if pdbe_publication_info:

                if len(pdbe_publication_info) > 0:

                    pdbe_downloaded.append(pdb_code)
                    
                    bibjson, has_doi, article_doi = pdbe_to_bibjson(pdb_code, pdbe_publication_info)
                    
                    publication_info = {
                        'open_access':None,
                        'in_pmc':None,
                        'in_pmce':None,
                        'abstract':None
                    }
                    if bibjson:
                        publication_info['bibjson'] = bibjson

                        if has_doi:
                            publication_info = fetch_pmce_data(article_doi, publication_info)                            
                            
                            publication_info['bibjson']['url'] = doi.get_real_url_from_doi(article_doi)
                            
                            write_facet(pdb_code, 'publication_details', publication_info)
                        else:
                            print ('NO DOI')
                            no_doi.append(pdb_code)
                            write_facet(pdb_code, 'publication_details', publication_info)
        
                        print(publication_info)
                    else:
                        print ('NO BIBJSON')
                        no_bibjson.append(pdb_code)

                else:
                    print ('EMPTY PDBe data')
                    pdbe_failed.append(pdb_code)
            else:
                print ('FAILED PDBe API call')
                pdbe_failed.append(pdb_code)
        else:
            print ('ALREADY STORED')
            already_stored.append(pdb_code)



        progress.update(task, advance=1)

print (f'{len(pdb_codes)} total records')
print (f'{len(already_stored)} already stored')
print (f'{len(pdbe_failed)} pdbe api failed')
print (f'{len(pdbe_downloaded)} pdbe api succeeded')
print (f'{len(pmce_failed)} pmce api failed')
print (f'{len(no_doi)} no DOI')
print (f'{len(no_bibjson)} no bibjson')


print ('PDBe API failures')
print (pdbe_failed)

print('PMCe API failures')
print (pmce_failed)

print ('No DOI')
print (no_doi)