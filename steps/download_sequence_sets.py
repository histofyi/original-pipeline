from typing import Dict
import os

from functions.files import write_file


sequence_collections = {
    'imgt_hla': {
        'url':'https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/hla_prot.fasta',
        'filename':'tmp/sequences/imgt_hla.fasta'
    },
    'ipd_mhc': {
        'url':'https://raw.githubusercontent.com/ANHIG/IPDMHC/Latest/MHC_prot.fasta',
        'filename':'tmp/sequences/ipd_mhc.fasta'
    }
}


def fetch_sequence_sets(sequence_collection:Dict) -> bool:
    """
    This function fetches bulk sequence collections from IPD

    Args:
        sequence_collection (Dict) - the url and filename for a sequence collection

    Returns:
        bool - whether or not the download is successful (based on amount of data returned)
    """

    filename = sequence_collection['filename']
    current_filesize = os.path.getsize(filename)
    url = sequence_collection['url']
    old_filename = f'{filename}_old'
    success = False
    if os.path.exists(filename):
        os.rename(filename, old_filename)
    download_string = f'curl {url} --output {filename}'
    os.system(download_string)
    if os.path.exists(filename):
        low = current_filesize - current_filesize*.1
        new_filesize = os.path.getsize(filename)
        if new_filesize > low:
            success = True
            os.remove(old_filename)
    return success

errors = []

for sequence_collection in sequence_collections:
    success = fetch_sequence_sets(sequence_collections[sequence_collection])
    if not success:
        errors.append(sequence_collection)
print (errors)








