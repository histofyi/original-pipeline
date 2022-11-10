from typing import Dict, List, Tuple
from fuzzywuzzy import fuzz


def pick_best_title(core:Dict, facet:Dict) -> str:
    best_title = None
    try:
        publication_title = (core['publication']['bibjson']['title'])
    except:
        publication_title = None
    if publication_title:
        score = fuzz.ratio(publication_title.lower(), facet['pdb_title_capitalized'].lower())
        if score > 80:
            best_title = publication_title
        else:
            best_title = facet['pdb_title_capitalized']
    else:
        best_title = facet['pdb_title_capitalized']
    return best_title


def clean_title(core:Dict, facet:Dict) -> str:
    cleaned_title = core['pdb_title']
    if 'peptide_sequence' in core:
        if core['peptide_sequence'] is not None:
            if core['peptide_sequence'].lower() in cleaned_title:
                cleaned_title = cleaned_title.replace(core['peptide_sequence'].lower(), core['peptide_sequence'].upper())
    capitalised_words = ['MHC', 'Class I','MR1', 'MAIT', 'CD1','LIR','KIR','H2-D','H2-L','H2-K', 'TCR','T cell']
    for capitalised_word in capitalised_words:
        if capitalised_word.lower() in cleaned_title.lower():
            cleaned_title = cleaned_title.replace(capitalised_word.lower(), capitalised_word)
    return cleaned_title



def build_titles(core:Dict) -> Tuple[str,str]:
    display_title = ''
    page_title = ''

    complex_type = core['complex_type']

    if 'peptide_sequence' in core:
        peptide = core['peptide_sequence']
    else:
        peptide = None

    if 'alpha' in core['allele']:
        alpha_chain = core['allele']['alpha']['name']
    else:
        alpha_chain = None

    alpha_beta = 'Alpha/Beta T cell receptor'
    gamma_delta = 'Gamma/Delta T cell receptor'
    pre_tcr_beta = 'pre TCR beta chain'

    cd8a = 'CD8a'

    plc = 'Peptide Loading Complex'
    tapasin = 'Tapasin'
    tapbpr = 'peptide editor TAPBPR'


    nk = 'Natural Killer'
    nkg2a = 'NKG2A'
    nkg2d = 'NKG2D'
    kir = 'KIR NK receptor'
    kir3 = 'KIR-3 NK receptor'
    lir1 = 'LIR-1 NK receptor'
    lir2 = 'LIR-2 NK receptor'
    lirb1 = 'LIRB-1 NK receptor'
    ly49a = 'Ly49a NK receptor'
    ly49c = 'Ly49c NK receptor'
    cd94 = 'CD94'


    chimeric = 'chimeric'
    single_chain = 'single chain'
    antibody = 'antibody'
       
    non_classical_class_i = 'Non-classical MHC Class I molecule'
    
    cd1a = 'CD1a'
    cd1b = 'CD1b'
    cd1c = 'CD1c'
    cd1d = 'CD1d'
    mr1 = 'MR1'
    mica = 'MICA'
    micb = 'MICB'
    fcrn = 'Feonatal Fc receptor (FcRn)'
    
    serum_albumin = 'serum albumin'
    
    # Class I with peptide
    if complex_type == 'class_i_with_peptide':
        display_title = f'{alpha_chain} binding "{peptide}"'
    elif complex_type == 'class_i_with_peptide_fragments':
        display_title = f'{alpha_chain} binding peptide fragments'

    # Class I with TCR/antibody
    elif complex_type == 'class_i_with_peptide_and_alpha_beta_tcr':
        display_title = f'{alpha_chain} presenting "{peptide}" to {alpha_beta}'
    elif complex_type == 'class_i_with_peptide_and_alpha_beta_tcr':
        display_title = f'{alpha_chain} presenting "{peptide}" to {gamma_delta}' 
    elif complex_type == 'class_i_with_peptide_and_antibody':
        display_title = f'{alpha_chain} binding "{peptide}" with {antibody}'
    
    # Class I with CD8
    elif complex_type == 'class_i_with_peptide_and_cd8a':
        display_title = f'{alpha_chain} binding "{peptide}" with {cd8a}'
    elif complex_type == 'class_i_with_cd8a':
        display_title = f'{alpha_chain}  with {cd8a}'
    elif complex_type == 'class_i_with_peptide_and_cd8_dimer':
        display_title = f'{alpha_chain} binding "{peptide}" with CD8 dimer'

     # Class I with peptide loading/editing
    elif complex_type == 'class_i_with_tapbpr':
        display_title = f'{alpha_chain} with {tapbpr}' 
    elif complex_type == 'class_i_with_tapasin':
        display_title = f'{alpha_chain} with {tapasin}' 
    elif complex_type == 'class_i_plc':
        display_title = f'{alpha_chain} with {plc}'
    elif complex_type == 'class_i_with_erp57_and_tapasin':
        display_title = f'{alpha_chain} with erp57 and {tapasin}' 

    # Class I with NK receptors
    elif complex_type == 'class_i_with_peptide_and_kir':
        display_title = f'{alpha_chain} binding "{peptide}" with {kir}'
    elif complex_type == 'class_i_with_peptide_and_kir3':
        display_title = f'{alpha_chain} binding "{peptide}" with {kir3}'
        
    elif complex_type == 'class_i_with_peptide_and_lir1':
        display_title = f'{alpha_chain} binding "{peptide}" with {lir1}'
    elif complex_type == 'class_i_with_peptide_and_lir2':
        display_title = f'{alpha_chain} binding "{peptide}" with {lir2}'
    elif complex_type == 'class_i_with_peptide_and_lirb1':
        display_title = f'{alpha_chain} binding "{peptide}" with {lirb1}'

    elif complex_type == 'class_i_with_peptide_and_ly49a':
        display_title = f'{alpha_chain} binding "{peptide}" with {ly49a}'            
    elif complex_type == 'class_i_with_peptide_and_ly49c':
        display_title = f'{alpha_chain} binding "{peptide}" with {ly49c}'
    elif complex_type == 'class_i_with_peptide_and_cd94_and_nkg2a':
        display_title = f'{alpha_chain} binding "{peptide}" with {cd94} and {nkg2a}'

    # Class I without distinct peptide etc
    elif complex_type == 'class_i_possibly_without_peptide':
        display_title = f'{alpha_chain} possibly without "peptide"'

    # Class I with viral proteins
    elif complex_type == 'class_i_with_peptide_and_us2':
        display_title = f'{alpha_chain} binding "{peptide}" with human cytomegalovirus protein US2'
    # TODO combine these two complex types
    elif complex_type == 'class_i_with_peptide_and_e3' or complex_type == 'class_i_with_peptide_and_e319k':
        display_title = f'{alpha_chain} binding "{peptide}" with Human adenovirus type 4 E3-19K '
    elif complex_type == 'class_i_with_peptide_and_cowpox_cpxv203':
        display_title = f'{alpha_chain} binding "{peptide}" with Cowpox CPXV203'

    # Truncated Class I
    elif complex_type == 'truncated_class_i_with_peptide':
        display_title = f'Truncated {alpha_chain} binding "{peptide}"'
    elif complex_type == 'truncated_class_i_with_peptide_and_alpha_beta_tcr':
        display_title = f'Truncated {alpha_chain} presenting "{peptide}" to {alpha_beta}'
    elif complex_type == 'truncated_class_i_with_peptide_and_single_chain_tcr_construct':
        display_title = f'Truncated {alpha_chain} presenting "{peptide}" to {single_chain} {alpha_beta} construct'
    elif complex_type == 'truncated_class_i_with_peptide_and_pre_tcr_beta':
        display_title = f'Truncated {alpha_chain} binding "{peptide}" with {pre_tcr_beta}'

    # Single chain Class I
    elif complex_type == 'single_chain_class_i_construct':
        display_title = f'Single chain construct of {alpha_chain} binding "{peptide}"'
    elif complex_type == 'single_chain_class_i_construct_with_gamma_delta_tcr':
        display_title = f'Single chain construct of {alpha_chain} presenting "{peptide}" to {gamma_delta}'

    # CD1
    elif complex_type == 'cd1a':
        display_title = f'{non_classical_class_i} {cd1a}'
    elif complex_type == 'cd1a_with_nkt_alpha_beta_tcr':
        display_title = f'{non_classical_class_i} {cd1a} with {nk} {alpha_beta}'
    elif complex_type == 'cd1a_with_nkt_gamma_delta_tcr':
        display_title = f'{non_classical_class_i} {cd1a} with {nk} {gamma_delta}'
    elif complex_type == 'cd1b':
        display_title = f'{non_classical_class_i} {cd1b}'
    elif complex_type == 'cd1b_with_nkt_alpha_beta_tcr':
        display_title = f'{non_classical_class_i} {cd1b} with {nk} {alpha_beta}'
        print (display_title)
    elif complex_type == 'cd1c':
        display_title = f'{non_classical_class_i} {cd1c}'
    elif complex_type == 'cd1d':
        display_title = f'{non_classical_class_i} {cd1d}'
    elif complex_type == 'cd1d_with_b2m_and_peptide':
        display_title = f'{non_classical_class_i} {cd1d} binding "{peptide}"'
    elif complex_type == 'cd1d_with_nkt_alpha_beta_tcr':
        display_title = f'{non_classical_class_i} {cd1d} with {nk} {alpha_beta}'
    elif complex_type == 'cd1d_with_chimeric_alpha_beta_tcr':
        display_title = f'{non_classical_class_i} {cd1d} with {chimeric} {alpha_beta}'
    elif complex_type == 'cd1d_with_antibody':
        display_title = f'{non_classical_class_i} {cd1d} with {antibody}'
    elif complex_type == 'cd1d_with_single_chain_trc':
        display_title = f'{non_classical_class_i} {cd1d} with {single_chain} gamma TCR construct'    

    # MR1
    elif complex_type == 'mr1':
        display_title = f'{non_classical_class_i} {mr1}'
    elif complex_type == 'mr1_with_alpha_beta_tcr':
        display_title = f'{non_classical_class_i} {mr1} with {alpha_beta}'
    elif complex_type == 'mr1_with_gamma_delta_tcr':
        display_title = f'{non_classical_class_i} {mr1} with {gamma_delta}'
    elif complex_type == 'mr1_with_cd8a':
        display_title = f'{non_classical_class_i} {mr1} with {cd8a}'
    elif complex_type == 'mr1_with_tapbpr':
        display_title = f'{non_classical_class_i} {mr1} with {tapbpr}'
       
    # MICA/MICB
    elif complex_type == 'mica':
        display_title = f'{non_classical_class_i} {mica}'
    elif complex_type == 'mica_with_nkg2d':
        display_title = f'{non_classical_class_i} {mica} with {nkg2d}'
    elif complex_type == 'micb':
        display_title = f'{non_classical_class_i} {micb}'

    # FcRn
    elif complex_type == 'fcrn':
        display_title = f'{non_classical_class_i} {fcrn}'
    elif complex_type == 'fcrn_with_antibody':
        display_title = f'{non_classical_class_i} {fcrn} with {antibody}'
    elif complex_type == 'fcrn_with_peptide_and_beta2m':
        display_title = f'{non_classical_class_i} {fcrn} with "{peptide}"'
    elif complex_type == 'fcrn_with_beta2m_and_serum_albumin':
        display_title = f'{non_classical_class_i} {fcrn} with {serum_albumin}'
    elif complex_type == 'fcrn_with_beta2m_and_ig_heavy_chain_gamma':
        display_title = f'{non_classical_class_i} {fcrn} with Ig heavy chain gamma'
    elif complex_type == 'fcrn_with_ig_gamma_and_serum_albumin':
        display_title = f'{non_classical_class_i} {fcrn} with Ig heavy chain gamma and {serum_albumin}'     
    elif complex_type == 'fcrn_with_somatotropin_and_serum_albumin':
        display_title = f'{non_classical_class_i} {fcrn} with somatoropin and {serum_albumin}' 
    elif complex_type == 'fcrn_with_echovirus_proteins':
        display_title = f'{non_classical_class_i} {fcrn} in complex with echovirus 18 proteins' 

    # Unusual ones
    elif complex_type == 'zag':
        display_title = f'{non_classical_class_i} Zinc Alpha-2 Glycoprotein'
    elif complex_type == 'zag_with_pip':
        display_title = f'{non_classical_class_i} Zinc Alpha-2 Glycoprotein complexed with Prolactin inducible protein (PIP)'

    elif complex_type == 'h2-t22':
        display_title = f'Mouse {non_classical_class_i} H2-T22'
    elif complex_type == 'h2_22_with_gamma_delta_tcr':
        display_title = f'Mouse {non_classical_class_i} H2-T22 with {gamma_delta}'
    elif complex_type == 'hfe2':
        display_title = f'{non_classical_class_i} haemochromatosis protein (HFE)'
    elif complex_type == 'hfe2_with_beta2m_and_transferrin_receptor':
        display_title = f'{non_classical_class_i} haemochromatosis protein (HFE) complexed with the transferrin receptor'




    else:
        print (core['pdb_code'])
        print (complex_type)

    if display_title == '':
        print (complex_type)
      
    page_title = f'{core["pdb_code"].upper()} | {display_title}'
    print (page_title)


    if core['resolution'] is not None and len(display_title) > 0:
        resolution = core['resolution']
        display_title = f'{display_title} at {resolution}&#8491; resolution'
    




    
    
    
    return display_title, page_title