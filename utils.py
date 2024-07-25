import tqdm
import requests
import parsel

def get_pdb_info(pdb:str=None,url:str=None)->str:
    """
    查询RCSB PDB数据，获取PDB的resolution.
    Example:
    get_pdb_info(pdb='1A22') 
    """
    results = {'PDB': pdb, 'Method':None, 'Resolution':None, 'Release Date':None, 'PubMed ID':None}
    if url==None and pdb!=None:
        url = f"https://www.rcsb.org/structure/{pdb}"
    try:
        response = requests.get(url)
    except:
        return results
    # 检查请求是否成功
    if response.status_code == 200:
        # 通过 xpath 定位解析网页内容
        
        # method
        selector = parsel.Selector(text=response.text)
        xpath_format = '//li[@id="exp_header_0_method"]/text()'
        get_results = selector.xpath(xpath_format).getall()
        if len(get_results)==0:
            results['Method'] = None
        else:
            assert len(set(get_results))==1
            results['Method'] = get_results[0]
            
        # resolution
        selector = parsel.Selector(text=response.text)
        xpath_format = '//li[@id="exp_header_0_diffraction_resolution"]/text()'
        get_results = selector.xpath(xpath_format).getall()
        if len(get_results)==0:
            results['Resolution'] = None
        else:
            assert len(set(get_results))==1
            results['Resolution'] = get_results[0]
            
        # Release Date
        selector = parsel.Selector(text=response.text)
        xpath_format = '//li[@id="header_deposited-released-dates"]/text()'
        get_results = selector.xpath(xpath_format).getall()
        if len(get_results)==0:
            results['Release Date'] = None
        else:
            get_results = [date.replace('&nbsp','') for date in get_results]
            assert len(set(get_results))==2
            results['Release Date'] = get_results[1]
            
        # PubMed ID
        selector = parsel.Selector(text=response.text)
        xpath_format = '//li[@id="pubmedLinks"]/a/text()'
        get_results = selector.xpath(xpath_format).getall()
        if len(get_results)==0:
            results['PubMed ID'] = None
        else:
            assert len(set(get_results))==1
            results['PubMed ID'] = get_results[0]
            
    return results


from datetime import datetime
import requests

# +
# def convert_to_iso_date(date_str:str, dst_datetime_format:str='%Y-%m-%d'):
#     formats = ['%m/%d/%y', '%b/%d/%y', '%Y-%m-%d', '%Y-%b-%d', '%Y %m', '%Y %b', '%Y %m %d', '%Y %b %d', '%Y']  # 添加更多可能的日期格式
#     for fmt in formats:
#         try:
#             date_obj = datetime.strptime(date_str, fmt)
#             return date_obj.strftime(dst_datetime_format)
#         except ValueError:
#             pass
#     return 'Invalid date format'

from dateutil import parser
# 定义一个函数来解析日期，返回年份
def parse_year(date_str:str):
    try:
        return parser.parse(date_str).year
    except (parser.ParserError, TypeError):
        return date_str

# 定义一个函数来解析日期,返回年-月-日
def parse_date(date_str:str):
    try:
        return parser.parse(date_str)
    except (parser.ParserError, TypeError):
        return date_str

# +
# def get_pubdate(pubmed_id:str, dst_datetime_format:str='%Y-%m-%d'):
#     """
#     根据输入的pumbed_id， 查询PubMed数据库，获取其发表日期。
#     Reference: https://www.ncbi.nlm.nih.gov/pmc/tools/get-metadata/
#     """
#     url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id={pubmed_id}&retmode=json"
#     try:
#         response = requests.get(url)
#         date_str = response.json()['result'][pubmed_id]['pubdate']
#         formatted_date = convert_to_iso_date(date_str, dst_datetime_format)
#         return {'PMID':pubmed_id, 'publish date':formatted_date}
#     except:
#         return {'PMID':pubmed_id, 'publish date':None}
    
def get_pubdate(pubmed_id:str, dst_datetime_format:str='%Y-%m-%d'):
    """
    根据输入的pumbed_id， 查询PubMed数据库，获取其发表日期。
    Reference: https://www.ncbi.nlm.nih.gov/pmc/tools/get-metadata/
    """
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id={pubmed_id}&retmode=json"
    try:
        response = requests.get(url)
        date_str = response.json()['result'][pubmed_id]['pubdate']
        return {'PMID':pubmed_id, 'publish date':date_str}
    except:
        return {'PMID':pubmed_id, 'publish date':None}


# -

import requests
import parsel
def get_PubMedID_by_title(title:str):
    """
    根据输入的标题， 查询PubMed数据库，获取其PubMed ID。
    Reference: https://www.ncbi.nlm.nih.gov/pmc/tools/get-metadata/
    """
    url = 'https://pubmed.ncbi.nlm.nih.gov/?term={}'.format(title)
    response = requests.get(url)
    selector = parsel.Selector(text=response.text)
    xpath_format = "//meta[@name='citation_pmid']/@content"
    pmid = selector.xpath(xpath_format).getall()
    
    if len(set(pmid))==1:
        pmid = pmid[0]
    elif len(set(pmid))==0:
        pmid = None
    else:
        # assert len(set(pmid))>1, "Error: more than one pmid is obtain"
        pass
    return {'title': title, 'PMID':pmid}



from Bio.PDB import PDBParser, is_aa
def check_chains(pdb_file, ligand_chains, receptor_chains):
    """
    检查一个pdb中的ligand_chains, receptor_chains是否氨基酸链。
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)
    
    chain_dict = {}
    
    for model in structure:
        for chain in model:
            chain_id = chain.get_id()
            if chain_id in ligand_chains + receptor_chains:
                #import pdb; pdb.set_trace()
                #is_aa_chain = all(is_aa(residue) or residue.resname=='HOH' for residue in chain)
                is_aa_chain = any(is_aa(residue) for residue in chain)
                chain_dict[chain_id] = is_aa_chain
    return chain_dict


import os
import requests
def download_pdb(pdb_code, save_dir):
    url = f"https://files.rcsb.org/download/{pdb_code}.pdb"
    response = requests.get(url)
    
    if response.status_code == 200:
        file_path = os.path.join(save_dir, f"{pdb_code}.pdb")
        
        with open(file_path, "wb") as file:
            file.write(response.content)
        
        return file_path
    else:
        return None
