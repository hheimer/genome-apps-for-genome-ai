import os
import sys
import json
import git
from datetime import datetime

PATH_TO_RESULT = 'results/result.json'
PATH_TO_TOOLS = 'tools/'


def timestamp_genome_app_result(result, path_to_app):
    """
    Add a timestamp key to the Genome App result.json to indicate its been run.
    :param result: dict; result dictionary
    :param path_to_app: str; path to Genome App folder
    """

    result['date_run'] = str(datetime.now().time())
    path_to_result = path_to_app + PATH_TO_RESULT
    sys.path.insert(0, path_to_result)
    with open(path_to_result, 'w') as outfile:
        json.dump(result, outfile, indent=True)


def read_genome_app_result(path_to_app):
    """
    Read a Genome App result json.
    :param path_to_app: str; path to Genome App folder
    """

    path_to_result = path_to_app + PATH_TO_RESULT
    with open(path_to_result) as fp:
        result = json.load(fp)
    return result


def check_genome_app_result(path_to_app):
    """
    Check if Genome App has been run.
    :param path_to_app: str; path to Genome App folder
    """
    
    if os.path.isfile(path_to_app + 'results/result.json'):
        result = read_genome_app_result(path_to_app)
        if 'date_run' not in result:
            to_run = True
        else:
            to_run = False
    else:
        to_run = True
    return to_run


def run_genome_app(path_to_app):
    """
    Check if Genome App has been run before and if not, runs it.
    :param path_to_app: str; path to Genome App folder
    """
    to_run = check_genome_app_result(path_to_app)
    if to_run == True:
        path_to_run = path_to_app + PATH_TO_TOOLS
        sys.path.insert(0, path_to_run)
        from run import run_app
        run_app(path_to_app)
        result = read_genome_app_result(path_to_app)
        result = timestamp_genome_app_result(result, path_to_app)
    result = read_genome_app_result(path_to_app)
    return result


def list_genome_apps(path_to_catalog):
    """
    Open and read catalog of Genome Apps.
    :param path_to_catalog: str; path to catalog.json
    """

    with open(path_to_catalog) as fp:
        catalog = json.load(fp)
    return catalog


# ==============================================================================
#  Functions for adding a Genome App
# ==============================================================================
def add_genome_app_to_catalog(path_to_app, path_to_catalog):
    """
    Adds Genome App to the Genome App catalog.
    :param path_to_app: str; path to Genome App folder
    :param path_to_catalog: str; path to catalog.json
    """

    # Get Genome App name
    app_name = os.path.basename(os.path.normpath(path_to_app))

    # Clean Genome App name
    edited_app_name = app_name.replace('-', ' ')
    edited_app_name = edited_app_name.replace('_', ' ')
    edited_app_name = edited_app_name.replace('.', ' ')
    edited_app_name = edited_app_name.title()

    # Add Genome App name to catalog.json
    with open(path_to_catalog) as fp:
        catalog = json.load(fp)

    ids = (list(catalog.keys()))
    ids = [int(n) for n in ids]
    new_id = max(ids) + 1

    app_to_add = {"name": edited_app_name}
    catalog[str(new_id)] = app_to_add

    # Save catalog.json
    with open(path_to_catalog, 'w') as outfile:
        json.dump(catalog, outfile, indent=True)


# TODO: Make this function work with Genome App store
def add_genome_app(path_to_app, path_to_catalog):
    """
    Add a Genome App.
    :param path_to_app: str; path to Genome App folder
    :param path_to_catalog: str; path to catalog.json
    """

    # Add Genome App repo to genome_apps folder
    os.chdir('../genome_apps')

    # Add github repo Genome App to genome_apps folder
    if 'https' in path_to_app:
        git clone -- recursive(path_to_app)
        # move repo into genome_apps folder

    # Add Genome App from USB to genome_apps folder
    else:
        # move gneome app folder from USB to genome_apps folder

    # Add Genome App name to catalog.json
    add_genome_app_to_catalog(path_to_app, path_to_catalog)

    # Link customer vcf and tbi to data folder of Genome App
    os.symlink('../../genome-data/variant/738.concat_snp_indel.extract_chr.recode.ann_snpeff.ann_clinvar.rename_chr.vcf.gz', '../app_name/data/')
    os.symlink('../../genome-data/variant/738.concat_snp_indel.extract_chr.recode.ann_snpeff.ann_clinvar.rename_chr.vcf.gz.tbi', '../app_name/data/')
