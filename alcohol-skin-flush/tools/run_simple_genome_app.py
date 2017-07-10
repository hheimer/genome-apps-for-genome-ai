from os.path import basename, dirname, join, realpath

from pandas import DataFrame

from ga import get_ga_file_path, read_ga, write_ga
from vcf import get_vcf_variants_by_tabix

# Set paths
GENOME_APP_DIRECTORY_PATH = dirname(dirname(realpath(__file__)))
GENOME_APP_NAME = basename(GENOME_APP_DIRECTORY_PATH)
VCF_FILE_PATH = join(GENOME_APP_DIRECTORY_PATH, 'input/dna.vcf.gz')


def run_simple_genome_app():
    """
    Runs this Simple Genome App.
    return: None
    """

    # Read .input.ga
    headers, header_d, input_ga_df = read_ga(
        get_ga_file_path(GENOME_APP_DIRECTORY_PATH, 'input'))

    # Analyze with .input.ga
    matches = []
    for i, row in input_ga_df.iterrows():

        feature, feature_type, region, state = row[:4]
        state = str(state)

        for f, ft, r, s in zip(
                feature.split(';'),
                feature_type.split(';'), region.split(';'), state.split(';')):

            variants = get_vcf_variants_by_tabix(VCF_FILE_PATH, query_str=r)

            # Check for matches in the .VCF
            if len(variants):

                if ft == 'variant':

                    sample_genotype = variants[0]['sample'][0]['genotype']

                    if '|' in s:
                        if (s.split('|')) == sample_genotype or s.split(
                                '|') == sample_genotype[::-1]:
                            matches.append(i)

                    else:
                        if s in sample_genotype:
                            matches.append(i)

                elif ft == 'gene':

                    found_fields = {'effect': [], 'impact': []}

                    for variant in variants:
                        found_fields['effect'].append(
                            variant['ANN'][0]['effect'])
                        found_fields['impact'].append(
                            variant['ANN'][0]['impact'])

                    if s.upper() in ['HIGH', 'MODERATE', 'LOW'
                                     ] and s in found_fields['impact']:
                        matches.append(i)

                    elif s in found_fields['effect']:
                        matches.append(i)

                    elif 'MODERATE' in found_fields['impact']:
                        matches.append(i)

    if len(matches):  # Keep matching results
        output_ga_df = input_ga_df.ix[matches]

    else:  # Make a DataFrame with the default values
        output_ga_df = DataFrame(columns=input_ga_df.columns)
        output_ga_df['RESULT'] = [header_d['RESULT']['default']]
        for key in header_d:
            if 'default' in header_d[key].keys():
                output_ga_df[key] = header_d[key]['default']

    # Write .output.ga
    write_ga(headers, output_ga_df,
             get_ga_file_path(GENOME_APP_DIRECTORY_PATH, 'output'))

    print(
        'This Genome App was run and the output was saved as a table in /output/{}.output.ga.'.
        format(GENOME_APP_NAME))
