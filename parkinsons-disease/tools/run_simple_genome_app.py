from os.path import basename, dirname, join, realpath

from pandas import DataFrame

from g2p import get_g2p_file_path, read_g2p, write_g2p
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

    # Read .input.g2p
    headers, header_d, input_g2p_df = read_g2p(
        get_g2p_file_path(GENOME_APP_DIRECTORY_PATH, 'input'))

    # Analyze with .input.g2p
    matches = []
    for i, row in input_g2p_df.iterrows():

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

                    found_impact = set()

                    for variant in variants:
                        found_impact.update([variant['ANN'][x]['impact'] for x in variant['ANN']])

                    if s.upper() in ['HIGH', 'MODERATE', 'LOW'] and s.upper() in found_impact:
                        matches.append(i)

                    elif s.upper() not in ['HIGH', 'MODERATE', 'LOW'] and 'MODERATE' in found_impact:
                        matches.append(i)

    if len(matches):  # Keep matching results
        output_g2p_df = input_g2p_df.ix[matches]

    else:  # Make a DataFrame with the default values
        output_g2p_df = DataFrame(columns=input_g2p_df.columns)
        output_g2p_df['RESULT'] = [header_d['RESULT']['default']]
        for key in header_d:
            if 'default' in header_d[key].keys():
                output_g2p_df[key] = header_d[key]['default']

    # Write .output.g2p
    write_g2p(headers, output_g2p_df,
             get_g2p_file_path(GENOME_APP_DIRECTORY_PATH, 'output'))

    print(
        'This Genome App was run and the output was saved as a table in /output/{}.output.g2p.'.
        format(GENOME_APP_NAME))

    print(output_g2p_df)
