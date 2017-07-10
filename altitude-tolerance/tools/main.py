import re
import json
import tabix

VCF_FIELD_CASTER = {
    'POS': int,
    'QUAL': float,
    'GT': lambda x: re.split('[|/]', x),
    'AD': lambda x: x.split(','),
    'VQSLOD': float,
    'CLNSIG': lambda x: max([int(s) for s in re.split('[,|]', x)]),
}


def parse_variant(vcf_row, n_anns=1):
    """

    :param vcf_row: iterable;
    :param n_anns: int;
    :return: dict; won't have the field if the its value is missing
    """

    variant = {
        'CHROM': _cast_vcf_field('CHROM', vcf_row[0]),
        'POS': _cast_vcf_field('POS', vcf_row[1]),
        'REF': _cast_vcf_field('REF', vcf_row[3]),
        'FILTER': _cast_vcf_field('FILTER', vcf_row[6]),
    }
    ref = variant['REF']

    # ID
    rsid = vcf_row[2]
    if rsid and rsid != '.':
        variant['ID'] = _cast_vcf_field('ID', rsid)

    # ALT
    alt = vcf_row[4]
    if alt and alt != '.':
        variant['ALT'] = _cast_vcf_field('ALT', alt)

    # QUAL
    qual = vcf_row[5]
    if qual and qual != '.':
        variant['QUAL'] = _cast_vcf_field('QUAL', qual)

    # Variant type
    if alt:
        vt = _get_variant_type(ref, alt)
        variant['variant_type'] = vt

    # Samples
    variant['samples'] = []
    format_ = vcf_row[8].split(':')
    for i, s in enumerate(vcf_row[9:]):
        s_d = {'sample_id': i + 1}

        # Sample
        for k, v in zip(format_, s.split(':')):
            s_d[k] = _cast_vcf_field(k, v)

        # Genotype
        if 'ALT' in variant:
            ref_alts = [variant['REF']] + variant['ALT'].split(',')
            s_d['genotype'] = [ref_alts[int(gt)] for gt in s_d['GT']]
        else:
            s_d['genotype'] = [variant['REF']] * 2

        # Allelic frequency
        if 'DP' in s_d and int(s_d['DP']):
            s_d['allelic_frequency'] = [round(int(ad) / int(s_d['DP']), 3) for ad in s_d['AD']]

        variant['samples'].append(s_d)

    info_split = vcf_row[7].split(';')
    for i_s in info_split:
        if i_s.startswith('ANN='):
            anns = {}
            for i, a in enumerate(i_s.split(',')[:n_anns]):
                a_split = a.split('|')

                anns[i] = {
                    'effect': a_split[1],
                    'putative_impact': a_split[2],
                    'gene_name': a_split[3],
                    'gene_id': a_split[4],
                    'feature_type': a_split[5],
                    'feature_id': a_split[6],
                    'transcript_biotype': a_split[7],
                    'rank': a_split[8],
                    'hgvsc': a_split[9],
                    'hgvsp': a_split[10],
                    'cdna_position': a_split[11],
                    'cds_position': a_split[12],
                    'protein_position': a_split[13],
                    'distance_to_feature': a_split[14],
                }
            variant['ANN'] = anns
        else:
            try:
                k, v = i_s.split('=')
                if v and v != '.':
                    # TODO: decode properly
                    variant[k] = _cast_vcf_field(k, v)
            except ValueError:
                pass
                # print('INFO error: {} (not key=value)'.format(i_s))

    return variant


def _cast_vcf_field(k, v, caster=VCF_FIELD_CASTER):
    """

    :param k: str;
    :param v: str;
    :param caster: dict;
    :return:
    """

    if k in caster:
        return caster[k](v)
    else:
        return v


def _get_variant_type(ref, alt):
    """

    :param ref: str;
    :param alt: str;
    :return: str;
    """

    if len(ref) == len(alt):
        if len(ref) == 1:
            vt = 'SNP'
        elif len(ref) == 2:
            vt = 'DNP'
        elif len(ref) == 3:
            vt = 'TNP'
        else:  # 4 <= len(ref)
            vt = 'ONP'

    elif len(ref) < len(alt):
        vt = 'INS'

    else:  # len(alt) < len(ref)
        vt = 'DEL'

    return vt


def _get_variant_start_and_end_positions(pos, ref, alt):
    """

    :param ref: str;
    :param alt: str;
    :return: (str, str);
    """

    if len(ref) == len(alt):
        s, e = pos, pos + len(alt) - 1

    elif len(ref) < len(alt):
        s, e = pos, pos + 1

    else:  # len(alt) < len(ref)
        s, e = pos + 1, pos + len(ref) - len(alt)

    return s, e


def get_variants_by_tabix(sample_vcf, contig=None, start=None, end=None, query_str=None, reference_vcf=None):
    """

    :param sample_vcf: str or pytabix handler;
    :param contig: str;
    :param start: int;
    :param end: int;
    :param query_str: int;
    :param reference_vcf: str or pytabix handler;
    :return: list; list of dict
    """

    if isinstance(sample_vcf, str):  # Open sample VCF
        sample_vcf = tabix.open(sample_vcf)

    if query_str:
        records = sample_vcf.querys(query_str)
    else:
        records = sample_vcf.query(contig, start, end)

    if reference_vcf and len(list(records)) == 0:  # If sample does not have the record, query reference if given

        if isinstance(reference_vcf, str):  # Open reference VCF
            reference_vcf = tabix.open(reference_vcf)

        records = reference_vcf.query(contig, start - 1, end)

    return [parse_variant(r) for r in records]


def _get_variants_from_regions(filepath_vcf, app):

    # Get variants
    variants = {}
    for r in app.get('regions'):
        for v in get_variants_by_tabix(filepath_vcf, query_str=r):
            if 'ID' in v:
                variants[v.pop('ID')] = v
            else:
                i = 1
                id_ = 'No ID ({})'.format(i)
                while id_ in variants:
                    i += 1
                    id_ = 'No ID ({})'.format(i)
                variants[id_] = v

    return variants


def _get_genes_from_variants(variants):

    # Get genes from variants
    genes = {}
    putative_impact = {'HIGH':4, 'MODERATE':3, 'LOW':2, 'MODIFIER':1}
    for id_, v in variants.items():
        a = v.get('ANN').get(0)
        if a['gene_name'] in genes:
            pi_score_existing = putative_impact[genes[a['gene_name']]['putative_impact']]
            pi_score_new = putative_impact[a['putative_impact']]
            if pi_score_new > pi_score_existing:
                genes[a.pop('gene_name')] = a
        else:
            genes[a.pop('gene_name')] = a

    return genes


def _inject_weight(number_to_inject, text):
    words = text.split()
    words = [w.replace('WWW', str(number_to_inject)) for w in words]
    words = ' '.join(words)
    return words


def _make_result_string(string, _list1, _list2):
    l = _list1 + _list2
    l.insert(0, string)
    l = ' '.join(l)
    return l


def _determine_result(app):
    result = {}
    result_based_info = {}
    n = 0
    variants = app['variants']
    genes = app['genes']
    for r, d in app.get('results', {}).items():
        logic = d.get('logic')
        analysis = d.get('analysis')

        variants_sub_results = []

        for item in analysis.keys():
            if 'variant' in item and len(analysis[item]['features']) > 0:
                variants_logic = analysis[item]['logic']

                if variants_logic == 'or':
                    variant_or_count = 0
                    for rs in analysis[item]['features']:
                        if rs in variants:
                            if variant_or_count == 0:
                                test_genotypes = analysis[item]['features'][rs]['genotypes']
                                customer_genotype = variants[rs]['samples'][0]['genotype']
                                for g in test_genotypes:
                                    if len(g) == 1:
                                        if g == customer_genotype[0] or customer_genotype[1]:
                                            variants_sub_results.append(analysis[item]['sub_result'])
                                            variant_or_count += 1
                                    else:
                                        if g == customer_genotype or g[::-1] == customer_genotype:
                                            variants_sub_results.append(analysis[item]['sub_result'])
                                            variant_or_count += 1

                if variants_logic == 'and':
                    count = 0
                    for rs in analysis[item]['features']:
                        if rs in variants:
                            test_genotypes = analysis[item]['features'][rs]['genotypes']
                            customer_genotype = variants[rs]['samples'][0]['genotype']
                            for g in test_genotypes:
                                if len(g) == 1:
                                    if g == customer_genotype[0] or customer_genotype[1]:
                                        count += 1
                                else:
                                    if g == customer_genotype or g[::-1] == customer_genotype:
                                        count += 1
                    if count == len(analysis[item]['features']):
                        variants_sub_results.append(analysis[item]['sub_result'])

                if variants_logic == 'add':
                    to_sum = []
                    for rs in analysis[item]['features']:
                        if rs in variants:
                            weight = analysis[item]['features'][rs]['weights']
                            test_genotypes = analysis[item]['features'][rs]['genotypes']
                            customer_genotype = variants[rs]['samples'][0]['genotype']
                            for g in test_genotypes:
                                if len(g) == 1:
                                    if g == customer_genotype[0] or customer_genotype[1]:
                                        to_sum.append(weight[test_genotypes.index(g)])
                                else:
                                    if g == customer_genotype or g[::-1] == customer_genotype:
                                        to_sum.append(weight[test_genotypes.index(g)])
                    if len(to_sum) > 0:
                        number_to_inject = sum(to_sum)
                        sub_result = _inject_weight(number_to_inject, analysis[item]['sub_result'])
                        variants_sub_results.append(sub_result)

                if variants_logic == 'multiply':
                    to_multiply = []
                    for rs in analysis[item]['features']:
                        if rs in variants:
                            weight = analysis[item]['features'][rs]['weights']
                            test_genotypes = analysis[item]['features'][rs]['genotypes']
                            customer_genotype = variants[rs]['samples'][0]['genotype']
                            for g in test_genotypes:
                                if len(g) == 1:
                                    if g == customer_genotype[0] or customer_genotype[1]:
                                        to_multiply.append(weight[test_genotypes.index(g)])
                                else:
                                    if g == customer_genotype or g[::-1] == customer_genotype:
                                        to_multiply.append(weight[test_genotypes.index(g)])
                    if len(to_multiply) > 0:
                        number_to_inject = 1
                        for n in to_multiply:
                            number_to_inject *= n
                        sub_result = _inject_weight(number_to_inject, analysis[item]['sub_result'])
                        variants_sub_results.append(sub_result)

        genes_sub_results = []

        for item in analysis.keys():
            if 'gene' in item and len(analysis[item]['features']) > 0:
                genes_logic = analysis[item]['logic']

                if genes_logic == 'or':
                    count = 0
                    for gene in analysis[item]['features']:
                        gene_features_logic = analysis[item]['features'][gene]['logic']
                        gene_features = analysis[item]['features'][gene]['fields']
                        if gene in genes:
                            if gene_features_logic == 'or':
                                if any(field for field in gene_features if gene_features[field] == genes[gene][field]):
                                    continue
                                else:
                                    count += 1
                            elif gene_features_logic == 'and':
                                for field in gene_features:
                                    if gene_features[field] != genes[gene][field]:
                                        count += 1
                    if count == 0:
                        genes_sub_results.append(analysis[item]['sub_result'])

                if genes_logic == 'and':
                    count = 0
                    for gene in analysis[item]['features']:
                        gene_features_logic = analysis[item]['features'][gene]['logic']
                        gene_features = analysis[item]['features'][gene]['fields']
                        if gene in genes:
                            if gene_features_logic == 'or':
                                if any(field for field in gene_features if gene_features[field] == genes[gene][field]):
                                    continue
                                else:
                                    count += 1
                            elif gene_features_logic == 'and':
                                for field in gene_features:
                                    if gene_features[field] != genes[gene][field]:
                                        count += 1
                    if count == 0:
                        genes_sub_results.append(analysis[item]['sub_result'])

        unwanted_keys = ['logic', 'analysis', 'result']

        if logic == 'or':
            if len(variants_sub_results) or len(genes_sub_results) is not 0:
                result[n] = _make_result_string(d['result'], variants_sub_results, genes_sub_results)
                for key, value in d.items():
                    if key not in unwanted_keys:
                        result_based_info[key] = d[key]

        if logic == 'and':
            if len(variants_sub_results) and len(genes_sub_results) is not 0:
                result[n] = _make_result_string(d['result'], variants_sub_results, genes_sub_results)
                for key, value in d.items():
                    if key not in unwanted_keys:
                        result_based_info[key] = d[key]

        n += 1
    if len(result.keys()) > 1:
        result = 'Inconclusive. This app returned multiple results.'
    elif len(result.keys()) == 1:
        for key in result.keys():
            result = result[key]

    return result, result_based_info


def run_simple_genome_app(path_to_app, filepath_vcf):

    """
    :param filepath_vcf: str;
    :param filepath_genome_app_json:
    :return:
    """
    path_to_app_json = path_to_app + 'data/app.json'
    with open(path_to_app_json) as fp:
        app = json.load(fp)

    app['variants'] = _get_variants_from_regions(filepath_vcf, app)
    app['genes'] = _get_genes_from_variants(app['variants'])

    app['result'], result_based_info = _determine_result(app)

    if len(app['result']) == 0:
        app['result'] = app['default_result']['result']
        for key in app['default_result']:
            if key != 'result':
                app[key] = app['default_result'][key]
    elif len(app['result']) > 1:
        if app['result'] != 'Inconclusive. This app returned multiple results.':
            for key in result_based_info:
                app[key] = result_based_info[key]

    del app['default_result']
    del app['results']
    del app['regions']

    path_to_result = path_to_app + 'results/result.json'
    with open(path_to_result, 'w') as outfile:
        json.dump(app, outfile, indent=True)

    return app
