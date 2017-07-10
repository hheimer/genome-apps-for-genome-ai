
CONTIG = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18',
          '19', '20', '21', '22', 'X', 'Y', 'MT']


def get_region_from_rsid(rsid_list):
    region_list = []

    if rsid_list:
        from biomart import BiomartServer
        server = BiomartServer("http://uswest.ensembl.org/biomart/")
        db = server.datasets['hsapiens_snp']
        response = db.search({'filters': {'variation_source': 'dbSNP',
                                          'snp_filter': rsid_list},
                              'attributes': ['refsnp_id', 'chr_name', 'chrom_start', 'chrom_end']})

        output = dict()
        for line in response.iter_lines():
            line = line.decode('utf-8').split("\t")
            try:
                if line[0] in output.keys():
                    print('Multiple entries for: {}'.format(line[0]))

                if line[1] in CONTIG:
                    output[line[0]] = {'chr': line[1],
                                       'start': line[2],
                                       'stop': line[3]}

            except:
                if line[0] not in output.keys():
                    print('Error getting chr pos for: {}'.format(line[0]))
                pass

        present = set(output.keys())
        if len(present) < len(rsid_list):
            missing = [x for x in rsid_list if x not in present]
            print("Missing following rsID: " + ', '.join(missing))

            response = db.search({'filters': {'variation_source': 'dbSNP',
                                              'snp_synonym_filter': missing},
                                  'attributes': ['refsnp_id', 'chr_name', 'chrom_start', 'chrom_end']})

            print('Converted the following rsID:')
            for i, line in enumerate(response.iter_lines()):
                line = line.decode('utf-8').split("\t")
                try:
                    if line[1] in CONTIG:
                        output[line[0]] = {'chr': line[1],
                                           'start': line[2],
                                           'stop': line[3]}
                    print(missing[i] + ' -> ' + line[0])
                except:
                    print("Error getting chr pos for: {}".format(line[0]))
                    pass

        region_list = [str(output[x]["chr"]) + ":" + str(output[x]["start"]) + "-" + str(output[x]["stop"]) for x in output
                       if output[x]["chr"] and output[x]["start"] and output[x]["stop"]]

    return region_list


def get_region_from_gene(gene_list):
    region_list = []

    if gene_list:
        from biomart import BiomartServer
        server = BiomartServer("http://uswest.ensembl.org/biomart/")
        db = server.datasets['hsapiens_gene_ensembl']
        response = db.search({'filters': {'hgnc_symbol': gene_list},
                              'attributes': ['hgnc_symbol', 'chromosome_name', 'start_position', 'end_position']})

        output = dict()
        for line in response.iter_lines():
            line = line.decode('utf-8').split("\t")
            try:
                if line[0] in output.keys():
                    print('Multiple entries for: {}'.format(line[0]))

                if line[1] in CONTIG:
                    output[line[0]] = {'chr': line[1],
                                       'start': line[2],
                                       'stop': line[3]}
            except:
                if line[0] not in output.keys():
                    print('Error getting chr pos for: {}'.format(line[0]))
                pass

        present = set(output.keys())
        if len(present) < len(gene_list):
            missing = [x for x in gene_list if x not in present]
            print("Missing following genes: " + ', '.join(missing))

        region_list = [str(output[x]["chr"]) + ":" + str(output[x]["start"]) + "-" + str(output[x]["stop"]) for x in
                       output
                       if output[x]["chr"] and output[x]["start"] and output[x]["stop"]]

    return region_list


def get_hgnc_symbol_from_gene(gene_list):
    hgnc_list = list()

    if gene_list:
        import gzip
        import numpy as np

        with gzip.open('../data/Homo_sapiens.gene_info.gz', 'rb') as f:
            for line in f:
                line = line.decode('utf-8').split('\t')
                genes = line[4].split('|')
                if any(x in genes for x in gene_list):
                    hgnc_list.append(line[2])
                    gene_list = np.setdiff1d(gene_list, genes + [line[2]])

            hgnc_list = gene_list.tolist() + hgnc_list

    return hgnc_list


def get_region_from_file(file):
    import pandas as pd
    import re

    data = pd.read_csv(file, header=None)
    data = data.iloc[:, 0]
    data_type = [1 if re.match(r"\brs\d*\b", x) else 0 for x in data]

    variant_regions = get_region_from_rsid([data[i] for i, x in enumerate(data_type) if x == 1])
    gene_regions = get_region_from_gene(get_hgnc_symbol_from_gene([data[i] for i, x in enumerate(data_type) if x == 0]))

    return variant_regions + gene_regions
