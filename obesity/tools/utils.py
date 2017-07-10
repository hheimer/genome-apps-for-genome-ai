
def get_region_from_rsid(rsid_list):
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
                print(line)

            output[line[0]] = {'chr': int(line[1]),
                               'start': int(line[2]),
                               'stop': int(line[3])}
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
                output[line[0]] = {'chr': int(line[1]),
                                   'start': int(line[2]),
                                   'stop': int(line[3])}
                print(missing[i] + ' -> ' + line[0])
            except:
                print("Error getting chr pos for: {}".format(line[0]))
                pass

    region_list = [str(output[x]["chr"]) + ":" + str(output[x]["start"]) + "-" + str(output[x]["stop"]) for x in output
                   if output[x]["chr"] and output[x]["start"] and output[x]["stop"]]

    return region_list


def get_region_from_gene(gene_list):
    from biomart import BiomartServer

    server = BiomartServer("http://uswest.ensembl.org/biomart/")
    db = server.datasets['hsapiens_gene_ensembl']
    response = db.search({'filters': {'source': 'ensembl_havana',
                                      'hgnc_symbol': gene_list},
                          'attributes': ['hgnc_symbol', 'chromosome_name', 'start_position', 'end_position']})

    output = dict()
    for line in response.iter_lines():
        line = line.decode('utf-8').split("\t")
        try:
            if line[0] in output.keys():
                print('Multiple entries for: {}'.format(line[0]))
                print(line)

            output[line[0]] = {'chr': int(line[1]),
                               'start': int(line[2]),
                               'stop': int(line[3])}
        except:
            if line[0] not in output.keys():
                print('Error getting chr pos for: {}'.format(line[0]))
            pass

    present = set(output.keys())
    if len(present) < len(gene_list):
        missing = [x for x in gene_list if x not in present]
        print("Missing following genes: " + ', '.join(missing))

    region_list = [str(output[x]["chr"]) + ":" + str(output[x]["start"]) + "-" + str(output[x]["stop"]) for x in output
                   if output[x]["chr"] and output[x]["start"] and output[x]["stop"]]

    return region_list
