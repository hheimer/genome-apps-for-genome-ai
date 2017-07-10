from os import listdir
from os.path import join, basename

from pandas import read_table


def get_ga_file_path(directory_path, type_):
    """
    Get .ga file path.
    :param directory_path: str; Genome App directory path
    :param type_: str; 'input' | 'output'
    :return: str; Genome App .ga file path
    """

    type_directory_path = join(directory_path, type_, '')

    try:
        file_path = join(type_directory_path, [
            fn for fn in listdir(type_directory_path)
            if fn.endswith('.{}.ga'.format(type_))][0])
    except IndexError:
        file_path = join(type_directory_path, '{}.{}.ga'.format(basename(directory_path), type_))

    return file_path


def read_ga(ga_file_path):
    """
    Read Genome App .ga by reading headers and table.
    :param ga_file_path: str; .ga file path
    :return: list & dict & DataFrame; .ga headers & header dict & table
    """

    # Read headers
    with open(ga_file_path) as f:

        headers = []
        header_d = {}

        for line in f:
            if line.startswith('#'):

                headers.append(line)

                line = line[1:].strip()

                ks, v = line.split('=')

                if '.' in ks:
                    leaf_d = header_d
                    ks = ks.split('.')
                    for i, k in enumerate(ks):
                        if k not in leaf_d:

                            if i != len(ks) - 1:
                                leaf_d[k] = {}
                            else:
                                leaf_d[k] = v

                        leaf_d = leaf_d[k]

                else:
                    header_d[ks] = v

    # Read table
    return headers, header_d, read_table(ga_file_path, comment='#')


def write_ga(headers, ga_df, file_path):
    """
    Write Genome App .ga.
    :param: headers: list; .ga header
    :param: ga_df: DataFrame; .ga table
    :param: str: .ga file path
    :return: None
    """

    with open(file_path, 'w') as f:

        # Write headers
        if len(headers):
            for h in headers:
                a, b = h.split('=')
                a_c, a_d = a.split('.')
                if a_d != 'default':
                    f.writelines(h)

        # Write table
        ga_df.to_csv(f, sep='\t', index=None)
