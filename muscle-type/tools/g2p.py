from os import listdir
from os.path import join, basename

from pandas import read_table


def get_g2p_file_path(directory_path, type_):
    """
    Get .g2p file path.
    :param directory_path: str; Genome App directory path
    :param type_: str; 'input' | 'output'
    :return: str; Genome App .g2p file path
    """

    type_directory_path = join(directory_path, type_, '')

    try:
        file_path = join(type_directory_path, [
            fn for fn in listdir(type_directory_path)
            if fn.endswith('.{}.g2p'.format(type_))][0])
    except IndexError:
        file_path = join(type_directory_path, '{}.{}.g2p'.format(basename(directory_path), type_))

    return file_path


def read_g2p(g2p_file_path):
    """
    Read Genome App .g2p by reading headers and table.
    :param g2p_file_path: str; .g2p file path
    :return: list & dict & DataFrame; .g2p headers & header dict & table
    """

    # Read headers
    with open(g2p_file_path) as f:

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
    return headers, header_d, read_table(g2p_file_path, comment='#')


def write_g2p(headers, g2p_df, file_path):
    """
    Write Genome App .g2p.
    :param: headers: list; .g2p header
    :param: g2p_df: DataFrame; .g2p table
    :param: str: .g2p file path
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
        g2p_df.to_csv(f, sep='\t', index=None)
