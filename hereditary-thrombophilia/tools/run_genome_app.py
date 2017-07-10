"""
Genome App entry point.
File for running Genome App from 1) Genome AI & 2) command line.
"""


def run_genome_app():
    """
    Required function for Guardiome to run this Genome App.
    :return: None
    """

    # Import run_simple_genome_app function from ./simple_genome_app.py
    from run_simple_genome_app import run_simple_genome_app

    # Run run_simple_genome_app on .VCF(.GZ)
    run_simple_genome_app()


# This Genome App can also be ran from comman line: $ python run_genome_app.py
if __name__ == '__main__':

    run_genome_app()
