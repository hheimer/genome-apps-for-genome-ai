from main import run_simple_genome_app


def run_app(path_to_app):
    path_to_customer_vcf = path_to_app + 'data/customer.vcf.gz'
    run_simple_genome_app(path_to_app, path_to_customer_vcf)
