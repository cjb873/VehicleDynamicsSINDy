import yaml


def load_parameters(file_name='params.yml'):

    params = None

    with open(file_name, 'r') as yaml_file:
        params = yaml.safe_load(yaml_file)

    return params
