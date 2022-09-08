import toml


def print_spacer():
    print (f'{"-"*68}\n')


def load_config():
    config = toml.load('config.toml')
    return config