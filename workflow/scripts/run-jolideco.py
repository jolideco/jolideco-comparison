from jolideco.core import MAPDeconvolver
from jolideco.models import FluxComponents


def read_datasets():
    pass


def run_jolideco(datasets):
    components = FluxComponents.from_dict()

    deco = MAPDeconvolver(**config)
    result = deco.run(datasets=datasets, components=components)
    result.write()


if __name__ == "__main__":
    datasets = read_datasets()
    run_jolideco(datasets)
