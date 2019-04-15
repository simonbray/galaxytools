import numpy as np
from rdkit import Chem
from rdkit.Chem import rdShapeHelpers
import argparse

# ligand_sdf = Chem.SDMolSupplier(ligand_path)
# if len([x for x in ligand_sdf]) != 1:
#     raise IOError


def get_params(options):
    # mol = Chem.MolFromMolFile('test-data/NUDT5A-x0114_2.mol')
    mol = Chem.MolFromMolFile(options.ligand_path)
    if not mol:
        raise IOError

    print('test')
    conf = mol.GetConformer()
    params = rdShapeHelpers.ComputeConfBox(conf)

    # mol = Chem.MolFromSdfFile(ligand)


    coords1 = np.array(params[0])
    coords2 = np.array(params[1])

    print(coords1)
    print(coords2)


    center = np.mean((coords1, coords2), axis=0)

    box_dims = np.abs(coords1 - coords2)

    print(center)
    print(box_dims)

    with open(options.output, 'w') as f:
        f.write(
            """
size_x =  {}
size_y =  {}
size_z =  {}
center_x =  {}
center_y =  {}
center_z =  {}
num_modes = 9999
energy_range = 9999
exhaustiveness = 10
cpu = 4
seed = 1
            """.format(box_dims[0], box_dims[1], box_dims[2], center[0], center[1], center[2])
        )



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--ligand', dest='ligand_path')
    parser.add_argument('--config', dest='output')

    options = parser.parse_args()
    get_params(options)