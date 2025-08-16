from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from PIL import Image, ImageDraw, ImageOps
import numpy as np
import os
from pathlib import Path
import gc

def make_cyclobutane_art(
    input_path='input.png',
    output_mol='output.mol',
    space=4.0 # Intermolecular distance (Å, recommended 3.0-5.0)
):
    RDLogger.DisableLog('rdApp.*')
    
    current = Path(__file__).parent
    input_img = current / input_path
    output_mol = current / output_mol

    if not input_img.exists():
        img = Image.new('L', (64, 64), 255) # 64x64
        draw = ImageDraw.Draw(img)
        draw.rectangle((16, 16, 48, 48), outline=0, width=1)
        img.save(input_img)

    def preprocess_image(img_path):
        img = Image.open(img_path).convert('L')
        img = img.resize((64, 64), Image.NEAREST) # 64x64
        img = ImageOps.autocontrast(img)
        return np.array(img) < np.mean(img)

    pixels = preprocess_image(input_img)

    template = Chem.MolFromSmiles('C1CCC1') # cyclobutane
    template = Chem.AddHs(template)
    AllChem.EmbedMolecule(template)
    AllChem.MMFFOptimizeMolecule(template)

    mol = Chem.RWMol()
    active_pixels = np.argwhere(pixels)
    
    pos_list = []
    for y, x in active_pixels:
        for atom in template.GetAtoms():
            p = template.GetConformer().GetAtomPosition(atom.GetIdx())
            new_p = (x*space + p.x, (63-y)*space + p.y, p.z)
            pos_list.append(new_p)
            mol.AddAtom(atom)

    conf = Chem.Conformer(len(pos_list))
    for i, p in enumerate(pos_list):
        conf.SetAtomPosition(i, p)

    per_ring = template.GetNumAtoms()
    for i in range(0, len(pos_list), per_ring):
        if i + per_ring - 1 < len(pos_list):
            for j in range(per_ring - 1):
                mol.AddBond(i+j, i+j+1, Chem.BondType.SINGLE)
            mol.AddBond(i+per_ring-1, i, Chem.BondType.SINGLE)

    mol.AddConformer(conf)
    os.makedirs(output_mol.parent, exist_ok=True)
    Chem.MolToMolFile(mol, str(output_mol))

    del mol, template, conf
    gc.collect()
    print(f"ok")

if __name__ == "__main__":
    make_cyclobutane_art()