from rdkit import Chem
from rdkit.Chem import AllChem
from collections import defaultdict

import matplotlib.pyplot as plt
import subprocess

def write_xyz(mol, file_name='temp.xyz'):
  number_of_atoms = mol.GetNumAtoms()
  symbols = [a.GetSymbol() for a in mol.GetAtoms()] 
  with open(file_name, "w") as file:
    file.write(str(number_of_atoms)+"\n")
    file.write("title\n")
    conf = mol.GetConformers()[0]
    for atom,symbol in enumerate(symbols):
      p = conf.GetAtomPosition(atom)
      line = " ".join((symbol,str(p.x),str(p.y),str(p.z),"\n"))
      file.write(line)

def show_mol(file_name, animate=False):
  xyz=open(file_name, 'r').read()
  p = py3Dmol.view(width=400,height=400)
  if animate:
    p.addModelsAsFrames(xyz,'xyz')
    p.animate({'loop': "forward",'reps': 5})
    #p.animate({'loop': 'backAndForth'})
  else:
    p.addModel(xyz,'xyz')
  p.setStyle({'stick':{}})
  p.setBackgroundColor('0xeeeeee')
  p.zoomTo()
  p.show()

def get_best_structure(mol,n_confs=10):
  new_mol = Chem.Mol(mol)

  AllChem.EmbedMultipleConfs(mol,numConfs=n_confs,useExpTorsionAnglePrefs=True,useBasicKnowledge=True)
  energies = AllChem.MMFFOptimizeMoleculeConfs(mol,maxIters=2000, nonBondedThresh=100.0)

  energies_list = [e[1] for e in energies]
  min_e_index = energies_list.index(min(energies_list))

  new_mol.AddConformer(mol.GetConformer(min_e_index))

  return new_mol

def shell(cmd, shell=False):
  #Written by Jimmy Kromann
  if shell:
    p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  else:
    cmd = cmd.split()
    p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

  output, err = p.communicate()
  return output

def run_xTB(mol, command, file_name='temp.xyz', conf_search=False):
  mol = Chem.AddHs(mol)
  if conf_search:
    mol = get_best_structure(mol)
  else:
    rdDistGeom.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol)
  write_xyz(mol)
  command = 'xtb ' + file_name + ' ' + command
  output = shell(command)
  output = str(output).replace('\\n','\n')
  
  return output

def get_energy(output):
  output = str(output)
  energy = float(output.split('TOTAL ENERGY')[1].split('Eh')[0])

  return energy

def run_rxn(reactant, smarts): 
  rxn = AllChem.ReactionFromSmarts(smarts)
  ps = rxn.RunReactants((reactant,))
  product = ps[0][0]
  Chem.SanitizeMol(product)
  
  return product

def reorder_product(product):
  reorder_inverse = [int(atom.GetProp('react_atom_idx')) for atom in product.GetAtoms()]
  reorder = len(reorder_inverse)*[0]

  for i in range(len(reorder_inverse)):
    reorder[reorder_inverse[i]] = i

  product = Chem.RenumberAtoms(product, reorder)

  return product

def label_atoms(mol):
  for i,atom in enumerate(mol.GetAtoms()):
    atom.SetAtomMapNum(i+1)
  
  return mol

def plot_PES():
  output = shell('grep energy xtbpath.xyz')
  energies = []
  for line in str(output).split('\\n'):
    if 'energy:' in line:
      energies.append(float(line.split('energy:')[1].split('xtb:')[0]))

  plt.plot(energies,'bo')
  plt.ylabel('Relative energy (kcal/mol)')
  plt.xlabel('Interpolation step')
  
def offset(mol, offset=5):
  smiles = Chem.MolToSmiles(mol, allHsExplicit=True)
  smilesA = smiles.split('.')[0]
  molA_idx = mol.GetSubstructMatch(Chem.MolFromSmiles(smilesA, sanitize=False))
  conf = mol.GetConformer()
  for i in range(mol.GetNumAtoms()):
    p = conf.GetAtomPosition(i)
    if i in molA_idx:
      conf.SetAtomPosition(i,(p.x,p.y,p.z+offset))

  return mol
