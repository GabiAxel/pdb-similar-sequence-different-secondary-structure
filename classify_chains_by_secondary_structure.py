import gzip
from argparse import ArgumentParser
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from enum import Enum
from pathlib import Path
from typing import List, Dict, Tuple, Set, DefaultDict, Any

from biotite.structure.io.pdbx import PDBxFile
from pydash import py_, invert_by, map_values
from tqdm import tqdm

MODELED_RESIDUES_THRESHOLD = 0.7
STRUCTURED_RESIDUES_THRESHOLD = 0.5


class SecondaryStructure(Enum):
    Helical = 'helical'
    Extended = 'extended'


def get_category_as_chain(cif: PDBxFile, category: str):
    dict_of_arrays = cif.get_category(category, expect_looped=True)
    if dict_of_arrays is None:
        return py_([])
    return py_([dict(zip(dict_of_arrays,i)) for i in zip(*dict_of_arrays.values())])


def get_chain_sequences(cif: PDBxFile) -> Dict[str, str]:
    return (get_category_as_chain(cif, 'entity_poly')
        .filter({'type': 'polypeptide(L)'})
        .map(lambda e: (e['pdbx_strand_id'].split(',')[0], e['pdbx_seq_one_letter_code_can']))
        .from_pairs()
        .value())


def get_chain_residues(cif: PDBxFile, chain_ids: List[str]) -> Dict[str, Set[int]]:
    return (get_category_as_chain(cif, 'pdbx_poly_seq_scheme')
        .filter(lambda r: r['pdb_strand_id'] in chain_ids)
        .group_by('pdb_strand_id')
        .to_pairs()
        .filter(lambda pair: py_(pair[1]).filter(lambda r: r['pdb_mon_id'] != '?').size().value() / len(pair[1]) >= MODELED_RESIDUES_THRESHOLD)
        .from_pairs()
        .map_values(lambda rs: set(map(lambda r: int(r['seq_id']), rs)))
        .value())


def get_secondary_structure(struct_conf_entry: Dict[str, Any]):
    if struct_conf_entry['conf_type_id'].startswith('HELX'):
        return SecondaryStructure.Helical
    if struct_conf_entry['conf_type_id'] == 'STRN':
        return SecondaryStructure.Extended
    return None


def get_structure_ranges(cif: PDBxFile, chains_ids: List[str]) -> DefaultDict[str, DefaultDict[SecondaryStructure, List[Tuple[int, int]]]]:
    struct_conf = get_category_as_chain(cif, 'struct_conf').map(lambda s: {**s, 'ss': get_secondary_structure(s)}).filter(lambda s: s['ss'] is not None)
    chain_structure_ranges: Dict[str, Dict[SecondaryStructure, Tuple[int, int]]] = defaultdict(lambda _: defaultdict(list))
    for chain_id in chains_ids:
        structure_ranges: Dict[SecondaryStructure, List[Tuple[int, int]]] = defaultdict(list)
        for s in struct_conf.filter({'beg_auth_asym_id': chain_id, 'end_auth_asym_id': chain_id}).value():
            structure_ranges[s['ss']].append((int(s['beg_label_seq_id']), int(s['end_label_seq_id'])))
        chain_structure_ranges[chain_id] = structure_ranges
    return chain_structure_ranges


def get_sheet_ranges(cif: PDBxFile, chains_ids: List[str]) -> Dict[str, List[Tuple[int, int]]]:
    struct_sheet_range = get_category_as_chain(cif, 'struct_sheet_range')
    chain_sheet_ranges: Dict[str, Tuple[int, int]] = dict()
    for chain_id in chains_ids:
        sheet_ranges: List[Tuple[str, str]] = []
        for s in struct_sheet_range.filter({'beg_auth_asym_id': chain_id, 'end_auth_asym_id': chain_id}).value():
            sheet_ranges.append((int(s['beg_label_seq_id']), int(s['end_label_seq_id'])))
        chain_sheet_ranges[chain_id] = sheet_ranges
    return chain_sheet_ranges


def classify_chain_structures(chain_residues: Dict[str, Set[int]], chain_structure_ranges: DefaultDict[str, DefaultDict[SecondaryStructure, List[Tuple[int, int]]]], sheet_structure_ranges: Dict[str, List[Tuple[int, int]]]):
    chain_classification: Dict[str, SecondaryStructure] = dict()
    for chain_id, sheet_ranges in sheet_structure_ranges.items():
        chain_structure_ranges[chain_id][SecondaryStructure.Extended].extend(sheet_ranges)
    for chain_id, residues in chain_residues.items():
        residue_count = len(residues)
        for structure_type, ranges in chain_structure_ranges[chain_id].items():
            residues_of_structure = len({r for r in residues for rr in ranges if rr[0] <= r <= rr[1]})
            if residues_of_structure / residue_count >= STRUCTURED_RESIDUES_THRESHOLD:
                chain_classification[chain_id] = structure_type
                break
    return invert_by(chain_classification)


def get_fasta(pdb_id: str, chain_sequences: Dict[str, str]):
    return ''.join([f'>{pdb_id}:{chain_id}\n{sequence}\n' for chain_id, sequence in chain_sequences.items()])


def process_cif_file(path: Path):
    with gzip.open(path, mode='rt') as file:
        cif = PDBxFile.read(file)
    pdb_id = cif['entry']['id']
    chain_sequences = get_chain_sequences(cif)
    chain_residues = get_chain_residues(cif, chain_sequences.keys())
    chain_structure_ranges = get_structure_ranges(cif, chain_sequences.keys())
    chain_sheet_ranges = get_sheet_ranges(cif, chain_sequences.keys())
    chain_classification = classify_chain_structures(chain_residues, chain_structure_ranges, chain_sheet_ranges)
    return map_values(chain_classification, lambda c: get_fasta(pdb_id, dict([cs for cs in chain_sequences.items() if cs[0] in c])))


def main():
    parser = ArgumentParser()
    parser.add_argument('-c', '--cif-parent-path', type=Path, required=True,
                        help='Input PDB mmCIF parent directory path')
    parser.add_argument('-f', '--fasta-dir', type=Path, required=True, help='Directory in which output FASTA files will be created')
    parser.add_argument('-t', '--threads', type=int, help='Maximum number of threads')
    args = parser.parse_args()

    paths = list(args.cif_parent_path.glob('*/*.cif.gz'))

    with ThreadPoolExecutor(args.threads) as executor, tqdm(total=len(paths), unit='', miniters=2) as progress:
        futures = [executor.submit(process_cif_file, path) for path in paths]
        for _ in as_completed(futures):
            progress.update()

    results = [future.result() for future in futures]

    for ss in SecondaryStructure:
        output_path = args.fasta_dif / f'pdb_{ss.value}.fasta'
        filtered = py_(results).map(lambda x: x.get(ss)).filter()
        count = filtered.size().value()
        content = filtered.join().value()
        with open(output_path, mode='wt') as file:
            file.write(content)
        print(f'✅ {output_path} : {count} sequences')


if __name__ == '__main__':
    main()