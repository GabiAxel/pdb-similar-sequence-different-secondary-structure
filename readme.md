## Prerequisites
* Python 3.8+
* Rsync
* [Diamond](https://github.com/bbuchfink/diamond)

## Procedure

Download PDB mmCIF. In the following command the target is `mmCIF` directory under the user home directory.
```shell
rsync -rlpt -v -z --delete --port=33444 rsync.rcsb.org::ftp_data/structures/divided/mmCIF/ ~/mmCIF
```

Install Python dependencies.
Run the following command from the directory containing the script.

```shell
pip install -r requirements.txt
```

Filter only polypeptide chains with at least 70% of their residues modeled and at least 50% of their residues making one secondary structure type.
Produce two FASTA files named `pdb_extended.fasta` and `pdb_helical.fasta`, each containing the chains classified as the same secondary structure.

| Short Parameter | Long Parameter     | Description                                            |
|-----------------|--------------------|--------------------------------------------------------|
| `-c`            | `--cif-parent-dir` | Path of mmCIF directory (downloaded in the first step) |
| `-f`            | `--fasta-dir`      | Path in which the FASTA files will be created          |
| `-t`            | `--threads`        | Number of threads (optional)                           |
| `-h`            | `--help`           | Show help message                                      |

```shell
python classify_chains_by_secondary_structure.py -c ~/mmCIF -f ~/
```

Build a Diamond database from the helical chains FASTA file produced in the previous step.

```shell
diamond makedb --in ~/pdb_helical.fasta -d ~/pdb_helical
```

Run a Diamond BLASTp for all extended chains against the helical chains database with minimum 80% identity and coverage.
The result tab-delimited file will contain the following columns:
1. Extended chain identifier
2. Extended chain coverage %
3. Helical chain identifier
4. Helical chain coverage %
5. Sequence identity %

```shell
diamond blastp -d pdb_helical -q pdb_extended.fasta -o pdb_blastp_80.tsv -f 6 qseqid qcovhsp sseqid scovhsp pident -k0 --id 80 --query-cover 80 --subject-cover 80
```