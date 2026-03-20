#!/usr/bin/env python3
"""
db2_to_binary_recursive.py - Recursively convert DB2 files and archives to db2bin format

This script is a standalone wrapper that:
1. Recursively finds `.db2`, `.db2.gz`, and `.db2.tgz` (or `.tar.gz`) files in an input directory.
2. Converts them to `.db2bin`, `.db2bin.gz`, and `.db2bin.tgz` respectively.
3. For archives, it processes the internals entirely in-memory using `tarfile` and `io.BytesIO`,
   avoiding any intermediate extraction to disk.
4. Replicates the relative directory structure if an output directory is provided.

Usage:
    python3 db2_to_binary_recursive.py /path/to/input/dir /path/to/output/dir

Author: Phong Lam, Uppsala University (Extended by AI)
Date: January 2026 (Updated March 2026)
"""

import struct
import gzip
import argparse
import tarfile
import io
try:   
    import zstandard
except ImportError:
    print("zstandard not installed, please install it using 'pip install zstandard'")
    exit(1)
from pathlib import Path

# Binary format constants
DB2BIN_MAGIC = b'DB2BIN'
DB2BIN_VERSION = 1


def parse_db2_stream(f_iter):
    """
    Parse an iterable of text DB2 lines (or decoded file stream) into structured data.
    Yields a dict for each molecule to keep memory footprint low.
    """
    data = {
        'refcod': '', 'protcod': '', 'smiles': '', 'molname': '', 'arbitrary': [],
        'total_atoms': 0, 'total_bonds': 0, 'total_coords': 0, 'total_confs': 0,
        'total_sets': 0, 'rigid_heavy_count': 0, 'total_clusters': 0, 'mlines': 0,
        'total_charge': 0.0, 'total_polar_solv': 0.0, 'total_apolar_solv': 0.0,
        'total_solv': 0.0, 'total_surfarea': 0.0,
        'atom_num': [], 'atom_name': [], 'atom_type': [], 'atom_vdwtype': [],
        'atom_color': [], 'atom_charge': [], 'atom_polsolv': [], 'atom_apolsolv': [],
        'atom_totalsolv': [], 'atom_surfarea': [], 'bond_num': [], 'bond_start': [],
        'bond_end': [], 'bond_type': [], 'coord_index': [], 'coords': [],
        'rigid_colors': [], 'rigid_coords': [], 'conf_coord': [], 'sets': [],
        'clusters': [],
    }

    line_idx = 0
    m_lines_remaining = 0
    atoms_seen = 0
    bonds_seen = 0
    coords_seen = 0
    rigid_seen = 0
    confs_seen = 0
    current_set = None
    set_lines_remaining = 0
    current_cluster = None
    cluster_match_remaining = 0
    for raw_line in f_iter:
        line = raw_line.strip()
        if isinstance(line, bytes):
            line = line.decode('utf-8', errors='replace')
            
        if not line:
            continue

        section = line[0] if line else ''

        if section == 'M':
            parts = line.split()
            if m_lines_remaining == 0:
                data['refcod'] = parts[1] if len(parts) > 1 else ''
                data['protcod'] = parts[2] if len(parts) > 2 else ''
                if len(parts) >= 11:
                    data['total_atoms'] = int(parts[3])
                    data['total_bonds'] = int(parts[4])
                    data['total_coords'] = int(parts[5])
                    data['total_confs'] = int(parts[6])
                    data['total_sets'] = int(parts[7])
                    data['rigid_heavy_count'] = int(parts[8])
                    data['mlines'] = int(parts[9])
                    data['total_clusters'] = int(parts[10])
                m_lines_remaining = data['mlines'] - 1
            elif m_lines_remaining == data['mlines'] - 1:
                if len(parts) >= 6:
                    data['total_charge'] = float(parts[1])
                    data['total_polar_solv'] = float(parts[2])
                    data['total_apolar_solv'] = float(parts[3])
                    data['total_solv'] = float(parts[4])
                    data['total_surfarea'] = float(parts[5])
                m_lines_remaining -= 1
            elif m_lines_remaining == data['mlines'] - 3:
                 data['molname'] = line[2:].strip() if len(line) > 2 else ''
                 m_lines_remaining -= 1
            elif m_lines_remaining == data['mlines'] - 2:
                data['smiles'] = line[2:].strip() if len(line) > 2 else ''
                m_lines_remaining -= 1
            else:
                arb_line = line[2:].strip() if len(line) > 2 else ''
                data['arbitrary'].append(arb_line)
                m_lines_remaining -= 1
        elif section == 'A':
            parts = line.split()
            if len(parts) >= 11:
                data['atom_num'].append(int(parts[1]))
                data['atom_name'].append(parts[2][:4])
                data['atom_type'].append(parts[3][:5])
                data['atom_vdwtype'].append(int(parts[4]))
                data['atom_color'].append(int(parts[5]))
                data['atom_charge'].append(float(parts[6]))
                data['atom_polsolv'].append(float(parts[7]))
                data['atom_apolsolv'].append(float(parts[8]))
                data['atom_totalsolv'].append(float(parts[9]))
                data['atom_surfarea'].append(float(parts[10]))
            atoms_seen += 1
        elif section == 'B':
            parts = line.split()
            if len(parts) >= 5:
                data['bond_num'].append(int(parts[1]))
                data['bond_start'].append(int(parts[2]))
                data['bond_end'].append(int(parts[3]))
                data['bond_type'].append(parts[4][:2])
            bonds_seen += 1
        elif section == 'X':
            parts = line.split()
            if len(parts) >= 6:
                data['coord_index'].append([int(parts[1]), int(parts[2]), int(parts[3])])
                data['coords'].append([float(parts[4]), float(parts[5]), float(parts[6])])
            coords_seen += 1
        elif section == 'R':
            parts = line.split()
            if len(parts) >= 5:
                data['rigid_colors'].append(int(parts[2]))
                data['rigid_coords'].append([
                    float(parts[3]), float(parts[4]), float(parts[5]) if len(parts) > 5 else float(parts[4])
                ])
            rigid_seen += 1
        elif section == 'C':
            parts = line.split()
            if len(parts) >= 3:
                data['conf_coord'].append([int(parts[2]), int(parts[3]) if len(parts) > 3 else int(parts[2])])
            confs_seen += 1
        elif section == 'S':
            parts = line.split()
            if len(parts) >= 5:
                if set_lines_remaining == 0:
                    set_num, total_lines, conf_count = int(parts[1]), int(parts[2]), int(parts[3])
                    current_set = {
                        'set_num': set_num, 'total_lines': total_lines, 'conf_count': conf_count,
                        'broken': int(parts[4]) if len(parts) > 4 else 0,
                        'energy': float(parts[6]) if len(parts) > 6 else 0.0,
                        'max_strain': float(parts[7]) if len(parts) > 7 else 0.0,
                        'confs': []
                    }
                    set_lines_remaining = total_lines
                    if set_lines_remaining <= 0:
                        data['sets'].append(current_set)
                        current_set = None
                else:
                    confs_on_line = int(parts[3])
                    if current_set is not None:
                        for i in range(confs_on_line):
                            if len(parts) > 4 + i:
                                current_set['confs'].append(int(parts[4 + i]))
                        set_lines_remaining -= 1
                        if set_lines_remaining <= 0:
                            data['sets'].append(current_set)
                            current_set = None
        elif section == 'D':
            parts = line.split()
            if cluster_match_remaining == 0 and len(parts) >= 6:
                current_cluster = {
                    'set_start': int(parts[2]), 'set_end': int(parts[3]),
                    'match_start': int(parts[5]), 'match_end': int(parts[6]) if len(parts) > 6 else int(parts[5]),
                    'add_match_colors': [], 'add_match_coords': []
                }
                cluster_match_remaining = int(parts[4])
                if cluster_match_remaining == 0:
                    data['clusters'].append(current_cluster)
                    current_cluster = None
            elif cluster_match_remaining > 0 and len(parts) >= 4:
                if current_cluster:
                    current_cluster['add_match_colors'].append(int(parts[2]))
                    current_cluster['add_match_coords'].append([
                        float(parts[3]), float(parts[4]) if len(parts) > 4 else 0.0,
                        float(parts[5]) if len(parts) > 5 else 0.0
                    ])
                cluster_match_remaining -= 1
                if cluster_match_remaining == 0 and current_cluster:
                    data['clusters'].append(current_cluster)
                    current_cluster = None
        elif section == 'E':
            yield data
            data = {
                'refcod': '', 'protcod': '', 'smiles': '', 'molname': '', 'arbitrary': [],
                'total_atoms': 0, 'total_bonds': 0, 'total_coords': 0, 'total_confs': 0,
                'total_sets': 0, 'rigid_heavy_count': 0, 'total_clusters': 0, 'mlines': 0,
                'total_charge': 0.0, 'total_polar_solv': 0.0, 'total_apolar_solv': 0.0,
                'total_solv': 0.0, 'total_surfarea': 0.0,
                'atom_num': [], 'atom_name': [], 'atom_type': [], 'atom_vdwtype': [],
                'atom_color': [], 'atom_charge': [], 'atom_polsolv': [], 'atom_apolsolv': [],
                'atom_totalsolv': [], 'atom_surfarea': [], 'bond_num': [], 'bond_start': [],
                'bond_end': [], 'bond_type': [], 'coord_index': [], 'coords': [],
                'rigid_colors': [], 'rigid_coords': [], 'conf_coord': [], 'sets': [],
                'clusters': [],
            }
            m_lines_remaining = 0
            atoms_seen = 0
            bonds_seen = 0
            coords_seen = 0
            rigid_seen = 0
            confs_seen = 0
            current_set = None
            current_cluster = None
            cluster_match_remaining = 0


def write_single_molecule(f, data):
    """Write a single molecule to an already-open binary file-like object."""
    f.write(DB2BIN_MAGIC)
    f.write(struct.pack('i', DB2BIN_VERSION))
    f.write(struct.pack('i', data['total_atoms']))
    f.write(struct.pack('i', data['total_bonds']))
    f.write(struct.pack('i', data['total_coords']))
    f.write(struct.pack('i', data['total_confs']))
    f.write(struct.pack('i', len(data['sets'])))
    f.write(struct.pack('i', data['rigid_heavy_count']))
    f.write(struct.pack('i', len(data['clusters'])))
    f.write(struct.pack('i', data['mlines']))
    
    f.write(struct.pack('f', data['total_charge']))
    f.write(struct.pack('f', data['total_polar_solv']))
    f.write(struct.pack('f', data['total_apolar_solv']))
    f.write(struct.pack('f', data['total_solv']))
    f.write(struct.pack('f', data['total_surfarea']))
    
    refcod = data['refcod'][:16].ljust(16, '\x00')
    f.write(refcod.encode('ascii'))
    protcod = data['protcod'][:9].ljust(9, ' ')
    f.write(protcod.encode('ascii'))
    smiles = data['smiles'][:78].ljust(78, ' ')
    f.write(smiles.encode('ascii', errors='replace'))
    molname = data['molname'][:78].ljust(78, ' ')
    f.write(molname.encode('ascii', errors='replace'))
    
    arb_count = len(data['arbitrary'])
    f.write(struct.pack('i', arb_count))
    for arb_line in data['arbitrary']:
        arb_bytes = arb_line[:78].ljust(78, ' ')
        f.write(arb_bytes.encode('ascii', errors='replace'))
    
    for val in data['atom_num']: f.write(struct.pack('i', val))
    for name in data['atom_name']: f.write(name[:4].ljust(4, ' ').encode('ascii'))
    for atype in data['atom_type']: f.write(atype[:5].ljust(5, ' ').encode('ascii'))
    for arr in ['atom_vdwtype', 'atom_color']:
        for val in data[arr]: f.write(struct.pack('i', val))
    for arr in ['atom_charge', 'atom_polsolv', 'atom_apolsolv', 'atom_totalsolv', 'atom_surfarea']:
        for val in data[arr]: f.write(struct.pack('f', val))
    
    for val in data['bond_num']: f.write(struct.pack('i', val))
    for val in data['bond_start']: f.write(struct.pack('i', val))
    for val in data['bond_end']: f.write(struct.pack('i', val))
    for btype in data['bond_type']: f.write(btype[:2].ljust(2, ' ').encode('ascii'))
    
    for i, (idx, coord) in enumerate(zip(data['coord_index'], data['coords'])):
        for val in idx: f.write(struct.pack('i', val))
        for val in coord: f.write(struct.pack('f', val))
    
    for i, (color, coord) in enumerate(zip(data['rigid_colors'], data['rigid_coords'])):
        f.write(struct.pack('i', color))
        for val in coord: f.write(struct.pack('f', val))
    
    for conf in data['conf_coord']:
        f.write(struct.pack('i', conf[0]))
        f.write(struct.pack('i', conf[1] if len(conf) > 1 else conf[0]))
    
    for s in data['sets']:
        f.write(struct.pack('i', len(s['confs'])))
        f.write(struct.pack('i', s['broken']))
        f.write(struct.pack('f', s['energy']))
        f.write(struct.pack('f', s['max_strain']))
        for conf in s['confs']: f.write(struct.pack('i', conf))
    
    for c in data['clusters']:
        f.write(struct.pack('i', c['set_start']))
        f.write(struct.pack('i', c['set_end']))
        f.write(struct.pack('i', c['match_start']))
        f.write(struct.pack('i', c['match_end']))
    for c in data['clusters']:
        for i, color in enumerate(c['add_match_colors']):
            f.write(struct.pack('i', color))
            for val in c['add_match_coords'][i]: f.write(struct.pack('f', val))


def convert_stream_to_binary_bytes(f_iter, compress=False):
    """Parse iterable stream, write to binary in-memory buffer, and return the bytes."""
    buf = io.BytesIO()
    
    if compress:
        with gzip.GzipFile(fileobj=buf, mode='wb') as gz:
            for data in parse_db2_stream(f_iter):
                write_single_molecule(gz, data)
    else:
        for data in parse_db2_stream(f_iter):
            write_single_molecule(buf, data)
            
    return buf.getvalue()


def process_tarball(input_tar_path, output_tar_path):
    """
    Open a tarball, iterate members, convert .db2 internally,
    and write to new tarball entirely in-memory.
    """
    print(f"Processing Archive: {input_tar_path} -> {output_tar_path}")
    
    # Create ZstdCompressor
    cctx = zstandard.ZstdCompressor(level=10)
    
    with tarfile.open(input_tar_path, 'r|gz') as in_tar, \
        open(output_tar_path, 'wb') as out_f, \
        cctx.stream_writer(out_f) as zstd_stream, \
        tarfile.open(fileobj=zstd_stream, mode='w|') as out_tar:
         
        for member in in_tar:
            if member.isfile() and (member.name.endswith('.db2') or member.name.endswith('.db2.gz')):
                # Skip if it is already a binary file (e.g. .db2bin.gz inside a .db2.tgz? edge case)
                if member.name.endswith(('.db2bin', '.db2bin.gz')):
                    continue
                f = in_tar.extractfile(member)
                if f is None:
                    continue
                
                # Read contents
                is_gz = member.name.endswith('.gz')
                try:
                    if is_gz:
                        with gzip.open(f, 'rt') as gz_f:
                            bin_data = convert_stream_to_binary_bytes(gz_f, compress=True)
                        new_name = member.name.replace('.db2.gz', '.db2bin.gz')
                    else:
                        bin_data = convert_stream_to_binary_bytes(f, compress=False)
                        new_name = member.name.replace('.db2', '.db2bin')
                    
                    # Add to new tar
                    new_member = tarfile.TarInfo(name=new_name)
                    new_member.size = len(bin_data)
                    new_member.mtime = member.mtime
                    
                    out_tar.addfile(new_member, io.BytesIO(bin_data))
                    # print(f"    - Converted in-memory: {member.name} -> {new_name}")
                    
                except Exception as e:
                    print(f"    - [ERROR] Failed converting inner file {member.name}: {e}")


def process_file(input_path, output_dir, relative_to=None):
    """
    Process a single file (db2, db2.gz, or db2.tgz), routing to proper handler.
    If `relative_to` is provided, re-create the directory structure inside `output_dir`.
    """
    input_path = Path(input_path)
    
    # Determine the output path preserving directory structure
    if relative_to:
        rel_path = input_path.relative_to(relative_to)
        out_path = Path(output_dir) / rel_path
    else:
        out_path = Path(output_dir) / input_path.name
        
    out_path.parent.mkdir(parents=True, exist_ok=True)
    
    filename = input_path.name.lower()
    
    try:
        if filename.endswith(('.tgz', '.tar.gz')):
            if filename.endswith('.tar.gz'):
                out_name = input_path.name[:-7] + ".db2bin.tar.zst" # .tar.gz -> .db2bin.tar.zst
            else:
                out_name = input_path.name[:-4] + ".db2bin.tar.zst" # .db2.tgz -> .db2bin.tar.zst
            final_out_path = out_path.with_name(out_name)
            process_tarball(input_path, final_out_path)
            
        elif filename.endswith('.db2.gz'):
            out_name = input_path.name.replace('.db2.gz', '.db2bin.gz')
            final_out_path = out_path.with_name(out_name)
            
            with gzip.open(input_path, 'rt') as f:
                bin_data = convert_stream_to_binary_bytes(f, compress=True)
            with open(final_out_path, 'wb') as out_f:
                out_f.write(bin_data)
                
            print(f"  Converted {input_path} -> {final_out_path}")
                
        elif filename.endswith('.db2'):
            out_name = input_path.name.replace('.db2', '.db2bin')
            final_out_path = out_path.with_name(out_name)
            
            with open(input_path, 'rt') as f:
                bin_data = convert_stream_to_binary_bytes(f, compress=False)
            with open(final_out_path, 'wb') as out_f:
                out_f.write(bin_data)
                
            print(f"  Converted {input_path} -> {final_out_path}")
    except BaseException as e:
        print(f"  [ERROR] Failed to process {input_path}: {e}")




def main():
    parser = argparse.ArgumentParser(
        description='Recursively convert text DB2 files / archives to binary db2bin format without intermediate extraction.'
    )
    parser.add_argument('input_dir', help='Input directory to search')
    parser.add_argument('output_dir', help='Output directory where results will be saved')
    
    args = parser.parse_args()
    
    input_dir = Path(args.input_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    
    if not input_dir.is_dir():
        print(f"Error: {input_dir} is not a directory.")
        return
        
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Scanning {input_dir} recursively for target files...")
    
    # Identify files
    targets = []
    targets.extend(input_dir.rglob("*.db2"))
    targets.extend(input_dir.rglob("*.db2.gz"))
    targets.extend(input_dir.rglob("*.db2.tgz"))
    targets.extend(input_dir.rglob("*.tar.gz"))
    
    # Filter out `.db2bin` so we don't process binary files by mistake
    targets = [t for t in targets if '.db2bin' not in t.name]
    
    if not targets:
        print("No .db2, .db2.gz, .db2.tgz, or .tar.gz files found.")
        return
        
    print(f"Found {len(targets)} files to process.")
    
    for t in targets:
        process_file(t, output_dir, relative_to=input_dir)

if __name__ == '__main__':
    main()
