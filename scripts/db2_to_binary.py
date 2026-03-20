#!/usr/bin/env python3
"""
db2_to_binary.py - Convert text DB2.gz files to binary DB2BIN.gz format

This script converts DOCK3 text-based DB2 files to the high-performance
binary format (db2bin.gz), which can be 5-20x faster to parse.

Usage:
    python3 db2_to_binary.py input.db2.gz output.db2bin.gz
    python3 db2_to_binary.py --dir /path/to/db2s/ /path/to/output/

Binary Format (db2bin.gz):
    - Header: magic "DB2BIN", version, counts, metadata
    - Arrays: contiguous binary arrays for fast I/O
    - All data is gzip compressed

Author: Phong Lam, Uppsala University
Date: January 2026
"""

import struct
import gzip
import argparse
from pathlib import Path


# Binary format constants
DB2BIN_MAGIC = b'DB2BIN'
DB2BIN_VERSION = 1


def parse_db2_text(input_path):
    """
    Parse a text DB2.gz file into structured data.
    
    Returns a dict with all the arrays needed for binary format.
    """
    data = {
        'refcod': '',
        'protcod': '',
        'smiles': '',
        'molname': '',
        'arbitrary': [],  # Lines 5+ of M section
        'total_atoms': 0,
        'total_bonds': 0,
        'total_coords': 0,
        'total_confs': 0,
        'total_sets': 0,
        'rigid_heavy_count': 0,
        'total_clusters': 0,
        'mlines': 0,
        'total_charge': 0.0,
        'total_polar_solv': 0.0,
        'total_apolar_solv': 0.0,
        'total_solv': 0.0,
        'total_surfarea': 0.0,
        # Arrays (populated during parsing)
        'atom_num': [],
        'atom_name': [],
        'atom_type': [],
        'atom_vdwtype': [],
        'atom_color': [],
        'atom_charge': [],
        'atom_polsolv': [],
        'atom_apolsolv': [],
        'atom_totalsolv': [],
        'atom_surfarea': [],
        'bond_num': [],
        'bond_start': [],
        'bond_end': [],
        'bond_type': [],
        'coord_index': [],  # 3 x total_coords
        'coords': [],       # 3 x total_coords
        'rigid_colors': [],
        'rigid_coords': [], # 3 x rigid_count
        'conf_coord': [],   # 2 x total_confs (start, end)
        'sets': [],         # List of dicts with set data
        'clusters': [],     # List of dicts with cluster data
    }
    
    open_func = gzip.open if str(input_path).endswith('.gz') else open
    
    with open_func(input_path, 'rt') as f:
        lines = f.readlines()
    
    line_idx = 0
    m_lines_remaining = 0
    current_section = None
    atoms_seen = 0
    bonds_seen = 0
    coords_seen = 0
    rigid_seen = 0
    confs_seen = 0
    current_set = None
    set_lines_remaining = 0
    set_confs_remaining = 0
    current_cluster = None
    cluster_match_remaining = 0
    molecules = []
    
    while line_idx < len(lines):
        line = lines[line_idx].strip()
        if not line:
            line_idx += 1
            continue
            
        section = line[0] if line else ''
        
        if section == 'M':
            parts = line.split()
            if m_lines_remaining == 0:
                # First M line: header counts
                # Format: M refcod protcod atoms bonds coords confs sets rigid mlines clusters
                # Example: M TS_01F_0tR fake 37 39 349 78 31 8 5 0
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
                # Second M line: charge and solvation
                # Format: M charge polar_solv apolar_solv total_solv surfarea
                # Example: M +0.0000 -13.330 +5.770 -7.560 321.870
                if len(parts) >= 6:
                    data['total_charge'] = float(parts[1])
                    data['total_polar_solv'] = float(parts[2])
                    data['total_apolar_solv'] = float(parts[3])
                    data['total_solv'] = float(parts[4])
                    data['total_surfarea'] = float(parts[5])
                m_lines_remaining -= 1
            elif m_lines_remaining == data['mlines'] - 2:
                # Third M line: smiles
                # Format: M smiles_string
                data['smiles'] = line[2:].strip() if len(line) > 2 else ''
                m_lines_remaining -= 1
            elif m_lines_remaining == data['mlines'] - 3:
                # Fourth M line: long name (molname)
                # Format: M long_name
                data['molname'] = line[2:].strip() if len(line) > 2 else ''
                m_lines_remaining -= 1
            else:
                # Lines 5+ are arbitrary lines (rig_frag_code, etc.)
                arb_line = line[2:].strip() if len(line) > 2 else ''
                data['arbitrary'].append(arb_line)
                m_lines_remaining -= 1
                
        elif section == 'A':
            # Atom line: A num name type vdwtype color charge polsolv apolsolv totalsolv surfarea
            # Example: A   1 C1   C.3    5  7   +0.0300     -0.040     +0.580     +0.540    11.000
            # Index:   0   1  2    3     4  5      6           7          8          9        10
            parts = line.split()
            if len(parts) >= 11:
                data['atom_num'].append(int(parts[1]))
                data['atom_name'].append(parts[2][:4])  # Max 4 chars
                data['atom_type'].append(parts[3][:5])  # Max 5 chars (sybyl type)
                data['atom_vdwtype'].append(int(parts[4]))
                data['atom_color'].append(int(parts[5]))
                data['atom_charge'].append(float(parts[6]))
                data['atom_polsolv'].append(float(parts[7]))
                data['atom_apolsolv'].append(float(parts[8]))
                data['atom_totalsolv'].append(float(parts[9]))
                data['atom_surfarea'].append(float(parts[10]))
            atoms_seen += 1
            
        elif section == 'B':
            # Bond line: B num start_atom end_atom type
            # Example: B   1   1   2 1   (single bond between atoms 1 and 2)
            # Bond types: 1=single, 2=double, 3=triple, ar=aromatic, am=amide
            parts = line.split()
            if len(parts) >= 5:
                data['bond_num'].append(int(parts[1]))
                data['bond_start'].append(int(parts[2]))
                data['bond_end'].append(int(parts[3]))
                data['bond_type'].append(parts[4][:2])  # Max 2 chars
            bonds_seen += 1
            
        elif section == 'X':
            # Coordinate line: X coord_num atom_num conf_num x y z
            parts = line.split()
            if len(parts) >= 6:
                data['coord_index'].append([
                    int(parts[1]),   # coord_num
                    int(parts[2]),   # atom_num  
                    int(parts[3])    # conf_num
                ])
                data['coords'].append([
                    float(parts[4]),  # x
                    float(parts[5]),  # y
                    float(parts[6])   # z
                ])
            coords_seen += 1
            
        elif section == 'R':
            # Rigid sphere: R num color x y z
            parts = line.split()
            if len(parts) >= 5:
                data['rigid_colors'].append(int(parts[2]))
                data['rigid_coords'].append([
                    float(parts[3]),
                    float(parts[4]),
                    float(parts[5]) if len(parts) > 5 else float(parts[4])
                ])
            rigid_seen += 1
            
        elif section == 'C':
            # Conformation: C conf_num start end
            parts = line.split()
            if len(parts) >= 3:
                data['conf_coord'].append([
                    int(parts[2]),  # start
                    int(parts[3]) if len(parts) > 3 else int(parts[2])  # end
                ])
            confs_seen += 1
            
        elif section == 'S':
            # Set line - multi-line format
            # First line:  S set_num total_lines conf_count broken hydro energy max_strain
            # Next lines:  S set_num line_num confs_on_line conf1 conf2 conf3 ...
            # Line numbers count DOWN from total_lines
            # Example: S      1      2  11 0 0      +0.000      +0.000  (header, 2 lines)
            #          S      1      1 8      1      2      3     14     15     23     34     35 (confs)
            #          S      1      2 3    145    146    312  (more confs, if line_num > 1)
            parts = line.split()
            if len(parts) >= 5:
                if set_lines_remaining == 0:
                    # First line of a new set (header) - has metadata
                    set_num = int(parts[1])
                    total_lines = int(parts[2])
                    conf_count = int(parts[3])
                    current_set = {
                        'set_num': set_num,
                        'total_lines': total_lines,
                        'conf_count': conf_count,
                        'broken': int(parts[4]) if len(parts) > 4 else 0,
                        'energy': float(parts[6]) if len(parts) > 6 else 0.0,
                        'max_strain': float(parts[7]) if len(parts) > 7 else 0.0,
                        'confs': []
                    }
                    # The number in parts[2] is the total lines (data lines only)
                    set_lines_remaining = total_lines
                    
                    if set_lines_remaining <= 0:
                        # Single-line set (header only) - unusual if lines>0
                        data['sets'].append(current_set)
                        current_set = None
                else:
                    # Continuation line with conf numbers  
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
            # Cluster or additional matching sphere
            parts = line.split()
            if cluster_match_remaining == 0 and len(parts) >= 6:
                # Cluster line: D cluster_num set_start set_end num_matches match_start match_end
                current_cluster = {
                    'set_start': int(parts[2]),
                    'set_end': int(parts[3]),
                    'match_start': int(parts[5]),
                    'match_end': int(parts[6]) if len(parts) > 6 else int(parts[5]),
                    'add_match_colors': [],
                    'add_match_coords': []
                }
                cluster_match_remaining = int(parts[4])
                if cluster_match_remaining == 0:
                    data['clusters'].append(current_cluster)
                    current_cluster = None
            elif cluster_match_remaining > 0 and len(parts) >= 4:
                # Additional matching sphere
                if current_cluster:
                    current_cluster['add_match_colors'].append(int(parts[2]))
                    current_cluster['add_match_coords'].append([
                        float(parts[3]),
                        float(parts[4]) if len(parts) > 4 else 0.0,
                        float(parts[5]) if len(parts) > 5 else 0.0
                    ])
                cluster_match_remaining -= 1
                if cluster_match_remaining == 0 and current_cluster:
                    data['clusters'].append(current_cluster)
                    current_cluster = None
                    
        elif section == 'E':
            # End of molecule - add to list and continue
            molecules.append(data)
            # Reset for next molecule
            data = {
                'refcod': '',
                'protcod': '',
                'smiles': '',
                'molname': '',
                'arbitrary': [],
                'total_atoms': 0,
                'total_bonds': 0,
                'total_coords': 0,
                'total_confs': 0,
                'total_sets': 0,
                'rigid_heavy_count': 0,
                'total_clusters': 0,
                'mlines': 0,
                'total_charge': 0.0,
                'total_polar_solv': 0.0,
                'total_apolar_solv': 0.0,
                'total_solv': 0.0,
                'total_surfarea': 0.0,
                'atom_num': [],
                'atom_name': [],
                'atom_type': [],
                'atom_vdwtype': [],
                'atom_color': [],
                'atom_charge': [],
                'atom_polsolv': [],
                'atom_apolsolv': [],
                'atom_totalsolv': [],
                'atom_surfarea': [],
                'bond_num': [],
                'bond_start': [],
                'bond_end': [],
                'bond_type': [],
                'coord_index': [],
                'coords': [],
                'rigid_colors': [],
                'rigid_coords': [],
                'conf_coord': [],
                'sets': [],
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
            
        line_idx += 1
    
    return molecules


def write_db2_binary(data, output_path):
    """
    Write parsed data to gzipped binary format.
    """
    open_func = gzip.open if str(output_path).endswith('.gz') else open
    
    with open_func(output_path, 'wb') as f:
        # Write magic (6 bytes)
        f.write(DB2BIN_MAGIC)
        
        # Write version (4 bytes int)
        f.write(struct.pack('i', DB2BIN_VERSION))
        
        # Write header integers (8 x 4 bytes)
        f.write(struct.pack('i', data['total_atoms']))
        f.write(struct.pack('i', data['total_bonds']))
        f.write(struct.pack('i', data['total_coords']))
        f.write(struct.pack('i', data['total_confs']))
        f.write(struct.pack('i', len(data['sets'])))  # actual sets parsed
        f.write(struct.pack('i', data['rigid_heavy_count']))
        f.write(struct.pack('i', len(data['clusters'])))  # actual clusters
        f.write(struct.pack('i', data['mlines']))
        
        # Write header reals (5 x 4 bytes)
        f.write(struct.pack('f', data['total_charge']))
        f.write(struct.pack('f', data['total_polar_solv']))
        f.write(struct.pack('f', data['total_apolar_solv']))
        f.write(struct.pack('f', data['total_solv']))
        f.write(struct.pack('f', data['total_surfarea']))
        
        # Write refcod (16 bytes, padded)
        refcod_bytes = data['refcod'].encode('ascii')[:16].ljust(16, b'\x00')
        f.write(refcod_bytes)
        
        # Write protcod (9 bytes, space-padded)
        protcod_bytes = data['protcod'].encode('ascii')[:9].ljust(9, b' ')
        f.write(protcod_bytes)
        
        # Write smiles (78 bytes, space-padded)
        smiles_bytes = data['smiles'].encode('ascii', errors='replace')[:78].ljust(78, b' ')
        f.write(smiles_bytes)
        
        # Write molname (78 bytes, space-padded)
        molname_bytes = data['molname'].encode('ascii', errors='replace')[:78].ljust(78, b' ')
        f.write(molname_bytes)
        
        # Write arbitrary lines count and data
        arb_count = len(data['arbitrary'])
        f.write(struct.pack('i', arb_count))
        for arb_line in data['arbitrary']:
            arb_bytes = arb_line.encode('ascii', errors='replace')[:78].ljust(78, b' ')
            f.write(arb_bytes)
        
        # Write atom arrays
        for val in data['atom_num']:
            f.write(struct.pack('i', val))
        for name in data['atom_name']:
            f.write(name[:4].ljust(4, ' ').encode('ascii'))
        for atype in data['atom_type']:
            f.write(atype[:5].ljust(5, ' ').encode('ascii'))
        for arr in ['atom_vdwtype', 'atom_color']:
            for val in data[arr]:
                f.write(struct.pack('i', val))
        for arr in ['atom_charge', 'atom_polsolv', 'atom_apolsolv', 
                    'atom_totalsolv', 'atom_surfarea']:
            for val in data[arr]:
                f.write(struct.pack('f', val))
        
        # Write bond arrays
        for val in data['bond_num']:
            f.write(struct.pack('i', val))
        for val in data['bond_start']:
            f.write(struct.pack('i', val))
        for val in data['bond_end']:
            f.write(struct.pack('i', val))
        for btype in data['bond_type']:
            f.write(btype[:2].ljust(2, ' ').encode('ascii'))
        
        # Write coord data - INTERLEAVED format: (3 ints + 3 floats) per coord
        for i, (idx, coord) in enumerate(zip(data['coord_index'], data['coords'])):
            # Write coord_index (3 ints: coord_num, atom_num, conf_num)
            for val in idx:
                f.write(struct.pack('i', val))
            # Write coords (3 floats: x, y, z)
            for val in coord:
                f.write(struct.pack('f', val))
        
        # Write rigid spheres - INTERLEAVED format: (1 int + 3 floats) per sphere
        for i, (color, coord) in enumerate(zip(data['rigid_colors'], data['rigid_coords'])):
            f.write(struct.pack('i', color))
            for val in coord:
                f.write(struct.pack('f', val))
        
        # Write conf_coord (2 x total_confs)
        for conf in data['conf_coord']:
            f.write(struct.pack('i', conf[0]))  # start
            f.write(struct.pack('i', conf[1] if len(conf) > 1 else conf[0]))  # end
        
        # Write set data
        for s in data['sets']:
            f.write(struct.pack('i', len(s['confs'])))  # conf_count
            f.write(struct.pack('i', s['broken']))
            f.write(struct.pack('f', s['energy']))
            f.write(struct.pack('f', s['max_strain']))
            for conf in s['confs']:
                f.write(struct.pack('i', conf))
        
        # Write cluster data
        for c in data['clusters']:
            f.write(struct.pack('i', c['set_start']))
            f.write(struct.pack('i', c['set_end']))
            f.write(struct.pack('i', c['match_start']))
            f.write(struct.pack('i', c['match_end']))
        
        # Write additional matching spheres
        for c in data['clusters']:
            for i, color in enumerate(c['add_match_colors']):
                f.write(struct.pack('i', color))
                for val in c['add_match_coords'][i]:
                    f.write(struct.pack('f', val))


def write_db2_binary_multi(molecules, output_path):
    """
    Write multiple molecules to a single gzipped binary file.
    Each molecule is written sequentially with its own magic header.
    """
    open_func = gzip.open if str(output_path).endswith('.gz') else open
    
    with open_func(output_path, 'wb') as f:
        for data in molecules:
            write_single_molecule(f, data)


def write_single_molecule(f, data):
    """Write a single molecule to an already-open file handle."""
    # Write magic (6 bytes)
    f.write(DB2BIN_MAGIC)
    
    # Write version (4 bytes int)
    f.write(struct.pack('i', DB2BIN_VERSION))
    
    # Write header integers (8 x 4 bytes)
    f.write(struct.pack('i', data['total_atoms']))
    f.write(struct.pack('i', data['total_bonds']))
    f.write(struct.pack('i', data['total_coords']))
    f.write(struct.pack('i', data['total_confs']))
    f.write(struct.pack('i', len(data['sets'])))  # use actual count
    f.write(struct.pack('i', data['rigid_heavy_count']))
    f.write(struct.pack('i', len(data['clusters'])))  # use actual count
    f.write(struct.pack('i', data['mlines']))
    
    # Write header reals (5 x 4 bytes)
    f.write(struct.pack('f', data['total_charge']))
    f.write(struct.pack('f', data['total_polar_solv']))
    f.write(struct.pack('f', data['total_apolar_solv']))
    f.write(struct.pack('f', data['total_solv']))
    f.write(struct.pack('f', data['total_surfarea']))
    
    # Write refcod (16 bytes, null-padded)
    refcod = data['refcod'][:16].ljust(16, '\x00')
    f.write(refcod.encode('ascii'))
    
    # Write protcod (9 bytes, space-padded)
    protcod = data['protcod'][:9].ljust(9, ' ')
    f.write(protcod.encode('ascii'))
    
    # Write smiles (78 bytes, space-padded)
    smiles = data['smiles'][:78].ljust(78, ' ')
    f.write(smiles.encode('ascii', errors='replace'))
    
    # Write molname (78 bytes, space-padded)
    molname = data['molname'][:78].ljust(78, ' ')
    f.write(molname.encode('ascii', errors='replace'))
    
    # Write arbitrary lines count and data
    arb_count = len(data['arbitrary'])
    f.write(struct.pack('i', arb_count))
    for arb_line in data['arbitrary']:
        arb_bytes = arb_line[:78].ljust(78, ' ')
        f.write(arb_bytes.encode('ascii', errors='replace'))
    
    # Write atom arrays
    for val in data['atom_num']:
        f.write(struct.pack('i', val))
    for name in data['atom_name']:
        f.write(name[:4].ljust(4, ' ').encode('ascii'))
    for atype in data['atom_type']:
        f.write(atype[:5].ljust(5, ' ').encode('ascii'))
    for arr in ['atom_vdwtype', 'atom_color']:
        for val in data[arr]:
            f.write(struct.pack('i', val))
    for arr in ['atom_charge', 'atom_polsolv', 'atom_apolsolv', 
                'atom_totalsolv', 'atom_surfarea']:
        for val in data[arr]:
            f.write(struct.pack('f', val))
    
    # Write bond arrays
    for val in data['bond_num']:
        f.write(struct.pack('i', val))
    for val in data['bond_start']:
        f.write(struct.pack('i', val))
    for val in data['bond_end']:
        f.write(struct.pack('i', val))
    for btype in data['bond_type']:
        f.write(btype[:2].ljust(2, ' ').encode('ascii'))
    
    # Write coord data - INTERLEAVED format
    for i, (idx, coord) in enumerate(zip(data['coord_index'], data['coords'])):
        for val in idx:
            f.write(struct.pack('i', val))
        for val in coord:
            f.write(struct.pack('f', val))
    
    # Write rigid spheres - INTERLEAVED format
    for i, (color, coord) in enumerate(zip(data['rigid_colors'], data['rigid_coords'])):
        f.write(struct.pack('i', color))
        for val in coord:
            f.write(struct.pack('f', val))
    
    # Write conf_coord
    for conf in data['conf_coord']:
        f.write(struct.pack('i', conf[0]))
        f.write(struct.pack('i', conf[1] if len(conf) > 1 else conf[0]))
    
    # Write set data
    for s in data['sets']:
        f.write(struct.pack('i', len(s['confs'])))
        f.write(struct.pack('i', s['broken']))
        f.write(struct.pack('f', s['energy']))
        f.write(struct.pack('f', s['max_strain']))
        for conf in s['confs']:
            f.write(struct.pack('i', conf))
    
    # Write cluster data
    for c in data['clusters']:
        f.write(struct.pack('i', c['set_start']))
        f.write(struct.pack('i', c['set_end']))
        f.write(struct.pack('i', c['match_start']))
        f.write(struct.pack('i', c['match_end']))
    
    for c in data['clusters']:
        for i, color in enumerate(c['add_match_colors']):
            f.write(struct.pack('i', color))
            for val in c['add_match_coords'][i]:
                f.write(struct.pack('f', val))


def convert_file(input_path, output_path):
    """Convert a single DB2 file to binary format."""
    print(f"Converting: {input_path}")
    molecules = parse_db2_text(input_path)
    write_db2_binary_multi(molecules, output_path)
    print(f"  -> {output_path}")
    print(f"     {len(molecules)} molecules")


def main():
    parser = argparse.ArgumentParser(
        description='Convert text DB2.gz files to binary db2bin.gz format'
    )
    parser.add_argument('input', help='Input DB2.gz file or directory')
    parser.add_argument('output', nargs='?', help='Output file or directory')
    parser.add_argument('--dir', action='store_true', 
                        help='Process entire directory')
    args = parser.parse_args()
    
    input_path = Path(args.input)
    
    if args.dir or input_path.is_dir():
        # Directory mode
        output_dir = Path(args.output) if args.output else input_path
        output_dir.mkdir(parents=True, exist_ok=True)
        
        for db2_file in input_path.glob('*.db2.gz'):
            out_name = db2_file.name.replace('.db2.gz', '.db2bin.gz')
            convert_file(db2_file, output_dir / out_name)
        for db2_file in input_path.glob('*.db2'):
            out_name = db2_file.name.replace('.db2', '.db2bin')
            convert_file(db2_file, output_dir / out_name)
    else:
        # Single file mode
        if args.output:
            output_path = Path(args.output)
        else:
            output_path = input_path.with_suffix('').with_suffix('.db2bin.gz')
        
        convert_file(input_path, output_path)
    
    print("\nConversion complete!")


if __name__ == '__main__':
    main()
