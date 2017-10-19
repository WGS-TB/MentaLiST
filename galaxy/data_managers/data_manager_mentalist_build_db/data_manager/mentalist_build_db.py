#!/usr/bin/env python
from __future__ import print_function

import argparse
import os
import subprocess
import sys
import errno
from json import dumps, loads

DEFAULT_DATA_TABLE_NAMES = ["mentalist_databases"]


def mentalist_build_db( data_manager_dict, database_name, kmer_size, profile, fasta_files, params, target_directory, data_table_names=DEFAULT_DATA_TABLE_NAMES ):
    args = [ 'mentalist', 'build_db', '--db', database_name, '-k', str(kmer_size)]
    if profile:
        args += ['--profile', profile]
    print(args)
    args += ['--fasta_files'] + fasta_files
    proc = subprocess.Popen( args=args, shell=False, cwd=target_directory )
    return_code = proc.wait()
    if return_code:
        print("Error building database.", file=sys.stderr)
        sys.exit( return_code )
    data_table_entry = dict( value=database_name, dbkey=database_name, name=database_name, path=database_name )
    for data_table_name in data_table_names:
        _add_data_table_entry( data_manager_dict, data_table_name, data_table_entry )


def _add_data_table_entry( data_manager_dict, data_table_name, data_table_entry ):
    data_manager_dict['data_tables'] = data_manager_dict.get( 'data_tables', {} )
    data_manager_dict['data_tables'][ data_table_name ] = data_manager_dict['data_tables'].get( data_table_name, [] )
    data_manager_dict['data_tables'][ data_table_name ].append( data_table_entry )
    return data_manager_dict


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('params')
    parser.add_argument( '-d', '--db', dest='database_name', default=None, help='Database Name' )
    parser.add_argument( '-f', '--fasta_files', dest='fasta_files', nargs='+', default=None, help='FASTA Filenames' )
    parser.add_argument( '-k', '--kmer_size', dest='kmer_size', type=int, default=None, help='kmer Size' )
    parser.add_argument( '-p', '--profile', dest='profile', type=int, default=None, help='Profile' )
    args = parser.parse_args()

    params = loads( open( args.params ).read() )
    target_directory = params[ 'output_data' ][0]['extra_files_path']

    try:
        os.mkdir( target_directory )
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir( target_directory ):
            pass
        else:
            raise

    data_manager_dict = {}

    # build the index
    mentalist_build_db( data_manager_dict, args.database_name, args.kmer_size, args.profile, args.fasta_files, params, target_directory, DEFAULT_DATA_TABLE_NAMES )

    # save info to json file
    open( args.params, 'wb' ).write( dumps( data_manager_dict ) )


if __name__ == "__main__":
    main()
