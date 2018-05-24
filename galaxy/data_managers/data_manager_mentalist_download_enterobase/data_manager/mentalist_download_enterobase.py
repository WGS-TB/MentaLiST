#!/usr/bin/env python

from __future__ import print_function

import argparse
import datetime
import errno
import os
import string
import subprocess
import sys

from json import dumps, loads


DEFAULT_DATA_TABLE_NAMES = ["mentalist_databases"]


def mentalist_download_enterobase( data_manager_dict, kmer_size, scheme, type, params, target_directory, data_table_names=DEFAULT_DATA_TABLE_NAMES ):
    char_to_full_organism_name = {
        'E': 'Escherichia/Shigella',
        'S': 'Salmonella',
        'Y': 'Yersinia'
    }
    translation_table = string.maketrans(string.punctuation, ("_" * 32))
    base_path = char_to_full_organism_name[scheme].lower().replace(" ", "_").translate(translation_table) + "_enterobase"
    today = datetime.date.today().isoformat()
    scheme_files_path = base_path + "_scheme_" + today
    database_path = base_path + "_k" + str(kmer_size) + "_" + today
    database_name = base_path + "_k" + str(kmer_size) + "_" + today + ".jld"
    display_name = char_to_full_organism_name[scheme] + " k=" + str(kmer_size) + " (Enterobase) " + today
    args = [ 'mentalist', 'download_enterobase', '-s', scheme, '-t', type, '-k', str(kmer_size), '--db', database_name, '-o', scheme_files_path]
    proc = subprocess.Popen( args=args, shell=False, cwd=target_directory )
    return_code = proc.wait()
    if return_code:
        print("Error building database.", file=sys.stderr)
        sys.exit( return_code )
    data_table_entry = dict( value=database_path, dbkey='Enterobase', name=display_name, path=database_name )
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
    parser.add_argument( '-s', '--scheme', dest='scheme', default=None, help="Scheme: ('E'=Escherichia/Shigella, 'S'=Salmonella, 'Y'=Yersinia)")
    parser.add_argument( '-k', '--kmer_size', dest='kmer_size', type=int, default=None, help='kmer Size' )
    parser.add_argument( '-t', '--type', dest='type', default=None, help="Type: ('cg'=cgMLST, 'wg'=wgMLST')")
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
    mentalist_download_enterobase( data_manager_dict, args.kmer_size, args.scheme, args.type, params, target_directory, DEFAULT_DATA_TABLE_NAMES )

    # save info to json file
    open( args.params, 'wb' ).write( dumps( data_manager_dict ) )


if __name__ == "__main__":
    main()
