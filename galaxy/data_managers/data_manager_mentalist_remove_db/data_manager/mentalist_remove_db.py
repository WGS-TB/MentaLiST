#!/usr/bin/env python
from __future__ import print_function

import argparse
import os
import subprocess
import sys
import errno
from json import dumps, loads

DEFAULT_DATA_TABLE_NAMES = ["mentalist_databases"]


def mentalist_remove_db( data_manager_dict, database_name, params, data_table_names=DEFAULT_DATA_TABLE_NAMES ):
    data_table_entry = dict( value=database_name, dbkey=database_name, name=database_name, path=database_name )
    for data_table_name in data_table_names:
        _remove_data_table_entry( data_manager_dict, data_table_name, data_table_entry )


def _remove_data_table_entry( data_manager_dict, data_table_name, data_table_entry ):
    data_manager_dict['data_tables'] = data_manager_dict.get( 'data_tables', {} )
    data_manager_dict['data_tables'][ data_table_name ] = data_manager_dict['data_tables'].get( data_table_name, [] )
    data_manager_dict['data_tables'][ data_table_name ].append( data_table_entry )
    return data_manager_dict


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('params')
    parser.add_argument( '-d', '--db', dest='database_name', default=None, help='Database Name' )
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

    # save info to json file
    open( args.params, 'wb' ).write( dumps( data_manager_dict ) )


if __name__ == "__main__":
    main()
