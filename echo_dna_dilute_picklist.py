#!/usr/bin/env python
"""
decho_dna_dilute_picklist.py 

Takes a csv input file with sample name and DNA concentration in separate cols
and creates an Echo pick list to dilute the DNA to the specified concentration
and volume. 

Input csv format:
Three columns, with a header row, in the following order:

Sample,Conc,Well
Bob,1,A1
Alice,5,C1
Dave,0.1,A3
Stacey,10,C3
"""

from __future__ import print_function
import sys
import os
import argparse
import unittest
from io import StringIO

parser = argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-i', '--in_csv', 
    type=str,
    help='path to csv with input concentrations')

parser.add_argument('-o', '--out_csv', 
    type=str,
    help='path to output csv for pick list (default: StdOut)')

parser.add_argument('-v', '--volume', 
    type=float, default=5.0,
    help='desired ending volume in destination plate (default: %(default)s)')

parser.add_argument('-c', '--concentration', 
    type=float, default=0.2,
    help='desired ending DNA concentration, in ng/uL (default: %(default)s)')

parser.add_argument('-w', '--water_plate_name', 
    type=str, default='Water',
    help='name to use for diluent source plate (default: %(default)s)')

parser.add_argument('-d', '--dna_plate_name', 
    type=str, default='DNA',
    help='name to use for DNA source plate (default: %(default)s)')

parser.add_argument('--destination_plate_name', 
    type=str, default='DilutedDNA',
    help='name to use for destination plate (default: %(default)s)')

parser.add_argument('--output_plate_format', 
    type=str, default='384',
    help='plate type for output (default: %(default)s)')

parser.add_argument('--input_plate_format', 
    type=str, default='384',
    help='plate type for input (default: %(default)s)')

parser.add_argument('--col_offset', 
    type=int, default=0,
    help='offset this number of columns in dest plate (default: %(default)s)')

parser.add_argument('-t', '--test', 
    action='store_true',
    help='run unittest')


class TestDilute(unittest.TestCase):
    def test_load_sample_concs(self):
        input_csv = ('Sample,Conc,Well\n'
                     'Bob,1.0,A1\n'
                     'Alice,5,B1\n'
                     'Dave,0.1,A2\n'
                     'Stacey,10,B2')

        input_f = StringIO(input_csv)

        exp_samples = ['Bob','Alice','Dave','Stacey']
        exp_concs = [1.0, 5.0, 0.1, 10.0]
        exp_rows = ['A','B','A','B']
        exp_cols = [1, 1, 2, 2]

        obs_samples, obs_concs, obs_rows, obs_cols = load_sample_concs(input_f,
                                                                       header=True)

        self.assertEqual(exp_samples, obs_samples)
        self.assertEqual(exp_concs, obs_concs)
        self.assertEqual(exp_rows, obs_rows)
        self.assertEqual(exp_cols, obs_cols)
    
    def test_step_round(self):
        down_input = 7.6
        up_input = 8.8

        down_exp = 7.5
        up_exp = 10

        self.assertEqual(down_exp, step_round(down_input, 2.5))
        self.assertEqual(up_exp, step_round(up_input, 2.5))
    
    def test_calc_dilution(self):

        low_input = 0.1
        high_input = 10

        obs_low_sample, obs_low_water = calc_dilution(low_input, 0.2, 5)
        exp_low_sample = 5.0
        exp_low_water = 0

        self.assertEqual(obs_low_sample,exp_low_sample)
        self.assertEqual(obs_low_water,exp_low_water)


        obs_high_sample, obs_high_water = calc_dilution(high_input, 0.2, 5)
        exp_high_sample = 0.1
        exp_high_water = 4.9

        self.assertEqual(obs_high_sample,exp_high_sample)
        self.assertEqual(obs_high_water,exp_high_water)

    def test_find_output_well(self):

        # 384 to 96

        input_row = 'C'
        input_col = 3

        obs_output_row, obs_output_col = find_output_well(input_row,
                                                          input_col,
                                                          input_fmt='384',
                                                          output_fmt='96',
                                                          condensed_output=False,
                                                          output_quadrant='NW')

        exp_output_row = 'B'
        exp_output_col = 2

        self.assertEqual(obs_output_row, exp_output_row)
        self.assertEqual(obs_output_col, exp_output_col)

        # 384 to 384

        obs_output_row, obs_output_col = find_output_well(input_row,
                                                          input_col,
                                                          input_fmt='384',
                                                          output_fmt='384',
                                                          condensed_output=False,
                                                          output_quadrant='NW')

        exp_output_row = 'C'
        exp_output_col = 3

        self.assertEqual(obs_output_row, exp_output_row)
        self.assertEqual(obs_output_col, exp_output_col)

        # 384 to 384, offset 6 cols

        obs_output_row, obs_output_col = find_output_well(input_row,
                                                          input_col,
                                                          input_fmt='384',
                                                          output_fmt='384',
                                                          condensed_output=False,
                                                          output_quadrant='NW',
                                                          col_offset = 6)

        exp_output_row = 'C'
        exp_output_col = 9

        self.assertEqual(obs_output_row, exp_output_row)
        self.assertEqual(obs_output_col, exp_output_col)

        # 96 to 384

        obs_output_row, obs_output_col = find_output_well(input_row,
                                                          input_col,
                                                          input_fmt='96',
                                                          output_fmt='384',
                                                          condensed_output=False,
                                                          output_quadrant='NW')

        exp_output_row = 'E'
        exp_output_col = 5

        self.assertEqual(obs_output_row, exp_output_row)
        self.assertEqual(obs_output_col, exp_output_col)

        # 384 to 384, condensed, NW

        obs_output_row, obs_output_col = find_output_well(input_row,
                                                          input_col,
                                                          input_fmt='384',
                                                          output_fmt='384',
                                                          condensed_output=True,
                                                          output_quadrant='NW')

        exp_output_row = 'B'
        exp_output_col = 2

        self.assertEqual(obs_output_row, exp_output_row)
        self.assertEqual(obs_output_col, exp_output_col)

        # 384 to 384, condensed, NW

        obs_output_row, obs_output_col = find_output_well(input_row,
                                                          input_col,
                                                          input_fmt='384',
                                                          output_fmt='384',
                                                          condensed_output=True,
                                                          output_quadrant='NE')

        exp_output_row = 'B'
        exp_output_col = 14

        self.assertEqual(obs_output_row, exp_output_row)
        self.assertEqual(obs_output_col, exp_output_col)

        # 384 to 384, condensed, SE

        obs_output_row, obs_output_col = find_output_well(input_row,
                                                          input_col,
                                                          input_fmt='384',
                                                          output_fmt='384',
                                                          condensed_output=True,
                                                          output_quadrant='SE')

        exp_output_row = 'J'
        exp_output_col = 14

        self.assertEqual(obs_output_row, exp_output_row)
        self.assertEqual(obs_output_col, exp_output_col)

        # 384 to 384, condensed, SW

        obs_output_row, obs_output_col = find_output_well(input_row,
                                                          input_col,
                                                          input_fmt='384',
                                                          output_fmt='384',
                                                          condensed_output=True,
                                                          output_quadrant='SW')

        exp_output_row = 'J'
        exp_output_col = 2

        self.assertEqual(obs_output_row, exp_output_row)
        self.assertEqual(obs_output_col, exp_output_col)

    def test_format_output(self):

        exp_water_output = 'Water,,,C3,,Stacey,,,DilutedDNA,,B2,4900.0\n'
        exp_sample_output = 'Sample,,,C3,,Stacey,,,DilutedDNA,,B2,100.0\n'

        obs_water_output = format_output(source_plate_name='Water',
                                         source_well='C3',
                                         sample_name='Stacey',
                                         destination_plate_name='DilutedDNA', 
                                         destination_well='B2',
                                         volume=4.9,
                                         multiplier=1000)

        obs_sample_output = format_output(source_plate_name='Sample',
                                         source_well='C3',
                                         sample_name='Stacey',
                                         destination_plate_name='DilutedDNA', 
                                         destination_well='B2',
                                         volume=0.1,
                                         multiplier=1000)

        self.assertEqual(obs_water_output, exp_water_output)
        self.assertEqual(obs_sample_output, exp_sample_output)


def run_unittests():
    TestDilute('test_load_sample_concs').test_load_sample_concs()
    TestDilute('test_step_round').test_step_round()
    TestDilute('test_calc_dilution').test_calc_dilution()
    TestDilute('test_find_output_well').test_find_output_well()
    TestDilute('test_format_output').test_format_output()

def format_output(source_plate_name,
                  source_well,
                  sample_name,
                  destination_plate_name, 
                  destination_well,
                  volume,
                  multiplier=1000):

    out_string = '{},{},{},{},{},{},{},{},{},{},{},{:.1f}\n'.format(source_plate_name,
          '',
          '',
          source_well,
          '',
          sample_name,
          '',
          '',
          destination_plate_name,
          '',
          destination_well,
          volume * multiplier)

    return(out_string)


def find_output_well(input_row,
                     input_col,
                     input_fmt='384',
                     output_fmt='96',
                     condensed_output=False,
                     output_quadrant=1,
                     row_offset=0,
                     col_offset=0):
    """
    Finds output well given input well and specified plate formats
    """
    input_row_ord = ord(input_row.upper()) - 65
    input_col = int(input_col)

    if input_fmt == '384' and output_fmt == '96':
        output_row = int(input_row_ord / 2 + 65)
        output_col = int((input_col - 1)/ 2 + 1)

    elif input_fmt == '96' and output_fmt == '384':
        output_row = int(input_row_ord * 2 + 65)
        output_col = int((input_col - 1) * 2 + 1)

    elif input_fmt == '384' and output_fmt == '384':
        if condensed_output:
            if output_quadrant is 'NW':
                output_row = int(input_row_ord / 2 + 65)
                output_col = int((input_col - 1)/ 2 + 1)
        
            if output_quadrant is 'NE':
                output_row = int(input_row_ord / 2 + 65)
                output_col = int((input_col - 1)/ 2 + 1 + 12)
        
            if output_quadrant is 'SE':
                output_row = int(input_row_ord / 2 + 65 + 8)
                output_col = int((input_col - 1)/ 2 + 1 + 12)
        
            if output_quadrant is 'SW':
                output_row = int(input_row_ord / 2 + 65 + 8)
                output_col = int((input_col - 1)/ 2 + 1)
        
        else:
            output_row = ord(input_row)
            output_col = input_col

    else:
        raise NotImplementedError('This conversion type not implemented')

    return(chr(output_row + row_offset), output_col + col_offset)


def load_sample_concs(csv_f, sep=',', header=False):
    """
    Loads sample concentrations from CSV file
    """
    samples = []
    concs = []
    rows = []
    cols = []

    for line in csv_f:
        if header:
            header=False
            continue
        sample, conc, well = line.strip().split(sep)
        samples.append(str(sample))
        concs.append(float(conc))
        rows.append(str(well[0]))
        cols.append(int(well[1:]))

    return(samples,concs,rows,cols)


def step_round(number, step):
    return(round(number / step) * step)


def calc_dilution(input_conc, output_conc, output_vol, increment=0.0025, max_avail=True):
    """
    Calculates dilution volumes 
    """
    input_conc = float(input_conc)
    output_conc = float(output_conc)
    output_vol = float(output_vol)

    if input_conc >= output_conc:
        sample_vol = step_round((output_conc/input_conc * output_vol), increment)
        water_vol = output_vol - sample_vol
    elif max_avail:
        sample_vol = output_vol
        water_vol = 0
    else:
        sample_vol = None
        water_vol = None

    return(sample_vol, water_vol)


def main():
    args = parser.parse_args()

    in_csv = args.in_csv
    out_csv = args.out_csv
    volume = args.volume
    concentration = args.concentration
    water_plate_name = args.water_plate_name
    dna_plate_name = args.dna_plate_name
    col_offset = args.col_offset
    output_plate_format = str(args.output_plate_format)
    input_plate_format = str(args.input_plate_format)
    destination_plate_name = args.destination_plate_name
    test = args.test
    
    if test:
        run_unittests()
        return(0)

    with open(in_csv) as input_f:
        samples, concs, rows, cols = load_sample_concs(input_f, header=True)

    header = ['Source Plate Name',
              'Source Plate Barcode',
              'Source Plate Type',
              'Source Well',
              'Sample ID',
              'Sample Name',
              'Sample Group',
              'Sample Comment',
              'Destination Plate Name',
              'Destination Plate Barcode',
              'Destination Well',
              'Transfer Volume']

    header= ','.join(header)

    water_str = ''
    sample_str = ''

    if input_plate_format == '96':
        print("Warning! Echo can not handle 96 well input.\n"
              "You will need to manually fix the input well values")


    for i in range(len(samples)):

        sample_vol, water_vol = calc_dilution(concs[i],
                                              concentration,
                                              volume,
                                              increment=0.0025,
                                              max_avail=True)

        output_row, output_col = find_output_well(rows[i],
                                                  cols[i],
                                                  input_fmt=input_plate_format,
                                                  output_fmt=output_plate_format,
                                                  col_offset = col_offset)

        source_well = str(rows[i]) + str(cols[i])
        destination_well = str(output_row) + str(output_col)

        water_str += format_output(water_plate_name,
                                   source_well,
                                   samples[i],
                                   destination_plate_name, 
                                   destination_well,
                                   water_vol,
                                   multiplier=1000)

        sample_str += format_output(dna_plate_name,
                                   source_well,
                                   samples[i],
                                   destination_plate_name, 
                                   destination_well,
                                   sample_vol,
                                   multiplier=1000)

    if out_csv is not None:
        with open(out_csv, 'w') as out_f:
            out_f.write(header + '\n')
            out_f.write(water_str + sample_str)

    else:
        print(header)
        print(water_str + sample_str)


if __name__ == "__main__":
    main()



