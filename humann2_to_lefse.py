#!/usr/bin/env python
"""
Takes a QIIME sample metadata mapping file and a Humann2 output table and makes
a lefse-compatable file using the supplied class and subclass categories. 

usage:

python humann2_to_lefse.py -i input_genetable.tsv -t output_genetable.tsv \
-m mapping.txt -c class -s subclass
"""

import argparse

parser = argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-i', '--input_fp', 
    type=str,
    help='path to Humann2 table')

parser.add_argument('-o', '--output_fp', 
    type=str,
    help='path to output lefse-formatted table')

parser.add_argument('-m', '--metadata_fp', 
    type=str, 
    help='path to qiime metadata table')

parser.add_argument('-c', '--class_cat', 
    type=str,
    help='metadata category to use as class')

parser.add_argument('-s', '--subclass_cat', 
    type=str,
    help='metadata category to use as subclass')


def parse_mapping_file(lines, strip_quotes=True, suppress_stripping=False):
    """Parser for map file that relates samples to metadata.
    Format: header line with fields
            optionally other comment lines starting with #
            tab-delimited fields
    Result: list of lists of fields, incl. headers.
    """
    if hasattr(lines, "upper"):
        # Try opening if a string was passed
        try:
            lines = open(lines, 'U')
        except IOError:
            raise IOError("A string was passed that doesn't refer "
                                  "to an accessible filepath.")

    if strip_quotes:
        if suppress_stripping:
            # remove quotes but not spaces
            strip_f = lambda x: x.replace('"', '')
        else:
            # remove quotes and spaces
            strip_f = lambda x: x.replace('"', '').strip()
    else:
        if suppress_stripping:
            # don't remove quotes or spaces
            strip_f = lambda x: x
        else:
            # remove spaces but not quotes
            strip_f = lambda x: x.strip()

    # Create lists to store the results
    mapping_data = []
    header = []
    comments = []

    # Begin iterating over lines
    for line in lines:
        line = strip_f(line)
        if not line or (suppress_stripping and not line.strip()):
            # skip blank lines when not stripping lines
            continue

        if line.startswith('#'):
            line = line[1:]
            if not header:
                header = line.strip().split('\t')
            else:
                comments.append(line)
        else:
            # Will add empty string to empty fields
            tmp_line = map(strip_f, line.split('\t'))
            if len(tmp_line) < len(header):
                tmp_line.extend([''] * (len(header) - len(tmp_line)))
            mapping_data.append(tmp_line)
    if not header:
        raise IOError("No header line was found in mapping file.")
    if not mapping_data:
        raise IOError("No data found in mapping file.")

    return mapping_data, header, comments


def parse_mapping_file_to_dict(*args, **kwargs):
    """Parser for map file that relates samples to metadata.
    input format: header line with fields
            optionally other comment lines starting with #
            tab-delimited fields
    calls parse_mapping_file, then processes the result into a 2d dict, assuming
    the first field is the sample id
    e.g.: {'sample1':{'age':'3','sex':'male'},'sample2':...
    returns the dict, and a list of comment lines
    """
    mapping_data, header, comments = parse_mapping_file(*args, **kwargs)
    return mapping_file_to_dict(mapping_data, header), comments


def mapping_file_to_dict(mapping_data, header):
    """processes mapping data in list of lists format into a 2 deep dict"""
    map_dict = {}
    for i in range(len(mapping_data)):
        sam = mapping_data[i]
        map_dict[sam[0]] = {}
        for j in range(len(header)):
            if j == 0:
                continue  # sampleID field
            map_dict[sam[0]][header[j]] = sam[j]
    return map_dict

def read_humann2_genetable_generator(f):
    header = f.readline().strip().split("\t")
    for l in f:
        line = l.strip().split('\t')
        gene_entry = line[0].split('|')
        gene = gene_entry[0]
        tax = None
        if len(gene_entry) > 1:
            tax = gene_entry[1]
        
        #abund = dict(zip(header[1:],line[1:]))
        yield gene, header[1:], line[1:], tax

def main():
    args = parser.parse_args()

    input_fp = args.input_fp
    output_fp = args.output_fp
    metadata_fp = args.metadata_fp
    class_cat = args.class_cat
    subclass_cat = args.subclass_cat


    md_dict = parse_mapping_file_to_dict(metadata_fp)[0]

    out_f = open(output_fp, 'w')

    with open(input_fp, 'r') as f:

        header = f.readline().strip().split("\t")

        h2_ids = header[1:]

        keep = []

        for i in h2_ids:
            if i in md_dict:
                keep.append(ids.index(i))
            else:
                print('Warning: %s not in metadata map; being dropped',
                      file=sys.stderr)

        ids = [h2_ids[j] for j in keep]

    h2g = read_humann2_genetable_generator(open(input_fp, 'r'))

    class_vals = [md_dict[x][class_cat] for x in ids]
        out_f.write('{0}\t{1}\n'.format(class_cat,class_vals.join('\t')))

    if subclass_cat:
        subclass_vals = [md_dict[x][subclass_cat] for x in ids]
        out_f.write('{0}\t{1}\n'.format(class_cat,class_vals.join('\t')))

    out_f.write('id\t{1}\n'.format(ids.join('\t')))

    for gene, samples, values, tax in h2g:

        out_f.write('{0}\t{1}\n'.format(gene,
                                        [values[s] for s in keep].join('\t'))

    out_f.close()   


if __name__ == "__main__":
    main()