# -*- coding: utf-8 -*-
"""
@author: Rachel Kelemen
@deputy wizard: Joshua Brown
Finalized April 2019

This module contains a function which calculates the enrichment of each sequence in each
selection and writes the results to a csv file.
"""

import csv
import numpy as np
from statistics import mean, stdev

def calculate_enrichments(sequence_count_dict, read_counts, biosamples, sel_1_name,
                          sel_2_name, sel_2_runs, output_format, full_sequence_format,
                          rand_expected_pairs, working_folder_name, lib_name, min_lib_count):
    """Calculates the enrichment of each sequence in eaach selection relative to the
    library and writes the results to results.csv.

    Parameters
    ---
    sequence_count_dict : dict
    biosample : list
    sel_1_name : str
    sel_2_name : str
    sel_2_runs : list
    output_format : str
    full_sequence_format : str
    rand_expected_pairs : list
    working_folder_name : str
    lib_name : str

    Returns
    ---
    results_table : list
        For each sequence the entry looks like
        [[randomized_sequence], []]

    Results table is written to working_folder\results\lib_name results.csv
    """

    # Setup
    # Set up column indices by name because it's clearer than referring to them by number.
    (Randomized_sequence,
    Raw_counts, Fraction_of_total,
    Fold_enrichment, Avg_fold_enrichment) = (0, 1, 2, 3, 4)

    # Column headers
    header_list = ['barcode',
                   'Raw counts']

    
    # Sequences will be added to one of two tables, sequence_table or
    # low_abundance_table.

    # If there are any sequences which are present in selection runs but not present in
    # the library, this will cause divide by 0 errors later in the analysis.
    # Therefore, only certain parts of the analysis - calculating paired bases and % of
    # total reads - can be performed for these sequences.

    # Additionally, sequences which are only observed a very small number of times in the
    # library could lead to inaccurate enrichment factors - if a sequence was seen twice
    # instead of once in the library, for instance, the enrichment factor would be off
    # by a factor of two.

    # sequence_table holds sequences which were found >= min_lib_count times in the
    # library biosample.
    sequence_table = []
    # Enrichment factors will be calculated for these sequences.

    # low_abundance_table holds sequences which were identified < min_lib_count times
    # in the library biosample.
    low_abundance_table = []
    # Enrichment factors will still be calculated for sequences which were identified
    # > 0 times in the library biosample, but these values should be used cautiously.
    # Track how many sequences were not observed at all
    not_observed = 0

    # Lastly, we'll need to find the highest fold enrichment value for each selection
    # before we can calculate the enrichment factors.
    max_fold_enrichments = np.zeros(len(biosamples))


    # Read sequence counts from sequence_count_dict.
    for sequence in sequence_count_dict:
        # Set up table row sequence_data
        sequence_data = []
        for header in header_list:
            sequence_data.append([])

        # 0. Randomized sequence
        # Format sequences for output
        # Replace T's with U's, but keep the original as key_sequence for sequence_count_dict
        key_sequence = sequence
        #sequence = sequence.replace('T', 'U')
        if output_format:
            # output_format is represented as 'NNN/NNN' where each N is to be replaced
            # with a base from sequence and all other characters are preserved.
            # Iterate through characters in output_format to generate formatted_sequence.
            formatted_sequence = ''
            # i represents the current position in sequence, starting at the beginning.
            i = 0
            for c in output_format:
                if c == 'N':
                    formatted_sequence += sequence[i]
                    i += 1
                else:
                    formatted_sequence += c
        else:
            formatted_sequence = sequence
        sequence_data[Randomized_sequence].append(formatted_sequence)


        # 3. Raw counts and
        # 4. Fraction of total
        for raw_count, total in zip(sequence_count_dict[key_sequence], read_counts):
            sequence_data[Raw_counts].append(raw_count)
            #sequence_data[Fraction_of_total].append(raw_count/total)

        # Sequences found in library biosample
        # Calculate the change in abundance before and after each selection,
        # relative to the library.
        if sequence_data[Raw_counts][0] > 0:
            # 6. Fold enrichments
            """library_abundance = sequence_data[Fraction_of_total][0]
            for selection_abundance, i in zip(sequence_data[Fraction_of_total],
                                              np.arange(len(max_fold_enrichments), dtype=int)):
                fold_enrichment = selection_abundance / library_abundance
                sequence_data[Fold_enrichment].append(fold_enrichment)
                # Update max_fold_enrichments unless the value comes from a low abundance
                # sequence.
                if sequence_data[Raw_counts][0] >= min_lib_count:
                    if fold_enrichment > max_fold_enrichments[i]:
                        max_fold_enrichments[i] = fold_enrichment"""

            if sequence_data[Raw_counts][0] >= min_lib_count:
                sequence_table.append(sequence_data)
            else:
                low_abundance_table.append(sequence_data)
                not_observed += 1

        # Sequences not found in library biosample
        # Enrichment cannot be calculated because the library abundance is 0.
        # Calculate the value if the library count was 1
        else:
            # What would library_abundance be if the library count was 1?
            library_abundance = 1 / read_counts[0]
            """for selection_abundance, i in zip(sequence_data[Fraction_of_total],
                                              np.arange(len(max_fold_enrichments))):
                fold_enrichment = selection_abundance / library_abundance
                sequence_data[Fold_enrichment].append(fold_enrichment)
                # Do not update fold enrichments."""
            low_abundance_table.append(sequence_data)


    # Print results
    print('/nStats for the library run:')
    print('Of {:,} sequences, {:,} were observed >= {} times'.format(
            len(sequence_count_dict), len(sequence_table), min_lib_count))
    print('{:,} sequences were observed < {} times, with {} not observed'.format(
            len(low_abundance_table), min_lib_count, not_observed))

    print('/nMost enriched sequence: ', sequence_table[0][Randomized_sequence])
    print('Least enriched sequence: ', sequence_table[-1][Randomized_sequence])

    print('/nData for the top sequence:')
    for header, value in zip(header_list, sequence_table[0]):
        print(header, value)


    # Write the results to csv files.
    # Tables are already sorted by average enrichment factor.
    # Table structure reminder:
    '''
    COLUMN   HEADER                     CONTENTS
    0        Randomized sequence        [formatted randomized sequence]
    1        Raw counts                 [lib_count, sel_1_count, sel_2_count...]
    2        Fraction of total          [lib_fraction, sel_1_fraction, sel_2_fraction...]
    3        Fold enrichment            [sel_1_fold_enrichment, sel_2_fold_enrichment...]
    4        Average enrichment factor  [average enrichment factor]
    '''

    for table, file_name_end in zip([sequence_table, low_abundance_table], [' Results.csv',
                               ' Low abundance sequences.csv']):
        output_file_address = working_folder_name + '//Results//' + lib_name + file_name_end
        with open(output_file_address, 'w', newline='') as output_file:
            # Clear any previous data
            output_file.truncate(0)

            # Write
            output_writer = csv.writer(output_file)
            # Write headers and subheaders
            header_row = []
            subheader_row = []
            # Randomized sequence
            # contain 1 value.
            for header in header_list[Randomized_sequence:Raw_counts]:
                header_row.append(header)
                subheader_row.append('_')
            # Raw counts, Fraction of total, Fold enrichment, Enrichment factors, and
            # Enrichment ranks have one value for each biosample.
            for header in header_list[Raw_counts:]:
                header_row.append(header)
                for l in range(len(biosamples)-1):
                    header_row.append(header)
                for biosample in biosamples:
                    subheader_row.append(biosample)
            # Average enrichment factor and Rank 3 values if there are distinct sel_1 and
            # sel_2 conditions, [Sel_1, Sel_2, and All selections].
            #if sel_2_runs:
            #    for header in header_list[Avg_fold_enrichment:]:
            #        header_row.append(header)
            #        header_row.append('')
            #        header_row.append('')
            #    for selection_name in [sel_1_name, sel_2_name, 'All selections']:
            #        subheader_row.append(selection_name)
            # Otherwise there is just one value for all selections.
            #else:
            #    for header in header_list[Avg_fold_enrichment:]:
            #        header_row.append(header)
            #        subheader_row.append('All selections')
            output_writer.writerow(header_row)
            output_writer.writerow(subheader_row)

            # Write data
            for sequence_data in table:
                data_row = []
                for column in sequence_data:
                    for value in column:
                        data_row.append(value)
                output_writer.writerow(data_row)

        print('/nData written to')
        print(output_file_address)


    return sequence_table, low_abundance_table
