from collections import Counter

from epitopes import iedb
from flask import Flask
from flask import redirect, request, render_template, g
import pandas as pd
from werkzeug.contrib.cache import SimpleCache


IEDB_COLS = [
    'Epitope Linear Sequence',
    'MHC Allele Name',
    'Host Organism Name',
    'Epitope Source Organism Name',
    'T Cell ID',
    'Reference ID',
    'PubMed ID',
    'Author',
    'Epitope Modification',
    'Epitope Modified Residues',
    'Epitope Starting Position',
    'Epitope Ending Position',
    'Epitope Source Molecule Name',
    'Method/Technique',
    'Assay Group',
    'Qualitative Measure'
]



cache = SimpleCache()

def iedb_tcell_entries():
    df = cache.get("df")
    if df is None:
        print "Loading IEDB T-Cell Entries"
        df = iedb.load_tcell()
        df = df[IEDB_COLS]
        cache.set("df", df)
    return df

def iedb_tcell_groupby_seq():
    df = cache.get("df_grouped")
    if df is None:
        print "Loading IEDB T-cell Groups"
        df = iedb.load_tcell_values()
        cache.set("df_grouped", df)
    return df

app = Flask(__name__)

@app.route('/search/<query>')
def search(query):
    df_grouped = iedb_tcell_groupby_seq()

    min_assay_count = request.args.get('min_assay_count', 1)
    min_positive_fraction = request.args.get('min_positive_fraction', 0.5)
    offset_from_start = request.args.get('offset_from_start', 0)
    offset_from_end = request.args.get('offset_from_end', 0)

    # TODO: fuzzy match using BLOSUM matrix
    def matches(seq):
        return query in seq[offset_from_start:len(seq)-offset_from_end]

    query_mask = df_grouped.index.map(matches)
    query_mask &= df_grouped['count'] >= min_assay_count

    result_groups = df_grouped[query_mask]
    n = len(result_groups)
    pos_group_mask = result_groups['value'] > min_positive_fraction
    n_pos = pos_group_mask.sum()
    pos_result_groups = result_groups[pos_group_mask]

    df_entries = iedb_tcell_entries()
    assay_results = df_entries['Qualitative Measure']
    pos_entries_mask = assay_results.str.startswith("Positive")
    pos_entries = df_entries[pos_entries_mask]

    # merge peptides matching search criteria with
    # all Positive T-cell assay entries
    merged = pd.merge(
        pos_result_groups,
        df_entries,
        left_index=True,
        right_on=['Epitope Linear Sequence'])

    merged_groups = merged.groupby('Epitope Linear Sequence')

    unique_organisms = merged_groups['Epitope Source Organism Name'].unique()
    organism_counter = Counter()
    for k, organism_list in unique_organisms.iteritems():

        for organism in organism_list:
            organism_counter[organism] += 1

    # Below the organism counts, we're going to display a table
    # of all the peptide/MHC entry combinations and their
    # counts in IEDB
    peptide_mhc_groups = merged.groupby(
        ['Epitope Linear Sequence', 'MHC Allele Name'])
    peptide_mhc_entries = []
    for (seq, allele), g in peptide_mhc_groups:
        d = {
            'Epitope Linear Sequence' : seq,
            'MHC Allele Name' : allele,
        }
        tcell_response = g['Qualitative Measure']
        n = len(tcell_response)
        n_pos = tcell_response.str.startswith("Positive").sum()
        d['Positive Count'] = n_pos
        d['Negative Count'] = n - n_pos
        d['Fraction'] = float(n_pos) / n
        organisms = g['Epitope Source Organism Name']
        organism_counts = organisms.value_counts().iteritems()
        d['Organisms'] = list(organism_counts)
        peptide_mhc_entries.append(d)
    def peptide_mhc_sort_key(d):
        fraction = d['Fraction']
        count = d['Positive Count']
        seq = d['Epitope Linear Sequence']
        return (fraction, count, seq)
    peptide_mhc_entries.sort(key=peptide_mhc_sort_key, reverse=True)
    return render_template(
        'search.html',
        query=query,
        organism_counts=organism_counter.most_common(),
        peptide_mhc_entries=peptide_mhc_entries)


if __name__ == '__main__':
    app.run(debug=True)