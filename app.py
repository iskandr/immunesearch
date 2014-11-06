from collections import Counter
from os import environ

from epitopes import iedb
from flask import Flask
from flask import redirect, request, render_template, g
import pandas as pd
from werkzeug.contrib.cache import SimpleCache

DEBUG = bool(environ.get("DEBUG", False))

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
        df = iedb.load_tcell(nrows=1000 if DEBUG else None)
        df = df[IEDB_COLS]
        cache.set("df", df)
    return df

def iedb_tcell_groupby_seq():
    df = cache.get("df_grouped")
    if df is None:
        print "Loading IEDB T-cell Groups"
        df = iedb.load_tcell_values(nrows=1000 if DEBUG else None)
        cache.set("df_grouped", df)
    return df

app = Flask(__name__)


def get_arg(name, value):
    converter = type(value)
    v = converter(request.args.get(name, value))
    try:
        return converter(v)
    except:
        return value

@app.route('/')
def search():
    query = request.args.get("query")
    print "Query:", query
    if query is None:
        return render_template('search.html')

    df_grouped = iedb_tcell_groupby_seq()

    min_assay_count = get_arg('minAssayCount', 1)
    min_positive_fraction = get_arg('minPositiveFraction', 0.5)
    offset_from_start = get_arg('offsetFromStart', 0)
    offset_from_end = get_arg('offsetFromEndend', 0)

    # TODO: fuzzy match using BLOSUM matrix
    def matches(seq):
        return query in seq[offset_from_start:len(seq)-offset_from_end]

    query_mask = df_grouped.index.map(matches)
    query_mask &= df_grouped['count'] >= min_assay_count

    result_groups = df_grouped[query_mask]

    n = len(result_groups)

    pos_group_mask = result_groups['value'] > min_positive_fraction
    total_n_pos = pos_group_mask.sum()
    total_n_neg = n - total_n_pos
    total_percent_pos =  float(total_n_pos) / n


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
        d['Positive Fraction'] = float(n_pos) / n
        organisms = g['Epitope Source Organism Name']
        organism_counts = organisms.value_counts().iteritems()
        d['Organisms'] = list(organism_counts)
        peptide_mhc_entries.append(d)

    return render_template(
        'results.html',
        n_pos=total_n_pos,
        n_neg=total_n_neg,
        total_percent_positive=total_percent_pos,
        min_positive_fraction=min_positive_fraction,
        min_assay_count=min_assay_count,
        offset_from_start=offset_from_start,
        offset_from_end=offset_from_end,
        query=query,
        organism_counts=organism_counter.most_common(),
        peptide_mhc_entries=peptide_mhc_entries)


if __name__ == '__main__':
    print "Debug:", DEBUG
    app.run(debug=DEBUG)