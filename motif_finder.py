import requests
import pandas as pd

def fetch_sequences_from_ncbi(accessions):
    """
    Fetch DNA sequences from the NCBI nucleotide database using accession IDs.

    Parameters
    ----------
    accessions : list of str
        A list of NCBI accession numbers.

    Returns
    -------
    dict
        Dictionary mapping accession IDs to their corresponding DNA sequences
        as strings.

    Notes
    -----
    If an accession cannot be retrieved, it is skipped.
    """
    results = {}
    for acc in accessions:
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        params = {
            "db": "nucleotide",
            "id": acc,  
            "rettype": "fasta",
            "retmode": "text"
        }

        try:
            response = requests.get(url,params)
            response.raise_for_status()
        except requests.exceptions.RequestException as e:
            print(f"Something went wrong with {acc}: {e}")
            continue
            
        response = response.text
        lines = response.splitlines()
        sequence = "".join(lines[1:])
        results[acc] = sequence
    return results 


def find_motif_positions(seq, motif):
    """
    Find all occurrences of a DNA motif within a sequence.

    Parameters
    ----------
    seq : str
        DNA sequence to be searched.
    motif : str
        DNA motif to search for.

    Returns
    -------
    list of int
        Zero-based positions where the motif starts within the sequence.
    """
    positions= []
    start = 0
    while True:
        index = seq.find(motif,start)
        if index == -1:
            break
        positions.append(index)
        start = index + 1
    return positions    


def reverse_complement(seq):
    """
    Return the reverse complement of a DNA sequence.

    Parameters
    ----------
    seq : str
        DNA sequence.

    Returns
    -------
    str
        Reverse complement of the input sequence.
    """
    table=str.maketrans("ATGC","TACG")
    return seq.translate(table)[::-1]


def find_motif_in_sequence_list(motif, accessions, save=False):
    """
    Search for a DNA motif in multiple sequences retrieved from NCBI.

    The motif is searched in both the forward strand and the reverse complement.
    Motif positions are reported as zero-based indices.

    Parameters
    ----------
    motif : str
        DNA motif to search for.
    accessions : list of str
        List of NCBI accession numbers.
    save : bool, optional
        If True, results are saved to a CSV file (default is False).

    Returns
    -------
    pandas.DataFrame
        DataFrame containing motif counts, positions, sequence length,
        and motif density for each accession.
    """
    results = {}
    sequences = fetch_sequences_from_ncbi(accessions)
    for acc, seq  in sequences.items():

        if not seq:
            print(f"No sequence retrieved for {acc}")
            continue
            
        positions_forward = find_motif_positions(seq, motif)
        reverse_seq = reverse_complement(seq)
        positions_reverse = find_motif_positions(reverse_seq, motif)
            
        results[acc] = {
            "count_forward" : len(positions_forward),
            "positions_forward" : positions_forward,
            "count_reverse_complement" : len(positions_reverse),
            "positions_reverse_complement" : positions_reverse,
            "length" : len(seq),
            "density_forward" : len(positions_forward)/len(seq),
            "density_reverse_complement": len(positions_reverse)/len(seq)
        }
        
    df = pd.DataFrame.from_dict(results, orient= "index")

    if save:
        df.to_csv("motif_results.csv")
        
    return df 