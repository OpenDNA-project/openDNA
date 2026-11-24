"""`gene` features: a gene feature that covers exactly the bases of a CDS is hidden by the display's "CDS only" switch (the CDS says it all).
Nothing groups features under a gene: the gene is an ordinary feature."""

GENE = "gene"


def span_of(features):
    starts = [s for f in features for s, _ in f.spans]
    ends = [e for f in features for _, e in f.spans]
    return min(starts), max(ends)


def redundant_genes(features):
    """Indices of the gene features at exactly the bases (and strand) of a CDS."""
    cds = {(tuple(f.spans), f.strand) for f in features if f.kind == "CDS"}
    return {f.index for f in features if f.kind == GENE and (tuple(f.spans), f.strand) in cds}
