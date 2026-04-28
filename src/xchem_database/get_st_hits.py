def get_st_hits(structure):
    hits = {}
    for model in structure:
        for chain in model:
            for res in chain.first_conformer():
                if res.name in ["LIG", "XXX"]:
                    hits[f"{chain.name}_{res.seqid.num}"] = res
    return hits