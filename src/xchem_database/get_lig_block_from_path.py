from rich import print as rprint
import gemmi

def get_lig_block_from_path(path):
    cif = gemmi.cif.read(str(path.resolve()))

    key = "comp_LIG"
    try:
        cif['comp_LIG']
    except:
        try:
            key = "comp_XXX"
            cif[key]
        except:
            try:
                key = "comp_UNL"
                cif[key]
            except:
                rprint(path)
    return cif[key]