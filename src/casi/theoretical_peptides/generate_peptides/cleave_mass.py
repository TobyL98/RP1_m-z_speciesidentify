"""
cleave_mass.py

Amending code that cleaves peptides
from pyteomics so that italso returns
start position, end position and
number of missed cleavages
"""

import re
import time

import pandas as pd  # type: ignore
from pyteomics import mass


def position_finder(seq: str,
                    rule: str) -> list:
    """
    Creates the position list where the digestion enzyme
    can cut based on the expasy rules

    args 
        seq (str): the peptide sequences
        rule (str): the digestion enzyme that the cleavage rules are based on
            in the expasy dictionary e.g., expasy["trypsin"]
    returns
        position_list (list): Ordered list of numbered cleavage positions (e.g., [1, 9, 24])
    """
    # create position list and assign the start position (0)
    position_list = []
    start_position = 0
    position_list.append(start_position)
    # find all cleave position in the peptide sequence
    for match in re.finditer(rule, seq):
        cleave_position = match.end()
        position_list.append(cleave_position)
    # assign end position as last positio in sequence
    end_position = len(seq)
    position_list.append(end_position)

    return position_list


def peptide_getter(seq: str,
                   position_list: list,
                   miss_number: int,
                   final_peptide_df: pd.DataFrame) -> pd.DataFrame:
    """
    Cleaves all the peptide fragments from a sequence using the
    position list. Accounts for a number of missed cleaves.
    For example, miss_number = 1 would use position 0-2, 1-3
    rather than 0-1, 1-2.

    Args:
        seq (str): the peptide sequence
        position_list (list): positions where a cleavage can occur
        miss_number (int): number of missed cleavages allowed

    Returns:
        final_peptide_df (pd.Dataframe): a dataframe containg the
        sequence of a cleaved peptide, the start and end position
        in the original peptide sequence and the number of missed
        clevages
    """

    num_cleave_sites = len(position_list)
    final_start_number = num_cleave_sites - (miss_number + 1)

    for start_number in range(0, num_cleave_sites):
        if start_number < final_start_number:
            start_position = position_list[start_number]
            end_number = start_number + (miss_number + 1)
            end_position = position_list[end_number]
            peptide_fragment = seq[start_position:end_position]
            # don't add any peptides with an unknown residue (X)
            if "X" not in peptide_fragment:
                single_peptide_list = [
                    peptide_fragment,
                    start_position + 1,
                    end_position,
                    miss_number,
                ]
                final_peptide_df.loc[len(final_peptide_df)] = single_peptide_list

    return final_peptide_df


def peptide_cleaver(seq: str,
                    position_list: list,
                    missed_cleavages: int) -> pd.DataFrame:
    """
    Iterates through the potential missed cleavages (0,1,2)
    Uses peptide_getter to obtain the peptides for the positions
    of a specific peptide clevage

    args
        seq (str): the COl1 peptide sequence
        positions_list (str): an ordered list of positions in the peptide
            sequence where cleavage can occur
        missed_cleavages (int): the number of missed cleavages allowed

    returns
        final_peptide_df (pd.Dataframe): a dataframe of cleaved peptides
            and their start, end positions and number of missed cleavages.
    """

    final_peptide_df = pd.DataFrame(
        [], columns=["seq", "seq_start", "seq_end", "missed_cleaves"]
    )
    for miss_number in range(0, missed_cleavages + 1):
        final_peptide_df = peptide_getter(
            seq, position_list, miss_number, final_peptide_df
        )
    final_peptide_df = final_peptide_df.sort_values(by=["seq_start"])

    return final_peptide_df


def mass_calculator(pep_df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculates the intial monoisotopics mass of
    each peptide in the dataframe. Using pyteomics
    fast mass function to calculate the mass.
    https://pyteomics.readthedocs.io/en/latest/index.html

    args
      pep_df (pd.DataFrame): dataframe of peptides

    returns 
      pep_df (pd.DataFrame): dataframe with additonal mass column
    """

    pep_df["mass1"] = pep_df["seq"].apply(
        lambda seq: mass.fast_mass(sequence=seq, ion_type="M", charge=1)
    )
    # filter for m/z values seen in spectrums
    pep_df = pep_df.loc[(pep_df["mass1"] <= 3500.0) & (pep_df["mass1"] >= 800.0)]

    return pep_df


def possible_ptms(pep_df: pd.DataFrame,
                  ptm_pattern: str,
                  column_name : str) -> pd.DataFrame:
    """
    Calculates the possible hydroxylations and deamidations
    and creates a row for each possible combination

    args
      pep_df (pd.DataFrame): dataframe of peptides
      ptm_pattern (str): pattern for hydroxylations or deamidations
      column_name (str): eithr nhyd (hydroxylations) or ndeam (deamidations)

    returns
      alt_hyd_pep_df (pd.DataFrame): same as pep_df but with news rows 
      in the peptide dataframe for each possible PTM combination and the 
      new mass for that combination.
    """

    # count total possible ptm modification (e.g., hydroxylations or deamidations)
    pep_df[column_name] = pep_df["seq"].str.count(ptm_pattern)

    # create df with all possible PTMS
    alt_hyd_pep_df = pd.DataFrame(columns=pep_df.columns)
    for index, row in pep_df.iterrows():
        num_hydroxylations = int(row[column_name])
        for number in range(0, num_hydroxylations + 1):
            row[column_name] = number
            alt_hyd_pep_df.loc[len(alt_hyd_pep_df)] = row

    del pep_df
    return alt_hyd_pep_df


def ptm_mass(pep_df: pd.DataFrame) -> pd.DataFrame:
    """
    Adds on the mass of the PTM (hydoryxlation or deamidation)
    to the predicted m/z value. PTM masses are the monisotopic
    mass from UNIMOD: https://www.unimod.org/modifications_list.php?
    
    args
        pep_df (pd.DataFrame): dataframe with peptide sequence, mass, 
        number of hydroxylations and number of deamidations
    returns 
        pep_df (pd.DataFrame): dataframe with peptide sequence, 
        mass (with PTM mass added), number of hydroxylations and 
        number of deamidations
    """

    pep_df["mass1"] = pep_df.apply(
        lambda pep_row: pep_row["mass1"]
        + (pep_row["nhyd"] * 15.99415)
        + (pep_row["ndeam"] * 0.984016),
        axis=1,
    )
    return pep_df


def cleave_and_mass(seq: str,
                    rule: str="trypsin",
                    missed_cleavages: int=0) -> pd.DataFrame:
    """
    Function takes a polypeptide sequences as an input
    and in silico digests it using an enzyme such as trypsin.
    Then calculates the mass accounting for possible hydroxylations and deamdiations
    as post-translational modifications.

    args   
        seq (string): collagen peptide sequence 
        rule (string): enzyme to use for cutting e.g., trypsin 
        missed_cleavages (int): number of missed cleavages allowed
    
    returns 
        peptide_df (pd.DataFrame): dataframe with peptide sequence, mass, 
        number of hydroxylations and number of deamidations
    """

    # defining the cleavage rules correctly
    expasy_rules = rules()
    # create position list
    position_list = position_finder(seq, expasy_rules[rule])
    # generates the potential peptides
    peptide_df = peptide_cleaver(seq, position_list, missed_cleavages)
    # calculates the initial mass of the peptides
    peptide_df = mass_calculator(peptide_df)
    # hydroxylations
    peptide_df = possible_ptms(peptide_df, "[PKM]", "nhyd")
    # deamidations
    peptide_df = possible_ptms(peptide_df, "[NQ]", "ndeam")
    # change masses to account for PTMs
    peptide_df = ptm_mass(peptide_df)
    
    return peptide_df


def rules():
    """Possible expasy rules that can be used.
    This dict contains regular expressions for cleavage rules of the most
    popular proteolytic enzymes. The rules were taken from the
    `PeptideCutter tool
    <http://ca.expasy.org/tools/peptidecutter/peptidecutter_enzymes.html>`_
    at Expasy.

    .. note::
        'trypsin_exception' can be used as `exception` argument when calling
        :py:func:`cleave` with 'trypsin' `rule`::

            >>> parser.cleave('PEPTIDKDE', parser.expasy_rules['trypsin'])
            {'DE', 'PEPTIDK'}
            >>> parser.cleave('PEPTIDKDE', parser.expasy_rules['trypsin'], \
    exception=parser.expasy_rules['trypsin_exception'])
            {'PEPTIDKDE'}
    """

    expasy_rules = {
        "arg-c": r"R",
        "asp-n": r"\w(?=D)",
        "bnps-skatole": r"W",
        "caspase 1": r"(?<=[FWYL]\w[HAT])D(?=[^PEDQKR])",
        "caspase 2": r"(?<=DVA)D(?=[^PEDQKR])",
        "caspase 3": r"(?<=DMQ)D(?=[^PEDQKR])",
        "caspase 4": r"(?<=LEV)D(?=[^PEDQKR])",
        "caspase 5": r"(?<=[LW]EH)D",
        "caspase 6": r"(?<=VE[HI])D(?=[^PEDQKR])",
        "caspase 7": r"(?<=DEV)D(?=[^PEDQKR])",
        "caspase 8": r"(?<=[IL]ET)D(?=[^PEDQKR])",
        "caspase 9": r"(?<=LEH)D",
        "caspase 10": r"(?<=IEA)D",
        "chymotrypsin high specificity": r"([FY](?=[^P]))|(W(?=[^MP]))",
        "chymotrypsin low specificity": r"([FLY](?=[^P]))|(W(?=[^MP]))|(M(?=[^PY]))|(H(?=[^DMPW]))",
        "clostripain": r"R",
        "cnbr": r"M",
        "enterokinase": r"(?<=[DE]{3})K",
        "factor xa": r"(?<=[AFGILTVM][DE]G)R",
        "formic acid": r"D",
        "glutamyl endopeptidase": r"E",
        "granzyme b": r"(?<=IEP)D",
        "hydroxylamine": r"N(?=G)",
        "iodosobenzoic acid": r"W",
        "lysc": r"K",
        "ntcb": r"\w(?=C)",
        "pepsin ph1.3": r"((?<=[^HKR][^P])[^R](?=[FL][^P]))|"
        r"((?<=[^HKR][^P])[FL](?=\w[^P]))",
        "pepsin ph2.0": r"((?<=[^HKR][^P])[^R](?=[FLWY][^P]))|"
        r"((?<=[^HKR][^P])[FLWY](?=\w[^P]))",
        "proline endopeptidase": r"(?<=[HKR])P(?=[^P])",
        "proteinase k": r"[AEFILTVWY]",
        "staphylococcal peptidase i": r"(?<=[^E])E",
        "thermolysin": r"[^DE](?=[AFILMV][^P])",
        "thrombin": r"((?<=G)R(?=G))|" r"((?<=[AFGILTVM][AFGILTVWA]P)R(?=[^DE][^DE]))",
        "trypsin": r"([KR](?=[^P]))|((?<=W)K(?=P))|((?<=M)R(?=P))",
        "trypsin_exception": r"((?<=[CD])K(?=D))|((?<=C)K(?=[HY]))|((?<=C)R(?=K))|((?<=R)R(?=[HR]))",
    }
    return expasy_rules


def main():
    """Runs the main function"""
    start = time.time()
    # input test sequence
    seq = "QMSYGYDEKSAGVSVPGPMGPSGPRGLPGPPGAPGPQGFQGPPGEPGEPGASGPMGPRGPPGPPGKNGDDGEAGKPGRPGERGPPGPQGARGLPGTAGLPGMKGHRGFSGLDGAKGDTGPAGPKGEPGSPGENGAPGQMGPRGLPGERGRPGPPGSAGARGNDGAVGAAGPPGPTGPTGPPGFPGAAGAKGEAGPQGARGSEGPQGVRGEPGPPGPAGAAGPAGNPGADGQPGAKGANGAPGIAGAPGFPGARGPSGPQGPSGAPGPKGNSGEPGAPGNKGDTGAKGEPGPAGVQGPPGPAGEEGKRGARGEPGPSGLPGPPGERGGPGSRGFPGADGVAGPKGPAGERGSPGPAGPKGSPGEAGRPGEAGLPGAKGLTGSPGSPGPDGKTGPPGPAGQDGRPGPAGPPGARGQAGVMGFPGPKGTAGEPGKAGERGVPGPPGAVGPAGKDGEAGAQGAPGPAGPAGERGEQGPAGSPGFQGLPGPAGPPGEAGKPGEQGVPGDLGAPGPSGARGERGFPGERGVQGPPGPAGPRGNNGAPGNDGAKGDTGAPGAPGSQGAPGLQGMPGERGAAGLPGPKGDRGDAGPKGADGSPGKDGVRGLTGPIGPPGPAGAPGDKGETGPSGPAGPTGARGAPGDRGEPGPPGPAGFAGPPGADGQPGAKGEPGDTGVKGDAGPPGPAGPAGPPGPIGNVGAPGPKGSRGAAGPPGATGFPGAAGRVGPPGPSGNAGPPGPPGPVGKEGGKGPRGETGPAGRPGEVGPPGPPGPAGEKGSPGADGPAGSPGTPGPQGIAGQRGVVGLPGQRGERGFPGLPGPSGEPGKQGPSGASGERGPPGPMGPPGLAGPPGESGREGSPGAEGSPGRDGAPGAKGDRGETGPAGPPGAPGAPGAPGPVGPAGKNGDRGETGPAGPAGPIGPAGARGPAGPQGPRGDKGETGEQGDRGIKGHRGFSGLQGPPGSPGSPGEQGPSGASGPAGPRGPPGSAGSPGKDGLNGLPGPIGPPGPRGRTGDSGPAGPPGPPGPPGPPGPPSGGYDFSFLPQPPQEKSQDGGRYYRARQYSDKGVSAGPGPMGLMGPRGPPGAVGAPGPQGFQGPAGEPGEPGQTGPAGSRGPAGPPGKAGEDGHPGKPGRPGERGVVGPQGARGFPGTPGLPGFKGIRGHNGLDGLKGQPGAQGVKGEPGAPGENGTPGQAGARGLPGERGRVGAPGPAGARGSDGSVGPVGPAGPIGSAGPPGFPGAPGPKGELGPVGNPGPAGPAGPRGEAGLPGLSGPVGPPGNPGANGLTGAKGATGLPGVAGAPGLPGPRGIPGPVGAAGATGPRGLVGEPGPAGSKGETGNKGEPGSAGAQGPPGPSGEEGKRGSPGEPGSAGPAGPPGLRGSPGSRGLPGADGRAGVMGPPGNRGSTGPAGVRGPNGDAGRPGEPGLMGPRGLPGSPGNVGPAGKEGPVGLPGIDGRPGPIGPAGPRGEAGNIGFPGPKGPSGDPGKPGEKGHPGLAGARGAPGPDGNNGAQGPPGPQGVQGGKGEQGPAGPPGFQGLPGPSGTAGEVGKPGERGLPGEFGLPGPAGPRGERGPPGESGAAGPSGPIGIRGPSGAPGPDGNKGEAGAVGAPGSAGASGPGGLPGERGAAGIPGGKGEKGETGLRGEIGNPGRDGARGAPGAIGAPGPAGASGDRGEAGAAGPSGPAGPRGSPGERGEVGPAGPNGFAGPAGSAGQPGAKGEKGTKGPKGENGIVGPTGPVGAAGPSGPNGPPGPAGSRGDGGPPGMTGFPGAAGRTGPPGPSGITGPPGPPGAAGKEGIRGPRGDQGPVGRTGEIGASGPPGFAGEKGPSGEPGTTGPPGTAGPQGLLGAPGILGLPGSRGERGQPGIAGALGEPGPLGIAGPPGARGPPGAVGSPGVNGAPGEAGRDGNPGSDGPPGRDGQPGHKGERGYPGNIGPTGAAGAPGPHGSVGPAGKHGNRGEPGPAGSVGPVGAVGPRGPSGPQGIRGDKGEPGDKGARGLPGLKGHNGLQGLPGLAGLHGDQGAPGPVGPAGPRGPAGPSGPIGKDGRSGHPGPVGPAGVRGSQGSQGPAGPPGPPGPPGPPGVSGGGYDFGFEGGFYRA"
    cleave_and_mass(seq,"trypsin", 1)
    end = time.time()
    elapsed_time = end - start
    print(f"Elapsed time = {elapsed_time}")


if __name__ == "__main__":
    main()
