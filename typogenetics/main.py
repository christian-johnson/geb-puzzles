import click
import numpy as np

AMINO_DICT = {
    "AA": {"amino": "punctuation", "direction": "N/A"},
    "AC": {"amino": "cut", "direction": "s"},
    "AG": {"amino": "del", "direction": "s"},
    "AT": {"amino": "swi", "direction": "r"},
    "CA": {"amino": "mvr", "direction": "s"},
    "CC": {"amino": "mvl", "direction": "s"},
    "CG": {"amino": "cop", "direction": "r"},
    "CT": {"amino": "off", "direction": "l"},
    "GA": {"amino": "ina", "direction": "s"},
    "GC": {"amino": "inc", "direction": "r"},
    "GG": {"amino": "ing", "direction": "r"},
    "GT": {"amino": "int", "direction": "l"},
    "TA": {"amino": "rpy", "direction": "r"},
    "TC": {"amino": "rpu", "direction": "l"},
    "TG": {"amino": "lpy", "direction": "l"},
    "TT": {"amino": "lpu", "direction": "l"},
}


def counterpart(B):
    if B == "A":
        return "T"
    if B == "T":
        return "A"
    if B == "C":
        return "G"
    if B == "G":
        return "C"
    if B == " ":
        return " "


def find_preferential_binding(dirs):
    """
    dirs: a list of directions, e.g. ["r","l","r"...]
    """
    net_right_turns = (dirs.count("r") - dirs.count("l")) % 4
    if net_right_turns == 0:
        return "A"
    if net_right_turns == 1:
        return "G"
    if net_right_turns == 2:
        return "T"
    return "C"


class Strand:
    def __init__(self, string):
        self.string = string

    def derive_codons(self):
        """
        Optional step to find the codons on this strand, for a ribosome to use
        Not needed during enzyme operations
        """
        self.codons = [
            self.string[2 * i: 2 * i + 2]
            for i in range(round(len(self.string) / 2 + 0.1))
        ]


class Enzyme:
    def __init__(self, seq, directions):
        self.seq = seq
        self.directions = directions
        self.preferential_binding = find_preferential_binding(self.directions)


class Operation:
    """
    An 'Operation' is the action of a single enzyme on
    a single strand. The outcome is one or more strands, possibly
    including the initial one.
    """

    def __init__(self, enzyme, strand, verbose=False):
        self.verbose = verbose
        self.enzyme = enzyme
        self.strand = strand
        # initially, this is nothing
        self.second_strand = Strand(len(self.strand.string) * " ")
        self.copy_mode = False
        # Store the outcome of this enzyme/strand set of operations
        self.final_strings = []
        # A placeholder array to store some extra strings
        self.strings = []
        self.operate()

    def operate(self):
        """
        Bind the enzyme & begin operating on the strand
        """
        # seq_i is the sequence number of the Enzyme, i.e. which operation we are on
        self.seq_i = 0

        # bound is the position of the enzyme on the Strand it's operating on
        self.bind()

        if self.verbose:
            print("")
            print("String: ", self.strand.string)
            print("Enzyme sequence: ", self.enzyme.seq)
            print("initial binding preference: ",
                  self.enzyme.preferential_binding)
            print("initial binding position: ", self.bound)
            print("")
        while (
            self.bound > -1
            and self.bound < len(self.strand.string) - 1
            and self.seq_i < len(self.enzyme.seq)
        ):
            self.iter()
            self.seq_i += 1

        # Strings created from cuts
        for string in self.strings:
            substrings = string.split(" ")
            for substring in substrings:
                if " " not in substring and len(substring):
                    self.final_strings.append(substring)

        # Ending strings
        self.primary_strings = "".join(self.strand.string).split(" ")
        self.secondary_strings = "".join(self.second_strand.string).split(" ")
        # Split up on empty sequences
        if self.verbose:
            print("self.primary_strings = ", self.primary_strings)
            print("self.second_string = ", self.second_strand.string)
        for entry in self.primary_strings:
            if len(entry) and " " not in entry:
                if self.verbose:
                    print("adding string: |" + entry + "|")
                self.final_strings.append(entry)
        for entry in self.secondary_strings:
            if len(entry) and " " not in entry:
                if self.verbose:
                    print("adding string: |" + entry + "|")
                self.final_strings.append(entry)

    def move_right(self):
        self.bound += 1
        if (
            self.bound < len(self.strand.string) - 1
            and self.copy_mode
            and self.second_strand.string[self.bound] == " "
        ):
            self.second_strand.string = (
                self.second_strand.string[: self.bound]
                + counterpart(self.strand.string[self.bound])
                + self.second_strand.string[self.bound + 1:]
            )

    def move_left(self):
        self.bound -= 1
        if (
            self.bound > -1
            and self.copy_mode
            and self.second_strand.string[self.bound] == " "
        ):
            self.second_strand.string = (
                self.second_strand.string[: self.bound]
                + counterpart(self.strand.string[self.bound])
                + self.second_strand.string[self.bound + 1:]
            )

    def bind(self):
        """
        Check whether the protein will bind at all to the strand.
        If so, set the position
        """
        self.bound = self.strand.string.find(self.enzyme.preferential_binding)

    def iter(self):
        """
        Perform the cut/insert/copy/etc operations
        """
        # First look up which operation we want
        action = self.enzyme.seq[self.seq_i]

        if self.verbose:
            print("Current action:", action)
            print("Secondary string: |" + self.second_strand.string + "|")
            print("Primary string:   |" + self.strand.string + "|")
            print(
                "                   "
                + " " * self.bound
                + "^"
                + " " * (len(self.strand.string) - self.bound)
            )

        if action == "punctuation":
            # end the operation
            self.bound = -1

        # These are easy
        if action == "off":
            self.copy_mode = False
        if action == "cop":
            self.copy_mode = True

        # insert & delete bases
        if action == "del":
            self.strand.string = (
                self.strand.string[: self.bound] +
                self.strand.string[self.bound + 1:]
            )
            self.strand.string += " "
            # This explicitly doesn't apply to the second strand

        if action == "ina":
            self.strand.string = (
                self.strand.string[: self.bound + 1]
                + "A"
                + self.strand.string[self.bound + 1:]
            )
            if self.copy_mode:
                self.second_strand.string = (
                    self.second_strand.string[: self.bound + 1]
                    + counterpart("A")
                    + self.second_strand.string[self.bound + 1:]
                )

            else:
                self.second_strand.string = (
                    self.second_strand.string[: self.bound + 1]
                    + " "
                    + self.second_strand.string[self.bound + 1:]
                )
        if action == "ing":
            self.strand.string = (
                self.strand.string[: self.bound + 1]
                + "G"
                + self.strand.string[self.bound + 1:]
            )
            if self.copy_mode:
                self.second_strand.string = (
                    self.second_strand.string[: self.bound + 1]
                    + counterpart("G")
                    + self.second_strand.string[self.bound + 1:]
                )
            else:
                self.second_strand.string = (
                    self.second_strand.string[: self.bound + 1]
                    + " "
                    + self.second_strand.string[self.bound + 1:]
                )

        if action == "int":
            self.strand.string = (
                self.strand.string[: self.bound + 1]
                + "T"
                + self.strand.string[self.bound + 1:]
            )
            if self.copy_mode:
                self.second_strand.string = (
                    self.second_strand.string[: self.bound + 1]
                    + counterpart("T")
                    + self.second_strand.string[self.bound + 1:]
                )
            else:
                self.second_strand.string = (
                    self.second_strand.string[: self.bound + 1]
                    + " "
                    + self.second_strand.string[self.bound + 1:]
                )
        if action == "inc":
            self.strand.string = (
                self.strand.string[: self.bound + 1]
                + "C"
                + self.strand.string[self.bound + 1:]
            )
            if self.copy_mode:
                self.second_strand.string = (
                    self.second_strand.string[: self.bound + 1]
                    + counterpart("C")
                    + self.second_strand.string[self.bound + 1:]
                )
            else:
                self.second_strand.string = (
                    self.second_strand.string[: self.bound + 1]
                    + " "
                    + self.second_strand.string[self.bound + 1:]
                )
        if action == "cut":
            # Cut the strand into two to the right of the current location
            # (Note this is ambiguous in the instructions in GEB, but matches
            # the example given)
            # This applies to both strands if there is a second one

            # save the cut off strings
            self.strings.append(self.strand.string[self.bound + 1:])
            self.strings.append(
                self.second_strand.string[self.bound + 1:][::-1])
            # Now cut
            self.strand.string = self.strand.string[: self.bound + 1]
            # also cut the second strand
            self.second_strand.string = self.second_strand.string[: self.bound + 1]
            # Now that it's been cut, we're at the last position
            self.bound = len(self.strand.string) - 1

        if action == "swi":
            # only switch if there's a base there, otherwise quit
            if self.second_strand.string[self.bound] == " ":
                self.bound = -1
            else:
                # swap which ones are the "strand" and "second strand"
                self.second_strand, self.strand = self.strand, self.second_strand
                # we need to switch left and right as well
                self.strand.string = self.strand.string[::-1]
                self.second_strand.string = self.second_strand.string[::-1]
                # and figure out our new 'bound' location
                self.bound = len(self.strand.string) - self.bound - 1

        if action == "mvr":
            self.move_right()

        if action == "mvl":
            self.move_left()

        if action == "rpy":
            self.move_right()

            # Move to the nearest pyrimidine (T or C) to the right
            if self.bound >= len(self.strand.string) - 1:
                pass
            elif (
                "T" in self.strand.string[self.bound + 1:]
                or "C" in self.strand.string[self.bound + 1:]
            ):
                while (
                    self.strand.string[self.bound] != "T"
                    and self.strand.string[self.bound] != "C"
                    and self.strand.string[self.bound] != " "
                ):
                    self.move_right()
                if self.strand.string[self.bound] == " ":
                    self.bound = -1
        if action == "rpu":
            # Move to the nearest purine (A or G) to the right
            self.move_right()
            if self.bound >= len(self.strand.string) - 1:
                pass
            elif (
                "A" in self.strand.string[self.bound + 1:]
                or "G" in self.strand.string[self.bound + 1:]
            ):
                while (
                    self.strand.string[self.bound] != "A"
                    and self.strand.string[self.bound] != "G"
                    and self.strand.string[self.bound] != " "
                ):
                    self.move_right()
                if self.strand.string[self.bound] == " ":
                    self.bound = -1
        if action == "lpy":
            # Move to the nearest pyrimidine (T or C) to the left
            if (
                "T" in self.strand.string[: self.bound]
                or "C" in self.strand.string[: self.bound]
            ):
                self.move_left()
                while (
                    self.strand.string[self.bound] != "T"
                    and self.strand.string[self.bound] != "C"
                    and self.strand.string[self.bound] != " "
                ):
                    self.move_left()
                if self.strand.string[self.bound] == " ":
                    self.bound = -1
        if action == "lpu":
            self.move_left()
            # Move to the nearest purine (A or G) to the left
            if (
                "A" in self.strand.string[: self.bound]
                or "G" in self.strand.string[: self.bound]
            ):
                while (
                    self.strand.string[self.bound] != "A"
                    and self.strand.string[self.bound] != "G"
                    and self.strand.string[self.bound] != " "
                ):
                    self.move_left()
                if self.strand.string[self.bound] == " ":
                    self.bound = -1


def create_enzyme_s(strand):
    """
    Perform the ribosome operation, i.e. given a DNA strand,
    produce one or more enzymes from it.
    """
    enzymes = []
    seq = []
    dirs = []
    for entry in strand.codons:
        if len(entry) == 1:
            break
        amino = AMINO_DICT[entry]["amino"]
        dir = AMINO_DICT[entry]["direction"]
        if amino == "punctuation":
            enzymes.append(Enzyme(seq, dirs))
            seq = []
            dirs = []
        else:
            seq.append(amino)
            dirs.append(dir)
    if seq != []:
        enzymes.append(Enzyme(seq, dirs))
    return enzymes


def create_random_strand(min_length=10, max_length=36):
    strand = ""
    choices = ["A", "T", "G", "C"]
    length = np.random.randint(min_length, max_length)
    for i in range(length):
        strand += np.random.choice(choices)

    return strand


def identify_duplicates(the_list):
    seen = set()
    duplicates = []
    for i in the_list:
        if i in seen:
            duplicates.append(i)
        else:
            seen.add(i)
    return set(duplicates)


def run_a_set_of_generations(starting_string, max_iterations, verbose):
    """
    Run the typogenetic "Central Dogma" with a starting-string
    """
    # Start with your starting strand
    if starting_string is None:
        starting_string = create_random_strand(min_length=12)

    strand = Strand(starting_string)
    strand.derive_codons()
    strands = [strand]

    strings_per_generation = {}
    broken = False
    print(f"Starting string: |{starting_string}|")
    for i in range(max_iterations):
        if verbose:
            print("Generation: ", i)
        # Apply ribosome action to that strand
        enzymes = []
        for strand in strands:
            enzymes.extend(create_enzyme_s(strand))

        if verbose:
            print(f"Created {len(enzymes)} enzymes from {
                  len(strands)} strands")
        # Act upon the strand with the resulting enzyme(s)
        # This nested loop resolves the ambiguity of which enzyme acts on which strand
        for enzyme in enzymes:
            temp_strands = []
            for strand in strands:
                operation = Operation(enzyme, strand, verbose)
                temp_strands.extend(operation.final_strings)
                strands = [Strand(s) for s in temp_strands]

        strings_per_generation[i] = temp_strands
        strands = [Strand(s) for s in temp_strands]

        if verbose:
            print(f"End of Generation {i}. {len(strands)} strands created:")
        for strand in strands:
            strand.derive_codons()

        if len(strands) > 20 or len(enzymes) > 10:
            print("Too many strands detected")
            broken = True
            break

    for i in range(max_iterations):
        if verbose:
            print(f"Generation {i}: ", strings_per_generation[i])

    # eval
    loop_found = False
    if not broken:
        for i in range(max_iterations)[::-1]:
            # We need to check for a duplicate - and whether that duplicate
            # has existed in a previous generation
            dupes = identify_duplicates(strings_per_generation[i])
            if len(dupes):
                for dupe in dupes:
                    # search backwards
                    for j in [0]:
                        if dupe == starting_string:
                            print("Loop found! with dupe: |" + dupe + "|")
                            print("Starting string was: |" +
                                  starting_string + "|")
                            loop_found = True
    if not loop_found:
        print("No loops found")
    return loop_found


@click.command()
@click.option("--starting-string", default=None, help="Starting typogenetic string.")
@click.option(
    "--max-iterations", default=10, help="Maximum number of iterations to run."
)
@click.option(
    "--verbose",
    default=False,
    help="Whether to display the full sequence of iterations.",
)
def run_iterations(starting_string, max_iterations, verbose):
    run_a_set_of_generations(starting_string, max_iterations, verbose)


@click.command()
@click.option(
    "--max-iterations", default=10, help="Maximum number of iterations to run."
)
@click.option("--n", default=1000, help="How many starting strings to try.")
def run_many_iterations(max_iterations, n):
    for k in range(n):
        loop_found = run_a_set_of_generations(
            None, max_iterations=max_iterations, verbose=False
        )
        if loop_found:
            break


if __name__ == "__main__":
    run_many_iterations()
