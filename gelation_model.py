import random as rnd


class Polymer:
    """Polymer class with a name and number of hydroxyl groups."""
    def __init__(self, name, hydroxyls):
        if not isinstance(name, str):
            raise TypeError("Only strings are accepted")
        if not isinstance(hydroxyls, int):
            raise TypeError("Only integers are accepted")
        self.name = name
        self.hydroxyls = hydroxyls


class Kneader:
    """
    Kneader class that models the crosslinking reaction in a stochastic fashion between polymers.

    The kneader class accepts polymers as feed using the specify_polymer and del_polymer methods to add or remove polymers. The crosslinking reaction is modelled using the crosslink method.
    """
    def __init__(self):
        self.polymers = dict()  # polymers and amounts
        self.blob_count = []  # each group of crosslinked polymers
        self.chain_count = []  # number of chains of each polymer in each blob
        self.hydroxyl_count = []  # number of hydroxyl groups of each polymer in each blob
        self.polymer_order = []  # order of polymer being crosslinked

    def _check_polymer(self, polymer):
        """Check polymer class."""
        if not isinstance(polymer, Polymer):
            raise TypeError("Only polymers are accepted")
        return

    def specify_polymer(self, polymer, amount):
        """Add an amount of a type of polymer."""
        self._check_polymer(polymer)
        self.polymers[polymer.name] = {
            "hydroxyls": polymer.hydroxyls,
            "amount": amount,
        }
        return

    def del_polymer(self, polymer):
        """Remove a polymer."""
        self._check_polymer(polymer)
        del(self.polymers[polymer.name])
        return

    def crosslink(self, crosslink_count):
        """Model the crosslinking reaction stochastically."""
        self._init_crosslink()
        for each in crosslink_count:
            for _ in range(int(each[2])):  # No. of crosslinks
                old_blob = -1
                for i in range(2):  # 2 polymers getting crosslinked
                    blob_idx, poly_idx = self._rand_blob(each[i])
                    old_blob, empty_pos = self._move_blob(blob_idx, poly_idx, old_blob)
                    if empty_pos:
                        self._collate_blobs(old_blob)
        return

    def _init_crosslink(self):
        """Initialise variables."""
        self.blob_count = []
        self.chain_count = []
        self.hydroxyl_count = []
        self.polymer_order = []
        self.polymer_count = len(self.polymers)
        for idx, poly in enumerate(self.polymers):
            blob = self.polymers[poly]["amount"]
            chain = [0]*self.polymer_count
            chain[idx] += 1
            hyd = [0]*self.polymer_count
            hyd[idx] += self.polymers[poly]["hydroxyls"]
            self.blob_count.append(blob)
            self.chain_count.append(chain)
            self.hydroxyl_count.append(hyd)
            self.polymer_order.append(poly)
        return

    def _rand_blob(self, polymer):
        """Return a random polymer chain chosen by randomly choosing a hydroxyl group."""
        poly_idx = self.polymer_order.index(polymer.name)
        hydroxyl_tot = sum(
            x[0]*x[1][poly_idx] for x in zip(
                self.blob_count,
                self.hydroxyl_count,
            )
        )
        hydroxyl_idx = int(rnd.getrandbits(100) % hydroxyl_tot + 1)
        blob_idx = self._find_blob(
            poly_idx,
            hydroxyl_idx,
            len(self.blob_count),
        )
        return blob_idx, poly_idx

    def _find_blob(self, poly_idx, hyd_idx, length):
        """Binary search for the index of a polymer based on the index of its hydroxyl group."""
        left = 0
        right = length
        while left < right:
            mid = left + (right-left)//2
            left_sum = sum(
                x[0]*x[1][poly_idx] for x in zip(
                    self.blob_count[:mid],
                    self.hydroxyl_count[:mid],
                )
            )
            if hyd_idx > left_sum:
                left = mid+1
            else:
                right = mid
        return right-1

    def _move_blob(self, blob_idx, poly_idx, old_blob):
        """Crosslink the two polymers."""
        empty_pos = False
        if self.blob_count[blob_idx] == 1 and (
                blob_idx == old_blob or old_blob == -1
        ):
            # Intramolecular cross-linking
            hyd = [0]*self.polymer_count
            hyd[poly_idx] -= 1
            self.hydroxyl_count[blob_idx] = [
                sum(x) for x in zip(
                    self.hydroxyl_count[blob_idx],
                    hyd,
                )
            ]
            return blob_idx, empty_pos
        if old_blob == -1 and self.blob_count[blob_idx] != 1:
            # Create ghost blob
            self.blob_count.append(1)
            self.chain_count.append([0]*self.polymer_count)
            self.hydroxyl_count.append([0]*self.polymer_count)
            old_blob = len(self.blob_count)-1
        self.chain_count[old_blob] = [
            sum(x) for x in zip(
                self.chain_count[blob_idx],
                self.chain_count[old_blob],
            )
        ]
        hyd = [0]*self.polymer_count
        hyd[poly_idx] -= 1
        self.hydroxyl_count[old_blob] = [
            sum(x) for x in zip(
                self.hydroxyl_count[blob_idx],
                self.hydroxyl_count[old_blob],
                hyd,
            )
        ]
        self.blob_count[blob_idx] -= 1
        if not self.blob_count[blob_idx]:  # clean up empty blob
            del(self.blob_count[blob_idx])
            del(self.chain_count[blob_idx])
            del(self.hydroxyl_count[blob_idx])
            empty_pos = True
            if old_blob > blob_idx:
                old_blob -= 1
        return old_blob, empty_pos

    def _collate_blobs(self, old_blob):
        """Order the list of polymer chains according to the number of each polymer types."""
        equal_list = list(
            id for id, x in enumerate(self.hydroxyl_count) if self.hydroxyl_count[id] == self.hydroxyl_count[old_blob] and id != old_blob
        )
        if not equal_list:
            return
        idx = equal_list[0]
        self.blob_count[idx] += 1
        del(self.blob_count[old_blob])
        del(self.chain_count[old_blob])
        del(self.hydroxyl_count[old_blob])
        return


# Example - Crosslinking of carboxymethyl cellulose (CMC) with starch (S)
points = [  # [reaction time, degree of crosslinking]
    [0.5, 3.5520782E-6],
    [1, 1.5063633E-5],
    [1.5, 3.453642E-5],
    [2, 6.197044E-5],
    [2.5, 9.7365686E-5],
    [3, 1.4072216E-4],
    [3.5, 1.9203988E-4],
    [4, 2.5081958E-4],
    [4.5, 3.1681114E-4],
    [5, 3.9050024E-4],
    [5.5, 4.7188686E-4],
    [6, 5.6097104E-4],
    [6.5, 6.577527E-4],
    [7, 7.622319E-4],
    [7.5, 8.744086E-4],
    [8, 9.936807E-4],
    [8.5, 0.0011203686],
    [9, 0.0012544971],
    [9.5, 0.0013960659],
]
for i in range(len(points)):
    frac_CMC = 0.160  # fraction of CMC-CMC crosslinks
    frac_CS = 0.288  # fraction of CMC-S crosslinks
    frac_SS = 0.552  # fraction of S-S crosslinks
    DC_total = points[i][1]
    CMC_hyd = 1000 * 2  # no. of hydroxyl groups per CMC chain
    Starch_hyd = 960 * 3  # no. of hydroxyl groups per starch chain
    total_mass = 100000  # total number of polymer chains
    wt_CMC = 0.7  # mass fraction of CMC
    mw_CMC = 240.19  # molecular weight of CMC
    mw_S = 162.14  # molecular weight of starch
    CMC_amt = total_mass * wt_CMC / mw_CMC
    Starch_amt = total_mass * (1-wt_CMC) / mw_S
    DC_CMC = DC_total * frac_CMC
    DC_CS = DC_total * frac_CS
    DC_starch = DC_total * frac_SS

    b = []  # list record of biggest blobs
    for j in range(10):  # number of stochastic simulations
        CMC = Polymer("CMC", CMC_hyd)
        Starch = Polymer("Starch", Starch_hyd)
        R301 = Kneader()
        R301.specify_polymer(CMC, CMC_amt)
        R301.specify_polymer(Starch, Starch_amt)
        R301.crosslink(
            [
                [CMC, CMC, int(CMC_hyd*CMC_amt*DC_CMC)],
                [CMC, Starch, int((CMC_hyd*CMC_amt + Starch_hyd*Starch_amt)*DC_CS)],
                [Starch, Starch, int(Starch_hyd*Starch_amt*DC_starch)],
            ]
        )
        size = max(sum(x) for x in R301.chain_count)  # biggest blob size
        b.append(list(x for x in R301.chain_count if sum(x) == size))  # add all blobs of that size
    number = sum(len(x) for x in b)  # number of blobs
    c = sum(sum(x[0] for x in y) for y in b)
    s = sum(sum(x[1] for x in y) for y in b)
    print(f"Data for {points[i][0]}: CMC - {c/number}, S - {s/number}")  # average CMC/S per largest blob
