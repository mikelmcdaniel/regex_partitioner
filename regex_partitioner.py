"""Module for partitioning a space defined by a regular expression.

For an NFA, maximum length, target ratio, tolerance ratio, lower bound and
upper bound, prints a list of  sequences, that are not necessarily accepted
by the NFA, that partitions the space of all sequences lexicographically
between the lower and upper bounds that the NFA accepts with a maximum length
such the target ratio +- the tolerance ratio of the sequences are
lexicographically less than result.
"""
import argparse
import collections
import copy
import fractions
import itertools
import re
import string
import sys

from typing import Dict, Iterable, Iterator, FrozenSet, List, Optional, Set, Text, Tuple


def sre_parse(regex_str: Text) -> Tuple:
    return re.sre_parse.parse(regex_str)  # type: ignore


def num_seqs_with_max_len(alphabet_size: int, max_len: int) -> int:
    assert max_len >= 0
    return (alphabet_size ** (max_len + 1) - 1) // (alphabet_size - 1)


def num_to_seq(num: int, alphabet: List[Text], max_len: int) -> Iterator[Text]:
    assert max_len >= 0
    num = int(num)
    len_seqs_less = max_len - 1
    while num:
        num_seqs_less = num_seqs_with_max_len(len(alphabet), len_seqs_less)
        yield alphabet[(num - 1) // num_seqs_less]
        num = (num - 1) % num_seqs_less
        len_seqs_less -= 1


def seq_to_num(seq: Iterable[Text], inverse_alphabet: Dict[Text, int], max_len: int) -> int:
    assert max_len >= 0
    len_seq_remaining = max_len - 1
    num = 0
    for c in itertools.islice(seq, 0, max_len):
        num_seqs_less = num_seqs_with_max_len(len(inverse_alphabet), len_seq_remaining)
        num += 1 + inverse_alphabet[c] * num_seqs_less
        len_seq_remaining -= 1
    return num


class NFA(object):
    """Class representing a Nondeterministic Finite Automata."""

    def __init__(self):
        self.nodes: Dict[int, Dict[Optional[Text], Set[int]]] = collections.defaultdict(
            lambda: collections.defaultdict(set)
        )
        self.start_nodes: Set[int] = set()
        self.accept_nodes: Set[int] = set()
        self.alphabet: List[Text] = []  # ['a', 'b', 'c', ..., 'z']
        self.inverse_alphabet: Dict[Text, int] = {}  # {'a': 0, 'b': 1, ..., 'z': 25}
        self._m = {}  # a cache for self.next_nodes(...)

    def add_node(self) -> int:
        new_node_id = len(self.nodes)
        # self.nodes is a default dictionary, so accessing it implicitly
        # creates an entry which can be iterated over.
        self.nodes[new_node_id]  # pylint: disable=pointless-statement
        return new_node_id

    def add_transition(self, from_node: int, to_node: int, transition: Optional[Text] = None) -> None:
        self.nodes[from_node][transition].add(to_node)
        if self._m:
            self._m = {}
        if transition is not None and transition not in self.inverse_alphabet:
            self.inverse_alphabet[transition] = len(self.alphabet)
            self.alphabet.append(transition)

    def possible_transitions(self, nodes: Iterable[int]) -> Set[Optional[Text]]:
        ts: Set[Optional[Text]] = set()
        for node in nodes:
            ts.update(self.nodes[node])
        ts.discard(None)
        return ts

    def _next_nodes_slow(self, cur_nodes: Iterable[int], sequence_element: Optional[Text]) -> FrozenSet[int]:
        """Returns the next nodes from the current nodes, given a transition.

        Args:
          cur_nodes: iterable, the current nodes
          sequence_element: object, the next item in the sequence, for a transition

        Returns:
          frozen_set: the next nodes
        """
        next_nodes: Set[int] = set()
        for cur_node in cur_nodes:
            next_nodes.update(self.nodes[cur_node].get(sequence_element, ()))
        # Follow epsilon transitions
        last_updates = next_nodes
        while last_updates:
            updates = set()
            for next_node in last_updates:
                if next_node in self.nodes:
                    for update_node in self.nodes[next_node][None]:
                        if update_node not in next_nodes:
                            updates.add(update_node)
            if not updates:
                break
            next_nodes.update(updates)
            last_updates = updates
        return frozenset(next_nodes)

    def next_nodes(self, cur_nodes: Iterable[int], sequence_element: Optional[Text]) -> FrozenSet[int]:
        cur_nodes = frozenset(cur_nodes)
        if (cur_nodes, sequence_element) not in self._m:
            self._m[cur_nodes, sequence_element] = self._next_nodes_slow(cur_nodes, sequence_element)
        return self._m[cur_nodes, sequence_element]

    def accepts(self, sequence: Iterable[Text]) -> bool:
        cur_nodes: FrozenSet[int] = self.start_nodes
        for item in sequence:
            cur_nodes = self.next_nodes(cur_nodes, item)
        return not self.accept_nodes.isdisjoint(cur_nodes)

    def prev_accepted(self, sequence: Iterable[Text], max_len: int) -> Optional[Text]:
        """Returns the largest string up to max_len characters accepted by this NFA less than sequence."""
        desired_num_accepted = self.num_accepts(max_len, sequence) + 1
        hi = seq_to_num(sequence, self.inverse_alphabet, max_len)

        # We don't have a lower bound on where the next_accepted string is, so we look at an exponentially increasing
        # gap above low until we find at least one accepted string.
        diff = 1
        while True:
            lo = hi - diff
            if lo <= 0:
                lo = 0
                break
            lo_seq = "".join(num_to_seq(lo, self.alphabet, max_len))
            lo_num_accepted = self.num_accepts(max_len, lo_seq)
            if lo_num_accepted > desired_num_accepted:
                break
            diff *= 2

        while lo <= hi:
            mid = (lo + hi) // 2
            mid_seq = "".join(num_to_seq(mid, self.alphabet, max_len))
            mid_num_accepted = self.num_accepts(max_len, mid_seq)
            if mid_num_accepted < desired_num_accepted:
                hi = mid - 1
            elif mid_num_accepted == desired_num_accepted and self.accepts(mid_seq):
                return mid_seq
            else:  # elif mid_num_accepted > desired_num_accepted:
                lo = mid + 1
        return None  # There is no previous sequence!

    def next_accepted(self, sequence: Iterable[Text], max_len: int) -> Optional[Text]:
        """Returns the smallest string up to max_len characters accepted by this NFA greater than sequence."""
        max_hi = num_seqs_with_max_len(len(self.alphabet), max_len)
        desired_num_accepted = self.num_accepts(max_len, sequence) - self.accepts(sequence)
        lo = seq_to_num(sequence, self.inverse_alphabet, max_len) + 1

        # We don't have an upper bound on where the next_accepted string is, so we look at an exponentially increasing
        # gap above low until we find at least one accepted string.
        diff = 1
        while True:
            hi = lo + diff
            if hi >= max_hi:
                hi = max_hi - 1
                break
            hi_seq = "".join(num_to_seq(hi, self.alphabet, max_len))
            hi_num_accepted = self.num_accepts(max_len, hi_seq)
            if hi_num_accepted < desired_num_accepted:
                break
            lo = hi
            diff *= 2

        # Now we know that the string we're looking for is [lo, hi]. The upper bound is inclusive!
        while lo <= hi:
            mid = (lo + hi) // 2
            mid_seq = "".join(num_to_seq(mid, self.alphabet, max_len))
            mid_num_accepted = self.num_accepts(max_len, mid_seq)
            if mid_num_accepted < desired_num_accepted:
                hi = mid - 1
            elif mid_num_accepted == desired_num_accepted and self.accepts(mid_seq):
                return mid_seq
            else:  # elif mid_num_accepted > desired_num_accepted:
                lo = mid + 1
        return None  # There is no next sequence!

    def _sum_tables(self, table: Dict[FrozenSet[int], int]) -> int:
        return sum(count for nodes, count in table.items() if any(node in self.accept_nodes for node in nodes))

    def num_accepts(self, max_len: int, bound: Iterable[Text] = ()) -> int:
        """Returns the number of sequences accepted by this NFA.

        Return the number of sequences with length less than or equal to max_len
        that are lexicographically greater than or equal to bound.

        Args:
          max_len: int, maximum length of sequences to consider
          bound: sequence below which, no sequences are considered

        Returns:
          int: number of sequences accepted considering max_len and bound
        """
        eq1: Dict[FrozenSet[int], int] = collections.defaultdict(int)
        eq2: Dict[FrozenSet[int], int] = collections.defaultdict(int)
        gt1: Dict[FrozenSet[int], int] = collections.defaultdict(int)
        gt2: Dict[FrozenSet[int], int] = collections.defaultdict(int)
        eq1[frozenset(self.start_nodes)] = 1
        result = 0
        c: Optional[Text] = None
        for _, c in zip(range(max_len), itertools.chain(bound, itertools.repeat(None))):
            for nodes, count in gt1.items():
                for element in self.possible_transitions(nodes):
                    next_nodes = frozenset(self.next_nodes(nodes, element))
                    gt2[next_nodes] += count
            for nodes, count in eq1.items():
                for element in self.possible_transitions(nodes):
                    next_nodes = frozenset(self.next_nodes(nodes, element))
                    if c is None or (element is not None and element > c):
                        gt2[next_nodes] += count
                    elif element == c:
                        eq2[next_nodes] += count
            result += self._sum_tables(gt2)
            eq1, eq2 = eq2, collections.defaultdict(int)
            gt1, gt2 = gt2, collections.defaultdict(int)
        try:
            len_bound: int = len(bound)
        except TypeError:
            len_bound = sum(1 for _ in bound)
        return result + (len_bound <= max_len and self.accepts(bound))

    def ensure_disjoint(self, other_nfa: "NFA") -> None:
        offset = len(other_nfa.nodes)
        self.start_nodes = set(n + offset for n in self.start_nodes)
        self.accept_nodes = set(n + offset for n in self.accept_nodes)
        nodes: Dict[int, Dict[Optional[Text], Set[int]]] = collections.defaultdict(lambda: collections.defaultdict(set))
        for node, transitions in self.nodes.items():
            nodes[node + offset] = collections.defaultdict(set)
            for transition, next_nodes in transitions.items():
                nodes[node + offset][transition] = set(n + offset for n in next_nodes)
        self.nodes = nodes

    def update_nodes(self, other_nfa: "NFA") -> None:
        for k, v in other_nfa.nodes.items():
            self.nodes[k] = v
        self.alphabet = sorted(set(self.alphabet + other_nfa.alphabet))
        self.inverse_alphabet = dict((v, k) for k, v in enumerate(self.alphabet))

    def build_dot_str(self) -> Text:
        """Returns a string representation of this NFA as a GraphViz *.dot file."""
        s = []
        s.append("digraph {")
        for node in self.nodes:
            label = str(node)
            if node in self.start_nodes:
                label += "S"
            if node in self.accept_nodes:
                label += "A"
            s.append(f'  "{node}" [label="{label}"];')
        s.append("")
        for from_node, transitions in self.nodes.items():
            for transition, to_nodes in transitions.items():
                if not transition:
                    transition = "&epsilon;"
                for to_node in to_nodes:
                    s.append(f'  "{from_node}" -> "{to_node}" [label="{transition}"];')
        s.append("}")
        return "\n".join(s)


class RE(object):
    """Base Class for representing regular expressions."""

    def match(self, unused_sequence: Text) -> Optional[int]:
        return None

    def __repr__(self) -> Text:
        return f"{type(self).__name__}({self})"

    def as_nfa(self) -> NFA:
        return NFA()


class RESequence(RE):
    """Class representing regular expression of a literal (e.g. "abc")."""

    def __init__(self, sequence: Text):
        self.sequence: Text = sequence

    def match(self, sequence: Text) -> Optional[int]:
        if sequence.startswith(self.sequence):
            return len(self.sequence)
        return None

    def __repr__(self) -> Text:
        return f"{type(self).__name__}({self.sequence!r})"

    def __str__(self) -> Text:
        return str(self.sequence)

    def as_nfa(self) -> NFA:
        """Return the NFA representation of this regular expression."""
        nfa = NFA()
        prev_node = nfa.add_node()
        next_node = prev_node
        nfa.start_nodes.add(prev_node)
        for item in self.sequence:
            next_node = nfa.add_node()
            nfa.add_transition(prev_node, next_node, item)
            prev_node = next_node
        nfa.accept_nodes.add(next_node)
        return nfa


class REConcat(RE):
    """Class representing concatenated regular expressions (e.g. "abc")."""

    def __init__(self, sub_trees: Iterable[RE]):
        self.sub_trees = list(sub_trees)

    def match(self, sequence: Text) -> Optional[int]:
        sub_trees = iter(self.sub_trees)
        first_match: Optional[int] = None
        try:
            first_match = next(sub_trees).match(sequence)
        except StopIteration:
            return None
        for sub_tree in sub_trees:
            m = sub_tree.match(sequence)
            if m is None:
                return None
            sequence = sequence[m:]
        return first_match

    def __repr__(self) -> Text:
        return f"{type(self).__name__}({self.sub_trees!r})"

    def __str__(self) -> Text:
        return f'({"".join(str(st) for st in self.sub_trees)})'

    def as_nfa(self) -> NFA:
        if not self.sub_trees:
            # The concatenation of nothing is equal to the empty sequence.
            nfa = NFA()
            node = nfa.add_node()
            nfa.start_nodes.add(node)
            nfa.accept_nodes.add(node)
            return nfa
        else:
            nfa = self.sub_trees[0].as_nfa()
            for sub_tree in itertools.islice(self.sub_trees, 1, len(self.sub_trees)):
                sub_nfa = sub_tree.as_nfa()
                sub_nfa.ensure_disjoint(nfa)
                nfa.update_nodes(sub_nfa)
                connecting_node = nfa.add_node()
                for node in nfa.accept_nodes:
                    nfa.add_transition(node, connecting_node, None)
                for node in sub_nfa.start_nodes:
                    nfa.add_transition(connecting_node, node, None)
                nfa.accept_nodes = sub_nfa.accept_nodes
            return nfa


class REOr(RE):
    """Class representing a union regular expression (e.g. "a|b")."""

    def __init__(self, sub_trees: Iterable[RE]):
        self.sub_trees = list(sub_trees)

    def match(self, sequence: Text) -> Optional[int]:
        m = None
        for sub_tree in self.sub_trees:
            m = sub_tree.match(sequence)
            if m is not None:
                break
        return m

    def __repr__(self) -> Text:
        return f"{type(self).__name__}({self.sub_trees!r})"

    def __str__(self) -> Text:
        return f'({"|".join(str(st) for st in self.sub_trees)})'

    def as_nfa(self) -> NFA:
        assert self.sub_trees
        nfa = self.sub_trees[0].as_nfa()
        for sub_tree in itertools.islice(self.sub_trees, 1, len(self.sub_trees)):
            sub_nfa = sub_tree.as_nfa()
            sub_nfa.ensure_disjoint(nfa)
            nfa.update_nodes(sub_nfa)
            nfa.start_nodes.update(sub_nfa.start_nodes)
            nfa.accept_nodes.update(sub_nfa.accept_nodes)
        return nfa


class RERepeat(RE):
    """Class representing a repeating regular expression (e.g. "a*")."""

    def __init__(self, sub_tree: RE):
        self.sub_tree = sub_tree

    def match(self, sequence: Text) -> Optional[int]:
        total = 0
        m = self.sub_tree.match(sequence)
        while m is not None:
            total += m
            sequence = sequence[m:]
            m = self.sub_tree.match(sequence)
        return total

    def __repr__(self) -> Text:
        return f"{type(self).__name__}({self.sub_tree!r})"

    def __str__(self) -> Text:
        return f"({self.sub_tree})*"

    def as_nfa(self) -> NFA:
        nfa = self.sub_tree.as_nfa()
        connecting_node = nfa.add_node()
        for accept_node in nfa.accept_nodes:
            nfa.add_transition(accept_node, connecting_node, None)
        for start_node in nfa.start_nodes:
            nfa.add_transition(connecting_node, start_node, None)
        nfa.accept_nodes.update(nfa.start_nodes)
        return nfa


def find_partition_seq(
    nfa: NFA,
    max_len: int,
    target_ratio=fractions.Fraction(1, 2),
    low: Iterable[Text] = (),
    high: Optional[Iterable[Text]] = None,
    tolerance_ratio: float = 0.0,
) -> Tuple[Text, ...]:
    """Return a sequence that partitions the space of accepted strings.

    Returns a sequence, results, that is not necessarily accepted by nfa that
    partitions the space of all sequences lexicographically between lo and hi that
    nfa accepts of length max_len or less such target_ratio +- tolerance_ratio of
    the sequences are lexicographically less than result.

    Args:
      nfa: NFA, an arbitrary NFA.
      max_len: int, maximum length of sequences to consider
      target_ratio: number, ratio of sequences to be less than the partition
      lo: iterable, lower bound sequence, below which sequences are not considered
      hi: iterable, upper bound sequence, above which sequences are not considered
      tolerance_ratio: number, ratio of tolerance within the target_ratio

    Returns:
      iterable: sequence that partitions the NFA space according to constraints
    """

    max_letter = max(nfa.alphabet)
    lo: int = seq_to_num(low, nfa.inverse_alphabet, max_len)
    hi: int = seq_to_num((max_letter for _ in range(max_len)) if high is None else high, nfa.inverse_alphabet, max_len)
    lo_num_accepts = nfa.num_accepts(max_len, num_to_seq(lo, nfa.alphabet, max_len))
    hi_num_accepts = nfa.num_accepts(max_len, num_to_seq(hi, nfa.alphabet, max_len))
    total = lo_num_accepts - hi_num_accepts
    assert total >= 0
    target = lo_num_accepts - int(total * target_ratio)
    tolerance = int(total * tolerance_ratio)
    for j in itertools.cycle((0, 1)):
        if j == 0 and hi_num_accepts != lo_num_accepts:
            mid: int = lo + (hi - lo) * (lo_num_accepts - target) // (lo_num_accepts - hi_num_accepts)
            assert mid >= lo
            assert mid <= hi
        else:
            mid = (lo + hi) // 2
        mid_str = num_to_seq(mid, nfa.alphabet, max_len)
        mid_num_accepts = nfa.num_accepts(max_len, mid_str)
        if lo >= hi or mid_num_accepts == target or abs(lo_num_accepts - hi_num_accepts) <= tolerance:
            break
        elif mid_num_accepts < target:
            hi = mid + 1
            hi_num_accepts = mid_num_accepts
        elif mid_num_accepts > target:
            lo = mid
            lo_num_accepts = mid_num_accepts
    return tuple(num_to_seq(mid, nfa.alphabet, max_len))


def find_partition_seqs(
    nfa: NFA,
    max_len: int,
    num_partitions: int = 1,
    lo: Iterable[Text] = (),
    hi: Optional[Iterable[Text]] = None,
    tolerance_ratio: float = 0.0,
) -> Iterable[Iterable[Text]]:
    return tuple(
        find_partition_seq(nfa, max_len, fractions.Fraction(j, num_partitions), lo, hi, tolerance_ratio)
        for j in range(1, num_partitions)
    )


MAX_REPEAT = sre_parse("a*")[0][1][1]
ALPHABET: FrozenSet[Text] = frozenset(string.printable)


def sre_to_re(sre, alphabet: FrozenSet[Text] = ALPHABET) -> RE:
    if isinstance(sre, tuple):
        op, args = sre
        op = str(op).lower()
        if op == "branch":  # or
            none, sres = args
            assert none is None  # I don't know how to interpret a non-None value.
            return REOr(map(sre_to_re, sres))
        elif op == "in":  # or
            sres = args
            return REOr(map(sre_to_re, sres))
        elif op == "category":  # or (\d, \W)
            if args == "category_digit":
                return REOr([RESequence(c) for c in string.digits if c in alphabet])
            elif args == "category_not_digit":
                return REOr([RESequence(d) for d in alphabet if d not in string.digits])
            elif args == "category_word":
                return REOr([RESequence(c) for c in string.ascii_letters if c in alphabet])
            elif args == "category_not_word":
                return REOr([RESequence(d) for d in alphabet if d not in string.ascii_letters])
            else:
                raise Exception(f'Unknown category type "{args!r}".')
        elif op == "assert":  # (?=REGEX)
            _, sres = args
            assert len(sres) == 2
            return sre_to_re(sres[1])
        elif op == "any":  # .
            return REOr([RESequence(c) for c in alphabet])
        elif op == "max_repeat" or op == "min_repeat":  # {m,n}, *, +
            min_repeat, max_repeat, sre = args
            res: List[RE] = []
            main_res = sre_to_re(sre)
            if max_repeat >= MAX_REPEAT:
                for _ in range(min_repeat):
                    res.append(copy.deepcopy(main_res))
                res.append(RERepeat(copy.deepcopy(main_res)))
            else:
                for _ in range(min_repeat):
                    res.append(copy.deepcopy(main_res))
                for _ in range(max_repeat - min_repeat):
                    res.append(REOr([copy.deepcopy(main_res), RESequence("")]))
            if len(res) > 1:
                return REConcat(res)
            else:
                return res[0]
            return RERepeat(res)
        elif op == "literal":  # A single character.
            return RESequence(chr(args))
        elif op == "not_literal":  # Any character,  except this single character.
            literal = chr(args)
            return REOr([RESequence(c) for c in alphabet if c != literal])
        elif op == "subpattern":
            _, _, _, sres = args
            return sre_to_re(sres)
        elif op == "range":
            min_literal, max_literal = args
            return REOr([RESequence(chr(c)) for c in range(min_literal, max_literal + 1)])
        elif op == "at":  # Anchor characters (e.g. "^" and "$")
            return RESequence("")
        else:
            raise Exception(f"Unknown op, args pair: {op!r}, {args!r}")
    else:  # concat
        if len(sre) == 1:
            return sre_to_re(sre[0])
        else:
            return REConcat(map(sre_to_re, sre))


def regex_str_to_re(regex_str: Text) -> RE:
    sre = sre_parse(regex_str)
    return sre_to_re(sre)


def main(argv: List[Text]) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("regex", type=str)
    parser.add_argument("-p", "--partitions", default=4, dest="num_partitions", type=int)
    parser.add_argument("-l", "--lower_bound", default="", dest="lower_bound", type=str)
    parser.add_argument("-u", "--upper_bound", default=None, dest="upper_bound")
    parser.add_argument("-t", "--tolerance_ratio", default=0.001, dest="tolerance_ratio", type=float)
    parser.add_argument("-m", "--max_str_len", default=10, dest="max_len", type=int)
    args = parser.parse_args(argv[1:])

    num_partitions = args.num_partitions
    regex_str = args.regex
    lo = args.lower_bound
    hi = args.upper_bound
    tolerance_ratio = args.tolerance_ratio
    max_len = args.max_len

    nfa = regex_str_to_re(regex_str).as_nfa()

    total_num_accepts = nfa.num_accepts(max_len)
    print(f"The regular expression {regex_str!r} accepts {total_num_accepts} strings with length <= {max_len}.")
    if lo or hi:
        print(
            f"Between {lo!r} and {hi!r}, {regex_str!r} accepts "
            f"{nfa.num_accepts(max_len, lo) - nfa.num_accepts(max_len, hi)} with len <= {max_len}."
        )
    print()
    print(f"These strings can be evenly split, within a {100 * tolerance_ratio:f}% tolerance by the strings:")
    for raw_seq in find_partition_seqs(nfa, max_len, num_partitions, lo, hi, tolerance_ratio):
        # For even nicer output, we use nfa.prev_accepted(...) to get a partition that the regex accepts.
        seq = "".join(raw_seq)
        accepted_seq = "".join(raw_seq if nfa.accepts(raw_seq) else nfa.next_accepted(raw_seq, max_len) or raw_seq)
        num_accepts = nfa.num_accepts(max_len, seq)
        print(
            f"  {regex_str!r} accepts {num_accepts} ({100 * num_accepts / total_num_accepts}%) "
            f"of strings after {accepted_seq!r} (or {seq!r})."
        )


if __name__ == "__main__":
    main(sys.argv)
