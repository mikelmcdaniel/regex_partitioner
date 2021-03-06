"""Module for partitioning a space defined by a regular expression.

For an NFA, maximum length, target ratio, tolerance ratio, lower bound and
upper bound, prints a list of  sequences, that are not necessarily accepted
by the NFA, that partitions the space of all sequences lexicographically
between the lower and upper bounds that the NFA accepts with a maximum length
such the target ratio +- the tolerance ratio of the sequences are
lexicographically less than result.
"""
import collections
import copy
import fractions
import itertools
import re
import string
import sys


def NumSeqsWithMaxLen(alphabet_size, max_len):
  assert max_len >= 0
  return (alphabet_size**(max_len + 1) - 1) / (alphabet_size - 1)


def NumToSeq(num, alphabet, max_len):
  assert max_len >= 0
  len_seqs_less = max_len - 1
  while num:
    num_seqs_less = NumSeqsWithMaxLen(len(alphabet), len_seqs_less)
    yield alphabet[(num - 1) / num_seqs_less]
    num = (num - 1) % num_seqs_less
    len_seqs_less -= 1


def SeqToNum(seq, inverse_alphabet, max_len):
  assert max_len >= 0
  len_seq_remaining = max_len - 1
  num = 0
  for c in itertools.islice(seq, 0, max_len):
    num_seqs_less = NumSeqsWithMaxLen(
        len(inverse_alphabet), len_seq_remaining)
    num += 1 + inverse_alphabet[c] * num_seqs_less
    len_seq_remaining -= 1
  return num


def TestSeqNum():
  for alphabet in ('ab', 'abc', 'abcd', 'abcde', 'abcdef'):
    inverse_alphabet = dict((v, k) for k, v in enumerate(alphabet))
    def AllSeqs(alphabet, max_len):
      def AllSeqsHelper(alphabet, max_len, seq):
        yield seq
        if max_len > 0:
          for _ in AllSeqsHelper(alphabet, max_len - 1, seq):
            for letter in alphabet:
              seq.append(letter)
              yield seq
              seq.pop()
      for s in AllSeqsHelper(alphabet, max_len, []):
        yield ''.join(s)

    for max_len in xrange(7):
      all_sequences = sorted(AllSeqs(alphabet, max_len))
      for expected_num, expected_str in enumerate(all_sequences):
        actual_num = SeqToNum(expected_str, inverse_alphabet, max_len)
        actual_str = ''.join(NumToSeq(expected_num, alphabet, max_len))
        try:
          assert expected_num == actual_num
          assert expected_str == actual_str
        except AssertionError:
          print 'Expected {} == {}'.format(expected_num, actual_num)
          print 'Expected {} == {}'.format(expected_str, actual_str)
          assert False


class NFA(object):
  """Class representing a Nondeterministic Finite Automata."""

  def __init__(self):
    # self.nodes[char or None] -> next_nodes
    self.nodes = collections.defaultdict(lambda: collections.defaultdict(set))
    self.start_nodes = set()
    self.accept_nodes = set()
    self.alphabet = []  # ['a', 'b', 'c', ..., 'z']
    self.inverse_alphabet = {}  # {'a': 0, 'b': 1, ..., 'z': 25}
    self._m = {}

  def AddNode(self):
    new_node_id = len(self.nodes)
    # self.nodes is a default dictionary, so accessing it implicitly
    # creates an entry which can be iterated over.
    self.nodes[new_node_id]  # pylint: disable=pointless-statement
    return new_node_id

  def AddTransition(self, from_node, to_node, transition=None):
    self.nodes[from_node][transition].add(to_node)
    if self._m:
      self._m = {}
    if transition is not None and transition not in self.inverse_alphabet:
      self.inverse_alphabet[transition] = len(self.alphabet)
      self.alphabet.append(transition)

  def PossibleTransitions(self, nodes):
    ts = set()
    for node in nodes:
      ts.update(self.nodes[node])
    ts.discard(None)
    return ts

  def _NextNodesSlow(self, cur_nodes, sequence_element):
    """Returns the next nodes from the current nodes, given a transition.

    Args:
      cur_nodes: iterable, the current nodes
      sequence_element: object, the next item in the sequence, for a transition

    Returns:
      frozen_set: the next nodes
    """
    next_nodes = set()
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

  def NextNodes(self, cur_nodes, sequence_element):
    cur_nodes = frozenset(cur_nodes)
    if (cur_nodes, sequence_element) not in self._m:
      self._m[cur_nodes, sequence_element] = self._NextNodesSlow(
          cur_nodes, sequence_element)
    return self._m[cur_nodes, sequence_element]

  def Accepts(self, sequence):
    cur_nodes = self.start_nodes
    for item in sequence:
      cur_nodes = self.NextNodes(cur_nodes, item)
    return not self.accept_nodes.isdisjoint(cur_nodes)

  def _SumTables(self, table):
    return sum(count
               for nodes, count in table.iteritems()
               if any(node in self.accept_nodes for node in nodes))

  def NumAccepts(self, max_len, bound=()):
    """Returns the number of sequences accepted by this NFA.

    Return the number of sequences with length less than or equal to max_len
    that are lexicographically greater than or equal to bound.

    Args:
      max_len: int, maximum length of sequences to consider
      bound: sequence below which, no sequences are considered

    Returns:
      int: number of sequences accepted considering max_len and bound
    """

    bound = iter(bound)
    eq1, eq2 = collections.defaultdict(int), collections.defaultdict(int)
    gt1, gt2 = collections.defaultdict(int), collections.defaultdict(int)
    eq1[frozenset(self.start_nodes)] = 1
    result = 0
    try:
      c = bound.next()
      has_next = True
    except StopIteration:
      has_next = False
    for _ in xrange(max_len):
      if not has_next:
        result += self._SumTables(eq1)
      result += self._SumTables(gt1)
      for nodes, count in gt1.iteritems():
        for element in self.PossibleTransitions(nodes):
          next_nodes = frozenset(self.NextNodes(nodes, element))
          gt2[next_nodes] += count
      for nodes, count in eq1.iteritems():
        for element in self.PossibleTransitions(nodes):
          next_nodes = frozenset(self.NextNodes(nodes, element))
          if not has_next or element > c:
            gt2[next_nodes] += count
          elif element == c:
            eq2[next_nodes] += count
      eq1, eq2 = eq2, collections.defaultdict(int)
      gt1, gt2 = gt2, collections.defaultdict(int)
      try:
        c = bound.next()
        has_next = True
      except StopIteration:
        has_next = False
    if not has_next:
      result += self._SumTables(eq1)
    result += self._SumTables(gt1)
    return result

  def EnsureDisjoint(self, other_nfa):
    offset = len(other_nfa.nodes)
    self.start_nodes = set(n + offset for n in self.start_nodes)
    self.accept_nodes = set(n + offset for n in self.accept_nodes)
    nodes = collections.defaultdict(lambda: collections.defaultdict(set))
    for node, transitions in self.nodes.iteritems():
      nodes[node + offset] = collections.defaultdict(set)
      for transition, next_nodes in transitions.iteritems():
        nodes[node + offset][transition] = set(n + offset for n in next_nodes)
    self.nodes = nodes

  def UpdateNodes(self, other_nfa):
    for k, v in other_nfa.nodes.iteritems():
      self.nodes[k] = v
    self.alphabet = sorted(set(self.alphabet + other_nfa.alphabet))
    self.inverse_alphabet = dict((v, k) for k, v in enumerate(self.alphabet))

  def BuildDotStr(self):
    """Returns a string representation of this NFA as a GraphViz *.dot file."""
    s = []
    s.append('digraph {')
    for node in self.nodes:
      label = str(node)
      if node in self.start_nodes:
        label += 'S'
      if node in self.accept_nodes:
        label += 'A'
      s.append('  "{}" [label="{}"];'.format(node, label))
    s.append('')
    for from_node, transitions in self.nodes.iteritems():
      for transition, to_nodes in transitions.iteritems():
        if not transition:
          transition = '&epsilon;'
        for to_node in to_nodes:
          s.append('  "{}" -> "{}" [label="{}"];'.format(
              from_node, to_node, transition))
    s.append('}')
    return '\n'.join(s)


class RE(object):
  """Base Class for representing regular expressions."""

  def Match(self, unused_sequence):
    return None

  def __repr__(self):
    return '{}({})'.format(type(self).__name__, str(self))


class RESequence(RE):
  """Class representing regular expression of a literal (e.g. "abc")."""

  def __init__(self, sequence):
    self.sequence = sequence

  def Match(self, sequence):
    if sequence.startswith(self.sequence):
      return len(self.sequence)
    return None

  def __str__(self):
    return str(self.sequence)

  def ToNFA(self):
    """Return the NFA representation of this regular expression."""
    nfa = NFA()
    prev_node = nfa.AddNode()
    next_node = prev_node
    nfa.start_nodes.add(prev_node)
    for item in self.sequence:
      next_node = nfa.AddNode()
      nfa.AddTransition(prev_node, next_node, item)
      prev_node = next_node
    nfa.accept_nodes.add(next_node)
    return nfa


class REConcat(RE):
  """Class representing concatenated regular expressions (e.g. "abc")."""

  def __init__(self, sub_trees):
    self.sub_trees = list(sub_trees)

  def Match(self, sequence):
    for sub_tree in self.sub_trees:
      m = sub_tree.Match(sequence)
      if m is None:
        return None
      sequence = sequence[m:]

  def __str__(self):
    return '({})'.format(''.join(str(st) for st in self.sub_trees))

  def ToNFA(self):
    if not self.sub_trees:
      # The concatenation of nothing is equal to the empty sequence.
      nfa = NFA()
      node = nfa.AddNode()
      nfa.start_nodes.add(node)
      nfa.accept_nodes.add(node)
      return nfa
    else:
      nfa = self.sub_trees[0].ToNFA()
      for sub_tree in itertools.islice(self.sub_trees, 1, len(self.sub_trees)):
        sub_nfa = sub_tree.ToNFA()
        sub_nfa.EnsureDisjoint(nfa)
        nfa.UpdateNodes(sub_nfa)
        connecting_node = nfa.AddNode()
        for node in nfa.accept_nodes:
          nfa.AddTransition(node, connecting_node, None)
        for node in sub_nfa.start_nodes:
          nfa.AddTransition(connecting_node, node, None)
        nfa.accept_nodes = sub_nfa.accept_nodes
      return nfa


class REOr(RE):
  """Class representing a union regular expression (e.g. "a|b")."""

  def __init__(self, sub_trees):
    self.sub_trees = list(sub_trees)

  def Match(self, sequence):
    m = None
    for sub_tree in self.sub_trees:
      m = sub_tree.Match(sequence)
      if m is not None:
        break
    return m

  def __str__(self):
    return '({})'.format('|'.join(str(st) for st in self.sub_trees))

  def ToNFA(self):
    assert self.sub_trees
    nfa = self.sub_trees[0].ToNFA()
    for sub_tree in itertools.islice(self.sub_trees, 1, len(self.sub_trees)):
      sub_nfa = sub_tree.ToNFA()
      sub_nfa.EnsureDisjoint(nfa)
      nfa.UpdateNodes(sub_nfa)
      nfa.start_nodes.update(sub_nfa.start_nodes)
      nfa.accept_nodes.update(sub_nfa.accept_nodes)
    return nfa


class RERepeat(RE):
  """Class representing a repeating regular expression (e.g. "a*")."""

  def __init__(self, sub_tree):
    self.sub_tree = sub_tree

  def Match(self, sequence):
    total = 0
    m = self.sub_tree.Match(sequence)
    while m is not None:
      total += m
      sequence = sequence[m:]
      m = self.sub_tree.Match(sequence)
    return total

  def __str__(self):
    return '({})*'.format(self.sub_tree)

  def ToNFA(self):
    nfa = self.sub_tree.ToNFA()
    connecting_node = nfa.AddNode()
    for accept_node in nfa.accept_nodes:
      nfa.AddTransition(accept_node, connecting_node, None)
    for start_node in nfa.start_nodes:
      nfa.AddTransition(connecting_node, start_node, None)
    nfa.accept_nodes.update(nfa.start_nodes)
    return nfa


def FindPartitionSeq(
    nfa, max_len, target_ratio=fractions.Fraction(1, 2),
    lo=(), hi=None, tolerance_ratio=0):
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
  lo = SeqToNum(lo, nfa.inverse_alphabet, max_len)
  hi = SeqToNum(
      (max_letter for _ in xrange(max_len)) if hi is None else hi,
      nfa.inverse_alphabet, max_len)
  lo_num_accepts = nfa.NumAccepts(max_len, NumToSeq(lo, nfa.alphabet, max_len))
  hi_num_accepts = nfa.NumAccepts(max_len, NumToSeq(hi, nfa.alphabet, max_len))
  total = lo_num_accepts - hi_num_accepts
  assert total >= 0
  target = lo_num_accepts - int(total * target_ratio)
  tolerance = int(total * tolerance_ratio)
  for j in itertools.cycle((0, 1)):
    if j == 0 and hi_num_accepts != lo_num_accepts:
      mid = (lo + (hi - lo) * (lo_num_accepts - target) /
             (lo_num_accepts - hi_num_accepts))
      assert mid >= lo
      assert mid <= hi
    else:
      mid = (lo + hi) / 2
    mid_str = NumToSeq(mid, nfa.alphabet, max_len)
    mid_num_accepts = nfa.NumAccepts(max_len, mid_str)
    if (lo >= hi or mid_num_accepts == target
        or abs(lo_num_accepts - hi_num_accepts) <= tolerance):
      break
    elif mid_num_accepts < target:
      hi = mid + 1
      hi_num_accepts = mid_num_accepts
    elif mid_num_accepts > target:
      lo = mid
      lo_num_accepts = mid_num_accepts
  return tuple(NumToSeq(mid, nfa.alphabet, max_len))


def FindPartitionSeqs(
    nfa, max_len, num_partitions=1, lo=(), hi=None, tolerance_ratio=0):
  return tuple(
      FindPartitionSeq(
          nfa, max_len,
          fractions.Fraction(j, num_partitions), lo, hi, tolerance_ratio)
      for j in xrange(1, num_partitions))


MAX_REPEAT = re.sre_parse.parse('a*')[0][1][1]
ALPHABET = frozenset(string.printable)


def SreToRE(sre, alphabet=ALPHABET):
  if isinstance(sre, tuple):
    op, args = sre
    if op == 'branch':  # or
      none, sres = args
      assert none is None  # I don't know how to interpret a non-None value.
      res = map(SreToRE, sres)
      return REOr(res)
    elif op == 'in':  # or
      sres = args
      res = map(SreToRE, sres)
      return REOr(res)
    elif op == 'category':  # or (\d, \W)
      if args == 'category_digit':
        return REOr([RESequence(c) for c in string.digits if c in alphabet])
      elif args == 'category_not_digit':
        return REOr([RESequence(d) for d in alphabet if d not in string.digits])
      elif args == 'category_word':
        return REOr([RESequence(c) for c in string.letters if c in alphabet])
      elif args == 'category_not_word':
        return REOr([
            RESequence(d) for d in alphabet if d not in string.letters])
      else:
        raise Exception('Unknown category type "{}".'.format(args))
    elif op == 'assert':  # (?=REGEX)
      _, sres = args
      assert len(sres) == 2
      return SreToRE(sres[1])
    elif op == 'any':  # .
      return REOr([RESequence(c) for c in alphabet])
    elif op == 'max_repeat' or op == 'min_repeat':  # {m,n}, *, +
      min_repeat, max_repeat, sre = args
      res = []
      main_res = SreToRE(sre)
      if max_repeat >= MAX_REPEAT:
        for _ in xrange(min_repeat):
          res.append(copy.deepcopy(main_res))
        res.append(RERepeat(copy.deepcopy(main_res)))
      else:
        for _ in xrange(min_repeat):
          res.append(copy.deepcopy(main_res))
        for _ in xrange(max_repeat - min_repeat):
          res.append(REOr([copy.deepcopy(main_res), RESequence('')]))
      if len(res) > 1:
        return REConcat(res)
      else:
        return res[0]
      return RERepeat(res)
    elif op == 'literal':  # A single character.
      return RESequence(chr(args))
    elif op == 'not_literal':  # Any character,  except this single character.
      literal = chr(args)
      return REOr([RESequence(c) for c in alphabet if c != literal])
    elif op == 'subpattern':
      _, sres = args
      return SreToRE(sres)
    elif op == 'range':
      min_literal, max_literal = args
      return REOr([RESequence(chr(c))
                   for c in xrange(min_literal, max_literal + 1)])
    elif op == 'at':  # Anchor characters (e.g. "^" and "$")
      return RESequence('')
    else:
      raise Exception('Unknown op, args pair: {}, {}'.format(op, args))
  else:  # concat
    if len(sre) == 1:
      return SreToRE(sre[0])
    else:
      res = map(SreToRE, sre)
      return REConcat(res)


def RegexStrToRE(regex_str):
  sre = re.sre_parse.parse(regex_str)
  return SreToRE(sre)


def main(argv):
  num_partitions = 4 if len(argv) <= 1 else int(argv[1])
  regex_str = '(foo_[a-z]+|[a-f0-9]*)'
  regex_str = regex_str if len(argv) <= 2 else argv[2]
  lo = () if len(argv) <= 3 else argv[3]
  hi = None if len(argv) <= 4 or argv[4] == 'None' else argv[4]
  tolerance_ratio = 0.25 / num_partitions if len(argv) <= 5 else float(argv[5])
  max_len = 10 if len(argv) <= 6 else int(argv[6])

  t = RegexStrToRE(regex_str)
  nfa = t.ToNFA()

  print '\n'.join(
      ''.join(x) for x in FindPartitionSeqs(
          nfa, max_len, num_partitions, lo, hi, tolerance_ratio))


if __name__ == '__main__':
  main(sys.argv)
