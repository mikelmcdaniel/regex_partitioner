import itertools
import unittest

import regex_partitioner


class Test(unittest.TestCase):
    def test_seq_num(self):
        for alphabet in ("ab", "abc", "abcd", "abcde", "abcdef"):
            inverse_alphabet = dict((v, k) for k, v in enumerate(alphabet))

            def all_seqs(alphabet, max_len):
                def all_seqs_helper(alphabet, max_len, seq):
                    yield seq
                    if max_len > 0:
                        for _ in all_seqs_helper(alphabet, max_len - 1, seq):
                            for letter in alphabet:
                                seq.append(letter)
                                yield seq
                                seq.pop()

                for s in all_seqs_helper(alphabet, max_len, []):
                    yield "".join(s)

            for max_len in range(7):
                all_sequences = sorted(all_seqs(alphabet, max_len))
                for expected_num, expected_str in enumerate(all_sequences):
                    actual_num = regex_partitioner.seq_to_num(expected_str, inverse_alphabet, max_len)
                    actual_str = "".join(regex_partitioner.num_to_seq(expected_num, alphabet, max_len))
                    self.assertEqual(expected_num, actual_num)
                    self.assertEqual(expected_str, actual_str)

    def test_nfa(self):
        nfa = regex_partitioner.NFA()
        # Make an NFA that accepts base-10 numbers divisible by 3
        # Note: We consider the string "" to be the same as "0"
        zero = nfa.add_node()
        one = nfa.add_node()
        two = nfa.add_node()
        nfa.start_nodes.add(zero)
        nfa.accept_nodes.add(zero)
        nodes = [zero, one, two]
        for start_val, start_node in enumerate(nodes):
            for next_digit in range(10):
                new_val = (start_val * 10 + next_digit) % 3
                new_node = nodes[new_val]
                nfa.add_transition(start_node, new_node, str(next_digit))

        self.assertTrue(nfa.accepts(""))
        self.assertTrue(nfa.accepts("0"))
        self.assertTrue(nfa.accepts("36"))
        self.assertFalse(nfa.accepts("1"))
        self.assertFalse(nfa.accepts("2"))
        self.assertFalse(nfa.accepts("3001"))
        for n in range(0, 1000, 3):
            self.assertTrue(nfa.accepts(str(n)))
            self.assertFalse(nfa.accepts(str(n + 1)))
            self.assertFalse(nfa.accepts(str(n + 2)))

        self.assertEqual(nfa.num_accepts(1), 5)  # "", "0", "3", "6", "9
        self.assertEqual(nfa.num_accepts(1, "555"), 2)  # "6", "9"

        max_len = 3
        all_strs = set()
        for str_len in range(max_len + 1):
            all_strs.update("".join(s) for s in itertools.product("0123456789", repeat=str_len))
        div_by_3_strs = {s for s in all_strs if s == "" or int(s) % 3 == 0}
        for bound in sorted(all_strs):
            expected_num_accepts = sum(s >= bound for s in div_by_3_strs)
            self.assertEqual(expected_num_accepts, nfa.num_accepts(max_len, bound), f"{bound=}")


if __name__ == "__main__":
    unittest.main()
