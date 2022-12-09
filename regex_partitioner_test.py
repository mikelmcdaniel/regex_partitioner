import fractions
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

        self.assertEqual(None, nfa.prev_accepted("", 3))
        self.assertEqual("0", nfa.prev_accepted("00", 3))
        self.assertEqual("987", nfa.prev_accepted("99", 3))
        self.assertEqual("993", nfa.prev_accepted("996", 3))
        self.assertEqual("996", nfa.prev_accepted("997", 3))
        self.assertEqual("996", nfa.prev_accepted("998", 3))
        self.assertEqual("996", nfa.prev_accepted("999", 3))
        self.assertEqual("", nfa.prev_accepted("0", 3))
        self.assertEqual("03", nfa.prev_accepted("030", 3))

        self.assertEqual("0", nfa.next_accepted("", 3))
        self.assertEqual(None, nfa.next_accepted("999", 3))
        self.assertEqual("6", nfa.next_accepted("5", 1))
        self.assertEqual("501", nfa.next_accepted("5", 3))
        self.assertEqual("504", nfa.next_accepted("502", 3))
        self.assertEqual("9", nfa.next_accepted("6", 1))
        self.assertEqual("60", nfa.next_accepted("6", 3))
        self.assertEqual("600", nfa.next_accepted("60", 3))
        self.assertEqual("600", nfa.next_accepted("60", 3))

        div_by_3_strs_iter = iter(sorted(div_by_3_strs))
        prev = next(div_by_3_strs_iter)
        for cur in div_by_3_strs_iter:
            self.assertEqual(prev, nfa.prev_accepted(cur, 3))
            self.assertEqual(cur, nfa.next_accepted(prev, 3))
            prev = cur

    def test_nfa_num_accepts_early_exit(self):
        nfa = regex_partitioner.regex_str_to_re("apple|banana|coconut").as_nfa()
        self.assertEqual(3, nfa.num_accepts(max_len=2 ** 31))

    def test_regex_str_to_re(self):
        nfa = regex_partitioner.regex_str_to_re("(apple)*|foo[0-9]{2,4}").as_nfa()
        self.assertTrue(nfa.accepts(""))
        self.assertTrue(nfa.accepts("apple"))
        self.assertTrue(nfa.accepts("appleapple"))
        self.assertTrue(nfa.accepts("appleappleapple"))
        self.assertTrue(nfa.accepts("foo01"))
        self.assertTrue(nfa.accepts("foo3210"))

        self.assertFalse(nfa.accepts("app"))
        self.assertFalse(nfa.accepts("applex"))
        self.assertFalse(nfa.accepts("foo"))
        self.assertFalse(nfa.accepts("foo0"))
        self.assertFalse(nfa.accepts("foo1"))
        self.assertFalse(nfa.accepts("foo12345"))

    def test_find_partition_seq(self):
        nfa = regex_partitioner.regex_str_to_re("[0-9]{4}").as_nfa()
        self.assertTrue(nfa.accepts("0000"))
        self.assertTrue(nfa.accepts("1234"))
        self.assertTrue(nfa.accepts("9876"))
        self.assertFalse(nfa.accepts(""))
        self.assertFalse(nfa.accepts("01"))
        self.assertFalse(nfa.accepts("123"))
        self.assertFalse(nfa.accepts("12345"))
        self.assertFalse(nfa.accepts("-1"))
        self.assertFalse(nfa.accepts("1.2"))
        self.assertEqual(10000, nfa.num_accepts(4))
        self.assertEqual("", "".join(regex_partitioner.find_partition_seq(nfa, max_len=4, target_ratio=0)))
        self.assertEqual("0000", "".join(regex_partitioner.find_partition_seq(nfa, max_len=4, target_ratio=0.0001)))
        self.assertEqual("9999", "".join(regex_partitioner.find_partition_seq(nfa, max_len=4, target_ratio=1)))
        self.assertEqual("4999", "".join(regex_partitioner.find_partition_seq(nfa, max_len=4, target_ratio=0.5)))
        self.assertEqual("3332", "".join(regex_partitioner.find_partition_seq(nfa, max_len=4, target_ratio=1 / 3)))

        for j in [1, 10000, 5000, 3333, 29, 4001, 7549]:
            self.assertEqual(
                f"{j - 1:04d}",
                "".join(regex_partitioner.find_partition_seq(nfa, 4, fractions.Fraction(j, 10000))),
            )


if __name__ == "__main__":
    unittest.main()
