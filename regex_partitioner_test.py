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


if __name__ == "__main__":
    unittest.main()
