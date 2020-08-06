#!/usr/bin/python3
""" Boyer-Moore String Matching. """

from typing import List, Tuple
import pprint
import time


def z_array(string: str) -> List[int]:
    """ Use Z algorithm to preprocess s.

    `Z[i]` is length of longest substring of `P` that starts at `i` and matches a prefix
    of `P`.
    """

    assert len(string) > 1, "Length of string must be greater than 1!"
    z_arr = [len(string)] + [0] * (len(string) - 1)

    # Initial comparison of s[1:] with prefix
    for i in range(1, len(string)):
        if string[i] != string[i - 1]:
            break

        z_arr[1] += 1

    right, left = 0, 0
    if z_arr[1] > 0:
        right, left = z_arr[1], 1

    for k in range(2, len(string)):
        assert z_arr[k] == 0

        # Case 1
        if k > right:
            for i in range(k, len(string)):
                if string[i] == string[i - k]:
                    z_arr[k] += 1
                else:
                    break
            right, left = k + z_arr[k] - 1, k

        # Calculate length of beta
        elif (right - k + 1) > z_arr[k - left]:
            z_arr[k] = z_arr[k - left]

        # Compare characters just past r
        else:
            match = 0
            for i in range(right + 1, len(string)):
                if string[i] != string[i - k]:
                    break
                match += 1
            left, right = k, right + match
            z_arr[k] = right - k + 1

    return z_arr


def n_array(string: str) -> List[int]:
    """ Compile the N array from the Z array.

    `N[i]` is the length of longest suffix of `P[:i]` which is also a substring of `P`.
    """

    return z_array(string[::-1])[::-1]


def big_l_prime_array(pattern: str, n_arr: List[int]) -> List[int]:
    """ Compile L' array using p and N array.

    `L'[i]` = largest index `j < m` such that `N[j] = |P[i:]|`. """

    l_prime = [0] * len(pattern)

    for j in range(len(pattern) - 1):
        i = len(pattern) - n_arr[j]
        if i < len(pattern):
            l_prime[i] = j + 1

    return l_prime


def big_l_array(pattern: str, l_prime_arr: List[int]) -> List[int]:
    """ Compile L array using p and L' array.

    `L[i]` = largest index `j < m` such that `N[j] >= |P[i:]|`. """

    l_arr = [0] * len(pattern)
    l_arr[1] = l_prime_arr[1]

    for i in range(2, len(pattern)):
        l_arr[i] = max(l_arr[i - 1], l_prime_arr[i])

    return l_arr


def small_l_prime_array(n_arr: List[int]) -> List[int]:
    """ Compile l' array using N array.

    `l'[i]`  = largest `j <= m - i` such that `N[j] = j`. """

    small_l_prime_arr = [0] * len(n_arr)

    for i, _ in enumerate(n_arr):
        if n_arr[i] == i + 1:  # Prefix matching a suffix
            small_l_prime_arr[len(n_arr) - i - 1] = i + 1

    for i in range(len(n_arr) - 2, -1, -1):  # "Smear" them out to the left
        if small_l_prime_arr[i] == 0:
            small_l_prime_arr[i] = small_l_prime_arr[i + 1]

    return small_l_prime_arr


def good_suffix_table(pattern: str) -> Tuple[List[int], List[int], List[int]]:
    """ Return tables needed to apply good suffix rule. """

    n_arr = n_array(pattern)
    l_prime_arr = big_l_prime_array(pattern, n_arr)
    return l_prime_arr, big_l_array(pattern, l_prime_arr), small_l_prime_array(n_arr)


def dense_bad_char_table(pattern: str, amap: dict) -> List[List[int]]:
    """ Given pattern string and list with ordered alphabet characters, create and
    return a dense bad character table. Table is indexed by offset then by character.
    """

    table = []
    nxt = [0] * len(amap)

    for i, _ in enumerate(pattern):
        char = pattern[i]
        assert char in amap, f"{char} not found in alphabet!"
        table.append(nxt[:])
        nxt[amap[char]] = i + 1

    return table


class BoyerMoore:
    """ Encapsulates pattern and associated Boyer-Moore preprocessing. """

    def __init__(self, pattern: str, alphabet: str):
        """ Initialize data structures. """

        # Create a map from alphabet characters to integers
        self.amap = {alphabet[i]: i for i in range(len(alphabet))}

        # Create bad character rule table
        self.bad_char = dense_bad_char_table(pattern, self.amap)

        # Create good suffix rule table
        _, self.big_l, self.small_l_prime = good_suffix_table(pattern)

    def bad_character_rule(self, offset: int, char: str) -> int:
        """ Return number of skips given by bad character rule at offset. """

        assert char in self.amap, f"{char} not found in alphabet!"
        assert offset < len(self.bad_char), f"Invalid offset: {offset}"
        index = self.amap[char]

        return offset - (self.bad_char[offset][index] - 1)

    def good_suffix_rule(self, offset: int) -> int:
        """ Given a mismatch at offset, return amount to shift as determined by good
        suffix table. """

        length = len(self.big_l)
        assert offset < length, f"Invalid offset: {offset}"

        if offset == length - 1:
            return 0

        offset += 1  # offset points to the leftmost matching position of p

        if self.big_l[offset] > 0:
            return length - self.big_l[offset]

        return self.small_l_prime[offset]

    def match_skip(self):
        """ Return amount to shift in case where P matches T. """

        return len(self.small_l_prime) - self.small_l_prime[1]


def visualize(
    text: str, pattern: str, t_off: int, p_off: int, match: bool, sleep_time: float
):
    """ Print text and pattern with appropriate colors. """

    green_bold = "\033[1;38;2;0;255;0m"
    red_bold = "\033[1;38;2;255;0;0m"
    grey = "\033[38;2;127;127;127m"
    end_format = "\033[0m"
    arrow_head = "\U000025C0"
    arrow_body = "\U0001F89C"

    for i in range(len(pattern) - 1, p_off - 1, -1):
        if match or i > p_off:
            print(
                text[: t_off + i]
                + f"{green_bold}{pattern[i:]}{end_format}"
                + text[t_off + len(pattern) :],
                flush=True,
            )
            print(
                " " * t_off + f"{pattern[:i]}{green_bold}{pattern[i:]}{end_format}",
                flush=True,
            )

        else:
            print(
                text[: t_off + p_off]
                + f"{red_bold}{text[t_off + p_off]}{end_format}"
                + f"{green_bold}{pattern[p_off + 1:]}{end_format}"
                + text[t_off + len(pattern) :],
                flush=True,
            )
            print(
                " " * t_off
                + pattern[:p_off]
                + f"{red_bold}{pattern[p_off]}{end_format}"
                + f"{green_bold}{pattern[p_off + 1:]}{end_format}",
                flush=True,
            )

        print(
            " " * (t_off + i)
            + f"{grey}"
            + arrow_head
            + arrow_body * (len(pattern) - i - 1)
            + f"{end_format}",
            flush=True,
        )

        time.sleep(sleep_time)
        if i > p_off:
            print("\033[3F", end="", flush=True)

    print()


def boyer_moore(
    pattern: str, bm_obj: BoyerMoore, text: str, sleep_time: float,
) -> Tuple[List[int], int, int]:
    """ Do Boyer-Moore matching with visualization """

    i = 0
    occurrences = []
    alignments = 0
    comparisons = 0

    while i < len(text) - len(pattern) + 1:
        shift = 1
        mismatched = False
        alignments += 1

        for j in range(len(pattern) - 1, -1, -1):
            comparisons += 1
            skip_bc = 0
            skip_gs = 0

            if pattern[j] != text[i + j]:
                skip_bc = bm_obj.bad_character_rule(j, text[i + j])
                skip_gs = bm_obj.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break

        if not mismatched:
            occurrences.append(i)
            skip_gs = bm_obj.match_skip()
            shift = max(shift, skip_gs)

        visualize(text, pattern, i, j, not mismatched, sleep_time)
        print(f"Comparisons: {comparisons}")

        if i < len(text) - len(pattern):
            if skip_bc:
                print(f"Bad character shift: {skip_bc}")

            if skip_gs:
                print(f"Good suffix shift: {skip_gs}")

            print("Press Enter to continue ...", end="\r", flush=True)
            input()
            print("\033[F", end="\033[K")

        print()
        i += shift

    return (occurrences, alignments, comparisons)


def main():
    """ Driver function. """

    alphabet = "abcdefghijklmnopqrstuvwxyz "

    # Time (in seconds) to sleep between character comparisons during visualization
    sleep_time = 0.25

    text = input("Enter text   : ")
    pattern = input("Enter pattern: ")
    bm_obj = BoyerMoore(pattern, alphabet)
    print()

    # Write data structures to file
    fout = open("data_structures.txt", "w")

    print(f"TEXT: {text}", file=fout)
    print(f"PATTERN: {pattern}", file=fout)

    print("\nAlphabet map:", file=fout)
    pprint.pprint(bm_obj.amap, stream=fout, compact=True)

    print("\nBad character table:", file=fout)
    pprint.pprint(bm_obj.bad_char, stream=fout, compact=True, width=85)

    print("\nL array: ", file=fout)
    pprint.pprint(bm_obj.big_l, stream=fout, compact=True)

    print("\nl' array: ", file=fout)
    pprint.pprint(bm_obj.small_l_prime, stream=fout, compact=True)

    fout.close()

    (occurrences, alignments, comparisons) = boyer_moore(
        pattern, bm_obj, text, sleep_time
    )
    print(
        f"Text length: {len(text)}",
        f"Pattern length: {len(pattern)}",
        f"Occurrences: {occurrences}",
        f"Alignments: {alignments}",
        f"Comparisons: {comparisons}",
        sep="\n",
    )


if __name__ == "__main__":
    main()
