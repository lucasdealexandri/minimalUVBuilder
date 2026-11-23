#!/usr/bin/env python3

"""
Unit tests for U1Rep, SU2Rep and SU3Rep using unittest style.
Expanded to include U1Rep and some edge cases.
"""

import unittest
from uvbuilder.tools.reps import U1Rep, SU2Rep, SU3Rep


def _to_tuples_su2(lst):
    return sorted([(x.dynkin,) for x in lst])


def _to_tuples_su3(lst):
    return sorted([(x.upper, x.lower) for x in lst])


class TestU1Rep(unittest.TestCase):
    def test_constructor_and_equality(self):
        for value in [0, 1, -1, 2, -3, 7, -13, 42]:
            rep = U1Rep(value)
            self.assertIsInstance(rep, U1Rep)
            self.assertEqual(rep, U1Rep(value))
            self.assertNotEqual(rep, U1Rep(value + 1))

    def test_zero_identity(self):
        zero = U1Rep(0)
        for value in range(-10, 11):
            rep = U1Rep(value)
            self.assertEqual(rep + zero, U1Rep(value))
            self.assertEqual(zero + rep, U1Rep(value))
            self.assertEqual(rep - zero, U1Rep(value))
            self.assertEqual(rep * zero, U1Rep(0))
            self.assertEqual(zero * rep, U1Rep(0))

    def test_addition_and_subtraction(self):
        for a in range(-10, 11):
            for b in range(-10, 11):
                self.assertEqual(U1Rep(a) + U1Rep(b), U1Rep(a + b))
                self.assertEqual(U1Rep(a) - U1Rep(b), U1Rep(a - b))
                self.assertEqual(U1Rep(a) + U1Rep(b), U1Rep(b) + U1Rep(a))

    def test_multiplication(self):
        for a in range(-6, 7):
            for b in range(-6, 7):
                self.assertEqual(U1Rep(a) * U1Rep(b), U1Rep(a * b))
                self.assertEqual(U1Rep(a) * U1Rep(b), U1Rep(b) * U1Rep(a))

    def test_fractional_construction_and_operations(self):
        a = U1Rep(1, 2)
        b = U1Rep(1, 3)
        self.assertEqual(a + b, U1Rep(5, 6))
        self.assertEqual(a * b, U1Rep(1, 6))
        self.assertEqual(a - b, U1Rep(1, 6))

    def test_self_operations(self):
        for value in range(-8, 9):
            rep = U1Rep(value)
            self.assertEqual(rep - rep, U1Rep(0))
            self.assertEqual(rep + (-rep), U1Rep(0))
            self.assertEqual(U1Rep(1) * rep, U1Rep(value))
            self.assertEqual(rep * U1Rep(1), U1Rep(value))


class TestSU2Rep(unittest.TestCase):
    def test_constructor_and_equality(self):
        for value in range(0, 8):
            rep = SU2Rep(value)
            self.assertIsInstance(rep, SU2Rep)
            self.assertEqual(rep, SU2Rep(value))
            self.assertNotEqual(rep, SU2Rep(value + 1))

    def test_dim_and_i3_range(self):
        for n in range(0, 8):
            rep = SU2Rep(n)
            self.assertEqual(rep.dim, n + 1)
            # i3_range: dynkin+1 physical i3 values, -j ... j in steps of 1 (U1Rep(k, 2) = k/2)
            i3 = rep.i3_range()
            self.assertEqual(len(i3), n + 1)
            self.assertEqual(i3[0], U1Rep(-n, 2))
            self.assertEqual(i3[-1], U1Rep(n, 2))
            # reverse should be reversed
            self.assertEqual(rep.i3_range(reverse=True), tuple(reversed(i3)))

    def test_additional_i3_range_examples(self):
        self.assertEqual(SU2Rep(2).i3_range(), (U1Rep(-1), U1Rep(0), U1Rep(1)))
        self.assertEqual(SU2Rep(2).i3_range(reverse=True), (U1Rep(1), U1Rep(0), U1Rep(-1)))
        self.assertEqual(SU2Rep(1).i3_range(), (U1Rep(-1, 2), U1Rep(1, 2)))

    def test_multiplication_and_commutativity(self):
        # Known decompositions
        a = SU2Rep(1)
        b = SU2Rep(1)
        expected = [SU2Rep(0), SU2Rep(2)]
        self.assertEqual(_to_tuples_su2(a * b), _to_tuples_su2(expected))
        self.assertEqual(_to_tuples_su2(b * a), _to_tuples_su2(expected))
        # commutativity for a few values
        for x in range(0, 6):
            for y in range(0, 6):
                left = _to_tuples_su2(SU2Rep(x) * SU2Rep(y))
                right = _to_tuples_su2(SU2Rep(y) * SU2Rep(x))
                self.assertEqual(left, right)

    def test_larger_dynkin_products(self):
        # product of dynkin 5 and 3 gives dynkin values 2,4,6,8
        res = SU2Rep(5) * SU2Rep(3)
        self.assertEqual(_to_tuples_su2(res), [(2,), (4,), (6,), (8,)])


class TestSU3Rep(unittest.TestCase):
    def test_constructor_and_equality(self):
        choices = [(0, 0), (1, 0), (0, 1), (1, 1), (2, 0), (2, 1)]
        for u, l in choices:
            rep = SU3Rep(u, l)
            self.assertIsInstance(rep, SU3Rep)
            self.assertEqual(rep, SU3Rep(u, l))
            # inequality with swapped indices when different
            if u != l:
                self.assertNotEqual(rep, SU3Rep(l, u))

    def test_dim_and_inv_and_repr(self):
        self.assertEqual(SU3Rep(0, 0).dim, 1)
        self.assertEqual(SU3Rep(1, 0).dim, 3)
        self.assertEqual(SU3Rep(0, 1).dim, 3)
        # inv swaps upper and lower
        self.assertEqual(SU3Rep(0, 1).inv(), SU3Rep(1, 0))
        self.assertEqual(SU3Rep(1, 1).inv(), SU3Rep(1, 1))
        # double inv gives original
        self.assertEqual(SU3Rep(2, 1).inv().inv(), SU3Rep(2, 1))
        # repr formatting
        self.assertEqual(repr(SU3Rep(0, 1)), "3*")
        self.assertEqual(repr(SU3Rep(1, 0)), "3")

    def test_multiplication_and_commutativity(self):
        # Known decompositions (order-insensitive)
        a = SU3Rep(1, 0)
        expected_a_a = [SU3Rep(2, 0), SU3Rep(0, 1)]
        self.assertEqual(_to_tuples_su3(a * a), _to_tuples_su3(expected_a_a))
        self.assertEqual(_to_tuples_su3(a * a), _to_tuples_su3(a.__rmul__(a)))

        # 3* x 3 -> 8 + 1 (order-insensitive)
        three = SU3Rep(1, 0)
        threebar = SU3Rep(0, 1)
        expected = [SU3Rep(1, 1), SU3Rep(0, 0)]
        self.assertEqual(_to_tuples_su3(threebar * three), _to_tuples_su3(expected))
        self.assertEqual(_to_tuples_su3(three * threebar), _to_tuples_su3(expected))

        # commutativity for a few examples
        samples = [SU3Rep(0, 0), SU3Rep(1, 0), SU3Rep(0, 1), SU3Rep(1, 1), SU3Rep(2, 1)]
        for x in samples:
            for y in samples:
                left = _to_tuples_su3(x * y)
                right = _to_tuples_su3(y * x)
                self.assertEqual(left, right)


if __name__ == "__main__":
    unittest.main(argv=["first-arg-is-ignored"], exit=False)
