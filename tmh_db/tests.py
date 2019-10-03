import unittest
from scripts.populate_tmh import Sort
from scripts.populate_tmh import odd_or_even


# run `./manage.py test`

# Create your tests here.
class TestStuff(unittest.TestCase):

    def test_sort(self):
        aim_result=[["a",1],["b",2],["c",3]]

        tricky_result = [["a",1],["b",2],["c",3]]
        observed_result=Sort(tricky_result)
        self.assertEqual(aim_result, observed_result)

        tricky_result = [["b",2],["c",3],["a",1]]
        observed_result=Sort(tricky_result)
        self.assertEqual(aim_result, observed_result)

        tricky_result = [["b",1],["a",2],["c",3]]
        observed_result=Sort(tricky_result)
        self.assertNotEqual(aim_result, observed_result)

    def test_odd_or_even(self):
        #number=1
        self.assertEqual(odd_or_even(1), "odd")
        self.assertEqual(odd_or_even(2), "even")
        self.assertEqual(odd_or_even(0), "even")



# allows to be run by python test_build
if __name__ == '__main__':
    unittest.main()
