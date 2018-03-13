import sys
sys.path.append('../')

import unittest
import multigroup_utilities


class TestMultigroupUtilities(unittest.TestCase):

    def setUp(self):
        pass

    def test_energy_groups(self):
        '''
        Test that the structures have the correct number of groups in the correct order
        '''

        structures = ['wims69', 'wims56', 'wims172', 'lwr32', 'lwr28', 'phoenix25',
                      'scale44', 'scale56', 'scale238', 'scale252', 'shem281',
                      'shem361', 'shem407', 'ga193', 'ga537', 'hr6',
                      'hr16', 'casmo2', 'casmo3', 'casmo4', 'casmo7', 'casmo8',
                      'casmo9', 'casmo12', 'casmo14', 'casmo16', 'casmo18',
                      'casmo23', 'casmo25', 'casmo40', 'casmo70', 'eurlib100',
                      'ecco1968']

        for struc in structures:
            nG = int(''.join([l for l in struc if l.isdigit()]))
            eb = multigroup_utilities.energy_groups(struc)
            assert len(eb) == nG + 1, '{} has the wrong number of groups'.format(struc)
            assert all(eb[1:] < eb[:-1]), '{} has groups out of order'.format(struc)

    def test_energy_groups_lethargy(self):
        '''
        Test that the lethargy structure produces the correct groups
        '''
        nG = 10
        eb = multigroup_utilities.energy_groups('lethargy', groups=nG, upper=1e7)

        assert len(eb) == nG + 525, 'lethargy has the wrong number of groups'
        assert all(eb[1:] < eb[:-1]), 'lethargy has groups out of order'

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
