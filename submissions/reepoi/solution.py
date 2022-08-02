import csv
import json
import numpy as np


def solution():
    r"""
    This is a test solution.
    """
    parameters = 0.5 + np.zeros((32, 6))
    header = 'Index,photocurrent,saturation_current,resistance_series,resistance_shunt,n,cells_in_series'.split(',')

    for name in ['case1.csv', 'case2.csv']:
        with open(name, 'w') as file:
            writer = csv.writer(file, delimiter=',')
            writer.writerow(header)
            for idx, row in enumerate(parameters):
                writer.writerow([1 + idx] + list(row))

if __name__ == '__main__':
    solution()

