import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Simple.Plot.Lib.')
    parser.add_argument("-i", "--i", help="input file (.csv) (default: %(default)s).")
    args = parser.parse_args()

    # new hydro class object






read_file = pd.read_csv(r'{}'.format(args.i),delimiter = ',')

x = read_file.iloc[:,0]
u = read_file.iloc[:,1]

plt.plot(x,u)
plt.xlabel('x')
plt.ylabel('u')
plt.show()
plt.savefig('test_xu{}.png'.format(args.i),dpi = 300)