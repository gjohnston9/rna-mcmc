import argparse
import os


def num_leaves(n, k):
    return (1.0 / n) * binomial(n, k) * binomial(n, k - 1)


def root_deg(n, r):
    return (r / (n * 1.0)) * binomial(2*n - 1 - r, n - 1)


def min_height(n, h):
    """
    returns the number of trees of size n with height >= h
    """
    total = 0
    for k in range(1, ((n+1)/h) + 3):
        # print(total)
        total += binomial(2*n, n + 1 - k*h) - 2*binomial(2*n, n - k*h) + binomial(2*n, n - 1 - k*h)
    return total


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('n', type=int, help='value of n to use')
    args = parser.parse_args()

    n = args.n
    Cn = catalan_number(n)

    out_dir = 'expectations'

    min_heights = {h : min_height(n, h) for h in range(1, n+3)}
    def tree_height(n, h, min_heights=min_heights):
        """
        returns the number of trees of size n with height exactly h
        """
        return min_heights[h] - min_heights[h+1]

    for name, func in (
        ('num_leaves', num_leaves), ### non-zero for 1:n
        ('root_degree', root_deg), ### non-zero for 1:n
        ('height', tree_height)): ### non-zero for 2:n+1

        filename = '{}_n={}_cdf.txt'.format(name, n)
        filename = os.path.join(out_dir, filename)
        with open(filename, 'w') as f:
            f.write('0\n')
            total = 0
            for k in range(1, n+2):
                total += func(n, k) / Cn
                f.write('{}\n'.format(total.n(digits=14)))
