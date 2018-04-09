import argparse
import os


def contacts(n,d): # cd for uniform distribution
    if d % 2 == 1:
        return 0
    answer = (1.0 / (d/2 +1)) * (binomial(d, d/2) * binomial(2*n - d -1, n - d/2 -1))
    if d <= 10:
        print('{} -> {}'.format(d, answer))
    return answer


def num_leaves(n, k):
    return (1.0 / n) * binomial(n, k) * binomial(n, k - 1)


def root_deg(n, r):
    return (r / (n * 1.0)) * binomial(2*n - 1 - r, n - 1)


def min_height(n, h):
    """
    returns the number of trees of size n with height >= h
    """
    h += 1
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

    for name, func, div, k_min, k_max in (
        ('num_leaves', num_leaves, Cn, 1, n+1),
        ('root_degree', root_deg, Cn, 1, n+1),
        ('height', tree_height, Cn, 1, n+1),
        ('cd_sums', contacts, Cn * n, 0, 2 * n - 1)):

        filename = '{}_n={}_cdf.txt'.format(name, n)
        filename = os.path.join(out_dir, filename)
        with open(filename, 'w') as f:
            total = 0
            for k in range(k_min, k_max):
                total += func(n, k) / div
                f.write('{}\n'.format(total.n(digits=14)))
