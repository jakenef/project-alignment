import matplotlib.pyplot as plt
from _runtimes import runtimes


def compute_coefficient(observed_performance, theoretical_order):
    return [
        time / theoretical_order(n) for n, time in observed_performance
    ]


def main():
    def theoretical_big_o(n):
        return n

    coeffs = compute_coefficient(runtimes, theoretical_big_o)

    # slice this list to use a subset for your estimate
    used_coeffs = coeffs[0:]

    coeff = sum(used_coeffs) / len(used_coeffs)
    print(coeff)

    plt.bar(range(len(coeffs)), coeffs)
    xlim = plt.xlim()
    plt.plot(xlim, [coeff, coeff], ls=':', c='k')
    plt.xlim(xlim)
    plt.title(f'coeff={coeff}')
    plt.show()


if __name__ == '__main__':
    main()
