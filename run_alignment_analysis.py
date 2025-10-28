from pprint import pprint
from time import time

from pathlib import Path

from alignment import align


def _analyze_alignment(N: int, algorithm, **kwargs):
    def read_sequence(file: Path) -> str:
        return ''.join(file.read_text().splitlines())

    seq1 = read_sequence(Path('test_files/bovine_coronavirus.txt'))[:N]
    seq2 = read_sequence(Path('test_files/murine_hepatitus.txt'))[:N]
    start = time()
    algorithm(seq1, seq2, **kwargs)
    end = time()
    runtime = end - start
    return N, runtime


def _compute_average_runtimes(runtimes):
    groups = {}
    for n, runtime in runtimes:
        key = n
        if key not in groups:
            groups[key] = []
        groups[key].append((n, runtime))

    return [
        (
            size,
            round(sum(t for _, t in stats) / len(stats), 3)
        )
        for size, stats in groups.items()
    ]


def _print_markdown_table(ave_runtimes, headers):
    header_widths = [len(header) for header in headers]

    rows = [
        '| ' + ' | '.join(headers) + ' |',
        '| ' + ' | '.join('-' * len(header) for header in headers) + ' |'
    ]

    rows += (
        '| ' + ' | '.join(
            f'{field:<{width}}'
            for field, width in zip(row, header_widths)
        ) + ' |'
        for row in ave_runtimes
    )

    print('\n'.join(rows))


def main(sizes, algorithm, file_name="_runtimes.py", **kwargs):

    runtimes = []
    for size in sizes:
        print('Running with size', size)
        for iteration in range(10):
            n, runtime = _analyze_alignment(size, algorithm, **kwargs)
            runtimes.append((n, runtime))

    ave_runtimes = _compute_average_runtimes(runtimes)

    print()
    print('Copy this markdown table into your report:  ')
    print()

    print()
    _print_markdown_table(
        ave_runtimes,
        [' N     ', 'Time (sec)']
    )

    # Print runtimes to a file
    with open(file_name, 'w') as file:
        print('runtimes = ', end='', file=file)
        pprint(runtimes, file)

    print()
    print(f'{file_name} written')


if __name__ == '__main__':
    sizes = [500, 1000, 1500, 2000, 2500, 3000]

    main(sizes=sizes,
         algorithm=align,
         banded_width = 5)
