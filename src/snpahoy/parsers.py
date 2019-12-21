from typing import Iterable, List

from snpahoy.core import Position


def get_positions(stream: Iterable[str]) -> List[Position]:
    result: List[Position] = []
    for line in stream:
        chromosome, start, stop = line.split('\t')
        for coordinate in range(int(start), int(stop)):
            result.append(Position(chromosome=chromosome, coordinate=coordinate))
    return result
