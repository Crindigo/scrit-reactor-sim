import math

def magnitude(vector: list[float]) -> float:
    mag = sum([v * v for v in vector])
    return math.sqrt(mag)

def normalize(vector: list[float]) -> list[float]:
    mag = magnitude(vector)
    return vector if mag == 0 else [v / mag for v in vector]

def linear_normalize(vector: list[float]) -> list[float]:
    vsum = sum(vector)
    return [v / vsum for v in vector]

def multiply(matrix: list[list[float]], vector: list[float]) -> list[float]:
    result = []
    for i, row in enumerate(matrix):
        for j, col in enumerate(row):
            result[i] += col * vector[j]
    return result
